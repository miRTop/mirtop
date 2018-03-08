""" Read prost files"""

import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam
from collections import defaultdict

from mirtop.mirna import mapper
from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits, make_id
from mirtop.bam import filter

logger = mylog.getLogger(__name__)

def header():
    return ""

def read_file(fn, precursors, database, mirna_gtf):
    """
    read bam file and perform realignment of hits
    """
    reads = defaultdict(dict)
    sample = os.path.splitext(os.path.basename(fn))[0]
    map_mir = mapper.read_gtf_to_mirna(mirna_gtf)
    non_mirna = 0
    non_chromosome_mirna = 0
    outside_mirna = 0
    lines_read = 0
    with open(fn) as handle:
        handle.readline()
        for line in handle:
            lines_read += 1
            cols = line.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if len(cols) < 12:
                non_mirna += 1
                continue
            miRNA = cols[11]
            if not miRNA:
                if cols[13]:
                    miRNA = cols[13]
                elif cols[15]:
                    miRNA = cols[15]
                else:
                    continue
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            for loc in cols[5].split(";")[:1]:
                if loc.find("-") < 0:
                    non_chromosome_mirna += 1
                    continue
                chrom = loc.split(":")[0]
                start, end = loc.split(":")[1].split("-")
                preName, reference_start =  genomic2transcript(map_mir[miRNA], chrom, int(start))
                if not chrom:
                    non_chromosome_mirna += 1
                    continue
                # reference_start = int(cols[4]) - 1
                logger.debug("\nPROST::NEW::query: {query_sequence}\n"
                             "  precursor {chrom}\n"
                             "  name:  {query_name}\n"
                             "  start: {start}\n"
                             "  reference_start: {reference_start}\n"
                             "  mirna: {miRNA}".format(**locals()))
                Filter = "PASS"
                hit = "NA"
                isoformat = _make_variant(cols[19:])
                idu = make_id(query_sequence)
                strand = "."
                counts = cols[9]
                cigar = "NA"
                score = "."
                source = "isomiR" if isoformat != "NA" else "ref_miRNA"
                attrb = ("Read {query_sequence}; UID {idu}; Name {miRNA}; Parent {preName}; Variant {isoformat}; Cigar {cigar}; Expression {counts}; Filter {Filter}; Hits {hit};").format(**locals())
                res = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t{score}\t{strand}\t.\t{attrb}").format(**locals())
                if start not in reads[chrom]:
                    reads[chrom][start] = []
                reads[chrom][start].append([idu, chrom, counts, sample, res])

    logger.info("Lines loaded: %s" % lines_read)
    logger.info("Skipped lines because non miRNA in line: %s" % non_mirna)
    logger.info("Skipped lines because non chromosome in GTF: %s" % non_chromosome_mirna)
    logger.info("Skipped lines because outside precursor: %s" % outside_mirna)
    logger.info("Hits: %s" % len(reads))
    return reads

def genomic2transcript(code, chrom, pos):
    for ref in code:
        if _is_chrom(chrom, code[ref][0]):
            if _is_inside(pos, code[ref][1:3]):
                return [ref, _transcript(pos, code[ref][1:4])]
    return [None, None]

def _is_chrom(chrom, annotated):
    logger.debug("TRANSCRIPT::CHROM::read position %s and db position %s" % (chrom, annotated))
    if chrom == annotated:
        return True
    if chrom == annotated.replace("chr", ""):
        return True
    return False

def _is_inside(pos, annotated):
    logger.debug("TRANSCRIPT::INSIDE::read position %s and db position %s" % (pos, annotated))
    if pos > annotated[0] and pos < annotated[1]:
        return True
    return False

def _transcript(pos, annotated):
    logger.debug("TRANSCRIPT::TRANSCRIPT::read position %s and db position %s" % (pos, annotated))
    if annotated[2] == "+":
        return pos - annotated[0]
    elif annotated[2] == "-":
        return annotated[1] - pos
    raise ValueError("Strand information is incorrect %s" % annotated[3])

def _make_variant(cols):
    logger.debug("PROST::variant: %s" % cols)
    variant = []
    if cols[0] != "0":
        variant.append("iso_5p:%s" % cols[0])
    if cols[1] != "0":
        variant.append("iso_3p:%s" % cols[1])
    if cols[2] != "0":
        variant.append("iso_add:%s" % cols[2])
    if cols[3] == "True":
        variant.append("iso_snp_seed")
    if cols[4] == "True":
        variant.append("iso_snp_central_offset")
    if cols[5] == "True":
        variant.append("iso_snp_central")
    if cols[6] == "True":
        variant.append("iso_snp_supp")
    if cols[7] == "True":
        variant.append("iso_snp")
    if not variant:
        return "NA"
    return ",".join(variant)

