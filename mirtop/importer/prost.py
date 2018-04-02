""" Read prost! files"""

import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam
from collections import defaultdict

from mirtop.mirna import mapper
from mirtop.mirna.fasta import read_precursor
from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits, make_id, get_mature_sequence, align
from mirtop.bam import filter

logger = mylog.getLogger(__name__)

def header():
    return ""

def read_file(fn, hairpins, database, mirna_gtf):
    """
    read bam file and perform realignment of hits
    """
    reads = defaultdict(hits)
    sample = os.path.splitext(os.path.basename(fn))[0]
    genomics = mapper.read_gtf_to_mirna(mirna_gtf)
    matures = mapper.read_gtf_to_precursor(mirna_gtf)
    non_mirna = 0
    non_chromosome_mirna = 0
    outside_mirna = 0
    lines_read = 0
    ann, ann_type = _group_seqs_by_ann(fn)
    with open(fn) as handle:
        handle.readline()
        for line in handle:
            lines_read += 1
            cols = line.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if not ann[query_sequence]:
                non_mirna += 1
                continue
            miRNA = ann_type[ann[query_sequence]][1]
            preNames = ann_type[ann[query_sequence]][0]
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            reads[query_name].set_sequence(query_sequence)
            reads[query_name].counts = cols[9]
            for preName in preNames.split(","):
                if preName in reads[query_name].precursors:
                    continue
                if preName not in hairpins:
                    non_chromosome_mirna += 1
                    continue
                reference_start = _align_to_mature(query_sequence, hairpins[preName], matures[preName][miRNA])
                logger.debug("\nPROST!::NEW::query: {query_sequence}\n"
                             "  precursor {preName}\n"
                             "  name:  {query_name}\n"
                             "  reference_start: {reference_start}\n"
                             "  mirna: {miRNA}".format(**locals()))
                iso = isomir()
                iso.align = line
                iso.set_pos(reference_start, len(reads[query_name].sequence))
                logger.debug("PROST!:: start %s end %s" % (iso.start, iso.end))
                if len(hairpins[preName]) < reference_start + len(reads[query_name].sequence):
                    continue
                iso.subs, iso.add, iso.cigar = filter.tune(reads[query_name].sequence,
                                                           hairpins[preName],
                                                           reference_start, None)
                logger.debug("PROST!::After tune start %s end %s" % (iso.start, iso.end))
                if len(iso.subs) < 2:
                    reads[query_name].set_precursor(preName, iso)
    logger.info("Lines loaded: %s" % lines_read)
    logger.info("Skipped lines because non miRNA in line: %s" % non_mirna)
    logger.info("Skipped lines because non chromosome in GTF: %s" % non_chromosome_mirna)
    logger.info("Skipped lines because outside precursor: %s" % outside_mirna)
    logger.info("Hits: %s" % len(reads))
    return reads

def _group_seqs_by_ann(fn):
    """Read file once to group sequences to same miRNA sequence"""
    ann = dict()
    ann_type = defaultdict(list)
    with open(fn) as inh:
        inh.readline()
        for line in inh:
            cols = line.strip().split("\t")
            ann[cols[0]] = cols[4]
            mirna = cols[11] if cols[11] else cols[13]
            hairpin = cols[15]
            if not cols[4] in ann_type:
                ann_type[cols[4]] = ["", ""]
            if mirna:
                ann_type[cols[4]][1] = mirna
            if hairpin:
                ann_type[cols[4]][0] = hairpin
    return [ann, ann_type]

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

def _align_to_mature(seq, hairpin, mature):
    """Get alignment between seq and mature"""
    mirna = get_mature_sequence(hairpin, mature)
    hit =  align(seq, mirna)
    start = hit[0][:8].count("-") - 4 + int(mature[0])
    logger.debug("PROST::align:sequence to mature %s" % hit[0])
    logger.debug("PROST::align:start: %s -> %s" % (mature[0], start))
    return start

def _cigar_to_variants(seq, mature, cigar):
    """From mature based cigar get variants"""
    return None

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

