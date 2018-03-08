""" Read isomiR GFF files"""

import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam
from collections import defaultdict, Counter

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna import mapper
from mirtop.mirna.realign import isomir, hits, expand_cigar, make_id
from mirtop.gff.body import read_attributes
from mirtop.bam import filter
from mirtop.importer.prost import genomic2transcript

logger = mylog.getLogger(__name__)

def header():
    h = ""
    return h

def read_file(fn, database, gtf):
    """
    read bam file and perform realignment of hits
    """
    map_mir = mapper.read_gtf_to_mirna(gtf)
    reads = defaultdict(dict)
    reads_in = 0
    sample = os.path.splitext(os.path.basename(gtf))[0]
    hits = _get_hits(fn)
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            attr = read_attributes(line, "=")
            query_name = attr['TS']
            query_sequence = attr['TS'].replace("U", "T")
            start = int(cols[3])
            end = int(cols[4])
            isomirseq_iso = attr['ISO']
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            counts = attr["TC"]
            chrom = cols[0]
            # logger.debug("SEQBUSTER:: cigar {cigar}".format(**locals()))
            cigar = attr['CI'].replace("U", "T")
            idu = make_id(query_sequence)
            isoformat = cigar2variants(cigar, query_sequence, attr['ISO'])
            logger.debug("\nSOMIRSEA::NEW::query: {query_sequence}\n"
                         "  precursor {chrom}\n"
                         "  name: {query_name}\n"
                         "  idu: {idu}\n"
                         "  start: {start}\n"
                         "  cigar: {cigar}\n"
                         "  iso: {isoformat}\n"
                         "  variant: {isoformat}".format(**locals()))
            source = "isomiR" if isoformat != "NA" else "ref_miRNA"
            strand = "+"
            database = cols[1]
            mirName = attr['MIN'].split()[0]
            preName = attr['PIN'].split()[0]
            score = "."
            Filter = attr['FILTER']
            isotag = attr['ISO']
            tchrom, tstart =  genomic2transcript(map_mir[mirName], chrom, start)
            start = start if not tstart else tstart
            chrom = chrom if not tstart else tchrom
            end = start + len(query_sequence)
            hit = hits[idu]
            attrb = ("Read {query_sequence}; UID {idu}; Name {mirName}; Parent {preName}; Variant {isoformat}; Isocode {isotag}; Cigar {cigar}; Expression {counts}; Filter {Filter}; Hits {hit};").format(**locals())
            res = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t{score}\t{strand}\t.\t{attrb}").format(**locals())
            if start not in reads[chrom]:
                reads[chrom][start] = []
            if Filter == "Pass":
                reads_in += 1
                reads[chrom][start].append([idu, chrom, counts, sample, res])

    logger.info("Hits: %s" % reads_in)
    return reads

def _get_hits(fn):
    hits = Counter()
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            attr = read_attributes(line, "=")
            query_sequence = attr['TS'].replace("U", "T")
            if query_sequence and query_sequence.find("N") > -1:
                continue
            idu = make_id(query_sequence)
            hits[idu] += 1
    return hits

def cigar2variants(cigar, sequence, tag):
    """From cigar to Variants in GFF format"""
    pos = 0
    iso5p = 0
    logger.debug("\nISOMIRSEA:: expanded: %s" % expand_cigar(cigar))
    for l in expand_cigar(cigar):
        if l == "I":
            iso5p += 1
        elif l == "D":
            iso5p -= 1
        else:
            break
    iso3p = 0
    for l in reversed(expand_cigar(cigar)):
        if l == "I":
            iso3p += 1
        elif l == "D":
            iso3p -= 1
        else:
            break
    isosnp = []
    for l in expand_cigar(cigar):
        if l in ['A', 'T', 'C', 'G' ]:
            isosnp.append([pos, sequence[pos], l])
        if l in ['D']:
            continue
        pos += 1
    iso5p = "iso_5p:%s" % _fix(iso5p) if iso5p else ""
    if tag[-1] == "T" or iso3p < 0:
        iso3p = "iso_3p:%s" % _fix(iso3p) if iso3p else ""
    else:
        iso3p = "iso_add:%s" % _fix(iso3p) if iso3p else ""

    variant = ""
    for iso in [iso5p, iso3p, _define_snp(isosnp)]:
        if iso:
            variant += "%s," % iso

    variant = "NA;" if not variant else variant
    return variant[:-1]

def _define_snp(subs):
    value = ""
    logger.debug("\nISOMIRSEA:: subs %s" % subs)
    for sub in subs:
        if sub:
            if sub[0] > 1 and sub[0] < 8:
                value += "iso_snp_seed,"
            elif sub[0] == 8:
                value += "iso_snp_central_offset,"
            elif sub[0] > 8 and sub[0] < 13:
                value += "iso_snp_central,"
            elif sub[0] > 12 and sub[0] < 18:
                value += "iso_snp_central_supp,"
            else:
                value += "iso_snp,"
    return value[:-1]

def _fix(n):
    if n > 0:
        return "+%s" % n
    return n
