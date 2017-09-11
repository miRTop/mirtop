""" Read bam files"""

import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam
from collections import defaultdict

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits
from mirtop.bam import filter

logger = mylog.getLogger(__name__)

def _sam_to_bam(bam_fn):
    if not bam_fn.endswith("bam"):
        bam_out = "%s.bam" % os.path.splitext(bam_fn)[0]
        cmd = "samtools view -Sbh {bam_fn} -o {bam_out}"
        do.run(cmd.format(**locals()))
        return bam_out
    return bam_fn

def _bam_sort(bam_fn):
    bam_sort_by_n = op.splitext(bam_fn)[0] + "_sort.bam"
    if not file_exists(bam_sort_by_n):
        do.run(("samtools sort -n -o {bam_sort_by_n} {bam_fn}").format(**locals()))
    return bam_sort_by_n

def read_bam(bam_fn, precursors):
    """
    read bam file and perform realignment of hits
    """
    bam_fn = _sam_to_bam(bam_fn)
    bam_fn = _bam_sort(bam_fn)
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    reads = defaultdict(hits)
    for line in handle:
        if line.reference_id < 0:
            logger.debug("Sequence not mapped: %s" % line.reference_id)
            continue
        query_name = line.query_name
        if query_name not in reads:
            reads[query_name].sequence = line.query_sequence
        if line.is_reverse:
            logger.debug("Sequence is reverse: %s" % line.query_name)
            continue
        chrom = handle.getrname(line.reference_id)
        #  print "%s %s %s %s" % (line.query_name, line.reference_start, line.query_sequence, chrom)
        cigar = line.cigartuples
        iso = isomir()
        iso.align = line
        iso.start = line.reference_start
        if len(precursors[chrom]) < line.reference_start + len(reads[query_name].sequence):
            continue
        iso.subs, iso.add, iso.end = filter.tune(reads[query_name].sequence, precursors[chrom], line.reference_start, cigar)
        if len(iso.subs) < 2:
            reads[query_name].set_precursor(chrom, iso)

    reads = filter.clean_hits(reads)
    return reads

def _coord(sequence, start, mirna, precursor, iso):
    """
    Define t5 and t3 isomirs
    """
    insertion = 0
    if iso.subs:
        insertion = 1 if iso.subs[0][-1] == "-" else 0
    end = start + (iso.end - len(iso.add) - insertion) - 1
    logger.debug("coor:s:%s len:%s flen:%s end:%s mirna:%s iso:%s" % (start, len(sequence), iso.end, end, mirna, iso.format()))
    dif = abs(mirna[0] - start)
    if start < mirna[0]:
        iso.t5 = sequence[:dif].upper()
    elif start > mirna[0]:
        iso.t5 = precursor[mirna[0] - 1:mirna[0] - 1 + dif].lower()
    elif start == mirna[0]:
        iso.t5 = "NA"
    if dif > 4:
        logger.debug("start > 3 %s %s %s %s %s" % (start, len(sequence), dif, mirna, iso.format()))
        return None

    dif = abs(mirna[1] - end)
    if iso.add:
        iso.add = iso.add.replace("-", "")
        sequence = sequence[:-len(iso.add)]
    # if dif > 3:
    #    return None
    if end > mirna[1]:
        iso.t3 = sequence[-dif:].upper()
    elif end < mirna[1]:
        iso.t3 = precursor[mirna[1] - dif:mirna[1]].lower()
    elif end == mirna[1]:
        iso.t3 = "NA"
    if dif > 4:
        logger.debug("end > 3 %s %s %s %s %s" % (len(sequence), end, dif, mirna, iso.format()))
        return None
    # logger.debug("coor end:%s %s %s %s %s iso:%s" % (start, len(sequence), end, dif, mirna, iso.format()))
    return True

def annotate(reads, mirbase_ref, precursors):
    """
    Using SAM/BAM coordinates, mismatches and realign to annotate isomiRs

    reads: dict object that comes from read_bam fn
    mirbase_ref: dict object that comers from mirtop.mirna.read_mature
    precursors: dict object (key : fasta) that comes from mirtop.mirna.fasta.read_precursor
    """
    for r in reads:
        for p in reads[r].precursors:
            start = reads[r].precursors[p].start + 1  # convert to 1base
            end = start + len(reads[r].sequence)
            for mature in mirbase_ref[p]:
                mi = mirbase_ref[p][mature]
                is_iso = _coord(reads[r].sequence, start, mi, precursors[p], reads[r].precursors[p])
                logger.debug(("read:{s} pre:{p} start:{start} is_iso:{is_iso} mir:{mature} mir_p:{mi} mir_s:{mature_s}").format(s=reads[r].sequence, mature_s=precursors[p][mi[0]-1:mi[1]], **locals()))
                logger.debug("annotation:%s iso:%s" % (r, reads[r].precursors[p].format()))
                if is_iso:
                    reads[r].precursors[p].mirna = mature
                    break
    return reads
