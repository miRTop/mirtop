""" Read seqbuster files"""

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

def header():
    h = ("## CMD: seqbuster: http://seqcluster.readthedocs.io/mirna_annotation.html#mirna-isomir-annotation-with-java\n"
        "# iso_snp are not filtered yet. Use isomiRs R pacakge to correct for error sequencing\n")
    return h

def read_file(fn, precursors):
    """
    read bam file and perform realignment of hits
    """
    reads = defaultdict(hits)
    with open(fn) as handle:
        handle.readline()
        for line in handle:
            cols = line.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            if query_name not in reads:
                reads[query_name].set_sequence(query_sequence)
                reads[query_name].counts = int(col[9])
            for loc in cols[5].split(";"):
                chrom = loc.split(":")[0]
                start, end = loc.split(":")[1].split("-")
                # map start to reference start
                miRNA = cols[11]
                reference_start = int(cols[4]) - 1
                seqbuster_iso = ":".join(cols[6:10])
                logger.debug("\nPROST::NEW::query: {query_sequence}\n"
                             "  precursor {chrom}\n"
                             "  name:  {query_name}\n"
                             "  start: {reference_start}\n"
                             "  iso: {seqbuster_iso}".format(**locals()))
                # logger.debug("PROST:: cigar {cigar}".format(**locals()))
                iso = isomir()
                iso.align = line
                iso.set_pos(reference_start, len(reads[query_name].sequence))
                logger.debug("PROST:: start %s end %s" % (iso.start, iso.end))
                if len(precursors[chrom]) < reference_start + len(reads[query_name].sequence):
                    continue
                iso.subs, iso.add, iso.cigar = filter.tune(reads[query_name].sequence,
                                                           precursors[chrom],
                                                           reference_start, None)
                logger.debug("PROST::After tune start %s end %s" % (iso.start, iso.end))
                if len(iso.subs) < 2:
                    reads[query_name].set_precursor(chrom, iso)
    logger.info("Hits: %s" % len(reads))
    return reads
