""" Read sRNAbench files"""

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


def read_file(fn, precursors):
    """
    read srnabench file and perform realignment of hits
    """
    reads = defaultdict(hits)
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            if cols[3].find("mature") == -1:
                continue
            if query_name not in reads:
                reads[query_name].set_sequence(query_sequence)
                reads[query_name].counts = _get_freq(int(cols[1]))

            for hit in cols[4].split("$"):
                logger.debug("SRNABENCH::line hit: %s" % hit)
                hit_info = hit.split("#")
                pos_info = hit_info[3].split(",")
                reference_start = int(pos_info[1]) - 1
                chrom = pos_info[0]
                logger.debug("SRNABENCH::query: {query_sequence}\n"
                             "  precursor {chrom}\n"
                             "  name:  {query_name}\n"
                             "  start: {reference_start}\n"
                             "  hit: {hit}".format(**locals()))
                iso = isomir()
                iso.align = line
                iso.set_pos(reference_start, len(reads[query_name].sequence))
                logger.debug("SRNABENCH:: start %s end %s" % (iso.start, iso.end))
                if len(precursors[chrom]) < reference_start + len(reads[query_name].sequence):
                    continue
                iso.subs, iso.add, iso.cigar = filter.tune(reads[query_name].sequence,
                                                           precursors[chrom],
                                                           reference_start, None)
                logger.debug("SRNABENCH::After tune start %s end %s" % (iso.start, iso.end))
                if len(iso.subs) < 3:
                    reads[query_name].set_precursor(chrom, iso)
    logger.info("Hits: %s" % len(reads))
    return reads

def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = int(name.split("_x")[1])
    except:
        return 0
    return counts

