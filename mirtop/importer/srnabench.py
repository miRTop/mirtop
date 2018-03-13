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


def read_file(folder,  precursors):
    """
    read srnabench file and perform realignment of hits
    """
    n_out = 0
    n_nonmature = 0
    n_ns = 0
    n_in = 0
    n_non_precursor = 0
    reads_anno = os.path.join(folder, "reads.annotation")
    reads_iso = os.path.join(folder, "microRNAannotation.txt")
    reads = defaultdict(hits)
    source_iso = _read_iso(reads_iso)
    logger.info("Reads with isomiR information %s" % len(source_iso))
    with open(reads_anno) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                n_ns += 1
                continue
            if cols[3].find("mature") == -1:
                n_nonmature += 1
                continue
            if query_name not in reads:
                reads[query_name].set_sequence(query_sequence)
                reads[query_name].counts = int(cols[1])

            for hit in cols[4].split("$"):
                logger.debug("SRNABENCH::line hit: %s" % hit)
                hit_info = hit.split("#")
                pos_info = hit_info[3].split(",")
                reference_start = int(pos_info[1]) - 1
                chrom = pos_info[0]
                iso = isomir()
                iso.align = line
                if (query_sequence, hit_info[1]) in source_iso:
                    iso.external = source_iso[(query_sequence, hit_info[1])]
                external = iso.external
                logger.debug("SRNABENCH::query: {query_sequence}\n"
                             "  precursor {chrom}\n"
                             "  name:  {query_name}\n"
                             "  start: {reference_start}\n"
                             "  external: {external}\n"
                             "  hit: {hit}".format(**locals()))
                iso.set_pos(reference_start, len(reads[query_name].sequence))
                logger.debug("SRNABENCH:: start %s end %s" % (iso.start, iso.end))
                if len(precursors[chrom]) < reference_start + len(reads[query_name].sequence):
                    n_out += 1
                    continue
                iso.subs, iso.add, iso.cigar = filter.tune(reads[query_name].sequence,
                                                           precursors[chrom],
                                                           reference_start, None)
                logger.debug("SRNABENCH::After tune start %s end %s" % (iso.start, iso.end))
                n_in += 1
                reads[query_name].set_precursor(chrom, iso)
            if len(reads[query_name].precursors) == 0:
                n_non_precursor += 1
    logger.info("Loaded %s reads with %s hits" % (len(reads), n_in))
    logger.info("Reads without precursor information: %s" % n_non_precursor)
    logger.info("Hit Filtered by having > 3 changes: %s" % n_out)
    logger.info("Hit Filtered by being non-mature: %s" % n_nonmature)
    return reads

def _read_iso(fn):
    """
    Read definitions of isomiRs by srnabench
    """
    iso = dict()
    with open(fn) as inh:
        h = inh.readline()
        for line in inh:
            cols = line.strip().split("\t")
            label = cols[3].split("$")
            mirnas = cols[1].split("$")
            if len(mirnas) != len(label):
                label.extend([label[-1]] * (len(mirnas) - len(label)))
            anno = dict(zip(mirnas, label))
            logger.debug("TRANSLATE::%s with %s" % (mirnas, label))
            for m in anno:
                iso[(cols[0], m)] = _translate(anno[m], cols[4])
                logger.debug("TRANSLATE::code %s" % iso[(cols[0], m)])
    return iso

def _translate(label, description):
    iso = []
    if label == "exact":
        return "NA"
    if label == "mv":
        return "notsure"
    number_nts = label.split("|")[-1].split("#")[-1]
    if number_nts.find("-") < 0:
        number_nts = "+%s" % number_nts
    if label.find("lv3p") > -1:
        iso.append("iso_3p:%s" % number_nts)
    if label.find("lv5p") > -1:
        iso.append("iso_5p:%s" % number_nts)
    if label.find("nta") > -1:
        iso.append("iso_add:%s" % number_nts)
    if label.find("NucVar") > -1:
        for nt in description.split(","):
            logger.debug("TRANSLATE::change:%s" % description)
            if nt == "-":
                return "notsure"
            iso.extend(_iso_snp(int(nt.split(":")[0])))
    return ",".join(iso)

def _iso_snp(pos):
    iso = []
    if pos > 1 and pos < 8:
        iso.append("iso_snp_seed")
    elif pos == 8:
        iso.append("iso_snp_central_offset")
    elif pos > 8 and pos < 13:
        iso.append("iso_snp+central")
    elif pos > 12 and pos < 18:
        iso.append("iso_snp_central_supp")
    else:
        iso.append("iso_snp")
    return iso
