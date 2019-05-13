""" Read prost! files"""

import os
from collections import defaultdict

from mirtop.mirna import mapper
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits, get_mature_sequence, align
from mirtop.bam import filter

logger = mylog.getLogger(__name__)


def header():
    """
    Custom header for PROST! importer.

    Returns:
        *(str)*: PROST! header string.
    """

    return ""


def read_file(fn, hairpins, database, mirna_gtf):
    """
    Read PROST! file and convert to mirtop GFF format.

    Args:
        *fn(str)*: file name with PROST output information.

        *database(str)*: database name.

        *args(namedtuple)*: arguments from command line.
            See *mirtop.libs.parse.add_subparser_gff()*.

    Returns:
        *reads*: dictionary where keys are read_id and values are *mirtop.realign.hits*

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
            if query_name not in reads and not query_sequence:
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
                iso.subs, iso.add, iso.cigar = filter.tune(
                    reads[query_name].sequence,
                    hairpins[preName],
                    reference_start, None)
                logger.debug("PROST!::After tune start %s end %s" % (
                    iso.start, iso.end))
                if len(iso.subs) < 2:
                    reads[query_name].set_precursor(preName, iso)
    logger.info("Lines loaded: %s" % lines_read)
    logger.info("Skipped lines because non miRNA in line: %s" % non_mirna)
    logger.info("Skipped lines because non chromosome in GTF:"
                " %s" % non_chromosome_mirna)
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


def _align_to_mature(seq, hairpin, mature):
    """Get alignment between seq and mature"""
    mirna = get_mature_sequence(hairpin, mature)
    hit = align(seq, mirna)
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
        variant.append("iso_5p:%s" % -1 * int(cols[0]))
    if cols[1] != "0":
        variant.append("iso_3p:%s" % cols[1])
    if cols[2] != "0":
        variant.append("iso_add3p:%s" % int(cols[2]))
    if cols[3] == "True":
        variant.append("iso_snv_seed")
    if cols[4] == "True":
        variant.append("iso_snv_central_offset")
    if cols[5] == "True":
        variant.append("iso_snv_central")
    if cols[6] == "True":
        variant.append("iso_snv_supp")
    if cols[7] == "True":
        variant.append("iso_snv")
    if not variant:
        return "NA"
    return ",".join(variant)
