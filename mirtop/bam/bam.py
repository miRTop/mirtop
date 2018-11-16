""" Read bam files"""

import os.path as op
import os
import pysam
from collections import defaultdict

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits
from mirtop.bam import filter

logger = mylog.getLogger(__name__)


def read_bam(bam_fn, args, clean=True):
    """
    Read bam file and perform realignment of hits

    Args:
        *bam_fn*: a BAM file with alignments to the precursor

        *precursors*: dict with keys being precursor names and values
            being sequences. Come from mirtop.mirna.fasta.read_precursor().

        *clean*: Use mirtop.filter.clean_hits() to remove lower score hits.

    Returns:
        *reads (dict)*:
             keys are read_id and values are *mirtop.realign.hits*

    """
    precursors = args.precursors
    bam_fn = _sam_to_bam(bam_fn)
    bam_fn = _bam_sort(bam_fn)
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    reads = defaultdict(hits)
    indels_skip = 0
    for line in handle:
        if line.reference_id < 0:
            logger.debug("READ::Sequence not mapped: %s" % line.reference_id)
            continue
        query_name = line.query_name
        if query_name not in reads and not line.query_sequence:
            continue
        logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, line.query_sequence))
        if line.query_sequence and line.query_sequence.find("N") > -1:
            continue
        if query_name not in reads:
            reads[query_name].set_sequence(line.query_sequence)
            reads[query_name].counts = _get_freq(query_name)
        if line.is_reverse:
            logger.debug("READ::Sequence is reverse: %s" % line.query_name)
            continue
        chrom = handle.getrname(line.reference_id)
        cigar = line.cigartuples
        if line.cigarstring.find("I") > -1:
            indels_skip += 1
            continue
        iso = isomir()
        iso.align = line
        iso.set_pos(line.reference_start, len(reads[query_name].sequence))
        logger.debug("READ::From BAM start %s end %s at chrom %s" % (iso.start, iso.end, chrom))
        if len(precursors[chrom]) < line.reference_start + len(reads[query_name].sequence):
            logger.debug("READ::%s start + %s sequence size are bigger than"
                         " size precursor %s" % (
                                                 line.reference_id,
                                                 len(reads[query_name].sequence),
                                                 len(precursors[chrom])))
            continue
        iso.subs, iso.add, iso.cigar = filter.tune(
            reads[query_name].sequence, precursors[chrom],
            line.reference_start, cigar)
        logger.debug("READ::After tune start %s end %s" % (iso.start, iso.end))
        logger.debug("READ::iso add %s iso subs %s" % (iso.add, iso.subs))

        reads[query_name].set_precursor(chrom, iso)
    logger.info("Hits: %s" % len(reads))
    logger.info("Hits skipped due to contain indels %s" % indels_skip)
    if clean:
        reads = filter.clean_hits(reads)
        logger.info("Hits after clean: %s" % len(reads))
    return reads


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
        do.run(("samtools sort -n -o {bam_sort_by_n} {bam_fn}").format(
            **locals()))
    return bam_sort_by_n


def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = int(name.split("_x")[1])
    except:
        return 0
    return counts
