""" Read bam files"""

import os.path as op
import os
import pysam
from collections import defaultdict

import pybedtools

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna.realign import isomir, hits, reverse_complement
from mirtop.bam import filter
from mirtop.mirna.mapper import liftover_genomic_precursor, read_gtf_chr2mirna

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
    if args.genomic:
        chrom2mirna = read_gtf_chr2mirna(args.gtf)
    for line in handle:
        if line.reference_id < 0:
            logger.debug("READ::Sequence not mapped: %s" % line.reference_id)
            continue
        if not line.cigarstring:
            logger.debug("READ::Sequence malformed: %s" % line)
            continue
        query_name = line.query_name
        if query_name not in reads and not line.query_sequence:
            continue
        sequence = line.query_sequence if not line.is_reverse else reverse_complement(line.query_sequence)
        logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, sequence))
        if line.query_sequence and line.query_sequence.find("N") > -1:
            continue
        if query_name not in reads:
            reads[query_name].set_sequence(sequence)
            reads[query_name].counts = _get_freq(query_name)
        if line.is_reverse and not args.genomic:
            logger.debug("READ::Sequence is reverse: %s" % line.query_name)
            continue
        chrom = handle.getrname(line.reference_id)
        start = line.reference_start
        # If genomic endcode, liftover to precursor position
        if args.genomic:
            chrom, start = _liftover(chrom, start + 1, line.reference_end,
                                     line.is_reverse, chrom2mirna)
        if not start:
            continue
        cigar = line.cigartuples
        if line.cigarstring.find("I") > -1:
            indels_skip += 1
        iso = isomir()
        iso.align = line
        iso.set_pos(start, len(reads[query_name].sequence))
        logger.debug("READ::From BAM start %s end %s at chrom %s" % (iso.start, iso.end, chrom))
        if len(precursors[chrom]) < start + len(reads[query_name].sequence):
            logger.debug("READ::%s start + %s sequence size are bigger than"
                         " size precursor %s" % (
                                                 line.reference_id,
                                                 len(reads[query_name].sequence),
                                                 len(precursors[chrom])))
            continue
        iso.subs, iso.add, iso.cigar = filter.tune(
            reads[query_name].sequence, precursors[chrom],
            start, cigar)
        logger.debug("READ::After tune start %s end %s" % (iso.start, iso.end))
        logger.debug("READ::iso add %s iso subs %s" % (iso.add, iso.subs))

        reads[query_name].set_precursor(chrom, iso)
    logger.info("Hits: %s" % len(reads))
    logger.info("Hits with indels %s" % indels_skip)
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


def _liftover(chr, start, end, is_reverse, chrom2mirna):
    strand = "+" if not is_reverse else "-"
    bed = _bed_with_mirna_in_chrom(chr, chrom2mirna)
    hit = "\t".join([chr, str(start), str(end), strand])
    hit = pybedtools.BedTool(hit, from_string=True)
    logger.debug("BAM::liftover:hit: %s" % hit)
    overlap = hit.intersect(bed, wo=True)
    logger.debug("BAM::liftover:overlap: %s" % overlap)
    # check only one overlap
    logger.debug("BAM::length overlap:%s" % len(overlap))
    ## read start and strand
    ## miRNA Name
    for align in overlap:
        print(align[8])
        read = {'start': int(align[1]), 'end': int(align[2]),
                'strand': align[3]}
        genomic = {'chrom': align[4], 'start': int(align[5]),
                   'end': int(align[6]),
                   'strand': align[7]}
        hairpin = {'start': int(align[8].split(":")[-1]),
                   'chrom': align[8].split(":")[0]}
        hairpin_position = liftover_genomic_precursor(read,
                                                      genomic,
                                                      hairpin)
        logger.debug("BAM::lifted:: %s" % [hairpin['chrom'], hairpin_position])
        return [hairpin['chrom'], hairpin_position]
    return [None, None]


def _bed_with_mirna_in_chrom(chr, chrom2mirna):
    lines = []
    if chr in chrom2mirna:
        for position in chrom2mirna[chr]:
            line = "\t".join([chr, position[1], position[2],  # start/end
                              position[3],  # strand
                              "%s:%s" % (position[4], position[5])])  # hairpin
            lines.append(line)
    lines = "\n".join(lines)
    logger.debug("BAM::CHROM: miRNAs %s in this chrom %s" % (lines, chr))
    return pybedtools.BedTool(lines, from_string=True)  # or bed_file_directly
