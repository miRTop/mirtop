""" Read bam files"""
from __future__ import print_function

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
from mirtop.gff import body
from mirtop.mirna.annotate import annotate

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
    bam_fn = _sam_to_bam(bam_fn)
    bam_fn = _bam_sort(bam_fn)
    reads = defaultdict(hits)
    if args.genomic:
        logger.warning("This is under development and variants can be unexact.")
        bed_fn = os.path.join(args.out, os.path.basename(bam_fn) + ".bed")
        _bed(bam_fn, bed_fn)
        intersect = _intersect(bed_fn, args.gtf)
        reads = _read_lifted_bam(intersect, reads, args, clean)
    else:
        reads = _read_original_bam(bam_fn, reads, args, clean)
    return reads


def low_memory_bam(bam_fn, sample, out_handle, args):
    if args.genomic:
        raise ValueError("low-memory option is not compatible with genomic coordinates.")
    precursors = args.precursors
    bam_fn = _sam_to_bam(bam_fn)
    bam_fn = _bam_sort(bam_fn)
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    lines = []
    current = None
    for line in handle:
        if not current or current == line.query_name:
            lines.append(line)
            current = line.query_name
        else:
            reads = _read_lines(lines, precursors, handle, args)
            ann = annotate(reads, args.matures, args.precursors, quiet=True)
            gff_lines = body.create(ann, args.database, sample, args, quiet=True)
            body.write_body_on_handle(gff_lines, out_handle)
            current = line.query_name
            lines = []
            lines.append(line)
    reads = _read_lines(lines, precursors, handle, args)
    ann = annotate(reads, args.matures, args.precursors, quiet=True)
    gff_lines = body.create(ann, args.database, sample, args, quiet=True)
    body.write_body_on_handle(gff_lines, out_handle)


def low_memory_genomic_bam(bam_fn, sample, out_handle, args):
    logger.info("Reading BAM file in low memory mode.")
    logger.warning("This is under development and variants can be unexact.")
    precursors = args.precursors
    bam_fn = _sam_to_bam(bam_fn)
    bam_fn = _bam_sort(bam_fn)

    bed_fn = os.path.join(args.out, os.path.basename(bam_fn) + ".bed")
    _bed(bam_fn, bed_fn)
    intersect = _intersect(bed_fn, args.gtf)
    lines = []
    current = None
    for line in intersect:
        # print(line)
        if not current or current == line[3]:
            lines.append(line)
            current = line[3]
        else:
            reads = _read_lifted_lines(lines, precursors, args)
            ann = annotate(reads, args.matures, args.precursors, quiet=True)
            gff_lines = body.create(ann, args.database, sample, args, quiet=True)
            body.write_body_on_handle(gff_lines, out_handle)
            current = line[3]
            lines = []
            lines.append(line)
    reads = _read_lifted_lines(lines, precursors, args)
    ann = annotate(reads, args.matures, args.precursors, quiet=True)
    gff_lines = body.create(ann, args.database, sample, args, quiet=True)
    body.write_body_on_handle(gff_lines, out_handle)


def _analyze_line(line, reads, precursors, handle, args):
    if line.reference_id < 0:
        logger.debug("READ::Sequence not mapped: %s" % line.reference_id)
        return None
    if not line.cigarstring:
        logger.debug("READ::Sequence malformed: %s" % line)
        return None
    query_name = line.query_name
    if query_name not in reads and not line.query_sequence:
        return None
    sequence = line.query_sequence if not line.is_reverse else reverse_complement(line.query_sequence)
    logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, sequence))
    if line.query_sequence and line.query_sequence.find("N") > -1:
        return None
    if query_name not in reads:
        reads[query_name].set_sequence(sequence)
        reads[query_name].counts = _get_freq(query_name)
    if line.is_reverse and not args.genomic:
        logger.debug("READ::Sequence is reverse: %s" % line.query_name)
        return None
    chrom = handle.getrname(line.reference_id)
    start = line.reference_start
    # If genomic endcode, liftover to precursor position
    if not start:
        return None
    cigar = line.cigartuples
    # if line.cigarstring.find("I") > -1:
    #     indels_skip += 1
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
        return None
    iso.subs, iso.add, iso.cigar = filter.tune(
        reads[query_name].sequence, precursors[chrom],
        start, cigar)
    logger.debug("READ::After tune start %s end %s" % (iso.start, iso.end))
    logger.debug("READ::iso add %s iso subs %s" % (iso.add, iso.subs))
    reads[query_name].set_precursor(chrom, iso)
    return reads


def _read_lines(lines, precursors, handle, args, clean=True):
    reads = defaultdict(hits)
    for line in lines:
        reads = _analyze_line(line, reads, precursors, handle, args)
    if clean:
        reads = filter.clean_hits(reads)
    return reads


def _read_original_bam(bam_fn, reads, args, clean):
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    indels_skip = 0
    precursors = args.precursors
    for line in handle:
        reads = _analyze_line(line, reads, precursors, handle, args)
    logger.info("Hits: %s" % len(reads))
    logger.info("Hits with indels %s" % indels_skip)
    if clean:
        reads = filter.clean_hits(reads)
        logger.info("Hits after clean: %s" % len(reads))
    return reads


def _analyze_lifted_line(line, reads, precursors, bed=False):
    start_idx = 9 if bed else 15
    end_idx = 10 if bed else 16
    attr_idx = 14 if bed else 20
    if str(line).find("miRNA_primary_transcript") < 0: # only working with mirbase
        return reads
    query_name = line[3]
    sequence = line[4]
    logger.debug(("READ::line name:{0}").format(line))
    if sequence and sequence.find("N") > -1:
        return reads
    if query_name not in reads:
        reads[query_name].set_sequence(sequence)
        reads[query_name].counts = _get_freq(query_name)
        reads[query_name].sequence = sequence

    chrom = line[attr_idx].strip().split("Name=")[-1]
    start = line[1]
    end = line[2]
    strand = line[5]
    if not start:
        return reads
    if strand == "+":
        start = int(start) - int(line[start_idx]) + 1
    else:
        start = int(line[end_idx]) - int(end)
    iso = isomir()
    iso.align = line
    iso.set_pos(start, len(reads[query_name].sequence))
    logger.debug("READ::From BAM start %s end %s at chrom %s" % (iso.start, iso.end, chrom))
    if len(precursors[chrom]) < start + len(reads[query_name].sequence):
        logger.debug("READ::%s start + %s sequence size are bigger than"
                     " size precursor %s" % (
                                             chrom,
                                             len(reads[query_name].sequence),
                                             len(precursors[chrom])))
        return reads
    iso.subs, iso.add, iso.cigar = filter.tune(
        reads[query_name].sequence, precursors[chrom],
        start, None)
    logger.debug("READ::After tune start %s end %s" % (iso.start, iso.end))
    logger.debug("READ::iso add %s iso subs %s" % (iso.add, iso.subs))

    reads[query_name].set_precursor(chrom, iso)
    return reads


def _read_lifted_lines(lines, precursors, args, clean=True):
    reads = defaultdict(hits)
    for line in lines:
        reads = _analyze_lifted_line(line, reads, precursors, bed=True)
    if clean:
        reads = filter.clean_hits(reads)
    return reads


def _read_lifted_bam(handle, reads, args, clean):
    indels_skip = 0
    precursors = args.precursors
    for line in handle:
        reads = _analyze_lifted_line(line, reads, precursors, bed=True)
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


def _bed(bam_fn, bed_fn):
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    current = None
    with open(bed_fn, 'w') as outh:
        for line in handle:
            if line.reference_id < 0:
                logger.debug("READ::Sequence not mapped: %s" % line.reference_id)
                continue
            if not line.cigarstring:
                logger.debug("READ::Sequence malformed: %s" % line)
                continue
            query_name = line.query_name
            if (not current or query_name!=current) and not line.query_sequence:
                continue
            if not current or query_name!=current:
                sequence = line.query_sequence if not line.is_reverse else reverse_complement(line.query_sequence)
            logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, sequence))
            if line.query_sequence and line.query_sequence.find("N") > -1:
                continue
            chrom = handle.getrname(line.reference_id)
            start = line.reference_start
            end = start + len(sequence) - 1
            current = query_name
            strand = "+" if not line.is_reverse else "-"
            bed_line = "\t".join(map(str, [chrom, start, end, query_name, sequence, strand]))
            outh.write(bed_line + '\n')


def _intersect(bam, gtf):
    bampy = pybedtools.BedTool(bam)
    gtfpy = pybedtools.BedTool(gtf)
    return bampy.intersect(gtfpy, wo=True, bed=True, s=True)


# deprecated from now on
def _read_quick_lines(lines):
    reads = defaultdict(hits)
    for line in lines:
        _analyze_quick_line(line, reads)
    return reads


def _analyze_quick_line(line, reads):
    if line.reference_id < 0:
        logger.debug("READ::Sequence not mapped: %s" % line.reference_id)
        return reads
    if not line.cigarstring:
        logger.debug("READ::Sequence malformed: %s" % line)
        return reads
    query_name = line.query_name
    if query_name not in reads and not line.query_sequence:
        return reads
    sequence = line.query_sequence if not line.is_reverse else reverse_complement(line.query_sequence)
    logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, sequence))
    if line.query_sequence and line.query_sequence.find("N") > -1:
        return reads
    if query_name not in reads:
        reads[query_name].set_sequence(sequence)
        reads[query_name].counts = _get_freq(query_name)
    return reads


def _read_quick_bam(bam_fn, reads):
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    for line in handle:
        reads = _analyze_quick_line(line, reads)
    return reads
