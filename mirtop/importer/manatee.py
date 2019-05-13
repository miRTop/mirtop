""" Read Manatee files"""
from __future__ import print_function

from collections import defaultdict
import os

import mirtop.libs.logger as mylog
from mirtop.bam.bam import intersect
from mirtop.bam import filter
from mirtop.mirna.mapper import get_primary_transcript, guess_database
from mirtop.mirna.realign import isomir, reverse_complement, make_id, hits
from mirtop.gff.body import paste_columns, variant_with_nt
# from mirtop.mirna import mapper
from mirtop.gff.classgff import feature
from mirtop.gff import body
from mirtop.mirna.annotate import annotate

logger = mylog.getLogger(__name__)


def read_file(fn, database, args):
    """
    Read Manatee file and convert to mirtop GFF format.

    Args:
        *fn(str)*: file name with Manatee output information.

        *database(str)*: database name.

        *args(namedtuple)*: arguments from command line.
            See *mirtop.libs.parse.add_subparser_gff()*.

    Returns:
        *reads (nested dicts)*:gff_list has the format as
            defined in *mirtop.gff.body.read()*.

    """
    reads = defaultdict(dict)
    sample = os.path.splitext(os.path.basename(fn))[0]
    precursors = args.precursors
    bed_fn = os.path.join(args.out, os.path.basename(fn) + ".bed")
    sep = " " if args.out_format == "gtf" else "="
    seen = set()
    with open(fn, 'r') as handle:
        if not os.path.exists(bed_fn):
            _bed(handle, bed_fn)
    intersect_fn = intersect(bed_fn, args.gtf)
    for line in intersect_fn:
        data = _analyze_line(line, precursors, database, sample, sep, args)
        if data:
            start = data["start"]
            chrom = data["chrom"]
            key = "%s:%s" % (data['mirna'], data["name"])
            if start not in reads[chrom]:
                reads[chrom][start] = []
            if key not in seen:
                seen.add(key)
                reads[chrom][start].append(data["line"])
    return reads


def _analyze_line(line, precursors, database, sample, sep, args):
    start_idx = 10
    end_idx = 11
    attr_idx = 15
    query_name = line[3]
    sequence = line[4]
    if str(line).find(get_primary_transcript(guess_database(args))) < 0: # only working with mirbase
        return None

    logger.debug(("READ::line name:{0}").format(line))
    if sequence and sequence.find("N") > -1:
        return None

    chrom = line[attr_idx].strip().split("Name=")[-1]
    start = line[1]
    end = line[2]
    strand = line[5]
    counts = float(line[6])
    Filter = "Pass"
    reads = dict()
    if not start:
        return None
    if strand == "+":
        start = int(start) - int(line[start_idx]) + 1
    else:
        start = int(line[end_idx]) - int(end)
    iso = isomir()
    iso.align = line
    iso.set_pos(start, len(sequence))
    logger.debug("READ::From BAM start %s end %s at chrom %s" % (iso.start, iso.end, chrom))
    if len(precursors[chrom]) < start + len(sequence):
        logger.debug("READ::%s start + %s sequence size are bigger than"
                     " size precursor %s" % (
                                             chrom,
                                             len(sequence),
                                             len(precursors[chrom])))
    iso.subs, iso.add, iso.cigar = filter.tune(
        sequence, precursors[chrom],
        start, None)
    logger.debug("READ::After tune start %s end %s" % (iso.start, iso.end))
    logger.debug("READ::iso add %s iso subs %s" % (iso.add, iso.subs))

    idu = make_id(sequence)
    reads[query_name] = hits()
    reads[query_name].set_sequence(sequence)
    reads[query_name].counts = counts
    reads[query_name].sequence = sequence
    reads[query_name].set_precursor(chrom, iso)
    reads = annotate(reads, args.matures, args.precursors, quiet=True)
    gff_line = body.create(reads, args.database, sample, args, quiet=True)
    if start not in gff_line[chrom]:
        return None
    line = gff_line[chrom][start][0][4]
    logger.debug("READ::line:%s" % line)
    if args.add_extra:
        extra = variant_with_nt(line, args.precursors,
                                args.matures)
        line = "%s Changes %s;" % (line, extra)

    line = paste_columns(feature(line), sep=sep)
    return {'chrom': chrom,
            'start': start,
            'name': query_name,
            'mirna': reads[query_name].precursors[chrom].mirna,
            'line': [idu, chrom, counts, sample, line]}


def _bed(handle, bed_fn):
    with open(bed_fn, 'w') as outh:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.strip().split()
            if cols[2]=="*":
                logger.debug("READ::Sequence not mapped: %s" % cols[0])
                continue
            query_name = cols[0]
            query_sequence = cols[9]
            counts = cols[14]
            start = int(cols[3])
            strand = cols[1]
            chrom = cols[2]
            # is there no hits
            # is the sequence always matching the read, assuming YES now
            # if not current or query_name!=current:
            query_sequence = query_sequence if not strand=="-" else reverse_complement(query_sequence)
            # logger.debug(("READ::Read name:{0} and Read sequence:{1}").format(line.query_name, sequence))
            if query_sequence and query_sequence.find("N") > -1:
                continue
            end = start + len(query_sequence) - 1
            bed_line = "\t".join(list(map(str, [chrom, start, end, query_name,
                                                query_sequence, strand, counts])))
            outh.write(bed_line + '\n')
