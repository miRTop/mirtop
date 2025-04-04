"""mirGFF3 proxy converter"""
from __future__ import print_function

import os.path as op

from mirtop.mirna import fasta, mapper
from mirtop.bam.bam import read_bam
from mirtop.importer import seqbuster, srnabench, prost, isomirsea, manatee, optimir
from mirtop.mirna.annotate import annotate
from mirtop.gff import body, header, merge, read
from mirtop.mirna.mapper import read_gtf_to_mirna
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)


def reader(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    if args.low_memory:
        read.reader(args)
        return None
    samples = []
    if args.database is None:
        database = mapper.guess_database(args)
    else:
        database = args.database
    args.database = database
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    args.precursors = precursors
    matures = mapper.read_gtf_to_precursor(args.gtf, database)
    args.matures = matures
    # TODO check numbers of miRNA and precursors read
    # TODO print message if numbers mismatch
    out_dts = dict()
    if args.keep_name and len(args.files) > 1:
        logger.warning("--keep-name when running multiple samples\n"
                       "can generate wrong results if the\n"
                       "name read is different across sample\n"
                       "for the same sequence.")
    for fn in args.files:
        fn = op.normpath(fn)
        if args.format != "gff":
            sample = op.splitext(op.basename(fn))[0]
            samples.append(sample)
            fn_out = op.join(args.out, sample + ".%s" % args.out_format)
        if args.format == "BAM":
            reads = _read_bam(fn, args)
        elif args.format == "seqbuster":
            reads = seqbuster.read_file(fn, args)
        elif args.format == "srnabench":
            out_dts[fn] = srnabench.read_file(fn, args)
        elif args.format == "prost":
            reads = prost.read_file(fn, precursors, database, args.gtf)
        elif args.format == "isomirsea":
            out_dts[fn] = isomirsea.read_file(fn, args)
        elif args.format == "manatee":
            out_dts[fn] = manatee.read_file(fn, database, args)
        elif args.format == "optimir":
            out_dts[fn] = optimir.read_file(fn, args)
        elif args.format == "gff":
            samples.extend(header.read_samples(fn))
            out_dts[fn] = body.read(fn, args)
            continue
        if args.format not in ["isomirsea", "srnabench", "manatee", 'optimir']:
            ann = annotate(reads, matures, precursors)
            out_dts[fn] = body.create(ann, database, sample, args)
        h = header.create([sample], database, header.make_tools(args.format))
        _write(out_dts[fn], h, fn_out, args)
    # merge all reads for all samples into one dict
    if args.low_memory:
        return None
    merged = merge.merge(out_dts, samples)
    fn_merged_out = op.join(args.out, "%s.%s" % (args.prefix, args.out_format))
    _write(merged, header.create(samples, database, header.make_tools([args.format])), fn_merged_out, args)


def _write(lines, header, fn, args = None):
    out_handle = open(fn, 'w')
    print(header, file=out_handle)
    database = mapper.guess_database(args)
    mapper = read_gtf_to_mirna(args.gtf, database)
    for m in lines:
        for s in sorted(lines[m].keys()):
            for hit in lines[m][s]:
                # TODO: convert to genomic if args.out_genomic
                if args and args.out_genomic:
                    lifted = body.lift_to_genome(hit[4], mapper)
                    print(lifted, file=out_handle)
                else:
                    print(hit[4], file=out_handle)
    out_handle.close()


def _read_bam(bam_fn, precursors):
    if bam_fn.endswith("bam") or bam_fn.endswith("sam"):
        logger.info("Reading %s" % bam_fn)
        reads = read_bam(bam_fn, precursors)
    else:
        raise ValueError("Format not recognized."
                         " Only working with BAM/SAM files.")
    return reads
