"""mirGFF3 proxy converter"""
from __future__ import print_function

import os.path as op

from mirtop.mirna import fasta, mapper
from mirtop.bam.bam import low_memory_bam, low_memory_genomic_bam
from mirtop.importer import seqbuster
from mirtop.gff import header
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)


def reader(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    samples = []
    database = mapper.guess_database(args)
    args.database = database
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    args.precursors = precursors
    matures = mapper.read_gtf_to_precursor(args.gtf)
    args.matures = matures
    # TODO check numbers of miRNA and precursors read
    # TODO print message if numbers mismatch
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
        h = header.create([sample], args.database, "")
        out_handle = open(fn_out, 'w')
        print(h, file=out_handle)
        if args.format == "BAM":
            if args.genomic:
                low_memory_genomic_bam(fn, sample, out_handle, args)
            else:
                low_memory_bam(fn, sample, out_handle, args)
        elif args.format == "seqbuster":
            seqbuster.read_file_low_memory(fn, sample, args, out_handle)
        out_handle.close()
