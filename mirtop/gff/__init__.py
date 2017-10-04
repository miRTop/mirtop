import os.path as op

from mirtop.mirna import fasta, mapper
from mirtop.bam.bam import read_bam, annotate
from mirtop.gff import body, header
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def reader(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    database = "miRBase21" # read from GTF mirna file
    # hairpin, mirna = download_mirbase(args)
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    # check numnbers of miRNA and precursors read
    # print message if numbers mismatch
    out_dts = dict()
    for bam_fn in args.files:
        sample = op.splitext(op.basename(bam_fn))[0]
        fn_out = op.join(args.out, sample + ".gff")
        if bam_fn.endswith("bam") or bam_fn.endswith("sam"):
            logger.info("Reading %s" % bam_fn)
            reads = read_bam(bam_fn, precursors)
        else:
            raise ValueError("Format not recognized. Only working with BAM/SAM files.")

        ann = annotate(reads, matures, precursors)
        h = header.create(bam_fn, [sample], database)
        out_dts[bam_fn] = body.create(ann, database, sample, fn_out, h)
    # merge all reads for all samples into one dicts
    # from dict with all samples convert each in a gff line
