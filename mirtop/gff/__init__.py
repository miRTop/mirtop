from mirtop.bam.bam import read_bam, annotate
from mirtop.gff import body
def reader(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    hairpin, mirna = download_mirbase(args)
    precursors = read_precursor(args.hairpin, args.sps)
    gtf = read_gtf(args.gtf)
    # check numnbers of miRNA and precursors read
    # print message if numbers mismatch
    out_dts = dict()
    for bam_fn in args.files:
        sample = op.splitext(op.basename(bam_fn))[0]
        if bam_fn.endswith("bam") or bam_fn.endswith("sam"):
            logger.info("Reading %s" % bam_fn)
            reads = read_bam(bam_sort_by_n, precursors)
        else:
            raise ValueError("Format not recognized. Only working with BAM/SAM files.")

        ann = annotate(reads, matures, precursors)
        out_dts[bam_fn] = body.create(reads)
    # merge all reads for all samples into one dicts
    # from dict with all samples convert each in a gff line
