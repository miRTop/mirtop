from mirtop.exporter import seqbuster, isomirs, fasta, vcf


def export(args):
    """
    GFF3 to others formats
    """
    if args.format == "seqbuster":
        seqbuster.convert(args)
    if args.format == "isomir":
        isomirs.convert(args)
    elif args.format == "fasta":
        fasta.convert(args)
    elif args.format == "vcf":
        vcf.convert(args)
