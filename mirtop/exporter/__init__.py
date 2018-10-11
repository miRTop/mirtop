from mirtop.exporter import isomirs, fasta


def export(args):
    """
    GFF3 to others formats
    """
    if args.format == "seqbuster":
        isomirs.convert(args)
    elif args.format == "fasta":
        fasta.convert(args)
