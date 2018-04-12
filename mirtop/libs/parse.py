import argparse
import sys


def parse_cl(in_args):
    print in_args
    sub_cmds = {"gff": add_subparser_gff,
                "stats": add_subparser_stats,
                "compare": add_subparser_compare,
                "target": add_subparser_target,
                "join": add_subparser_join,
                "simulator": add_subparser_simulator,
                "counts": add_subparser_counts,
                "export": add_subparser_export
                }
    parser = argparse.ArgumentParser(description="small RNA analysis")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="mirtop supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        print "use %s" % sub_cmds.keys()
        sys.exit(0)
    args = parser.parse_args()

    assert sub_cmd is not None
    kwargs = {"args": args, sub_cmd: True}
    return kwargs


def _add_debug_option(parser):
    parser.add_argument("-d", "--debug", action="store_true",
                        dest="debug", help="max verbosity mode", default=False)
    parser.add_argument("-vd", "--print_debug", action="store_true",
                        help="print debug messageson terminal", default=False)
    return parser


def add_subparser_stats(subparsers):
    parser = subparsers.add_parser("stats", help="show general stats for each sample.")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_compare(subparsers):
    parser = subparsers.add_parser("compare", help="Compare two GFF files.")
    parser.add_argument("files", nargs="*", help="Files to compare."
                                                 "First will be used as reference.")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_gff(subparsers):
    parser = subparsers.add_parser("gff", help="realign miRNA BAM file")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps", required=1,
                        help="species")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf", help="gtf/gff file with precursor and mature position to genome.")
    parser.add_argument("--format", help="Input format, default BAM file.",
                        choices=['BAM', 'seqbuster', 'srnabench', 'prost', 'isomirsea'], default="BAM")
    parser.add_argument("--out-format", help="Supported formats: gff3 or gtf", default = "gtf")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_export(subparsers):
    parser = subparsers.add_parser("export", help="export GFF into other format")
    parser.add_argument("files", nargs="*", help="GFF files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps", required=1,
                        help="species")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf", help="gtf file with precursor and mature position to genome.")
    parser.add_argument("--format", help="Output format",
                        choices=['seqbuster'], default="seqbuster")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_join(subparsers):
    parser = subparsers.add_parser("join", help="join data")
    parser.add_argument("-f", "--fastq", dest="fastq", required=1,
                         help="fastq file"),
    parser.add_argument("-m", "--min", dest="minimum", default=1,
                        type=int,
                         help="Minimum number of counts required."
                              "Not recomended > 1. Could bias downstream"
                              "Analysis.")
    parser.add_argument("-o", "--out",
                         dest="out", help="output file", required=1)
    parser = _add_debug_option(parser)
    return parser

    parser = _add_debug_option(parser)
    return parser


def add_subparser_target(subparsers):
    parser = subparsers.add_parser("target", help="Annotate miRNA targets.")
    parser.add_argument("--input", required=1,
                        help="list of miRNAs in 1 column format")
    parser.add_argument("--sps", required=1,
                        help="species")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--annotation", required=1,
                        help="Folder with tarets annotation. If bcbio installed would be the srnaseq ffolder")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_simulator(subparsers):
    parser = subparsers.add_parser("simulator", help="simulate small RNAfrom fasta/bed file")
    parser.add_argument("--bed",
                        help="bed file with position of precursors <=200 nt")
    parser.add_argument("--fasta", help = "fasta with precursors.")
    parser.add_argument("--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("-r", "--reference", dest="ref",
                        help="reference fasta file with index"),
    parser = _add_debug_option(parser)
    return parser


def add_subparser_counts(subparsers):

    parser = subparsers.add_parser("counts", help="extract expression counts for each sample for mirna/variants")
    parser.add_argument("--gff",
                        help="/path/to/GFF/file/file.gff", required=True)
    parser.add_argument("--out", 
                        required=True,
                        help="/path/to/output/directory")
    parser = _add_debug_option(parser)
    return parser

