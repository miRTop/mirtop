import argparse
import sys


def parse_cl(in_args):
    print in_args
    sub_cmds = {"annotate": add_subparser_mirbuster,
                "target": add_subparser_target,
                "collapse": add_subparser_collapse,
                "simulator": add_subparser_simulator
                }
    parser = argparse.ArgumentParser(description="small RNA analysis")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="seqcluster supplemental commands")
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


def add_subparser_mirbuster(subparsers):
    parser = subparsers.add_parser("annotate", help="realign miRNA BAM file")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps", required=1,
                        help="species")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf", help="gtf file with precursor position to genome.")
    parser.add_argument("--mirna", help="miRNA.str")
    parser.add_argument("--miraligner", action="store_true",
                        help="align with JAVA version.", default=False)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_collapse(subparsers):
    parser = subparsers.add_parser("collapse", help="collapse data")
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
    parser = subparsers.add_parser("simulator", help="simulate small RNA  from bed file")
    parser.add_argument("--bed",
                        help="bed file with position of precursors <=200 nt")
    parser.add_argument("--fasta", help = "fasta with precursors.")
    parser.add_argument("--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("-r", "--reference", dest="ref",
                        help="reference fasta file with index"),
    parser = _add_debug_option(parser)
    return parser
