from __future__ import print_function

import argparse
import sys


def parse_cl(in_args):
    """Function to parse the subcommands arguments.
    """
    print(in_args)
    sub_cmds = {"gff": _add_subparser_gff,
                "stats": _add_subparser_stats,
                "compare": _add_subparser_compare,
                "target": _add_subparser_target,
                "simulator": _add_subparser_simulator,
                "counts": _add_subparser_counts,
                "export": _add_subparser_export,
                "validator": _add_subparser_validator,
                "spikein": _add_subparser_spikein,
                "update": _add_subparser_update
                }
    parser = argparse.ArgumentParser(description="small RNA analysis")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="mirtop supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        print("use %s" % sub_cmds.keys())
        sys.exit(0)
    args = parser.parse_args()
    if "files" in args:
        if not args.files:
            print("use %s -h to see help." % in_args[0])
            sys.exit(1)

    assert sub_cmd is not None
    kwargs = {"args": args, sub_cmd: True}
    return kwargs


def _add_debug_option(parser):
    parser.add_argument("-d", "--debug", action="store_true",
                        dest="debug", help="max verbosity mode", default=False)
    parser.add_argument("-vd", "--print_debug", action="store_true",
                        help="print debug messages on terminal", default=False)
    return parser


def _add_subparser_stats(subparsers):
    parser = subparsers.add_parser("stats", help="show general stats for each sample.")
    parser.add_argument("files", nargs="*", help="GFF/GTF files.")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_compare(subparsers):
    parser = subparsers.add_parser("compare", help="Compare two GFF files.")
    parser.add_argument("files", nargs="*", help="Files to compare."
                                                 "First will be used as reference.")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_gff(subparsers):
    parser = subparsers.add_parser("gff", help="realign miRNA BAM file")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps",
                        help="species")
    parser.add_argument("--keep-name", action="store_true",
                        default=False,
                        help="Use sequence name in the Attribute column.")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf",
                        help="GFF file with precursor and mature position to genome.")
    parser.add_argument("--format", help="Input format, default BAM file.",
                        choices=['BAM', 'seqbuster', 'srnabench',
                                 'prost', 'isomirsea', 'optimir',
                                 'manatee', 'gff'], default="BAM")
    parser.add_argument("--out-format", help="Supported formats: gff3 or gtf",
                        choices=["gff", "gft"], default="gff")
    parser.add_argument("--add-extra", help="Add extra attributes to gff",
                        action="store_true")
    parser.add_argument("--database", help="Custom database name",
                        default=None)
    parser.add_argument("--genomic", action="store_true", default=False,
                        help="BAM file is mapped against genome.")
    parser.add_argument("--low-memory", action="store_true", default=False,
                        help="Read File by chunks. Only supported for BAM files.")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_export(subparsers):
    parser = subparsers.add_parser("export", help="export GFF into other format")
    parser.add_argument("files", nargs="*", help="GFF files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps",
                        help="species")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf", help="gtf file with precursor and mature position to genome.")
    parser.add_argument("--format", help="Output format",
                        choices=['seqbuster', 'fasta', 'vcf'], default="seqbuster")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_target(subparsers):
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


def _add_subparser_simulator(subparsers):
    parser = subparsers.add_parser("simulator", help="simulate small RNAfrom fasta/bed file")
    parser.add_argument("--bed",
                        help="bed file with position of precursors <=200 nt")
    parser.add_argument("--fasta", help="fasta with precursors.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("-r", "--reference", dest="ref",
                        help="reference fasta file with index"),
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_counts(subparsers):

    parser = subparsers.add_parser("counts", help="extract expression counts for each sample for mirna/variants")
    parser.add_argument("--gff",
                        help="/path/to/GFF/file/file.gff", required = 1)
    parser.add_argument("-o", "--out",
                        required=True,
                        help="/path/to/output/directory")
    parser.add_argument("--add-extra", help="Add extra attributes to gff", action="store_true")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--gtf", help="gtf/gff file with precursor and mature position to genome.")
    parser.add_argument("--sps",
                        help="species")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_validator(subparsers):
    parser = subparsers.add_parser("validator", help="validate if the file has the correct format")
    parser.add_argument("files", nargs="*", help="GFF files")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_spikein(subparsers):
    parser = subparsers.add_parser("spikein",
                                   help="Work with spike-ins.")
    parser.add_argument("file", help="FASTA file with spikeins.")
    parser.add_argument("-o", "--out", dest="out",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser


def _add_subparser_update(subparsers):
    parser = subparsers.add_parser("update",
                                   help="update GFF to current format version.")
    parser.add_argument("files", nargs="*", help="GFF/GTF files.")
    parser.add_argument("-o", "--out", dest="out", default="tmp_mirtop",
                        help="folder of output files")
    parser = _add_debug_option(parser)
    return parser
