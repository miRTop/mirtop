from __future__ import print_function

import argparse
import sys

from mirtop import __version__

def parse_cl(in_args):
    """Function to parse the subcommands arguments.
    """
    sub_cmds = {"gff": _add_subparser_gff,
                "stats": _add_subparser_stats,
                "compare": _add_subparser_compare,
                "target": _add_subparser_target,
                "simulator": _add_subparser_simulator,
                "counts": _add_subparser_counts,
                "export": _add_subparser_export,
                "validate": _add_subparser_validator,
                "spikein": _add_subparser_spikein,
                "update": _add_subparser_update,
		"sql": _add_subparser_sql
                }
    parser = argparse.ArgumentParser(description="small RNA analysis")
    parser.add_argument("--version", action="store_true",help="show version.")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        print(in_args)
        subparsers = parser.add_subparsers(help="mirtop supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    elif (len(in_args) > 0):
        args = parser.parse_args()
        if args.version:
            print("mirtop %s" % __version__)
            sys.exit(0)
    else:
        print("use %s" % sub_cmds.keys())
        print("mirtop %s" % __version__)
        sys.exit(0)

    args = parser.parse_args()
    if "files" in args:
        if not args.files:
            print("use %s -h to see help." % in_args[0])
            print("mirtop %s" % __version__)
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
    parser.add_argument("--prefix", dest="prefix", required=0,
                        default="mirtop", help="prefix for output file")
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
    parser.add_argument("--out-genomic", action="store_true", default=False,
                        help="output in genomic coordinates.")
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
                        choices=['seqbuster', 'fasta', 'vcf', 'isomir'], default="isomir")
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
    parser = subparsers.add_parser("validate", help="validate if the file has the correct format")
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


def _add_subparser_sql(subparsers):
    parser = subparsers.add_parser("sql", help="SQL create or query from GFF.", formatter_class=argparse.RawTextHelpFormatter, )
    #parser.add_argument("--gff", help="GFF file with precursor and mature position to genome.")
    parser.add_argument('--db', metavar='', action='store', help='SQL Database name. (default: mirtop.db)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c','--create', help="Creates a SQLite database from GFF", action='store_true') 
    group.add_argument('-q', '--query', help="Query from a SQLite database",  action='store_true')

    group1 = parser.add_argument_group('SQL create usage mode')
    group1.add_argument("--gff", metavar='', help="GFF file with precursor and mature position to genome")
    group1.add_argument("-o", "--out", metavar='', dest="out", default="tmp_mirtop", help="Directory of output files")

    group2 = parser.add_argument_group('SQL query usage mode') 
    group2.add_argument("-t", "--table", metavar='', help="Specify table name to use")
    #group2.add_argument("-txti", "--txt-in", metavar='', help="Provide the list of miRNA's in the text file as input. NOTE: List of miRNA's should be separated by new line")
    group2.add_argument("-txto", "--txtout", metavar='', help="Writes the output of the query to a file speficied. Format (-fmt) is a tab-delimited text file by default")
    #group2.add_argument("-a", "--all", metavar='', help="Selects all the columns from the table")
    group2.add_argument("-col", "--columns", metavar='', help="Select specific columns from the table to display (Default: all columns), or use with -n option to return n-counts. For information of the available columns see 'show-schema' or 'show-columns'. NOTE: options -e select must be applied!.")
    group2.add_argument("-n", "--count", metavar='', help="Returns 'n' counts for the query. Options 'T' for True, if not 'F' (Default: -n F). NOTE: options -e select must be applied! and accepts only one column from -col option." )
    group2.add_argument("-miR", "--miRNA", metavar='', help="Specify the miRNA names to query. For multiple miRNAs use comma(,) as separator; or text file (.txt) separated with new line character")
    group2.add_argument("-pm", "--miRNA_prefix", metavar='', help="Specify the prefix name for miRNAs to query. Example: -pm hsa -miR let-7a-5p results into querying hsa-let-7a-5p")
    group2.add_argument("-var", "--variant", metavar='', help="""Specify one or more types of variants to query. Use comma(,) as separator 
    Choices supports the following:
        iso_5p                  - indicates the shift at the reference 5' miRNA
        iso_3p                  - indicates the shift at the reference 3' miRNA
        iso_add3p               - Number of non-template nucleotides added at 3p
        iso_add5p               - Number of non-template nucleotides added at 5p
        iso_snv_seed            - when affected nucleotides are between [2-7]
        iso_snv_central_offset  - when affected nucleotides is at position [8]
        iso_snv_central         - when affected nucleotides are between [9-12]
        iso_snv_central_supp    - when affected nucleotides are between [13-17]
        iso_snv                 - anything else
        """)
    group2.add_argument("-f", "--filter", metavar='', help="Specify Filter tag attribute. Options: Pass, Reject. (Default: None)")
    group2.add_argument("-l", "--limit", metavar='', help="Specify the number of rows to output. (Example: --limit 30, to limit the first 30 rows)")
    #  group2.add_argument("-imiR", "--isomiR", metavar='', help="Specify the miRNA name to query")
    #  group2.add_argument("-s", "--schema", metavar='', help="Show the schema of the select tables; (-s <table_name>)")
    group2.add_argument("-e", "--expr", metavar='', 
        help="""Expression is the query that you want to run; (-e \"<statement>\")
    Choices supports the following: 
       show-tables              - Displays tables in the database (default: mirtop.db) 
       show-schema              - Displays the table schema (requires -t)
       show-columns             - Displays available columns in the table
       describe-gff             - Prints out the header information from the GFF file
       isomirs-per-mirna        - Displays the count of isomiRs for miRNA (requires -miR)
       select                   - Allows specific query construction. 
                                  Example: mirtop sql --db tmp_mirtop/SRR333680_revised2.db -qe select -var iso_5p,iso_3p -miR hsa-let-7a-5p,hsa-let-7d-5p -l 30
                                  The above expression evaluates to selecting miRNAs in -miR with variants in -var and prints out first 30 rows in --limit 
    """)

    parser = _add_debug_option(parser)
    return parser

