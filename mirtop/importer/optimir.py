""" Read isomiR GFF files from optimir tool"""

from collections import defaultdict

import mirtop.libs.logger as mylog
from mirtop.gff.header import read_samples
from mirtop.gff.body import paste_columns
from mirtop.gff.body import variant_with_nt
from mirtop.gff.classgff import feature

logger = mylog.getLogger(__name__)


def read_file(fn, args):
    """
    Read OptimiR file and convert to mirtop GFF format.

    Args:
        *fn(str)*: file name with isomiR-SEA output information.

        *database(str)*: database name.

        *args(namedtuple)*: arguments from command line.
            See *mirtop.libs.parse.add_subparser_gff()*.

    Returns:
        *reads (nested dicts)*:gff_list has the format as
            defined in *mirtop.gff.body.read()*.

    """
    database = args.database
    gtf = args.gtf
    sep = " " if args.out_format == "gtf" else "="
    sample = read_samples(fn)
    reads = defaultdict(dict)
    logger.debug("OPTIMIR::SAMPLE::%s" % sample)
    with open(fn) as handle:
        for line in handle:
            gff = feature(line)
            fixed_line = line
            if gff.columns:
                if "Variant" not in gff.attributes:
                    gff.attributes["Variant"] = "NA"

                logger.debug("OPTIMIR::Chrom update from %s to %s" % (gff.columns["chrom"], gff.attributes["Parent"]))
                gff.columns["chrom"] = gff.attributes["Parent"].split(",")[0]
                fixed_line = gff.paste_columns(sep=sep)
                if args.add_extra:
                    extra = variant_with_nt(fixed_line, args.precursors, args.matures)
                    fixed_line = "%s Changes %s;" % (fixed_line, extra)

                fixed_line = paste_columns(feature(fixed_line), sep=sep)
                counts = gff.attributes["Expression"].split(",")
                chrom = gff.columns["chrom"]
                start = gff.columns["start"]
                if start not in reads[chrom]:
                    reads[chrom][start] = []
                reads[chrom][start].append([gff.attributes["UID"],
                                            gff.columns["chrom"],
                                            counts,
                                            sample, fixed_line])
    return reads
