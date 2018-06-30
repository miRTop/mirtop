"""Helpers to define the header fo the GFF file"""

import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)


def create(samples, database, custom, filter=None):
    """Create header for GFF file.

    Args:
        *samples (list)*: character list with names for samples

        *database (str)*: name of the database.

        *custom (str)*: extra lines.

        *filter (list)*: character list with filter definition.

    Returns:
        *header (str)*: header string.
    """
    header = ""
    header += _get_gff_version()
    header += _get_database(database)
    header += _get_samples(samples)
    header += custom
    return header

def _get_gff_version():
    return "## GFF3 adapted for miRNA sequencing data. VERSION 0.0.1\n"

def _get_samples(samples):
    return "## COLDATA: %s" % ",".join(samples)

def _get_database(database):
    return "## source-ontology: %s\n" % database


def _filter(filters):
    if not filters:
        return "## FILTER: PASS\n"
    return "## FILTER: %s" % ";\n".join(filters)


def read_samples(fn):
    """Read samples from the header of a GFF file.

    Args:
        *fn(str)*: GFF file to read.

    Returns:
        *(list)*: character list with sample names.
    """
    with open(fn) as inh:
        for line in inh:
            if line.startswith("## COLDATA"):
                return line.strip().split(": ")[1].strip().split(",")
    raise ValueError("%s doesn't contain COLDATA header." % fn)
