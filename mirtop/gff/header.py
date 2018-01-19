import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def create(samples, database, custom, filter = None):
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
