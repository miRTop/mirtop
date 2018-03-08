"""Main functions to simulates reads"""

import sys

from mirtop.libs.logger import initialize_logger
from mirtop.libs.parse import parse_cl
from mirtop.libs.simulator import simulate
from mirtop.gff import reader
from mirtop.gff.stats import stats
from mirtop.gff.compare import compare
from mirtop.gff.convert import convert_gff_counts
from mirtop.exporter import isomirs
import mirtop.libs.logger as mylog

import time


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out, kwargs['args'].debug, kwargs['args'].print_debug)
    logger = mylog.getLogger(__name__)
    start = time.time()
    if "gff" in kwargs:
        logger.info("Run annotation")
        reader(kwargs["args"])
    elif "stats" in kwargs:
        logger.info("Run stats.")
        stats(kwargs["args"])
    elif "compare" in kwargs:
        logger.info("Run compare.")
        compare(kwargs["args"])
    elif "simulator" in kwargs:
        logger.info("Run simulation")
        simulate(kwargs["args"])
    elif "counts" in kwargs:
        logger.info("Run convert of GFF to TSV containing expression")
        convert_gff_counts(kwargs["args"])
    elif "export" in kwargs:
        logger.info("Run export of GFF into other format.")
        isomirs.convert(kwargs["args"])
    elif "join" in kwargs["args"]:
        logger.info("Not yet ready: This will join multiple GFF files.")
    elif "check" in kwargs["args"]:
        logger.info("Not yet ready: This will check GFF files.")
    elif "query" in kwargs["args"]:
        logger.info("Not yet ready: This will allow queries to GFF files.")
    elif "convert" in kwargs["args"]:
        logger.info("Not yet ready: This will output tabular format from GFF files.")
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))
