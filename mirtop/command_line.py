"""Main functions to simulates reads"""

import sys

from mirtop.libs.logger import initialize_logger
from mirtop.libs.parse import parse_cl
from mirtop.libs.simulator import simulate
from mirtop.gff import reader
from mirtop.gff.stats import stats
from mirtop.gff.compare import compare
from mirtop.gff.convert import convert_gff_counts
from mirtop.exporter import export
from mirtop.gff import validator
from mirtop.libs import spikeins
from mirtop.gff import update
import mirtop.libs.logger as mylog

import time


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out, kwargs['args'].debug,
                      kwargs['args'].print_debug)
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
        export(kwargs["args"])
    elif "validator" in kwargs:
        logger.info("Run validator.")
        validator.check_multiple(kwargs["args"])
    elif "query" in kwargs:
        logger.info("Not yet ready: This will allow queries to GFF files.")
    elif "spikein" in kwargs:
        logger.info("Run spike-in tools")
        spikeins.convert(kwargs["args"])
    elif "update" in kwargs:
        logger.info("Run update tools")
        update.convert(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))
