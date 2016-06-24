from mirtop.libs.logger import initialize_logger
import sys
from mirtop.libs.parse import parse_cl
from mirtop import mirna
import mirtop.libs.logger as mylog
import time


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out, kwargs['args'].debug, kwargs['args'].print_debug)
    logger = mylog.getLogger(__name__)
    start = time.time()
    if "annotate" in kwargs:
        logger.info("Run annotation")
        mirna.miraligner(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))
