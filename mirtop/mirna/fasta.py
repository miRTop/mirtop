import os.path as op
import os
import re
from collections import defaultdict

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)

def read_precursor(precursor, sps):
    """
    Load precursor file for that species
    """
    hairpin = defaultdict(str)
    name = None
    with open(precursor) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                if name in hairpin:
                    hairpin[name] = hairpin[name] + "NNNNNNNNNNNN"
                if not sps or line.find(sps) > -1:
                    name = line.strip().replace(">", " ").split()[0]
                else:
                    name = None
                logger.debug(name)
            elif name:
                hairpin[name] += line.strip().replace("U", "T")
                logger.debug(hairpin[name])
        if name:
            hairpin[name] = hairpin[name] + "NNNNNNNNNNNN"
    return hairpin
