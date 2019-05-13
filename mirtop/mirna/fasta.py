"""Read precursor fasta file"""

from collections import defaultdict

import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def read_precursor(precursor, sps=None):
    """
    Load precursor file for that species

    Args:
        *precursor(str)*: file name with fasta sequences

        *sps(str)*: if any, select species to keep.
            It'll do a `header_sequence.find(sps)`.

    Returns:
        *hairpin(dict)*: keys are precursor names and
            values are precursor sequences.
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
                logger.debug("PRECURSOR::name %s" % name)
            elif name:
                hairpin[name] += line.strip().replace("U", "T")
                logger.debug("PRECURSOR::sequence %s" % hairpin[name])
        if name:
            hairpin[name] = hairpin[name] + "NNNNNNNNNNNN"
    return hairpin
