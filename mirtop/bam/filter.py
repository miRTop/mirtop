from collections import defaultdict
from mirtop.mirna.realign import hits, cigar_correction, make_cigar, align
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def tune(seq, precursor, start, cigar):
    """
    The actual fn that will realign the sequence to find the nt changes
    at 5', 3' sequence and nt variations.

    Args:
        *seq (str)*: sequence of the read.

        *precursor (str)*: sequence of the precursor.

        *start (int)*: start position of sequence on the precursor, +1.

        *cigar (str)*: similar to SAM CIGAR attribute.

    Returns:

        *list* with:

            subs (list): substitutions

            add (list): nt added to the end

            cigar (str): updated cigar
    """
    if cigar:
        seq, mature = cigar_correction(cigar, seq, precursor[start:])
    else:
        seq, mature, score, p, size = align(seq, precursor[start:start + len(seq)])
        cigar = make_cigar(seq, mature)
    if seq.startswith("-"):
        seq = seq[1:]
    if seq.endswith("-"):
        seq = seq[:-1]
    logger.debug("TUNE:: %s %s %s" % (cigar, seq, mature))

    error = set()
    for pos in range(0, len(seq)):
        if seq[pos] != mature[pos]:
            error.add(pos)

    subs, add = [], []

    prob = 0
    add_position = []
    for e in range(len(seq) - 1, len(seq) - 6, -1):
        if e in error:
            prob = 1
        if prob == 1:
            add.append(seq[e])
            add_position.append(e)
        if e not in error and prob == 0 and seq[e] in ["A", "T"]:
            add.append(seq[e])
            add_position.append(e)
            continue
        if e not in error:
            if add:
                add.pop()
                add_position.pop()
            if prob == 0:
                add = []
                add_position = []
            break

    for e in error:
        if e not in add_position:
            subs.append([e, seq[e], mature[e]])

    logger.debug("TUNE:: %s %s" % (subs, add))

    return subs, "".join(add), make_cigar(seq, mature)


def clean_hits(reads):
    """
    Select only best matches from a list of hits from the same read.

    Args:
        *reads*: dictionary as:

        >>> {'read_id': mirtop.realign.hits, ...}

    Returns:

        *reads*: same than input but with best hits only.
    """
    new_reads = defaultdict(hits)
    for r in reads:
        world = {}
        sc = 0
        for p in reads[r].precursors:
            world[p] = reads[r].precursors[p].get_score(len(reads[r].sequence))
            if sc < world[p]:
                sc = world[p]
        new_reads[r] = reads[r]
        for p in world:
            logger.debug("CLEAN::score %s %s %s" % (r, p, world[p]))
            if sc != world[p]:
                logger.debug("CLEAN::remove %s %s %s" % (r, p, world[p]))
                new_reads[r].remove_precursor(p)

    return new_reads
