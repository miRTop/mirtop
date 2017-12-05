"""
Produce stats from GFF3 format
"""

import os
import pandas as pd

from mirtop.gff import header
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

# Add check first 

def stats(args):
    """
    From a list of files produce stats
    """
    out = dict()
    for fn in args.files:
        if not os.path.exists(fn):
            raise IOError("%s doesn't exist" %s)
        out[fn] = _calc_stats(fn)
    print out
    # merge samples in fn with pandas

def _get_samples(fn):
    with open(fn) as inh:
        for line in inh:
            if line.startswith("## COLDATA"):
                return line.strip().split(": ")[1].strip().split(",")
    raise ValueError("%s doesn't contain COLDATA header." % fn)
def _calc_stats(fn):
    """
    Read files and parse into categories
    """
    samples = _get_samples(fn)
    print samples
    lines = []
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            logger.debug("## STATS: attribute %s" % cols[8])
            print cols[8]
            attr_v = [v.strip().split(" ")[1] for v in cols[8].strip().split(";")[:-1]]
            attr_k = [v.strip().split(" ")[0] for v in cols[8].strip().split(";")[:-1]]
            attr = dict(zip(attr_k, attr_v))
            lines.extend(_classify(cols[2], attr, samples))
            # table [category, sample, counts] (add lines)
    df = _summary(lines)
            # then pandas.spread()
    return df

def _classify(srna_type, attr, samples):
    """
    Parse the line and return one
    line for each category and sample.
    """
    # iso_5p, iso_3p, iso_add ...
    # FILTER :: exact/isomiR_type
    lines = []
    counts = dict(zip(samples, attr['Expression'].split(",")))
    for s in counts:
        lines.append([srna_type, s, counts[s]])
        if attr['Variant'].find("iso") == -1:
            continue
        for v in attr['Variant'].split(","):
            lines.append([v.split(":")[0], s, counts[s]])
    return lines

def _summary(lines):
    """
    Summarize long table according to thresholds
    """
    # summarize  > 0, > 5, >10, >50, >100 and bind rows
    labels = ["category", "sample", "counts"]
    df = pd.DataFrame.from_records(lines, columns=labels)
    df.counts = df.counts.astype(int)
    df = df.groupby(['category', 'sample'], as_index=False).sum()
    return df
