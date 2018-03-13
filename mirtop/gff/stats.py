"""
Produce stats from GFF3 format
"""

import os
import pandas as pd

from mirtop.gff import header
from mirtop.gff.body import read_attributes
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

# Add check first 

def stats(args):
    """
    From a list of files produce stats
    """
    out = list()
    for fn in args.files:
        if not os.path.exists(fn):
            raise IOError("%s doesn't exist" %s)
        logger.info("Reading: %s" % fn)
        out.append(_calc_stats(fn))
    df_final = pd.concat(out)
    outfn = os.path.join(args.out, "mirtop_stats.txt")
    if args.out != "tmp_mirtop":
        df_final.to_csv(outfn)
        logger.info("Stats saved at %s" % outfn)
    else:
        print df_final

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
    lines = []
    seen = set()
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            logger.debug("## STATS: attribute %s" % cols[8])
            attr = read_attributes(line)
            if "-".join([attr['Variant'], attr['Name']]) in seen:
                continue
            seen.add("-".join([attr['Variant'], attr['Name']]))
            lines.extend(_classify(cols[2], attr, samples))
    df = _summary(lines)
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
    df_sum = df.groupby(['category', 'sample'], as_index=False).sum()
    df_sum['category'] = ["%s_sum" % r for r in df_sum['category']]
    df_count = df.groupby(['category', 'sample'], as_index=False).count()
    df_count['category'] = ["%s_count" % r for r in df_count['category']]
    df_mean = df.groupby(['category', 'sample'], as_index=False).mean()
    df_mean['category'] = ["%s_mean" % r for r in df_mean['category']]
    df = pd.concat([df_sum, df_count, df_mean])
    return df
