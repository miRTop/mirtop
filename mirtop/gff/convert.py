import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)

def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = int(name.split("_x")[1])
    except:
        return 0
    return counts

def _tab_output(reads, out_file, sample):
    """Create info matrix"""
    seen = set()
    lines = []
    lines_pre = []
    seen_ann = {}
    dt = None
    with open(out_file, 'w') as out_handle:
        print >>out_handle, "seq\tname\tfreq\tchrom\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\tName"
        for r, read in reads.iteritems():
            hits = set()
            [hits.add(mature.mirna) for mature in read.precursors.values() if mature.mirna]
            hits = len(hits)
            for p, iso in read.precursors.iteritems():
                if len(iso.subs) > 3 or not iso.mirna:
                    continue
                if (r, iso.mirna) not in seen:
                    seen.add((r, iso.mirna))
                    chrom = iso.mirna
                    if not chrom:
                        chrom = p
                    count = _get_freq(r)
                    seq = reads[r].sequence
                    if iso.get_score(len(seq)) < 1:
                        continue
                    if iso.subs:
                        iso.subs = [] if "N" in iso.subs[0] else iso.subs
                    annotation = "%s.%s" % (chrom, iso.format_id(sep="."))
                    idname = ("{0}.{1}").format(chrom, iso.format_id(sep="."))
                    res = ("{seq}\t{r}\t{count}\t{chrom}\tNA\tNA\t{format}\tNA\tNA\tmiRNA\t{p}\t{hits}\t{idname}").format(format=iso.format().replace("NA", "0"), **locals())
                    if annotation in seen_ann and seq.find("N") < 0 and seen_ann[annotation].split("\t")[0].find("N") < 0:
                        logger.warning("Same isomir %s from different sequence: \n%s and \n%s" % (annotation, res, seen_ann[annotation]))
                    seen_ann[annotation] = res
                    lines.append([annotation, chrom, count, sample, hits])
                    lines_pre.append([annotation, chrom, p, count, sample, hits])
                    print >>out_handle, res

    if lines:
        dt = pd.DataFrame(lines)
        dt.columns = ["isomir", "chrom", "counts", "sample", "hits"]
        dt = dt[dt['hits']>0]
        dt = dt.loc[:, "isomir":"sample"]
        dt = dt.groupby(['isomir', 'chrom', 'sample'], as_index=False).sum()
        dt.to_csv(out_file + "_summary")
        dt_pre = pd.DataFrame(lines_pre)
        dt_pre.columns = ["isomir", "mature", "chrom", "counts", "sample", "hits"]
        dt_pre = dt_pre[dt_pre['hits']==1]
        dt_pre = dt_pre.loc[:, "isomir":"sample"]
        dt_pre = dt_pre.groupby(['isomir', 'chrom', 'mature', 'sample'], as_index=False).sum()
        return out_file, dt, dt_pre
    return None

def _merge(dts):
    """
    merge multiple samples in one matrix
    """
    df= pd.concat(dts)

    ma = df.pivot(index='isomir', columns='sample', values='counts')
    ma_mirna = ma
    ma = ma.fillna(0)
    ma_mirna['mirna'] = [m.split(":")[0] for m in ma.index.values]
    ma_mirna = ma_mirna.groupby(['mirna']).sum()
    ma_mirna = ma_mirna.fillna(0)
    return ma, ma_mirna

def _create_counts(out_dts, out_dir):
    """Summarize results into single files."""
    ma, ma_mirna = _merge(out_dts)
    out_ma = op.join(out_dir, "counts.tsv")
    out_ma_mirna = op.join(out_dir, "counts_mirna.tsv")
    ma.to_csv(out_ma, sep="\t")
    ma_mirna.to_csv(out_ma_mirna, sep="\t")
    return out_ma_mirna, out_ma

