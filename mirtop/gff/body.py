import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def create(reads, database, sample, fn):
    """Read https://github.com/miRTop/mirtop/issues/9"""
    seen = set()
    lines = []
    seen_ann = {}
    # print >>out_handle, "seq\tname\tfreq\tchrom\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\tName"
    out_handle = open(fn, 'w')
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
                # count = _get_freq(r)
                seq = reads[r].sequence
                if iso.get_score(len(seq)) < 1:
                    continue
                if iso.subs:
                    iso.subs = [] if "N" in iso.subs[0] else iso.subs
                annotation = "%s.%s" % (chrom, iso.format_id(sep="."))
                idseq = reads[r].idseq
                source = "ref_miRNA" if iso.is_iso() else "isomiR"
                strand = iso.strand
                start, end = iso.start, iso.end
                score = iso.map_score
                filter = iso.filter
                mirName = iso.mirna
                preName = p
                Variant = iso.formatGFF()
                Cigar = iso.cigar
                counts = read.counts
                Filter = iso.filter
                attrb = ("ID {idseq}; Name {mirName}; Parent {preName}; Variant {Variant}; Cigar {Cigar}; Expression {counts}; Filter {Filter};").format(**locals())
                res = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t{score}\t{strand}\t.\t{attrb}").format(**locals())
                if annotation in seen_ann and seq.find("N") < 0 and seen_ann[annotation].split("\t")[0].find("N") < 0:
                    logger.warning("Same isomir %s from different sequence: \n%s and \n%s" % (annotation, res, seen_ann[annotation]))
                seen_ann[annotation] = res
                lines.append([annotation, chrom, counts, sample])
                logger.debug("GFF::%s" % res)
                # lines_pre.append([annotation, chrom, p, count, sample, hits])
                print >>out_handle, res
    out_handle.close()
    return lines

def _merge(lines):
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
