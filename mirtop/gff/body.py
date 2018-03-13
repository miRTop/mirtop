from collections import defaultdict, OrderedDict

import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def create(reads, database, sample):
    """Read https://github.com/miRTop/mirtop/issues/9"""
    seen = set()
    lines = defaultdict(defaultdict)
    seen_ann = {}
    filter_precursor = 0
    filter_score = 0
    n_hits = 0
    n_reads = 0
    n_seen = 0
    for r, read in reads.iteritems():
        hits = set()
        [hits.add(mature.mirna) for mature in read.precursors.values() if mature.mirna]
        hits = len(hits)
        if len(read.precursors) > 0:
            n_reads += 1
        for p, iso in read.precursors.iteritems():
            if not iso.mirna:
                filter_precursor += 1
                continue
            if (r, iso.mirna) not in seen:
                seen.add((r, iso.mirna))
                chrom = p
                seq = reads[r].sequence
                if iso.get_score(len(seq)) < 1:
                    filter_score += 1
                    continue
                if iso.subs:
                    iso.subs = [] if "N" in iso.subs[0] else iso.subs
                idseq = reads[r].idseq
                annotation = "%s.%s" % (chrom, idseq)
                source = "ref_miRNA" if not iso.is_iso() else "isomiR"
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
                attrb = ("Read {r}; UID {idseq}; Name {mirName}; Parent {preName}; Variant {Variant}; Cigar {Cigar}; Expression {counts}; Filter {Filter}; Hits {hits};").format(**locals())
                res = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t{score}\t{strand}\t.\t{attrb}").format(**locals())
                if annotation in seen_ann and seq.find("N") < 0 and seen_ann[annotation].split("\t")[0].find("N") < 0:
                    logger.warning("Same isomir %s from different sequence: \n%s and \n%s" % (annotation, res, seen_ann[annotation]))
                seen_ann[annotation] = res
                logger.debug("GFF::external %s" % iso.external)
                if start not in lines[chrom]:
                    lines[chrom][start] = []
                lines[chrom][start].append([annotation, chrom, counts, sample, res])
                logger.debug("GFF::%s" % res)
                n_hits += 1
            else:
                n_seen += 1

    logger.info("GFF miRNAs: %s" % len(lines))
    logger.info("GFF hits %s by %s reads" % (n_hits, n_reads))
    logger.info("Filtered by being duplicated: %s" % n_seen)
    logger.info("Filtered by being outside miRNA positions: %s" % filter_precursor)
    logger.info("Filtered by being low score: %s" % filter_score)
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

def read_attributes(gff_line, sep = " "):
    gff_line = gff_line.strip().split("\t")[8]
    gff_dict = OrderedDict()
    for gff_item in gff_line.strip().split(";"):
        item_pair = gff_item.strip().split(sep)
        if len(item_pair) > 1:
            gff_dict[item_pair[0]] = item_pair[1]
    return gff_dict

