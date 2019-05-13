"""GFF reader and creator helpers"""
from __future__ import print_function

from collections import defaultdict, OrderedDict
from mirtop.mirna.realign import get_mature_sequence, align_from_variants, \
    read_id, variant_to_5p, variant_to_3p, variant_to_add, \
    is_sequence, make_id
from mirtop.gff.header import read_samples
from mirtop.gff.classgff import feature

import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)


def read(fn, args):
    """Read GTF/GFF file and load into annotate, chrom counts, sample, line"""
    samples = read_samples(fn)
    lines = defaultdict(dict)
    sep = " " if args.out_format == "gtf" else "="
    corrupted_uid = 0
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            line = paste_columns(feature(line), sep=sep)
            gff = feature(line)
            cols = gff.columns
            attr = gff.attributes
            if attr['UID'] and not read_id(attr['UID']):
                corrupted_uid += 1
                continue
            if 'UID' not in attr:
                msg = "UID not found."
                if 'Read' not in attr:
                    if not is_sequence(attr['Read']):
                        msg = msg + " Sequence not valid in Read attribute."
                    else:
                        attr['UID'] = make_id(attr['Read'])
                if 'sequence' not in attr:
                    msg = msg + " Sequence not found in sequence attribute."
                    if not is_sequence(attr['sequence']):
                        msg = msg + " Sequence not valid in sequence attribute."
                    else:
                        attr['UID'] = make_id(attr['Read'])
            if 'UID' not in attr:
                logger.warning("Line is not a valid GFF3 line: %s" %
                               line.strip())
                logger.warning(msg)

            if cols['start'] not in lines[cols['chrom']]:
                lines[cols['chrom']][cols['start']] = []
            uid = "%s-%s-%s" % (attr['UID'],
                                attr['Variant'],
                                attr['Name'])
            if args.keep_name:
                uid = "%s-%s" % (uid, attr['Read'])
            lines[cols['chrom']][cols['start']].append(
                [uid,
                 cols['chrom'],
                 attr['Expression'].strip().split(","),
                 samples,
                 line.strip()])
    logger.info("Lines skipped due to corrupted UID: %s" % corrupted_uid)
    return lines


def write_body_on_handle(lines, out_handle):
    for m in lines:
        for s in sorted(lines[m].keys()):
            for hit in lines[m][s]:
                print(hit[4], file=out_handle)


def create(reads, database, sample, args, quiet=False):
    """Read https://github.com/miRTop/mirtop/issues/9"""
    sep = " " if args.out_format == "gtf" else "="
    seen = set()
    lines = defaultdict(defaultdict)
    seen_ann = {}
    filter_precursor = 0
    filter_score = 0
    n_hits = 0
    n_reads = 0
    n_seen = 0
    if args.add_extra:
        precursors = args.precursors
        matures = args.matures
    for (r, read) in reads.items():
        hits = set()
        [hits.add(mature.mirna) for mature in read.precursors.values()
            if mature.mirna]
        hits = len(hits)
        if len(read.precursors) > 0:
            n_reads += 1
        for (ps, iso) in read.precursors.items():
            p = list(ps)[0]
            if not iso.mirna:
                filter_precursor += 1
                continue
            if (r, iso.mirna) not in seen:
                seen.add((r, iso.mirna))
                chrom = p
                seq = reads[r].sequence
                seq_name = seq if not args.keep_name else r
                if iso.get_score(len(seq)) < 1:
                    filter_score += 1
                    continue
                if iso.subs:
                    iso.subs = [] if "N" in iso.subs[0] else iso.subs
                idseq = reads[r].idseq
                source = "ref_miRNA" if not iso.is_iso() else "isomiR"
                strand = iso.strand
                start, end = iso.start, iso.end
                score = iso.map_score
                mirName = iso.mirna
                preName = p
                Variant = iso.formatGFF()
                Cigar = iso.cigar
                counts = read.counts
                Filter = iso.filter
                annotation = "%s.%s.%s" % (chrom, idseq, seq_name)
                # TODO:  This need to be moved to use the feature class
                # It needs a dict with all variable in keys
                attrb = ("Read {seq_name};UID {idseq};Name {mirName};"
                         "Parent {preName};"
                         "Variant {Variant};Cigar {Cigar};"
                         "Expression {counts};"
                         "Filter {Filter};Hits {hits};").format(**locals())
                line = ("{chrom}\t{database}\t{source}\t{start}\t{end}"
                        "\t{score}\t{strand}\t.\t{attrb}").format(**locals())
                logger.debug("GFF::%s" % line)
                if args.add_extra:
                    extra = variant_with_nt(line, precursors, matures)
                    line = "%s Changes %s;" % (line, extra)

                line = feature(line).paste_columns(sep)
                if annotation in seen_ann and seq.find("N") < 0 and (
                        seen_ann[annotation].split("\t")[0].find("N") < 0):
                    logger.warning(
                        "Same isomir %s from different sequence:"
                        " \n%s and \n%s" % (annotation, line,
                                            seen_ann[annotation]))
                seen_ann[annotation] = line
                logger.debug("GFF::external %s" % iso.external)
                if start not in lines[chrom]:
                    lines[chrom][start] = []
                lines[chrom][start].append([annotation, chrom,
                                            counts, sample, line])
                logger.debug("GFF::%s" % line)
                n_hits += 1
            else:
                n_seen += 1
    if not quiet:
        logger.info("GFF miRNAs: %s" % len(lines))
        logger.info("GFF hits %s by %s reads" % (n_hits, n_reads))
        logger.info("Filtered by being duplicated: %s" % n_seen)
        logger.info("Filtered by being outside miRNA positions:"
                    " %s" % filter_precursor)
        logger.info("Filtered by being low score: %s" % filter_score)
    return lines


def create_line(read, name, database, args):
    sep = " " if args.out_format == "gtf" else "="

    if args.add_extra:
        precursors = args.precursors
        matures = args.matures

    for (ps, iso) in read.precursors.items():
        p = list(ps)[0]
        if not iso.mirna:
            continue
        chrom = p
        seq = read.sequence
        seq_name = seq if not args.keep_name else name
        if iso.get_score(len(seq)) < 1:
            continue
        if iso.subs:
            iso.subs = [] if "N" in iso.subs[0] else iso.subs
        idseq = read.idseq
        source = "ref_miRNA" if not iso.is_iso() else "isomiR"
        strand = iso.strand
        start, end = iso.start, iso.end
        score = iso.map_score
        mirName = iso.mirna
        preName = p
        Variant = iso.formatGFF()
        Cigar = iso.cigar
        counts = read.counts
        Filter = iso.filter
        annotation = "%s.%s.%s" % (chrom, idseq, seq_name)
        # This get correctly formated with paste_columns below
        attrb = ("Read {seq_name};UID {idseq};Name {mirName};"
                 "Parent {preName};"
                 "Variant {Variant};Cigar {Cigar};"
                 "Expression {counts};"
                 "Filter {Filter};").format(**locals())
        line = ("{chrom}\t{database}\t{source}\t{start}\t{end}"
                "\t{score}\t{strand}\t.\t{attrb}").format(**locals())
        logger.debug("GFF::%s" % line)
        if args.add_extra:
            extra = variant_with_nt(line, precursors, matures)
            line = "%s Changes %s;" % (line, extra)

        line = feature(line).paste_columns(sep)
        return line


def guess_format(line):
    return "=" if line.find("Name=") > -1 else " "


def paste_columns(line, sep=" "):
    """
    Create GFF/GTF line from read_gff_line
    """
    cols = line.columns
    attr = line.attributes
    attr_paste = "; ".join(
        "%s%s%s" % (a, sep, attr[a]) for a in attr)
    return "\t".join([cols['chrom'], cols['source'], cols['type'],
                      cols['start'], cols['end'], cols['score'],
                      cols['strand'], cols['ext'], attr_paste])


def read_variant(attrb, sep=" "):
    """
    Read string in variants attribute.

    Args:
        *attrb(str)*: string in Variant attribute.

    Returns:
        *(gff_dict)*: dictionary with:
            >>> {'iso_3p': -3, ...}
    """
    gff_dict = OrderedDict()
    logger.debug("variant: %s" % attrb)
    for gff_item in attrb.strip().split(","):
        item_pair = gff_item.strip().split(":")
        if len(item_pair) > 1:
            gff_dict[item_pair[0].strip()] = int(item_pair[1].strip())
        else:
            gff_dict[item_pair[0].strip()] = True
    logger.debug("Keys found: %s" % gff_dict.keys())
    logger.debug("Values found: %s" % gff_dict.values())
    return gff_dict


def read_attributes(gff_line, sep=" "):
    gff_line = gff_line.strip().split("\t")[8]
    gff_dict = OrderedDict()
    for gff_item in gff_line.strip().split(";"):
        item_pair = gff_item.strip().split(sep)
        if len(item_pair) > 1:
            gff_dict[item_pair[0].strip()] = item_pair[1].strip()
    return gff_dict


def read_gff_line(line):
    """
    Read GFF/GTF line and return dictionary with fields
    """
    if line.startswith("#"):
        return line
    cols = line.strip().split("\t")
    sep = guess_format(line)
    if len(cols) != 9:
        raise ValueError("Line has less than 9 elements: %s" % line)
    fields = {'chrom': cols[0],
              'source': cols[1],
              'type': cols[2],
              'start': cols[3],
              'end': cols[4],
              'score': cols[5],
              'strand': cols[6],
              'ext': cols[7],
              'attrb': read_attributes(line, sep)}
    return fields


def variant_with_nt(line, precursors, matures):
    """
    Return nucleotides changes for each variant type
    using Variant attribute, precursor sequences and
    mature position.
    """
    gff = feature(line)
    attr = gff.attributes
    read = read_id(attr["UID"])
    attr["Parent"] = attr["Parent"].split(",")[0]
    if attr["Parent"] not in matures:
        logger.warning("Parent miRNA not found in database %s" % attr["Parent"])
        return ""
    if attr["Name"] not in matures[attr["Parent"]]:
        logger.warning("miRNA not found in database %s" % attr["Name"])
        return ""

    logger.debug("GFF::BODY::precursors %s" % precursors[attr["Parent"]])
    logger.debug("GFF:BODY::mature %s" % matures[attr["Parent"]][attr["Name"]])

    t5 = variant_to_5p(precursors[attr["Parent"]],
                       matures[attr["Parent"]][attr["Name"]],
                       attr["Variant"])
    t3 = variant_to_3p(precursors[attr["Parent"]],
                       matures[attr["Parent"]][attr["Name"]],
                       attr["Variant"])
    add = variant_to_add(read,
                         attr["Variant"])
    mature_sequence = get_mature_sequence(
        precursors[attr["Parent"]],
        matures[attr["Parent"]][attr["Name"]],
        nt=8)
    logger.debug("GFF::BODY::mature_sequence %s" % mature_sequence)
    mm = align_from_variants(read,
                             mature_sequence,
                             attr["Variant"])
    if mm == "Invalid":
        return mm
    if len(mm) > 0:
        mm = "".join(["".join([str(v) for v in m]) for m in mm])
    else:
        mm = "0"
    return "iso_5p:%s,iso_3p:%s,iso_add3p:%s,iso_snv:%s" % (t5, t3, add, mm)
