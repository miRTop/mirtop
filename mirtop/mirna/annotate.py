""" Read bam files"""
import copy

import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def _coord(sequence, start, mirna, precursor, iso):
    """
    Define t5 and t3 isomirs
    """
    insertion = 0
    deletion = 0
    if iso.subs:
        insertion = 1 if iso.subs[0][-1] == "-" else 0
    if iso.subs:
        deletion = 1 if iso.subs[0][1] == "-" else 0
    end = (iso.end - len(iso.add) - insertion + deletion)
    logger.debug("COOR:: s:%s len:%s end:%s fixedEnd:%s mirna:%s iso:%s" % (start, len(sequence), iso.end, end, mirna, iso.format()))
    dif = abs(mirna[0] - start)
    if start < mirna[0]:
        iso.t5 = sequence[:dif].upper()
    elif start > mirna[0]:
        iso.t5 = precursor[mirna[0]:mirna[0] + dif].lower()
    elif start == mirna[0]:
        iso.t5 = 0
    if dif > 4:
        logger.debug("COOR::start > 3 %s %s %s %s %s" % (start, len(sequence), dif, mirna, iso.format()))
        return None

    dif = abs(mirna[1] - end)
    if iso.add:
        iso.add = iso.add.replace("-", "")
        sequence = sequence[:-len(iso.add)]
    # if dif > 3:
    #    return None
    if end > mirna[1]:
        iso.t3 = sequence[-dif:].upper()
    elif end < mirna[1]:
        iso.t3 = precursor[mirna[1] + 1 - dif:(mirna[1] + 1)].lower()
    elif end == mirna[1]:
        iso.t3 = 0
    if dif > 4:
        logger.debug("COOR::end > 3 %s %s %s %s %s" % (len(sequence), end, dif, mirna, iso.format()))
        return None
    # logger.debug("coor end:%s %s %s %s %s iso:%s" % (start, len(sequence), end, dif, mirna, iso.format()))
    return True

def annotate(reads, mature_ref, precursors):
    """
    Using coordinates, mismatches and realign to annotate isomiRs

    reads: dict object that comes from read_bam fn
    mirbase_ref: dict object that comers from mirtop.mirna.read_mature
    precursors: dict object (key : fasta) that comes from mirtop.mirna.fasta.read_precursor

    Return: dict object with reasd as Keys and...
    """
    n_iso = 0
    for r in reads:
        for p in reads[r].precursors:
            start = reads[r].precursors[p].start
            end = reads[r].precursors[p].end
            for mature in mature_ref[p]:
                mi = mature_ref[p][mature]
                logger.debug(("\nANN::NEW::read:{s}\n pre:{p} start:{start} end: {end} "
                              "cigar: {cigar} "
                              "\n mir:{mature} mir_pos:{mi}\n mir_seqs:{mature_s}"
                              ).format(s=reads[r].sequence,
                                       mature_s = precursors[p][mi[0]:mi[1] + 1],
                                       cigar = reads[r].precursors[p].cigar,
                                       **locals()))
                iso_copy =  copy.deepcopy(reads[r].precursors[p])
                is_iso = _coord(reads[r].sequence, start, mi, precursors[p], iso_copy)
                logger.debug(("ANN::is_iso:{is_iso}").format(**locals()))
                logger.debug("ANN::annotation:%s iso:%s" % (r, reads[r].precursors[p].format()))
                logger.debug("ANN::annotation:%s Variant:%s" % (r, reads[r].precursors[p].formatGFF()))
                if is_iso:
                    n_iso += 1
                    reads[r].precursors[p] = iso_copy
                    reads[r].precursors[p].mirna = mature
                    # break
    logger.info("Valid hits (+/-3 reference miRNA): %s" % n_iso)
    return reads
