from Bio import pairwise2
from Bio.Seq import Seq
from collections import defaultdict

from mirtop.mirna.keys import *
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)

class hits:

    def __init__(self):
        self.sequence = ""
        self.idseq = ""
        self.precursors = defaultdict(isomir)
        self.score = []
        self.best_hits = [] # maybe sam object?
        self.counts = 0

    def set_sequence(self, seq):
        self.sequence = seq
        self.idseq = make_id(seq)

    def set_precursor(self, precursor, isomir):
        self.precursors[precursor] = isomir

    def remove_precursor(self, precursor):
        del self.precursors[precursor]

class isomir:

    def __init__(self):
        self.t5 = []
        self.t3 = []
        self.add = []
        self.subs = []
        self.external = ""
        self.align = None
        self.cigar = None
        self.filter = "Pass"
        self.map_score = 0
        self.end = None
        self.start = None
        self.mirna = None
        self.strand = "+"

    def set_pos(self, start, l, strand = "+"):
        self.strand = strand
        self.start = start
        self.end = start + l - 1
        if strand == "-":
            self.start = start + l - 1
            self.end = start

    def formatGFF(self):
        value = ""
        subs = self.subs
        if self.external != "notsure" and self.external != "":
            return self.external
        for sub in subs:
            if sub:
                if sub[0] > 1 and sub[0] < 8:
                    value += "iso_snp_seed,"
                elif sub[0] == 8:
                    value += "iso_snp_central_offset,"
                elif sub[0] > 8 and sub[0] < 13:
                    value += "iso_snp_central,"
                elif sub[0] > 12 and sub[0] < 18:
                    value += "iso_snp_central_supp,"
                else:
                    value += "iso_snp,"

        if self.add:
            value += "iso_add:+%s," % len(self.add)
        if self.t5:
            size = len(self.t5)
            direction = "+" if self.t5.isupper() else "-"
            value += "iso_5p:%s%s," % (direction, size)
        if self.t3:
            size = len(self.t3)
            direction = "+" if self.t3.isupper() else "-"
            value += "iso_3p:%s%s," % (direction, size)
        if not value:
            value = "NA;"
        return value[:-1]

    def format(self, sep="\t"):
        subs = "".join(["".join(map(str, mism)) for mism in self.subs])
        if not subs:
            subs = "0"
        add = "0" if not self.add else self.add
        return "%s%s%s%s%s%s%s" % (subs, sep, add, sep,
                                   self.t5, sep, self.t3)

    def format_id(self, sep="\t"):
        subs = ["".join(["".join([c[2], str(c[0]), c[1]]) for c in self.subs])]
        if not subs:
            subs = []
        add = [] if not self.add else ["e%s" % self.add]
        t5 = ["s%s" % self.t5] if self.t5 and self.t5 != "NA" else []
        t3 = ["%s" % self.t3] if self.t3 and self.t3 != "NA" else []
        full = t5 + subs + t3 + add
        return sep.join([f for f in full if f])

    def get_score(self, sc):
        for a in self.add:
            if a in ['A', 'T']:
                sc -= 0.25
            else:
                sc -= 0.75
        for e in self.subs:
            sc -= 1
        return sc

    def is_iso(self):
        if self.t5 or self.t3 or self.add or self.subs or self.external != "":
            return True
        return False

def read_id(idu):
    """
    Inspared in MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates
    """
    seq = ""
    for i in idu:
        if i == "1" or i == "2":
            return seq[:-int(i)]
        else:
            seq += CODE2NT[i]
    return seq

def make_id(seq):
    """
    Inspared in MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates
    """
    start = 0
    idName = ""
    for i in range(0, len(seq) + 1, 3):
        if i == 0:
            continue
        trint = seq[start:i]
        idName += NT2CODE[trint]
        start = i
    if len(seq) > i:
        dummy = "A" * (3 - (len(seq) - i))
        trint = seq[i:len(seq)]
        idName += NT2CODE["%s%s" % (trint, dummy)]
        idName += str(len(dummy))
    return idName

def align(x, y):
    """
    https://medium.com/towards-data-science/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f
    """
    aligned_x =  pairwise2.align.globalms(x, y, 1, -1, -1, -0.5)[0]
    aligned_x = list(aligned_x)
    n_x = aligned_x[0]
    if "N" in n_x:
        N_indices = [i for i, ltr in enumerate(n_x) if ltr == 'N']
        n_x = list(n_x)
        for N_index in N_indices:
            n_x[N_index] = y[N_index]
        n_x = ''.join(n_x)

    aligned_x[0] = n_x
    return tuple(aligned_x)

def _add_cigar_char(counter, cache):
    if counter == 1:
        return cache
    else:
        return str(counter) + cache

def make_cigar(seq, mature):
    """
    Function that will create CIGAR string from aligment
    between read and reference sequence.
    """
    cigar = ""
    for pos in range(0,len(seq)):
        if seq[pos] == mature[pos]:
            cigar += "M"
        elif seq[pos] != mature[pos] and seq[pos] != "-" and mature[pos] != "-":
            cigar += seq[pos]
        elif seq[pos] == "-":
            cigar += "D"
        elif mature[pos] == "-":
            cigar += "I"

    cache = ""
    counter = 1
    short = ""
    for c in cigar:
        if c != cache and cache != "" and cache == "M":
            short += _add_cigar_char(counter, cache)
            counter = 1
            cache = c
        if c != "M":
            short += c
        if c == cache and c == "M":
            counter += 1
        cache = c

    if cache == "M":
        short += _add_cigar_char(counter, cache)
    return short

def cigar_correction(cigarLine, query, target):
    """Read from cigar in BAM file to define mismatches"""
    query_pos = 0
    target_pos = 0
    query_fixed = []
    target_fixed = []
    for (cigarType, cigarLength) in cigarLine:
        if cigarType == 0: #match 
            query_fixed.append(query[query_pos:query_pos+cigarLength])
            target_fixed.append(target[target_pos:target_pos+cigarLength])
            query_pos = query_pos + cigarLength
            target_pos = target_pos + cigarLength
        elif cigarType == 1: #insertions
            query_fixed.append(query[query_pos:query_pos+cigarLength])
            target_fixed.append("".join(["-"] * cigarLength))
            query_pos = query_pos + cigarLength
        elif cigarType == 2: #deletion
            target_fixed.append(target[target_pos:target_pos+cigarLength])
            query_fixed.append("".join(["-"] * cigarLength))
            target_pos = target_pos + cigarLength
    return ["".join(query_fixed), "".join(target_fixed)]

def expand_cigar(cigar):
    """
    From short CIGAR version to long CIGAR version
    where each character is each nts in the sequence
    """
    cigar_long = ""
    n = 0
    for nt in cigar:
        if nt in ["D", "M", "I", "A", "T", "C", "G"]:
            if n > 0:
                cigar_long += nt * int(n)
            else:
                cigar_long += nt
            n = 0
        else:
            if n > 0:
                n = int("%s%s" % (n, nt))
            else:
                n = int(nt)
    return cigar_long

def cigar2snp(cigar, reference):
    """
    From a CIGAR string and reference sequence
    return position of mismatches (indels included) as
      [pos, seq_nt, ref_nt]
    """
    snp = []
    pos_seq = 0
    pos_ref = 0
    for nt in expand_cigar(cigar):
        if nt != "M":
            if nt == "I":
                snp.append([pos_seq, nt, "-"])
                pos_seq += 1
            elif nt == "D":
                snp.append([pos_seq, "-", reference[pos_ref]])
                pos_ref += 1
            else:
                snp.append([pos_seq, nt, reference[pos_ref]])
                pos_ref += 1
                pos_seq += 1
        else:
            pos_ref += 1
            pos_seq += 1
    return snp

def reverse_complement(seq):
    return Seq(seq).reverse_complement()

def get_mature_sequence(precursor, mature):
    """From precursor FASTA and mature positions
       Get mature sequence +- 4 flanking nts
    """
    return precursor[mature[0] - 4 :mature[1] + 5]

def align_from_variants(sequence, mature, variants):
    """Giving the sequence read,
       the mature from get_mature_sequence,
       and the variant GFF annotation:
            Get a list of substitutions
    """
    snps = []
    k = [v.split(":")[0] for v in variants.split(",") if v.find(":") > -1]
    v = [int(v.split(":")[1]) for v in variants.split(",") if v.find(":") > -1]
    var_dict = dict(zip(k, v))
    logger.debug("realign::align_from_variants::variants %s" % variants)
    snp = [v for v in variants.split(",") if v.find("snp") > -1]
    if "iso_5p" in k:
        fix_5p = 4 - var_dict["iso_5p"]
        mature = mature[fix_5p:]
    if "iso_add" in k:
        sequence = sequence[:-1 * var_dict["iso_add"]]
    if "iso_3p" in k and var_dict["iso_3p"] > 0:
        sequence = sequence[:-1 * var_dict["iso_3p"]]
    logger.debug("realign::align_from_variants::snp %s" % snp)
    for pos in range(0, len(sequence)):
        if sequence[pos] != mature[pos]:
            value = ""
            if pos > 1 and pos < 8:
                value = "iso_snp_seed"
            elif pos == 8:
                value = "iso_snp_central_offset"
            elif pos > 8 and pos < 13:
                value = "iso_snp_central"
            elif pos > 12 and pos < 18:
                value = "iso_snp_central_supp"
            else:
                value = "iso_snp"
            logger.debug("realign::align_from_variants::value %s" % value)
            if value in snp:
                snps.append([pos, sequence[pos], mature[pos]])
    logger.debug("realign::align_from_variants::snps %s" % snps)
    return snps
