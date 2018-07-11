from __future__ import print_function

from optparse import OptionParser
import sys
import os
import re
import random
import numpy

from mirtop.mirna import fasta
from mirtop.mirna import mapper
from mirtop.mirna import realign
from mirtop.gff import body, header
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def write_collapse_fastq(reads, out_fn):
    idx = 0
    with open(out_fn, 'a') as outh:
        for r in reads:
            idx += 1
            print(">name%s_x%s" % (idx, r[1]), file=outh)
            print(r[0], file=outh)

def write_fastq(reads, out_fn):
    idx = 0
    with open(out_fn, 'a') as outh:
        for r in reads:
            idx += 1
            print("@name_read:%s" % idx, file=outh)
            print(r, file=outh)
            print("+", file=outh)
            print("I" * len(r), file=outh)

def create_read(read, count, adapter="TGGAATTCTCGGGTGCCAAGGAACTC", size=36):
    reads = list()
    for i in  range(0, count):
        rest = size - len(read)
        part = adapter[:rest]
        reads.append(read + part)
    return reads

def variation(info, seq):
    randS = random.randint(info[0] - 2, info[0] + 2) + 1
    randE = random.randint(info[1] - 1, info[1] + 2) + 1
    if randS < 1:
        randS = 1
    if randE > len(seq):
        randE = info[1] - 1
    randSeq = seq[randS:randE]
    t5Lab = ""
    t5Lab = seq[randS:info[0]] if randS < info[0] else t5Lab
    t5Lab = seq[info[0]:randS].lower() if randS > info[0] else t5Lab
    t3Lab = ""
    t3Lab = seq[randE:info[1] + 1].lower() if randE < info[1] + 1 else t3Lab
    t3Lab = seq[info[1] + 1:randE] if randE > info[1] + 1 else t3Lab
    # mutation
    isMut = random.randint(0, 10)
    mutLab = []
    if isMut == 3:
        ntMut = random.randint(0, 3)
        posMut = random.randint(0, len(randSeq) - 1)
        if not randSeq[posMut] == nt[ntMut]:
            temp = list(randSeq)
            mutLab = [[posMut, nt[ntMut], randSeq[posMut]]]
            temp[posMut] = nt[ntMut]
            randSeq = "".join(temp)

    # addition
    isAdd = random.randint(0, 3)
    addTag = ""
    if isAdd == 2:
        posAdd = random.randint(1, 3)
        for numadd in range(posAdd):
            ntAdd = random.randint(0, 1)
            print([randSeq, seq[randS + len(randSeq)]])
            if nt[ntAdd] == seq[randS + len(randSeq)]:
                ntAdd = 1 if ntAdd == 0 else 0
            randSeq += nt[ntAdd]
            addTag += nt[ntAdd]
            print([randSeq, randE, info[1]])
    return [randSeq, randS, t5Lab, t3Lab, mutLab, addTag]

def create_iso(name, mir, seq, numsim, exp):
    data = dict()
    reads = dict()
    full_read = list()
    clean_read = list()
    seen = set()
    for mirna in mir[name]:
        info = mir[name][mirna]
        mirSeq = seq[info[0]:info[1] + 1]
        for rand in range(int(numsim)):
             # expression
            e = 1
            if exp:
                trial = random.randint(1, 100)
                p = random.randint(1, 50) / 50.0
                e = numpy.random.negative_binomial(trial, p, 1)[0]
            iso = realign.isomir()
            randSeq, iso.start, iso.t5, iso.t3, iso.subs, iso.add = variation(info, seq)
            if randSeq in seen:
                continue
            seen.add(randSeq)
            iso.end = iso.start + len(randSeq)
            aln = realign.align(randSeq, seq[iso.start:iso.end])
            iso.cigar = realign.make_cigar(aln[0], aln[1])
            iso.mirna = mirna
            query_name = "%s.%s.%s" % (mirna, iso.format_id("."), randSeq)
            reads[query_name] = realign.hits()
            reads[query_name].set_sequence(randSeq)
            reads[query_name].counts = e
            reads[query_name].set_precursor(name, iso)
            full_read.extend(create_read(randSeq, e))
            clean_read.append([randSeq, e,])
            # print([randSeq, mutLab, addTag, t5Lab, t3Lab, mirSeq])
            # data[randSeq] = [exp, iso] # create real object used in code to generate GFF
    write_fastq(full_read, full_fq)
    write_collapse_fastq(clean_read, clean_fq)
    gff = body.create(reads, "miRBase21", "sim1")
    return gff

def _write(lines, header, fn):
    out_handle = open(fn, 'w')
    print(header, file=out_handle)
    for m in lines:
        for s in sorted(lines[m].keys()):
            for hit in lines[m][s]:
                print(hit[4], file=out_handle)
    out_handle.close()

usagetxt = "usage: %prog  --fa precurso.fa --gtf miRNA.gtf -n 10"

parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("--fa",
                  help="", metavar="FILE")
parser.add_option("--gtf",
                                help="", metavar="FILE")
parser.add_option("-n", "--num", dest="numsim",
                                help="")
parser.add_option("-e", "--exp", dest="exp", action="store_true",
                                help="give expression", default=False)
parser.add_option("-p", "--prefix", help="output name")
parser.add_option("--seed", help="set up seed for reproducibility.", default = None)


(options, args) = parser.parse_args()

if options.seed:
    random.seed(options.seed)

full_fq = "%s_full.fq" % options.prefix
clean_fq = "%s_clean.fq" % options.prefix
out_gff = "%s.gff" % options.prefix
if os.path.exists(full_fq):
    os.remove(full_fq)
if os.path.exists(clean_fq):
    os.remove(clean_fq)

pre = fasta.read_precursor(options.fa, "")
mir = mapper.read_gtf_to_precursor(options.gtf)

nt = ['A', 'T', 'G', 'C']
gffs = dict()
h = header.create(["sampleX"], "miRBase1", "")
for precursor in pre:
    seq = pre[precursor]
    gffs.update(create_iso(precursor, mir, seq, options.numsim, options.exp))


_write(gffs, h, out_gff)
