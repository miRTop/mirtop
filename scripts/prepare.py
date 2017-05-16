import os
import sys
from collections import defaultdict

from argparse import ArgumentParser

def _read_pri(fn):
    pri = dict()
    with open(fn) as inh:
        for line in inh:
            if line.startswith(">") & line.strip().endswith("pri"):
                name = line.strip()[1:-4]
            else:
                pri[name] = line.strip()
    return pri

def _read_bed(fn):
    bed = defaultdict(dict)
    with open(fn) as inh:
        for line in inh:
            cols = line.strip().split("\t")
            if cols[3].find("pri") > 0:
                continue
            if cols[3].find("loop") > 0:
                continue
            if cols[3].find("seed") > 0:
                continue
            if cols[3].find("motif") > 0:
                continue
            if cols[3].find("co") > 0:
                continue
            bed[cols[3].split("_")[0]].update({cols[3]: [int(cols[1]), int(cols[2]), cols[5]]})
    return bed


def _download(url, outfn):
    if os.path.isfile(outfn):
        return outfn
    os.system('wget -O %s %s' % (outfn, url))
    return outfn

if __name__ == "__main__":
    parser = ArgumentParser(description="Prepare files from mirGeneDB to be used with seqbuster")
    parser.add_argument("--bed", help="bed file with position of all sequence", required=1)
    parser.add_argument("--precursor30", help="file or url with fasta of precursor + 30 nt", required=1)
    args = parser.parse_args()

    sps = os.path.basename(args.precursor30).split("-")[0]

    if os.path.isfile(args.bed):
        fnbed = args.bed
    else:
        fnbed = _download(args.bed, "%s.bed" % sps)
    if os.path.isfile(args.precursor30):
        fnfa = args.precursor30
    else:
        fnfa = _download(args.precursor30, "%s.fa" % sps)
    fa = _read_pri(fnfa)
    bed = _read_bed(fnbed)
    OUT = open("%s.miRNA.str" % sps, 'w')
    OUTP = open("%s.hairpin.fa" % sps, 'w')
    for mir in fa:
        if mir in bed:
            precursor = bed[mir][mir + "_pre"]
            print precursor
            mir5p = ""
            mir3p = ""
            for mature in bed[mir]:
                info = bed[mir][mature]
                # print info
                if mature.endswith("pre"):
                    continue
                if precursor[2] == "-":
                    start = int(precursor[1]) - int(info[1]) + 31
                    end = int(precursor[1]) - int(info[0]) + 30
                else:
                    start = int(info[0]) - int(precursor[0]) + 31
                    end = int(info[1]) - int(precursor[0]) + 30
                    # print [mature, start, end, fa[mir][start:end]]
                if mature.find("5p") > 0:
                    mir5p = "[%s:%s-%s]" % (mature, start, end)
                if mature.find("3p") > 0:
                    mir3p = "[%s:%s-%s]" % (mature, start, end)

            print >>OUT, ">%s (X) %s %s" % (mir, mir5p, mir3p)
            print >>OUTP, ">%s\n%s" % (mir, fa[mir])
OUT.close()
OUTP.close()
