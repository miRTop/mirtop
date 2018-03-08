import os.path as op

from mirtop.mirna import fasta, mapper
from mirtop.bam.bam import read_bam
from mirtop.importer import seqbuster, srnabench, prost, isomirsea
from mirtop.mirna.annotate import annotate
from mirtop.gff import body, header, merge
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def reader(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    samples = []
    database = mapper.guess_database(args.gtf)
    # hairpin, mirna = download_mirbase(args)
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    # check numnbers of miRNA and precursors read
    # print message if numbers mismatch
    out_dts = dict()
    for fn in args.files:
        sample = op.splitext(op.basename(fn))[0]
        samples.append(sample)
        fn_out = op.join(args.out, sample + ".%s" % args.out_format)
        if args.format == "BAM":
            reads = _read_bam(fn, precursors)
        elif args.format == "seqbuster":
            reads = seqbuster.read_file(fn, precursors)
            custom = seqbuster.header()
        elif args.format == "srnabench":
            reads = srnabench.read_file(fn, precursors)
        elif args.format == "prost":
            out_dts[fn] = prost.read_file(fn, precursors, database, args.gtf)
        elif args.format == "isomirsea":
            out_dts[fn] = isomirsea.read_file(fn, database, args.gtf)
        if args.format not in  ["isomirsea", "prost"]:
            ann = annotate(reads, matures, precursors)
            out_dts[fn] = body.create(ann, database, sample)
        h = header.create([sample], database, "")
        _write(out_dts[fn], h, fn_out)
    # merge all reads for all samples into one dict
    merged = merge.merge(out_dts)
    fn_merged_out = op.join(args.out, "mirtop.%s" % args.out_format)
    _write(merged, header.create(samples, database, ""), fn_merged_out)

def _write(lines, header, fn):
    out_handle = open(fn, 'w')
    print >>out_handle, header
    for m in lines:
        for s in sorted(lines[m].keys()):
            for hit in lines[m][s]:
                print >>out_handle, hit[4]
    out_handle.close()

def _read_bam(bam_fn, precursors):
    if bam_fn.endswith("bam") or bam_fn.endswith("sam"):
        logger.info("Reading %s" % bam_fn)
        reads = read_bam(bam_fn, precursors)
    else:
        raise ValueError("Format not recognized. Only working with BAM/SAM files.")
    return reads
