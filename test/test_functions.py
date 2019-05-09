"""
This directory is setup with configurations to run the main functional test.

Inspired in bcbio-nextgen code
"""
from __future__ import print_function
import os
import unittest
import argparse
import contextlib
import shutil

from nose.plugins.attrib import attr


@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    dirname = os.path.join("test", "test_automated_output")
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        yield dirname
    finally:
        os.chdir(orig_dir)


def annotate(fn, read_file, load=False, create=True, keep_name=False,
             gtf=None, genomic=None):
    args = argparse.Namespace()
    args.hairpin = "data/examples/annotate/hairpin.fa"
    args.sps = "hsa"
    args.gtf = "data/examples/annotate/hsa.gff3"
    args.out = "test/test_automated_output"

    if gtf:
        args.gtf = gtf
    args.genomic = genomic
    args.add_extra = True
    args.out_format = "gtf"
    args.keep_name = keep_name
    from mirtop.mirna import fasta, mapper
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    args.precursors = precursors
    args.matures = matures
    args.database = mapper.guess_database(args)
    from mirtop.mirna import annotate
    from mirtop.gff import body
    if not load:
        reads = read_file(fn, args)
    else:
        reads = read_file
    if create:
        ann = annotate.annotate(reads, matures, precursors)
        body = body.create(ann, "miRBase21", "Example", args)
    return body


class FunctionsTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    @attr(database=True)
    def test_database(self):
        from mirtop.mirna import mapper
        args = argparse.Namespace()
        args.gtf = "data/examples/annotate/hsa.gff3"
        db = mapper.guess_database(args)
        print("Database is %s" % db)
        if db != "miRBasev21":
            raise ValueError("%s not eq to miRBasev21" % db)

    @attr(read_hairpin=True)
    def test_read_hairpin(self):
        from mirtop.mirna import mapper, fasta
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = mapper.read_gtf_to_precursor(
            "data/examples/annotate/hsa.gff3")
        print(map_mir)
        if map_mir["hsa-let-7a-1"]["hsa-let-7a-5p"][0] != 5:
            raise ValueError("GFF is not loaded correctly.")
        fasta_precursor = fasta.read_precursor(
            "data/examples/annotate/hairpin.fa", "hsa")
        print(fasta_precursor)
        fasta_precursor2 = fasta.read_precursor(
            "data/examples/annotate/hairpin.fa", None)
        print(fasta_precursor2)
        if fasta_precursor != fasta_precursor2:
            raise ValueError("species value generates two different dicts.")
        # read data/aligments/let7-perfect.bam
        return True

    @attr(read_hairpin_mirgenedb=True)
    def test_read_hairpin_mirgenedb(self):
        from mirtop.mirna import mapper
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = mapper.read_gtf_to_precursor(
            "data/db/mirgenedb/hsa.gff")
        print(map_mir)

    @attr(read_mir2chr=True)
    def test_read_mir2chr(self):
        from mirtop.mirna import mapper
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = mapper.read_gtf_chr2mirna("data/examples/annotate/hsa.gff3")
        print(map_mir)
        # print(mapper.read_gtf_chr2mirna2("data/examples/annotate/hsa.gff3"))

    @attr(read_mir2genomic=True)
    def test_read_mir2genomic(self):
        from mirtop.mirna import mapper
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = mapper.read_gtf_to_mirna("data/examples/annotate/hsa.gff3")
        print(map_mir)

    @attr(read_line=True)
    def test_read_line(self):
        """Read GFF/GTF line"""
        from mirtop.gff.body import read_gff_line
        with open("data/examples/gff/2samples.gff") as inh:
            for line in inh:
                print(read_gff_line(line))

    @attr(code=True)
    def test_code(self):
        """testing code correction function"""
        from mirtop.mirna.realign import make_id, read_id

        def _convert(s, test, reverse=False):
            code = read_id(s) if reverse else make_id(s)
            if code != test:
                raise ValueError("%s didn't result on %s but in %s" %
                                 (s, test, code))

        _convert("AAACCCTTTGGG", "iso-12-B1NY4")
        _convert("AAACCCTTTGGGA", "iso-13-B1NYDX")
        _convert("AAACCCTTTGGGAT", "iso-14-B1NYI7")
        _convert("iso-12-B1NY4", "AAACCCTTTGGG", True)
        _convert("iso-13-B1NYDX", "AAACCCTTTGGGA", True)
        _convert("iso-14-B1NYI7", "AAACCCTTTGGGAT", True)

        # if make_id("AGTFCVS"):
        #     raise ValueError("This should be False. Not valid sequence.")
        # if read_id("asD(-"):
        #     raise ValueError("This should be False, Not valid code.")

    @attr(code_convert=True)
    def test_code_convert(self):
        """testing code correction function"""
        from mirtop.mirna.realign import make_id
        from mirtop.gff.update import read_uid_10

        if not make_id(read_uid_10("@#%$")) == "iso-12-B1NY4":
            raise ValueError("Update ID is not working.")
        if not make_id(read_uid_10("@#%$@2")) == "iso-13-B1NYDX":
            raise ValueError("Update ID is not working.")

    @attr(cigar=True)
    def test_cigar(self):
        """testing cigar correction function"""
        cigar = [[0, 14], [1, 1], [0, 5]]
        from mirtop.mirna.realign import cigar_correction, make_cigar, \
            cigar2snp, expand_cigar
        fixed = cigar_correction(cigar, "AAAAGCTGGGTTGAGGAGGA",
                                 "AAAAGCTGGGTTGAGAGGA")
        if not fixed[0] == "AAAAGCTGGGTTGAGGAGGA":
            raise ValueError("Sequence 1 is not right.")
        if not fixed[1] == "AAAAGCTGGGTTGA-GAGGA":
            raise ValueError("Sequence 2 is not right.")
        if not make_cigar("AAA-AAATAAA", "AGACAAA-AAA") == "MAMD3MI3M":
            raise ValueError("Cigar not eq to MAMD3MI3M: %s" %
                             make_cigar("AAA-AAATAAA", "AGACAAA-AAA"))
        # test expand cigar
        if not expand_cigar("3MA3M") == "MMMAMMM":
            raise ValueError("Cigar 3MA3M not eqaul to MMMAMMM but to %s" %
                             expand_cigar("3MA3M"))
        # test cigar to snp
        if not cigar2snp("3MA3M", "AAATCCC")[0] == [3, "A", "T"]:
            raise ValueError("3MA3M not equal AAATCCC but %s" %
                             cigar2snp("3MA3M", "AAATCCC"))

    @attr(sequence=True)
    def test_is_sequence(self):
        """testing if string is valid sequence"""
        from mirtop.mirna.realign import is_sequence
        if not is_sequence("ACTGC"):
            raise ValueError("ACTGC should return true.")
        if is_sequence("AC2TGC"):
            raise ValueError("AC2TGC should return false.")

    @attr(locala=True)
    def test_locala(self):
        """testing pairwise alignment"""
        from mirtop.mirna.realign import align
        print("\nExamples of perfect match, deletion, mutation")
        print(align("TGAGTAGTAGGTTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0])
        print(align("TGAGGTGTAGGTTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0])
        print(align("TGAGGTAGTAGGCTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0])
        print(align("TGANTAGTAGNTTGTATNGTT", "TGAGTAGTAGGTTGTATAGTTT")[0])
        print(align("TGANTAGTNGNTTGTATNGTT", "TGAGTATAGGCCTTGTATAGTT")[0])
        print(align("NCANAGTCCAAGNTCATN", "TCATAGTCCAAGGTCATG")[0])

    @attr(reverse=True)
    def test_reverse(self):
        """Test reverse complement function"""
        from mirtop.mirna.realign import reverse_complement
        print("Testing ATGC complement")
        if "GCAT" != reverse_complement("ATGC"):
            raise ValueError("ATGC complement is not: %s" %
                             reverse_complement("ATGC"))

    @attr(class_gff=True)
    def test_class(self):
        """Test class to read GFF line"""
        from mirtop.gff.classgff import feature
        gff = feature("hsa-let-7a-5p\tmiRBasev21\tisomiR\t4\t25\t0\t+\t.\t"
                      "Read hsa-let-7a-1_hsa-let-7a-5p_5:26_-1:-1_mut:"
                      "null_add:null_x861; UID bhJJ5WJL2;"
                      " Name hsa-let-7a-5p; Parent hsa-let-7a-1;"
                      " Variant iso_5p:+1,iso_3p:-1; Cigar 22M;"
                      " Expression 861; Filter Pass; Hits 1;")
        print(gff.columns)
        print(gff.attributes)

    @attr(merge=True)
    def test_merge(self):
        """Test merge functions"""
        from mirtop.gff import merge
        example_line = "hsa-let-7a-5p\tmiRBasev21\tisomiR\t4\t25"
        if merge._chrom(example_line) != "hsa-let-7a-5p":
            raise ValueError("Chrom should be hsa-let-7a-5p.")
        if merge._start(example_line) != "4":
            raise ValueError("Start should be 4.")
        expression = merge._convert_to_string({'s': 1, 'x': 2}, ['s', 'x'])
        print(merge._fix("hsa-let-7a-5p\tmiRBasev21\tisomiR\t4\t25\t0\t+\t.\t"
                         "Read hsa-let-7a-1_hsa-let-7a-5p_5:26_-1:-1_mut:"
                         "null_add:null_x861; UID bhJJ5WJL2;"
                         " Name hsa-let-7a-5p; Parent hsa-let-7a-1;"
                         " Variant iso_5p:+1,iso_3p:-1; Cigar 22M;"
                         " Expression 861; Filter Pass; Hits 1;", expression))
        if expression != "1,2":
            raise ValueError("This is wrong: %s" % expression)

    @attr(align_mature=True)
    def test_variant(self):
        """testing get mature sequence"""
        from mirtop.mirna import fasta, mapper
        from mirtop.mirna.realign import get_mature_sequence, \
            align_from_variants
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa",
                                          "hsa")
        matures = mapper.read_gtf_to_precursor(
            "data/examples/annotate/hsa.gff3")
        res = get_mature_sequence("GAAAATTTTTTTTTTTAAAAG", [5, 15])
        if res != "AAAATTTTTTTTTTTAAAA":
            raise ValueError("Results for GAAAATTTTTTTTTTTAAAAG was %s" % res)
        mature = get_mature_sequence(precursors["hsa-let-7a-1"],
                                     matures["hsa-let-7a-1"]["hsa-let-7a-5p"], nt = 8)
        if mature != "NNTGGGATGAGGTAGTAGGTTGTATAGTTTTAGGGT":
            raise ValueError("Results for hsa-let-7a-5p is %s" % mature)

        res = align_from_variants("AGGTAGTAGTTGTATAGTT", mature,
                                  "iso_5p:+2,iso_snv_central")
        if not res or res[0][0] != 10:
            raise ValueError("Wrong alignment for test 0 %s" % res)

        res = align_from_variants("AGGTAGTAGGTTGTATAGTT", mature,
                                  "iso_5p:+2")
        if res:
            raise ValueError("Wrong alignment for test 1 %s" % res)

        res = align_from_variants("GATGAGGTAGTAGGTTGTATAGTT", mature,
                                  "iso_5p:-2")
        if res:
            raise ValueError("Wrong alignment for test 2 %s" % res)

        res = align_from_variants("AGGTAGTAGGTTGTATAGTTTT", mature,
                                  "iso_5p:+2,iso_add3p:2")
        if res:
            raise ValueError("Wrong alignment for test 3 %s" % res)

        res = align_from_variants("AGGTAGTAGGTTGTATAGTTTT", mature,
                                  "iso_5p:+2,iso_3p:2")
        if res:
            raise ValueError("Wrong alignment for test 4 %s" % res)

        res = align_from_variants("AGGTAGTAGGTTGTATAG", mature,
                                  "iso_5p:+2,iso_3p:-2")
        if res:
            raise ValueError("Wrong alignment for test 5 %s" % res)

        res = align_from_variants("AGGTAGTAGGTTGTATAGAA", mature,
                                  "iso_5p:+2,iso_3p:-2,iso_add3p:2")
        if res:
            raise ValueError("Wrong alignment for test 6 %s" % res)

        res = align_from_variants("AGGTAGTAGGATGTATAGTT", mature,
                                  "iso_5p:+2,iso_snv_central")
        if not res:
            if res[0][0] != 10:
                raise ValueError("Wrong alignment for test 7 %s" % res)
        res = align_from_variants("AGGTAGTAGGATGTATAGAA", mature,
                                  "iso_5p:+2,iso_3p:-2,iso_add3p:2")
        if res:
            raise ValueError("Wrong alignment for test 8 %s" % res)

    @attr(alignment=True)
    def test_alignment(self):
        """testing alignments function"""
        from mirtop.bam import bam
        print("\nlast1D\n")
        print(annotate("data/aligments/let7-last1D.sam", bam.read_bam))
        # mirna TGAGGTAGTAGGTTGTATAGTT
        # seq   AGAGGTAGTAGGTTGTA
        print("\n1D\n")
        print(annotate("data/aligments/let7-1D.sam", bam.read_bam))
        # mirna TGAGGTAG-TAGGTTGTATAGTT
        # seq   TGAGGTAGGTAGGTTGTATAGTTA
        print("\nlast7M1I\n")
        print(annotate("data/aligments/let7-last7M1I.sam", bam.read_bam))
        # mirna TGAGGTAGTAGGTTGTATAGTT
        # seq   TGAGGTAGTAGGTTGTA-AGT
        print("\nmiddle1D\n")
        print(annotate("data/aligments/let7-middle1D.sam", bam.read_bam))
        # mirna TGAGGTAGTAGGTTGTATAGTT
        # seq   TGAGGTAGTAGGTTGTATAGTT
        print("\nperfect\n")
        print(annotate("data/aligments/let7-perfect.sam", bam.read_bam))
        # mirna TGAGGTAGTAGGTTGTATAGTT
        # seq   TGAGGTAGTAGGTTGTATAG (3tt 3TT)
        print("\ntriming\n")
        print(annotate("data/aligments/let7-triming.sam", bam.read_bam))

    @attr(alignment_genomic=True)
    def test_alignment_genomic(self):
        """testing alignments function"""
        from mirtop.bam import bam
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        # print(annotate("data/examples/annotate/hsa-let-7a-5ploss1_neg.sam",
        #                bam.read_bam,
        #                gtf="data/db/hsa.gff3", genomic=True))
        print("\ngenomic\n")
        with make_workdir():
            for example in ["hsa-let-7a-nm", "hsa-let-7a-5ploss1",
                            "hsa-let-7a-3ploss1", "hsa-let-7a-5ploss1_neg"]:
                print(annotate("data/examples/annotate/%s.sam" % example,
                               bam.read_bam,
                               gtf="data/db/mirbase/hsa.gff3", genomic=True))

    @attr(keep_name=True)
    def test_keep_name(self):
        from mirtop.bam import bam
        line = annotate("data/aligments/let7-perfect.sam",
                        bam.read_bam,
                        keep_name=True)
        print(line)
        if line["hsa-let-7a-1"][5][0][4].find("seq_perfect_x2") < 0:
            raise ValueError("Keep name failed: %s" % line)

    @attr(seqbuster=True)
    def test_seqbuster(self):
        """testing reading seqbuster files function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.importer import seqbuster
        print("\nperfect\n")
        annotate("data/examples/seqbuster/reads20.mirna", seqbuster.read_file)
        print("\naddition\n")
        annotate("data/examples/seqbuster/readsAdd.mirna", seqbuster.read_file)
        print("\nno frequency\n")
        annotate("data/examples/seqbuster/seqbuster_nofreq.mirna", seqbuster.read_file)

    @attr(srnabench=True)
    def test_srnabench(self):
        """testing reading srnabench files function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.importer import srnabench
        annotate("data/examples/srnabench", srnabench.read_file, create=False)

    @attr(optimir=True)
    def test_optimir(self):
        """testing reading optimir files function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.importer import optimir
        annotate("data/examples/optimir/synthetic_100_full.gff3", optimir.read_file, create=False)

    @attr(prost=True)
    def test_prost(self):
        """testing reading prost files function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.mirna import fasta
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa",
                                          "hsa")
        fn = "data/examples/prost/prost.example.txt"
        from mirtop.importer import prost
        reads = prost.read_file(
            fn, precursors, "miRBasev21", "data/examples/annotate/hsa.gff3")
        annotate("data/example/prost/prost.example.txt", reads, True)

    @attr(gff=True)
    def test_gff(self):
        """testing GFF function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.bam import bam
        bam_fn = "data/aligments/let7-perfect.sam"
        annotate(bam_fn, bam.read_bam)
        return True

    @attr(collapse=True)
    def test_collapse(self):
        """testing GFF function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.bam import bam
        bam_fn = "data/aligments/collapsing-isomirs.sam"
        annotate(bam_fn, bam.read_bam)
        return True

    @attr(counts=True)
    def test_counts(self):
        """testing convert_gff_counts in convert.py function"""
        from mirtop.libs import logger
        from mirtop.gff.convert import convert_gff_counts
        import argparse

        logger.initialize_logger("test counts", True, True)
        logger = logger.getLogger(__name__)

        args = argparse.Namespace()
        args.hairpin = "data/examples/annotate/hairpin.fa"
        args.sps = "hsa"
        args.gtf = "data/examples/annotate/hsa.gff3"
        args.gff = 'data/examples/synthetic/let7a-5p.gff'
        args.out = 'data/examples/synthetic'
        args.add_extra = True
        convert_gff_counts(args)
        os.remove(os.path.join(args.out, "let7a-5p.tsv"))

        return True

    @attr(stats=True)
    def test_stats(self):
        """testing stats function"""
        from mirtop.gff import stats
        from mirtop import version
        df = stats._calc_stats("data/examples/gff/correct_file.gff")
        stats._dump_log(df, version, None)
        print(df)

    @attr(variant=True)
    def test_string_variant(self):
        """testing parsing string variants"""
        from mirtop.gff import body
        gff = body.read_variant("iso_5p:-1,iso_add3p:2,iso_snp_central_supp")
        truthk = ["iso_5p", "iso_add3p", "iso_snp_central_supp"]
        truthv = [-1, 2, True]
        if len(gff) != 3:
            raise ValueError("Error size of output. Expectd 3.")
        if (truthk > list(gff.keys())) - (list(gff.keys()) > truthk):
            raise ValueError("Not found expected keys.")
        if not isinstance(gff["iso_snp_central_supp"], bool):
            raise ValueError("iso_snp_central_supp should be boolean.")
        if (truthv > list(gff.values())) - (list(gff.values()) > truthv):
            raise ValueError("Not found expected Values.")

    @attr(validator=True)
    def test_validator(self):
        """test validator functions"""
        from mirtop.gff.validator import _check_file
        _check_file("data/examples/gff/2samples.gff")
        _check_file("data/examples/gff/coldata_missing.gff")
        _check_file("data/examples/gff/3wrong_type.gff")

    @attr(spikeins=True)
    def test_spikeins(self):
        """Test spikeins reading and annotation"""
        from mirtop.libs import spikeins
        from mirtop.mirna.realign import get_mature_sequence
        load = spikeins.read_spikeins("data/examples/spikeins/spikeins.fa")
        print(load)
        load1 = load['spikein-1']
        mature_from_data = get_mature_sequence(load1['precursor'],
                                               load1['position'],
                                               exact=True)
        if mature_from_data != load1['mature']:
            raise ValueError("Sequences doesn't match \n%s\n%s" %
                             (mature_from_data, load1['mature']))

        file_fasta = "data/examples/spikeins/spikeins_pre.fasta"
        file_gff = "data/examples/spikeins/spikeins_pre.gff"
        spikeins.write_precursors(load, file_fasta)
        spikeins.write_gff(load, file_gff)

        from mirtop.mirna import mapper, fasta
        map_mir = mapper.read_gtf_to_mirna(file_gff)
        print(map_mir)
        fasta_precursor = fasta.read_precursor(file_fasta, None)
        print(fasta_precursor)

    @attr(export_fasta=True)
    def test_export_fasta(self):
        from mirtop.exporter.fasta import _process
        print("\n")
        _process("data/examples/gff/2samples.gff", None)

    @attr(update=True)
    def test_update(self):
        from mirtop.gff.update import update_file
        print("\n")
        update_file("data/examples/versions/version1.0.gff", None)
