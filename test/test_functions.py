"""This directory is setup with configurations to run the main functional test.

Inspired in bcbio-nextgen code
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import functools

from nose import SkipTest
from nose.plugins.attrib import attr
import yaml


class FunctionsTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    @attr(read=True)
    def test_read(self):
        from mirtop.mirna import mapper, fasta
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = mapper.read_gtf_to_precursor("data/examples/annotate/hsa.gff3")
        print map_mir
        if map_mir["hsa-let-7a-1"]["hsa-let-7a-5p"][0] != 5:
            raise ValueError("GFF is not loaded correctly.")
        fasta_precursor = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        # read data/aligments/let7-perfect.bam
        return True

    @attr(code=True)
    def test_code(self):
        """testing code correction function"""
        from mirtop.mirna.realign import make_id
        print make_id("AAACCCTTTGGG")
        print make_id("AAACCCTTTGGGA")
        print make_id("AAACCCTTTGGGAT")

    @attr(cigar=True)
    def test_cigar(self):
        """testing cigar correction function"""
        cigar = [[0, 14], [1, 1], [0, 5]]
        from mirtop.mirna.realign import cigar_correction, make_cigar
        fixed = cigar_correction(cigar, "AAAAGCTGGGTTGAGGAGGA", "AAAAGCTGGGTTGAGAGGA")
        if not fixed[0] == "AAAAGCTGGGTTGAGGAGGA":
            raise ValueError("Sequence 1 is not right.")
        if not fixed[1] == "AAAAGCTGGGTTGA-GAGGA":
            raise ValueError("Sequence 2 is not right.")
        if not make_cigar("AAA-AAATAAA", "AGACAAA-AAA") == "MGMD3MI3M":
            raise ValueError("Cigar not eq to MAMDMMMIMMM: %s" % make_cigar("AAA-AAATAAA", "AGACAAA-AAA"))

    @attr(locala=True)
    def test_locala(self):
        """testing pairwise alignment"""
        from mirtop.mirna.realign import align
        print "\nExamples of perfect match, deletion, mutation"
        print align("TGAGGTAGTAGGTTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0]
        print align("TGAGGTGTAGGTTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0]
        print align("TGAGGTAGTAGGCTGTATAGTT", "TGAGGTAGTAGGTTGTATAGTT")[0]

    @attr(alignment=True)
    def test_alignment(self):
        """testing alignments function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.mirna import fasta
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        matures = {}
        # matures = mirtop.mirna.read_mature("data/examples/annotate/mirnas.gff", "hsa")
        def annotate(fn, precursors, matures):
            from mirtop.bam import bam
            reads = bam.read_bam(fn, precursors)
            # ann = mirtop.bam.bam.annotate(reads, matures, precursors)
            return True
        print "\nlast1D\n"
        annotate("data/aligments/let7-last1D.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   AGAGGTAGTAGGTTGTA
        print "\n1D\n"
        annotate("data/aligments/let7-1D.sam", precursors, matures)
        #mirna TGAGGTAG-TAGGTTGTATAGTT
        #seq   TGAGGTAGGTAGGTTGTATAGTTA
        print "\nlast7M1I\n"
        annotate("data/aligments/let7-last7M1I.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTA-AGT
        print "\nmiddle1D\n"
        annotate("data/aligments/let7-middle1D.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTATAGTT
        print "\nperfect\n"
        annotate("data/aligments/let7-perfect.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTATAG (3tt 3TT)
        print "\ntriming\n"
        annotate ("data/aligments/let7-triming.sam", precursors, matures)

    @attr(gff=True)
    def test_gff(self):
        """testing GFF function"""
        from mirtop.libs import logger
        from mirtop.mirna import mapper, fasta
        from mirtop.gff import body, header
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        # depend on https://github.com/miRTop/mirtop/issues/6
        matures = mapper.read_gtf_to_precursor("data/examples/annotate/hsa.gff3")
        # matures = mirtop.mirna.read_mature("data/examples/annotate/mirnas.gff", "hsa")
        from mirtop.bam import bam
        bam_fn = "data/aligments/let7-perfect.sam"
        reads = bam.read_bam(bam_fn, precursors)
        ann = bam.annotate(reads, matures, precursors)
        fn = bam_fn + ".gff"
        h = header.create(bam_fn, ["example"], "miRBase21")
        print h
        gff = body.create(ann, "miRBase21", "example", fn, header)
        return True

    @attr(collapse=True)
    def test_collapse(self):
        """testing GFF function"""
        from mirtop.libs import logger
        from mirtop.mirna import mapper, fasta
        from mirtop.gff import body, header
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        # depend on https://github.com/miRTop/mirtop/issues/6
        matures = mapper.read_gtf_to_precursor("data/examples/annotate/hsa.gff3")
        # matures = mirtop.mirna.read_mature("data/examples/annotate/mirnas.gff", "hsa")
        from mirtop.bam import bam
        bam_fn = "data/aligments/collapsing-isomirs.sam"
        reads = bam.read_bam(bam_fn, precursors)
        ann = bam.annotate(reads, matures, precursors)
        fn = bam_fn + ".gff"
        h = header.create(bam_fn, ["example"], "miRBase21")
        gff = body.create(ann, "miRBase21", "example", fn, header)
        print gff
        return True


