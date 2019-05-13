"""This directory is setup with configurations to run the main functional test.

Inspired in bcbio-nextgen code
"""
from __future__ import print_function
import os
import subprocess
import unittest
import shutil
import contextlib
import functools

from nose import SkipTest
from nose.plugins.attrib import attr


@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    dirname = os.path.join(os.path.dirname(__file__), "test_automated_output")
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)


def expected_failure(test):
    """Small decorator to mark tests as expected failure.
    Useful for tests that are work-in-progress.
    """
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner


class AutomatedAnalysisTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        #       self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print(dirname)
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))

    @attr(simulate=True)
    def test_simulate(self):
        """Check simulated data"""
        mirna = "TGAGGTAGTAGGTTGTATAGTT"
        correct = 0
        n = 0
        with open("data/examples/simulation/res/reads.mirna") as inh:
            header = inh.readline()
            for line in inh:
                cols = line.strip().split()
                mut, add, t5, t3 = cols[6:10]
                seq = cols[0]
                if mut!="0":
                    pos = int(mut[:-2])
                    nt1 = mut[-2]
                    nt2 = mut[-1]
                    seql = list(seq)
                    seql[pos] = nt2
                    seq = "".join(seql)
                if t5!="0" and t5.islower():
                    seq = "%s%s" % (t5.upper(), seq)
                elif t5!="0" and t5.isupper():
                    seq = seq[len(t5):]
                if add!="0":
                    seq = seq[:-len(add)]
                if t3!="0" and t3.islower():
                    seq = "%s%s" % (seq, t3.upper())
                elif t3!="0" and t3.isupper():
                    seq = seq[:-len(t3)]
                if seq == mirna:
                    correct += 1
                else:
                    print("\nerror:\n%s\n%s" % (seq, mirna))
                n += 1
        print("rate %s/%s" % (correct, n))

    @attr(complete=True)
    @attr(annotate=True)
    @attr(bam=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_bam(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--sps", "hsa", "--add-extra",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/annotate/sim_isomir.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(low_memory=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_bam_chunk(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff", "--low-memory",
                      "--sps", "hsa", "--add-extra",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/annotate/sim_isomir.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(cmd_bam_genomic=True)
    @attr(complete=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_genomic_bam(self):
        """Run genomic bam analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--sps", "hsa", "--add-extra", "--genomic",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/db/mirbase/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/annotate/hsa-let-7a-nm.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(cmd_bam_genomic_low_memory=True)
    @attr(complete=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_genomic_bam_low_memory(self):
        """Run genomic bam analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff", "--genomic", "--low-memory",
                      "--sps", "hsa", "--add-extra",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/db/mirbase/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/annotate/hsa-let-7a-nm.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_seqbuster=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_seqbuster(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "seqbuster",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/seqbuster/reads.mirna"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_seqbuster_low_memory=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_seqbuster_low_memory(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff", "--low-memory",
                      "--format", "seqbuster",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/seqbuster/reads.mirna"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_isomirsea=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_isomirsea(self):
        """Run isomirsea analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "isomirsea",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/isomir-sea/tagMir-all.gff",
                      "-d", "-vd"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_srnabench=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_srnabench(self):
        """Run srnabench analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "srnabench",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/srnabench/",
                      "-d", "-vd"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_optimir=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_optimir(self):
        """Run optimir analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "optimir",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/optimir/synthetic_100_full.gff3",
                      "-d", "-vd"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_manatee=True)
    @attr(cmd=True)
    def test_srnaseq_annotation_manatee(self):
        """Run Manatee analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "manatee",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/examples/manatee/simulated.sam",
                      "-d", "-vd"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_stats=True)
    @attr(cmd=True)
    def test_srnaseq_stats(self):
        """Run stats analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "stats",
                      "-o", "test_out_mirs",
                      "../../data/examples/gff/correct_file.gff"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_stats=True)
    @attr(cmd=True)
    def test_srnaseq_stats(self):
        """Run stats analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "stats",
                      "../../data/examples/gff/correct_file.gff"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_merge=True)
    @attr(cmd=True)
    def test_merge_bam(self):
        """
        Run collapse two samples
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--sps", "hsa", "--add-extra",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs",
                      "../../data/merge/samples1.sam",
                      "../../data/merge/samples2.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_export_seqbuster=True)
    @attr(cmd=True)
    def test_export_seqbuster(self):
        """
        Run export command
        """
        with make_workdir():
            clcode = ["mirtop",
                      "export",
                      "-o", "test_out_mirs",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "../../data/examples/gff/correct_file.gff"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_export_vcf=True)
    @attr(cmd=True)
    def test_export_vcf(self):
        """
        Run export command
        """
        with make_workdir():
            clcode = ["mirtop",
                      "export",
                      "-o", "test_out_mirs",
                      "--format", "vcf",
                      "-d", "-vd",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "../../data/examples/gff/correct_file.gff"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_count=True)
    @attr(cmd=True)
    def test_count(self):
        """
        Run count command
        """
        with make_workdir():
            clcode = ["mirtop",
                      "count",
                      "-o", "test_out_mirs",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "../../data/examples/synthetic/let7a-5p.gtf"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_spikeins=True)
    @attr(cmd=True)
    def test_spikeins_cmd(self):
        """Run spikeins analysis
        """
        import platform
        with make_workdir():
            shutil.copy("../../data/examples/spikeins/spikeins.fa",
                        "spikeins.fa")
            clcode = ["mirtop",
                      "spikein",
                      "spikeins.fa",
                      "-o",
                      "test_out_spikeins"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

            if platform.system() == "Linux":
                clcode = ["razers3", "-dr", "0", "-i", "80", "-rr", "90",
                          "-f", "-o", "spikeins.sam",
                          "test_out_spikeins/spikeins_pre.fasta",
                          "../../data/examples/spikeins/test-spikeins.fa"]
                print(" ".join(clcode))
                subprocess.check_call(clcode)
            else:
                shutil.copy("../../data/examples/spikeins/spikeins.sam",
                            "spikeins.sam")
            clcode = ["mirtop",
                      "gff",
                      "--add-extra",
                      "--hairpin", "test_out_spikeins/spikeins_pre.fasta",
                      "--gtf", "test_out_spikeins/spikeins_pre.gff",
                      "-o", "test_out_mirs",
                      "spikeins.sam"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(cmd_update=True)
    @attr(cmd=True)
    def test_update_cmd(self):
        """Run update analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "update",
                      "-o", "test_out_mirs",
                      "../../data/examples/versions/version1.0.gff"]
            print("")
            print(" ".join(clcode))
            subprocess.check_call(clcode)
