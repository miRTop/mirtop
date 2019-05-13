# mirtop <img src="artwork/logo.png" width="80" height="80">

[![](https://travis-ci.org/miRTop/mirtop.svg?branch=master)](https://travis-ci.org/miRTop/mirtop#)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

Command line tool to annotate with a standard naming miRNAs e isomiRs.

This tool adapt the miRNA GFF3 format agreed on here: https://github.com/miRTop/mirGFF3

Chat
----

[Ask question, ideas](https://gitter.im/mirtop/Lobby#)
[Contributors to code](https://gitter.im/mirtop/devel)

Cite
---------

http://mirtop.github.io

Contributing
------------

Everybody is welcome to contribute, fork the `devel` branch and start working!

If you are interesting in miRNA or small RNA analysis, you can jump into the incubator issue pages to propose/ask or say hi:

https://github.com/miRTop/incubator/issues

About
-----

Join the team: https://orgmanager.miguelpiedrafita.com/join/15463928

Read more: http://mirtop.github.io

Installation
------------

### Bioconda

`conda install mirtop -c bioconda`

### PIP

`pip install mirtop`

### develop version

Thes best solution is to install conda to get an independent enviroment.

```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

bash Miniconda-latest-Linux-x86_64.sh -b -p ~/mirtop_env

export PATH=$PATH:~/mirtop_env

conda install -c bioconda pysam pybedtools pandas biopython samtools

git clone http://github.com/miRTop/mirtop

cd mirtop

python setup.py develop
```

Quick start
-----------

Read complete commands at: https://mirtop.readthedocs.org

```
git clone mirtop
cd mirtop/data/examples/annotate
mirtop gff --sps hsa --hairpin hairpin.fa --gtf hsa.gff3 -o test_out sim_isomir.bam
```

Output
------

The `mirtop gff` generates the GFF3 adapter format to capture miRNA variations. The output is explained [here](https://github.com/miRTop/incubator/blob/master/format/definition.md).

Contributors
------------

* [Lorena Pantano](https://github.com/lpantano) (Bioinformatic Core, Harvard Chan School, Boston, USA)
* [Shruthi Bhat Bandyadka](https://github.com/sbb25) (Partners Personalized Medicine, Cambridge MA, USA)
* [Iñaki Martínez de Ilarduya](http://www.germanstrias.org/technology-services/high-performance-computing/contact/)(HPC core, IGTP, Badalona, Spain)
* Rafael Alis
* [Victor Barrera](https://github.com/vbarrera) (Bioinformatic Core, Harvard Chan School, Boston, USA)
* [Steffen Möller](https://github.com/smoe) (University of Rostock)
* [Kieran O'Neill](https://github.com/oneillkza)
* [Roderic Espin](Universitat Oberta de Barcelona)

Citizens
--------

Here we cite any person who has contribute somehow to the project different than through code development and/or
bioinformatic concepts.

Gianvito Urgese,
Jan Oppelt(CEITEC Masaryk University, Brno, Czech Republic),
Thomas Desvignes,
Bastian,
Kieran O'Neill (BC Cancer),
Charles Reid (University of California Davis),
Radhika Khetani (Harvard Chan School of Public Health),
Shannan Ho Sui (Harvard Chan School of Public Health),
Simonas Juzenas(CAU),
Rafael Alis (Catholic University of Valencia),
Aida Arcas (Instituto de Neurociencias (CSIC-UMH)),
Yufei Lin (Harvard University),
Victor Barrera(Harvard Chan School of Public Health),
Marc Halushka (Johns Hopkins University)
