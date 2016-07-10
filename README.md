mirtop
---------

Command line tool to annotate with a standard naming miRNAs e isomiRs.

See more at: [isomiRs naming discussion](https://github.com/miRTop/incubator/blob/master/isomirs/isomir_naming.md)

Cite
---------

http://mirtop.github.io

Installation
------------
Thes best solution is to install conda to get an independent enviroment.

`wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh`

`bash Miniconda-latest-Linux-x86_64.sh -b -p ~/mirtop_env`

`export PATH=$PATH:~/mirtop_env`

`git clone http://github.com/miRTop/mirtop`

`cd mirtop`

`python setup.py develop`

Quick start
-----------

`cd data`
`mirtop annotate --sps hsa --hairpin examples/annotate/hairpin.fa --mirna examples/annotate/miRNA.str -o test_out_mirs_fasta examples/annotate/sim_isomir.fa`

Contributors
------------

* [Lorena Pantano](https://github.com/lpantano) (Bioinformatic Core, Harvard Chan School, Boston, USA)

