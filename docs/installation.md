# Installation

## Bioconda

coming soon

## PIP

`pip install mirtop`

## develop version

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