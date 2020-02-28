# CCP: Classification by clustering with Pfam
This tool searches the Pfam domain in the input FAST file and annotates species names.

## Recommended environment

Ubuntu 18.04

Memory > 60GB

## Install

Download and extract source files.

```
wget http://marine-meta.healthscience.sci.waseda.ac.jp/omd/ccp_v0.5.tar.gz
tar vxf ccp_v0.5.tar.gz
```

Insatall transeq (EMBOSS), Moose and BioPerl module.

```
sudo apt update
sudo apt -y install emboss libmoose-perl bioperl
```

If you are not an Ubuntu user, set PATH to transeq and install perl modules as follows.

```
sudo cpan install Moose
sudo cpan install Bio::SeqFeature::Generic
```

## How to use

```
/path/to/ccp/ccp.sh input.fasta
```

This FASTA file contains one or more contigs of a single bacteria.

## Example

```
./ccp.sh ecoli.fasta
     .
     .
     .
Calculating correlation coefficient...
Escherichia_coli_str._K-12_substr._MG1655_taxid511145   0.99631245221891
Escherichia_coli_K-12_taxid83333        0.995419660449845
Escherichia_coli_K-12_taxid83333        0.995419660449845
Escherichia_coli_KLY_taxid1435461       0.995328155844128
Escherichia_coli_BW25113_taxid679895    0.994095064570474
Escherichia_coli_DH1_taxid536056        0.993808667795345
Escherichia_coli_str._K-12_substr._MC4100_taxid1403831  0.992066376413065
Escherichia_coli_BW2952_taxid595496     0.991978772282405
Escherichia_coli_ER2796_taxid1245474    0.991231631038488
Escherichia_coli_str._K-12_substr._MG1655star_taxid879462       0.990685357989476
```
