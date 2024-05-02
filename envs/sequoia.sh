#!/bin/bash
. ~/.bashrc
mamba create -n sequoia 

mamba install \
    conda-forge::r-base \
    conda-forge::r-ggplot2 \
    conda-forge::r-ape \
    bioconda::r-seqinr \
    r::r-stringr \
    conda-forge::r-data.table \
    conda-forge::r-tidyr \
    conda-forge::r-dplyr \
    r::r-vgam \
    r::r-mass \
    conda-forge::r-devtools \
    bioconda::bioconductor-rsamtools \
    bioconda::bioconductor-genomicranges
mamba install conda-forge::r-optparse 
mamba install conda-forge/label/cf201901::r-biocmanager