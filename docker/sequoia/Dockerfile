# Use the base container phuongle2510/lcm_phylogeny:0.1.1
FROM phuongle2510/lcm_phylogeny:0.1.1

# Set environment variables (optional, in case needed)
ENV DEBIAN_FRONTEND=noninteractive

# Install required CRAN packages and Bioconductor packages
RUN Rscript -e " \
    install.packages(c('seqinr', 'stringr', 'tidyr', 'dplyr', 'data.table', 'VGAM', 'MASS'), repos='http://cran.r-project.org'); \
    if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('Rsamtools', 'GenomicRanges')) \
"

# Confirm installation by running Rscript with the installed packages
CMD ["Rscript", "-e", "library(seqinr); library(tidyverse); library(data.table); library(VGAM); library(MASS); library(Rsamtools); library(GenomicRanges); print('Required packages installed successfully!')"]
