# Building trees from low-input DNA

## Introduction

This [Nextflow](https://www.nextflow.io/) pipeline was written to generate
phylogenetic trees from low input WGS. It takes as input `CaVEMAN` (SNV) 
and `Pindel` (indel) VCF files generated from the alignment files from 
sequencing. The steps are as follows.

1. Filter the `CaVEMAN` SNV output.
    
    a. Run Mathis Sanders' [SangerLCMFiltering](https://github.com/MathijsSanders/SangerLCMFiltering). 
    
    b. Post-filter the VCF (FILTER = PASS & CLPM=0.00 ASRD>=0.87).

2. Filter the `Pindel` indel output.
    
    a. Filter the VCF (FILTER = PASS).

3. Run `cgpVAF` ((Sanger's internal module)[https://confluence.sanger.ac.uk/pages/viewpage.action?pageId=22710418])
    a.k.a. [vafCorrect](https://github.com/cancerit/vafCorrect). This re-checks
    variants in individual samples that were non-significant, but that were
    found to be significant in other samples of that donor. `cgpVAF` also
    facilitates unbiased `pileup` (SNVs) and `exonerate` (indels)-based VAF
    calculation for the union of variant sites in the set of related samples
    from each donor. 
    
4. Run [Sequoia](https://github.com/TimCoorens/Sequoia).
    

