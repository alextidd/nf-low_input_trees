# Building trees from low-input DNA

This [Nextflow](https://www.nextflow.io/) pipeline was written to generate
phylogenetic trees from low input WGS/WES. It takes as input `CaVEMAN` (SNV) 
and `Pindel` (indel) VCF files generated from the alignment files from 
sequencing. The steps are as follows.

1. Reflag the VCFs to remove some filters that have been found to be too 
stringent. The flags that are removed depend on the type of experiment (WGS vs
WES vs TGS) and the variant type (SNV vs indel).

2. Filter the variants.

  a. Filter the `CaVEMAN` SNV output. We run Mathis Sanders' 
  [`SangerLCMFiltering`](https://github.com/MathijsSanders/SangerLCMFiltering) 
  a.k.a. [`hairpin`](https://confluence.sanger.ac.uk/display/CAS/hairpin). This 
  filters out SNVs based on ASMD and CLPM values to remove cruciform DNA 
  artifacts that commonly form in low input DNA sequencing (FILTER = PASS & 
  CLPM=0.00 & ASRD>=0.87). 

  b. Filter the `Pindel` indel output (FILTER = PASS).

3. Run [`cgpVAF`](https://confluence.sanger.ac.uk/pages/viewpage.action?pageId=22710418)
a.k.a. [`vafCorrect`](https://github.com/cancerit/vafCorrect). This revisits
variants in individual samples that were non-significant, but that were found to
be significant in other samples of the same donor. `cgpVAF` also facilitates
unbiased `pileup` (SNVs) and `exonerate` (indels)-based VAF calculation for the
union of variant sites in the set of related samples from each donor. 

4. Run [Sequoia](https://github.com/TimCoorens/Sequoia) for tree building.