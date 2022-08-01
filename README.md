# RedaCT-Seq
Reduction of ac4C and NaBH4 sequencing

Code for processing and analysis of Reduction of ac4C and NaBH4 sequencing (RedaCT-Seq) data.
This protocol identifies RNA modifications (specifically, ac4C) by analyzing basecalling changes following reduction of ac4C with Sodium Borohydride (NaBH4).


<!-- GETTING STARTED -->
## Getting Started

The steps below will guide through the steps required to perform the analysis. It's recommended to perform this analysis in a cluster or high-performance computing environment

### Prerequisites

This is a list of components needed:
* [samtools](http://www.htslib.org/) v.1.11
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/) v.1.6
* [STAR aligner] (https://github.com/alexdobin/STAR) v.2.7.5a
* [mpileup2readcounts] (https://github.com/IARCbioinfo/mpileup2readcounts)
* [R] (https://www.r-project.org/) v.4.1.1 
* [R-studio] (https://www.rstudio.com/) v.2022.02.3 Build 492
R packages required:
* data.table v.1.14.2
Matt Dowle and Arun Srinivasan (2021)
https://CRAN.R-project.org/package=data.table
* dplyr v.1.0.8
Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2022). dplyr: A Grammar of Data Manipulation.
https://CRAN.R-project.org/package=dplyr 
* Genomation v1.24.0
Akalin A, Franke V, Vlahovicek K, Mason CE, Schubeler D, (2014) genomation: a toolkit to summarize, annotate and visualize genomic intervals. Bioinformatics. doi: 10.1093/bioinformatics/btu775
https://bioconductor.org/packages/release/bioc/html/genomation.html
* GenomicFeatures v1.44.2
Lawrence M, Huber W, Pages H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
* rtracklayer v.1.52.1
M. Lawrence, R. Gentleman, V. Carey: "rtracklayer: an {R} package for interfacing with genome browsers". Bioinformatics
  25:1841-1842.
https://bioconductor.org/packages/release/bioc/html/rtracklayer.html
* stringr v.1.4.0
Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String Operations.
https://CRAN.R-project.org/package=stringr




<!-- WORKFLOW -->
## Retrieve the data
The repository below includes sample pileup data, corresponding to the input of the parsing step #2 in the Workflow section
   ```sh
   https://figshare.com/s/d7ab88c65c69ec0e097d
   (Replace with DOI when complete)
   ```
## Workflow
1. Perform the pileup
   ```sh
   samtools mpileup -A -R -Q20 -C0 -d 100000 --ff UNMAP,SECONDARY,QCFAIL,DUP -f /data/Soberlab/indexes/STAR/hg19_UCSC/ref.fa WT.BH4.bam KO.BH4.chr19.bam WT.Ctrl.chr19.bam | sed 's/        /    *     */g' | mpileup2readcounts 0 -5 true 0 0 > mpileup_output/mpileup_output.txt ;
   ```
   To reduce file size, you could require a minimum depth of 10 in each sample, for example:
   ```sh
   samtools mpileup -A -R -Q20 -C0 -d 100000 --ff UNMAP,SECONDARY,QCFAIL,DUP -f /data/Soberlab/indexes/STAR/hg19_UCSC/ref.fa WT.BH4.bam KO.BH4.chr19.bam WT.Ctrl.chr19.bam | sed 's/        /    *     */g' | mpileup2readcounts 0 -5 true 0 0 | awk '$4 >= 10 && $15 >= 10 && $26 >= 10' > mpileup_output/mpileup_output_chr19_min10.txt ;
   ```
2. Parse the pileup results
  ```sh
    scripts/redact_parse_script.pl sampledata/mpileup_output_chr19_min10.txt 3 > sampledata/mpileup_output_chr19_min10_parsed.txt ;
   ```
3. Open results in R:  Downstream analysis is performed in [R] (https://www.r-project.org/).  Follow the example in the [Markdown file](RedaCT-Seq.md)


<!-- CONTACT -->
## Contact
Dave Sturgill - dave.sturgill@gmail.com
or
Daniel Arango - dany33co@gmail.com

