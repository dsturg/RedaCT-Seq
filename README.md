# RedaCT-Seq
Reduction of ac4C and NaBH4 sequencing

Code for processing and analysis of Reduction of ac4C and NaBH4 sequencing (RedaCT-Seq) data.
This protocol identifies RNA modifications (specifically, ac4C) by analyzing basecalling changes following reduction of ac4C with Sodium Borohydride (NaBH4).


<!-- GETTING STARTED -->
## Getting Started

The steps below will guide through the steps required to perform the analysis. It's recommended to perform this analysis in a cluster or high-performance computing environment

### Prerequisites

This is a list of components needed
* [samtools](http://www.htslib.org/) v.1.10 or 1.11
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/) v.4.0
* [STAR aligner] (https://github.com/alexdobin/STAR) v.2.7.5a
* [mpileup2readcounts] (https://github.com/IARCbioinfo/mpileup2readcounts)
* [R] (https://www.r-project.org/)
* [R-studio] (https://github.com/IARCbioinfo/mpileup2readcounts)



<!-- WORKFLOW -->
## Retrieve the data
   ```sh
   https://figshare.com/s/d7ab88c65c69ec0e097d
   ```
## Workflow
1. Perform the pileup
   ```sh
   samtools mpileup -A -R -Q20 -C0 -d 100000 --ff UNMAP,SECONDARY,QCFAIL,DUP -f /data/Soberlab/indexes/STAR/hg19_UCSC/ref.fa WT.BH4.bam KO.BH4.chr19.bam WT.Ctrl.chr19.bam | sed 's/        /    *     */g' | mpileup2readcounts 0 -5 true 0 0 > mpileup_output/mpileup_output.txt ;
   ```
2. Parse the pilup results
  ```sh
    perl redact_parse_script.pl mpileup_output/mpileup_output.txt 3 > mpileup_output/mpileup_output_parsed.txt
   ```
3. Open results in R




<!-- CONTACT -->
## Contact

Daniel Arango - dany33co@gmail.com
or Dave Sturgill
dave.sturgill@gmail.com
