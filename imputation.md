
1) two plink files from ASA, one is forward and the other one is raw, 
we need to use the forward one since it has been aligned to the same trand (forward strand)

2) in total, 743,722 SNPs in 188 samples, we need to do pre-imputation filtering

*** plink -- select biallelic snps only
(743,722 SNPs and 188 samples left)

*** plink -- remove calling rate <98% samples, e.g. one sample has less than 98% SNPs are present
(743,722 SNPs and 185 samples left)

*** plink -- remove SNPs with missing rate >5%
(733,059 SNPs and 185 samples left)

*** plink -- remove SNPs with MAF <1%, NOTE. this step filters out lots of loci
(504,414 SNPs and 185 samples left)

*** plink -- remove SNPs with HWE <1e-06
(503,336 SNPs and 185 samples left)

3) then, it is important to align the locus to reference genome, here we use 1000G

*** Database preparation, wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
*** Tool preparation, download HRC-1000G-check-bim.pl from https://www.chg.ox.ac.uk/~wrayner/tools/
*** PLINK preparation, plink -- calulate freq

*** perl HRC-1000G-check-bim.pl -b input.bim -f freq.freq -r HRC-1000G-check-bim -g -p EAS

*** lots of files will be genarted, e.g., Exclude-data.bia.mind.geno.maf.hwe-1000G.txt
*** this file contains bad SNPs, including #non-match with 1000G, #freq different from 1000G, #palindromic SNPs and # indels and ambigous SNPs
*** plink -- exclude the SNPs in this file
(491,989 SNPs and 185 samples left)

ONCE you run the perl command, there is a Run-plink.sh file, JUST bash Run-plink.sh, it contains the following steps
See https://genepi.github.io/michigan-imputationserver/prepare-your-data/

*** plink -- flip strand and plink -- force alleles and plink -- update and plink -- convert to vcf

4) then submit all vcf.gz and tabix files to server
   choose references and click run

Here shows an example of statistics:

******************************************************************************
Input Validation
22 valid VCF file(s) found.

Samples: 185
Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9
SNPs: 475198
Chunks: 292
Datatype: unphased
Build: hg19
Reference Panel: apps@topmed-r3@1.0.0 (hg38)
Population: all
Phasing: eagle
Mode: imputation


Quality Control
Uploaded data is hg19 and reference is hg38.

Lift Over

Calculating QC Statistics


Statistics:
Alternative allele frequency > 0.5 sites: 87,834
Reference Overlap: 87.36 %
Match: 413,953
Allele switch: 753
Strand flip: 0
Strand flip and allele switch: 0
A/T, C/G genotypes: 10
Filtered sites:
Filter flag set: 0
Invalid alleles: 0
Multiallelic sites: 0
Duplicated sites: 0
NonSNP sites: 0
Monomorphic sites: 0
Allele mismatch: 57
SNPs call rate < 90%: 0

Excluded sites in total: 810
Remaining sites in total: 413,963
See snps-excluded.txt for details
Typed only sites: 59,995
See typed-only.txt for details
************************************************************************************

5) after running, download all the log files and md5 vcf files

*** check md5, then QC
for i in {22..22}; do plink2 --vcf chr${i}.dose.vcf --make-pgen --out Pgen.chr${i}.dose; plink2 --pfile Pgen.chr${i}.dose --maf 0.01 --hwe 1e-06 --geno 0.05 --mind 0.05 --ex

