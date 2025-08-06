# QTL project

## 1. Genetic data

We first process genetics data, due to
1) imputation (See imputation.md)
2) repeated measurements, e.g., one individual has multiple biopsy RNAseq
3) calcualte kinship matrix

Lets prepare post-imputated and QCed plink ped and map files,for each individual, lets generate 6 duplicated IDs which is enough for all samples (one individual has up to 5 samples)
```
input_ped = "Plink.merged.dose.filtered.newID.biallelic.clean.ped"
output_ped = "Plink.merged.dose.filtered.newID.biallelic.clean.duplicated.ped"

with open(input_ped) as fin, open(output_ped, "w") as fout:
    for line in fin:
        cols = line.strip().split()
        fid, iid = cols[0], cols[1]
        rest = cols[2:]
        # write original
        fout.write("\t".join([fid, iid] + rest) + "\n")
        # write duplicates
        for i in range(1,7):
            new_iid = f"{iid}_{i}"
            fout.write("\t".join([fid, new_iid] + rest) + "\n")
~                                                                
```

Okay, now we have duplicated plink files, but, this is not end, we need to manually match RNAseq data to plink data. For example, in ASCEND cohort, paired genetics-RNAseq is 408, so we have to extract these 408 samples to build new plink files which is the final genetics data we use in analysis

```
plink \
  --bfile /groups/g5840105/home/share/QTLproject/QTL_geneticData/Genetic_ASCEND/PLINK_408/plink.408sample
  --indep-pairwise 50 5 0.2
  --out pruned_snps

plink \
  --bfile /groups/g5840105/home/share/QTLproject/QTL_geneticData/Genetic_ASCEND/PLINK_408/plink.408sample
  --extract pruned_snps.prune.in
  --make-bed
  --out plink.408sample.prune

gemma --bfile plink.408sample.pruned -gk 1 -o IBS.matrix

```

## 2. Data prepare and harmonize

Before analysis, we need to harmonize genotype data, phenotype (covariates) data and expression (e.g., gene/isoform/AS) data.

（1） identify lanten factors (lets take gene data for example)

```
# import cinical data and gene count
meta_ASCEND=read.table("five_cohorts/ASCEND/clinial_data/ASCEN_clinical.harmonize.txt",header = T,stringsAsFactors = F)
data_gene_ASCEND=read.table("five_cohorts/ASCEND/gene/txi_gene_count_ascend_513.tsv",header = T,stringsAsFactors = F,row.names = 1)
data_gene_ASCEND=data_gene_ASCEND[rowSums(data_gene_ASCEND>0)>(0.8*ncol(data_gene_ASCEND)),]
data_gene_ASCEND <- DGEList(counts = data_gene_ASCEND, genes = rownames(data_gene_ASCEND))
data_gene_ASCEND <- calcNormFactors(data_gene_ASCEND)
data_gene_ASCEND <- edgeR::cpm(data_gene_ASCEND, log = T, prior.count = 0.1) 
data_gene_ASCEND=as.data.frame(t(data_gene_ASCEND))

library(MOFA2)
# lets check and impute phenotype first
meta_ASCEND=meta_ASCEND[,c("RNA_ID","age","sex","BMI","inflammation","location")]
colnames(meta_ASCEND)[1]='ID'
meta_ASCEND=meta_ASCEND[meta_ASCEND$ID %in% rownames(data_gene_ASCEND),]
which(is.na(meta_ASCEND))

tmp.matrix=data_gene_ASCEND
tmp.pheno=meta_ASCEND
rownames(tmp.pheno)=tmp.pheno$ID
tmp.matrix=tmp.matrix[rownames(tmp.matrix) %in% rownames(tmp.pheno),]
tmp.pheno=tmp.pheno[rownames(tmp.pheno) %in% rownames(tmp.matrix),]
tmp.matrix <- tmp.matrix[match(rownames(tmp.pheno), rownames(tmp.matrix)), ]

stopifnot(rownames(tmp.pheno)==rownames(tmp.matrix))
tmp.matrix=as.matrix(tmp.matrix)

tmp.matrix_clean <- apply(tmp.matrix, 2, function(feature) {
  residuals(lm(
    feature ~ tmp.pheno$age + tmp.pheno$sex + tmp.pheno$BMI + tmp.pheno$inflammation + tmp.pheno$location,
    na.action = na.exclude
  ))
})
tmp.matrix_clean <- t(tmp.matrix_clean)

tmp.omics=list(RNA=tmp.matrix_clean)
mofa_object <- create_mofa(data = tmp.omics)

train_opts <- get_default_training_options(mofa_object)
train_opts$seed <- 666  
data_opts <- get_default_data_options(mofa_object)
data_opts$scale_views=F
model_opts <- get_default_model_options(mofa_object)
model_opts$num_factors=10

mofa_object <- prepare_mofa(
  mofa_object,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
mofa_model <- run_mofa(mofa_object, use_basilisk = TRUE)

write.table(tmp.matrix,file = "QTL_files/ASCEND/RNA_gene.ASCEND.txt",row.names = T,quote = F,sep = "\t")
write.table(tmp.pheno,file = "QTL_files/ASCEND/Covariates.gene.ASCEND.txt",row.names = T,quote = F,sep = "\t")
```
Okay, now we have log transformed **gene matrix** and **covairates** (age, sex, BMI, inflammation and location)
（2）Lets harmonize genetic data with gene matrix.*Note*, we have to generate coupling file contains IDs of gene matrix and geneitc individuals

```
# change IDs
gene_ASCEND=read.table("QTL_files/ASCEND/RNA_gene.ASCEND.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
cov_gene_ASCEND=read.table("QTL_files/ASCEND/Covariates.gene.ASCEND.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
meta_ASCEND=read.table("five_cohorts/ASCEND/clinial_data/ASCEN_clinical.harmonize.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t")
plink_ASCEND=read.table("QTL_files/ASCEND/plink.final.fam",header = F,stringsAsFactors = F)
coupling_ASCEND=read.table("QTL_files/ASCEND/coupling_ASCEND.txt",header = T,stringsAsFactors = F)

plink_ASCEND$sample=gsub("^.*_FKDN","FKDN",plink_ASCEND$V2)

# make duplicates
# Add a new column with suffixes for duplicates
coupling_ASCEND <- coupling_ASCEND %>%
  group_by(ASA_ID) %>%
  mutate(
    new_ASA_ID = if (all(is.na(ASA_ID))) {  # If all ASA_IDs in group are NA
      NA_character_
    } else if (n() == 1) {  # If only 1 RNA_ID per ASA_ID (and not NA)
      ASA_ID
    } else {  # If duplicates exist
      case_when(
        row_number() == 1 ~ ASA_ID,
        TRUE ~ paste0(ASA_ID, "_", row_number() - 1)
      )
    }
  ) %>%
  ungroup()

# note, 3 ASA samples are filtered out 
coupling_ASCEND=merge(coupling_ASCEND,plink_ASCEND,by.x="new_ASA_ID",by.y="sample",all=F)
coupling_ASCEND=na.omit(coupling_ASCEND)

coupling_ASCEND=coupling_ASCEND[,c(1:4,6)]
colnames(coupling_ASCEND)=c("new_ASA_ID","RNA_ID","Patient_ID","ASA_ID","Plink_ID")
write.table(coupling_ASCEND,file="QTL_files/ASCEND/QTL.linkID.file.txt",quote = F,sep = "\t",row.names = F)

```

```
new_ASA_ID	RNA_ID	Patient_ID	ASA_ID	Plink_ID
FKDN240483734-1A	R24118680	CD0809	FKDN240483734-1A	77_FKDN240483734-1A
FKDN240483734-1A_1	R24118681	CD0809	FKDN240483734-1A	77_FKDN240483734-1A_1
FKDN240483735-1A	R24115452	CD0611	FKDN240483735-1A	78_FKDN240483735-1A<img width="565" height="65" alt="image" src="https://github.com/user-attachments/assets/b9281ad6-bbc6-4829-ad6c-574eb5f5f725" />
```

Okay, now we have **QTL.linkID.file.txt** file, which is used for link PLINK fam file and gene matrix. Then lets take residuals from linear regression

(3) Here, we regress out age, sex, BMI, inflammation and location, and 10 factors to identify broader eQTLs

```
gene_ASCEND=read.table("QTL_files/ASCEND/RNA_gene.ASCEND.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
cov_gene_ASCEND=read.table("QTL_files/ASCEND/Covariates.gene.ASCEND.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
cov_gene_ASCEND$ID=NULL
stopifnot(rownames(gene_ASCEND)==rownames(cov_gene_ASCEND))

tmp.residuals <- apply(gene_ASCEND, 2, function(feature) {
  residuals(lm(
    feature ~ cov_gene_ASCEND$age + cov_gene_ASCEND$sex + cov_gene_ASCEND$BMI + cov_gene_ASCEND$inflammation + cov_gene_ASCEND$location +
      cov_gene_ASCEND$Factor1 +
      cov_gene_ASCEND$Factor2 +
      cov_gene_ASCEND$Factor3 +
      cov_gene_ASCEND$Factor4 +
      cov_gene_ASCEND$Factor5 +
      cov_gene_ASCEND$Factor6 +
      cov_gene_ASCEND$Factor7 +
      cov_gene_ASCEND$Factor8 +
      cov_gene_ASCEND$Factor9 +
      cov_gene_ASCEND$Factor10,
    na.action = na.exclude
  ))
})
tmp.residuals=as.data.frame(tmp.residuals)
write.table(tmp.residuals,file = "QTL_files/ASCEND/residuals.ASCEND.gene.txt",row.names = T,quote = F,sep = "\t")
```

Okay, now upload **QTL.linkID.file.txt**, **residuals.ASCEND.gene.txt** to cluster

(4) On cluster, lets harmonize genetics and gene matrix

residuals.ASCEND.gene.txt should be like:
```
ENSG00000000003.16      ENSG00000000005.6
R24118608       -0.18438795110267       0.330104350271361
R24118607       -0.243135003150785      -1.18534908752482
```
QTL.linkID.file.txt should be like:
```
new_ASA_ID      RNA_ID  Patient_ID      ASA_ID  Plink_ID
FKDN240483734-1A        R24118680       CD0809  FKDN240483734-1A        77_FKDN240483734-1A
FKDN240483734-1A_1      R24118681       CD0809  FKDN240483734-1A        77_FKDN240483734-1A_1
FKDN240483735-1A        R24115452       CD0611  FKDN240483735-1A        78_FKDN240483735-1A
```

```
Must re-order the tissue order matching to plink.fam file, otherwise all the analysis would be wrong!!!!

1) extract individuals from link file, e.g., Link.ASCEND.ID, ASCEND has paired 408 samples with RNA and ASA
*** plink -- keep
2) cp the fam, plink.408sample.fam,

3) run Re-order.R, note, all fam sample must be included in expression table
4) check sample orders
   awk '{print $2}' fam > final_sample1.txt
   head -1 expression.txt | tr '\t' '\n' | tail -n +2 > final_sample2.txt

5) diff final_sample1.txt final_sample2.txt
   this should return empty
```

```
library(tidyverse)

# Read files
fam_samples <- read.table("plink.408sample.fam", header = FALSE) %>%
  unite("sample_id", V2, sep = "") %>%
  pull(sample_id)

cat("============= fam samples is",length(fam_samples),"\n")

exp_data <- read.table("residuals.test.txt", header = TRUE, check.names = FALSE,row.names=1)
#exp_samples <- colnames(exp_data)[-1]  # Assuming first column is gene_id

coupling=read.table("Link.ASCEND.ID",header=T)
exp_data=exp_data[rownames(exp_data) %in% coupling$RNA_ID,]
coupling=coupling[coupling$RNA_ID %in% rownames(exp_data),]

stopifnot(nrow(exp_data)==nrow(coupling))

coupling=coupling[order(coupling$RNA_ID),]
exp_data=exp_data[order(rownames(exp_data)),]

stopifnot(coupling$RNA_ID==rownames(exp_data))

rownames(exp_data)=coupling$Plink_ID

exp_samples <- rownames(exp_data) 

cat("============= expression samples is",length(exp_samples),"\n")

# Identify mismatches
missing_in_fam <- setdiff(exp_samples, fam_samples)  # In expression but not genotypes
missing_in_exp <- setdiff(fam_samples, exp_samples)  # In genotypes but not expression

# Report mismatches
if (length(missing_in_fam) > 0) {
  writeLines("\n[WARNING] Samples in expression BUT NOT in .fam file:")
  print(missing_in_fam)
  writeLines("")
}

if (length(missing_in_exp) > 0) {
  writeLines("\n[WARNING] Samples in .fam file BUT NOT in expression:")
  print(missing_in_exp)
  writeLines("")
}

# Keep only common samples (strict mode)
common_samples <- intersect(fam_samples, exp_samples)
# Keep only common samples and reorder
exp_reordered <- exp_data[common_samples, , drop = FALSE]

# Transpose to make genes rows and samples columns
exp_final <- t(exp_reordered)


# Save results
write.table(exp_final, "expression_reordered.txt", sep = "\t", 
            quote = FALSE, row.names = T)

# Save list of mismatched samples (optional)
writeLines(c(missing_in_fam, missing_in_exp), "mismatched_samples.txt")

```

Okay, now we have **expression_reordered.txt**

## 3. gene-snp pairs data

For each gene, we have to know what SNPs we should use in QTL analysis

```

To perform e/s/iQTL analysis, we need to generate a file containing paired SNP-Gene;

1) generate a bed file to point out the gene positions from genecode.gtf file for each ensemble gene name
2) generate a bed file to point out all the SNP positions from PLINK (after imputation) file
3) use bedtools to select, e.g. SNPs around the TSS +/- 500,000 bp for each ensemble gene

bedtools window -a genecode.bed -b Imputation.bed -w 500000 > paired.SNP.Gene.bed

NOTE: the following is to define cis-gene window, then use this in bedtools command

awk 'BEGIN{OFS="\t"} 
$3 == "gene" {
  chr = $1;
  start = $4;
  end = $5;
  strand = $7;
  center = (strand == "+" ? start : end);
  
  # Extract gene_id and gene_name using regex
  match($0, /gene_id "([^"]+)"/, a);
  match($0, /gene_name "([^"]+)"/, b);
  gene_id = a[1];
  gene_name = b[1];
  
  print chr, center-1, center, gene_id, gene_name, strand
}' gencode.v38.annotation.gtf > genes_center.bed


the code above is how to extract position for a single gene, here, the gene center is defined as TSS region,
if the strand is +, the TSS would start from $4, namely the starting postion; if strand is -, the TSS would start from $5.

4) after extract paired SNPs from whole plink file, then we have to hand-make a new plink genotype file to containing all tissues from
   RNAseq data; because we only have 188 ASA, but more than 500 RNAseq, so, IMPORTANTLY!!!! make sure that tissue IDs follow the same orders
   for both plink.fam and gene/AS/isoform expression matrix
```

Okay, now we have **Paired.SNP.gene** file, looks like:
```
ENSG00000228463.10 rs200188737
ENSG00000228463.10 rs61769351
ENSG00000228463.10 rs745495619
ENSG00000228463.10 rs142559957
ENSG00000228463.10 rs74879860
```

## 4. Finally, lets set up QTL analysis using GEMMA

```
#!/bin/bash
#SBATCH --job-name=Job.aa
#SBATCH --error=Job.aa.err
#SBATCH --output=Job.aa.out
#SBATCH --mem=15gb
#SBATCH --cpus-per-task=6


# Before running, several things need to check

# 1) the expression matrix, make sure the row is gene, column is sample; the order of sample must be the same as plink.fam: expression_reordered.txt

# 2) IBS.matrix, should be re-genrated if plink file is changed

# 3) genotype, the one used to match expression matrix: plink files

# 4) Paired.snp.bed, be careful with hg19 or hg38

path_IBS="/groups/g5840105/home/share/QTLproject/pipeline_QTL/IBS_matrix/ASCEND/output/"
path_geno="/groups/g5840105/home/share/QTLproject/QTL_geneticData/Genetic_ASCEND/PLINK_408/"
path_pheno="/groups/g5840105/home/share/QTLproject/pipeline_QTL/Phenotype_matrix/ASCEND/Gene/"
path_tmp="/groups/g5840105/home/share/QTLproject/pipeline_QTL/TMP_files/"
path_eQTL="/groups/g5840105/home/share/QTLproject/pipeline_QTL/Result_eQTL/"

#gemma -bfile $path_geno/all_duplicates -gk 1 -o IBS.matrix
#mv output/IBS.matrix* $path_IBS/

source /opt/app/anaconda3/bin/activate
conda activate bcftools_env

log_file="$path_tmp/QTLlog.ASCEND.gene.log"
# Create log header
echo -e "Gene\tSNP_Count\tStatus" > "$log_file"

cat /groups/g5840105/home/share/QTLproject/pipeline_QTL/Phenotype_matrix/ASCEND/Gene/Batch/split_file_aa.txt | while read gene
do
  # Extract expression for the gene
  awk -v g=$gene 'NR==1 || $1==g' $path_pheno/expression_reordered.txt > $path_tmp/tmp.$gene.txt

echo -e "======================================================================================= start $gene ==="

# extract SNPs
LC_ALL=C grep -wF "$gene" ///groups/g5840105/home/share/QTLproject/pipeline_QTL/Gene_rsID/Pair_ASCEND/Paired.SNP.gene | awk '{print $2}' > "$path_tmp/tmp.$gene.SNP"

# Count SNPs
    snp_count=$(wc -l < "$path_tmp/tmp.$gene.SNP" | awk '{print $1}')
    
    if [ "$snp_count" -eq 0 ]; then
        echo "WARNING: No SNPs found for gene $gene - skipping"
        echo -e "$gene\t0\tSkipped (no SNPs)" >> "$log_file"
        # Clean up temporary files
        rm -f "$path_tmp/tmp.$gene.txt" "$path_tmp/tmp.$gene.SNP"
        continue
    fi
    
    echo "Found $snp_count SNPs for gene $gene"
    echo -e "$gene\t$snp_count\tProcessing" >> "$log_file"


# make tmp genotype
plink --bfile $path_geno/plink.408sample --extract $path_tmp/tmp.$gene.SNP --make-bed --out $path_tmp/tmp.$gene.genotype --allow-no-sex --threads 8

# check sample orders
if ! diff <(awk '{print $2}' "$path_tmp/tmp.$gene.genotype.fam") \
          <(head -1 "$path_tmp/tmp.$gene.txt" | tr '\t' '\n' | tail -n +2); then
  echo "ERROR: Sample mismatch detected between .fam and .txt files!"
  exit 1
fi

echo -e "======================================================================================= Checking sample orders before GEMMA === Ture ==="

# convert phenotype
awk 'NR==2 {for (i=1; i<=NF; i++) print $i}' $path_tmp/tmp.$gene.txt > $path_tmp/tmp.$gene.exp.txt
tail -n +2 $path_tmp/tmp.$gene.exp.txt > $path_tmp/tmp.$gene.txt
mv $path_tmp/tmp.$gene.txt $path_tmp/tmp.$gene.exp.txt

echo -e "======================================================================================= GEMMA ====="

# Run GEMMA LMM
  gemma -bfile $path_tmp/tmp.$gene.genotype \
    -p $path_tmp/tmp.$gene.exp.txt \
    -k $path_IBS/IBS.matrix.cXX.txt \
    -lmm 4 \
    -o eQTL_${gene}.result

mv output/eQTL_${gene}.result* $path_eQTL/
rm $path_tmp/tmp.$gene*

echo -e "======================================================================================= GEMMA $gene done ==="
done 
```

Okay, for each gene feature, we have two outputs

**eQTL_ENSG00000103061.13.result.assoc.txt**
```
chr     rs      ps      n_miss  allele1 allele0 af      beta    se      logl_H1 l_remle l_mle   p_wald  p_lrt   p_score
16      rs76230155      67812791        0       T       A       0.022   -2.375880e-02   3.362992e-02    2.696864e+02    1.362581e+00    1.354717e+00    4.802962e-01    4.793
16      rs35316276      67816797        0       T       C       0.109   -4.929886e-02   1.593310e-02    2.741563e+02    1.211374e+00    1.204070e+00    2.110566e-03    2.122
16      rs75334109      67818742        0       T       C       0.036   -2.028597e-02   2.694611e-02    2.697204e+02    1.363160e+00    1.355064e+00    4.519853e-01    4.509
```
**eQTL_ENSG00000103061.13.result.log.txt**
```
## GEMMA Version    = 0.98.3 (2020-11-28)
## Build profile    = 
## GCC version      = 7.5.0
## GSL Version      = 1.16
## OpenBlas         = OpenBLAS 0.3.12  - OpenBLAS 0.3.21 DYNAMIC_ARCH NO_AFFINITY SkylakeX MAX_THREADS=128
##   arch           = SkylakeX
##   threads        = 6
##   parallel type  = threaded
##
## Command Line Input = gemma -bfile /groups/g5840105/home/share/QTLproject/pipeline_QTL/TMP_files//tmp.ENSG00000103061.13.genotype -p /groups/g5840105/home/share/QTLproject
##
## Date = Tue Aug  5 15:46:35 2025
##
## Summary Statistics:
## number of total individuals = 408
## number of analyzed individuals = 408
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs/var = 1550
## number of analyzed SNPs/var = 1516
## REMLE log-likelihood in the null model = 268.2
## MLE log-likelihood in the null model = 269.436
## pve estimate in the null model = 0.157143
## se(pve) in the null model = 0.0649112
## vg estimate in the null model = 0.0183289
## ve estimate in the null model = 0.013448
## beta estimate in the null model =   -0.000127553
## se(beta) =   0.00574115
##
## Computation Time:
## total computation time = 0.0378333 min 
## computation time break down: 
##      time on eigen-decomposition = 0.00166667 min 
##      time on calculating UtX = 0.00733333 min 
##      time on optimization = 0.0183333 min 
```
Okay, we also have a log file, **QTLlog.ASCEND.gene.log**
```
Gene    SNP_Count       Status
ENSG00000000003.16      0       Skipped (no SNPs)
ENSG00000100842.13      2183    Processing
ENSG00000118960.13      2391    Processing
ENSG00000185658.14      3292    Processing
ENSG00000170006.13      0       Skipped (no SNPs)
ENSG00000224051.7       2273    Processing
ENSG00000280893.1       989     Processing
ENSG00000256128.8       0       Skipped (no SNPs)
```





















