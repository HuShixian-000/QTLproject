**QTL project**

*1. Data prepare and harmonize

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
Lets harmonize genetic data with gene matrix.*Note*, we have to generate coupling file contains IDs of gene matrix and geneitc individuals

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

Okay, now we have **QTL.linkID.file.txt** file



