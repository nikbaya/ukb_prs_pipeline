#!/usr/bin/env Rscript

# Following instructions from the following page:
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

# install.packages("bigsnpr")
library("bigsnpr")

######################
# Read in genotypes
######################

# Read in UKB imputed v3 ref panel of 5k EUR individuals
# Only need to run once if rds file already exists
# Old version: 
# snp_readBed("/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ref/ukb_wes_200k/ukb_ref_panel/data/plink/ukb_imp_v3_eur_ref_panel_5k_chr22.bed") 
# New version (10k variants):
# snp_readBed("/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ref/ukb_wes_200k/ukb_ref_panel/data/plink/tmp-ukb_imp_v3_eur_ref_panel_5k_chr22-ldpred2_10k_variants.bed") 

obj.bigSNP <- snp_attach("/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ref/ukb_wes_200k/ukb_ref_panel/data/plink/ukb_imp_v3_eur_ref_panel_5k_chr22.rds")
# obj.bigSNP <- snp_attach("/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ref/ukb_wes_200k/ukb_ref_panel/data/plink/tmp-ukb_imp_v3_eur_ref_panel_5k_chr22-ldpred2_10k_variants.rds")

# Rename columns
map = setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a0", "a1"))

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()

# Old version: Choose validation and test sets of individuals
# set.seed(1)
# n_samples = nrow(G)
# ind.val <- sample(n_samples, round(n_samples*(3/4))) # choose 3/4 of the samples to be in the validation set
# ind.test <- setdiff(rows_along(G), ind.val) # remaining 1/4 is assigned to test set

######################
# Read in sumstats
######################

# Read external summary statistics
# Ideal columns: 
# chr, pos, a0, a1, beta, beta_se

# Old attempt:
# sumstats <- bigreadr::fread2("/gpfs3/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/prs_cs/input_sumstats/t2d.tsv")
# > str(sumstats)
# 'data.frame':    13791467 obs. of  5 variables:
# $ SNP : chr  "rs6602381" "rs7899632" "rs61875309" "rs573719411" ...
# $ A1  : chr  "A" "A" "A" "C" ...
# $ A2  : chr  "G" "G" "C" "T" ...
# $ BETA: num  -0.000427 -0.0006 0.000332 -0.011798 0.003944 ...
# $ P   : num  0.3965 0.2337 0.5911 0.0981 0.6152 ...
#
# Need to add chr and pos columns
# variant_id_map = bigreadr::fread2("/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/neale_lab_round2/variant_id_map.chr22.tsv")
# Merge on "rsid" column
# sumstats = merge(sumstats, variant_id_map, by="rsid")
# sumstats = na.omit(sumstats) # remove rows with NAs

# New attempt
ss = bigreadr::fread2("/gpfs3/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/neale_lab_round2/2443.gwas.imputed_v3.both_sexes.chr22.tsv")
colnames(ss) = c("variant", "minor_allele", "maf", "exp_case_mac", 
		 "low_confidence", "n_complete_samples", "ac", "ytx", 
		 "beta", "beta_se", "tstat", "pval")
ss = na.omit(ss)
ss = ss[,c("variant", "n_complete_samples", "beta", "beta_se")]
library(stringr)
cols = data.frame(str_split_fixed(ss$variant, ":", 4))
ss$chr = as.numeric(as.character(cols$X1))
ss$pos = as.numeric(as.character(cols$X2))
ss$a0 = as.character(cols$X3)
ss$a1 = as.character(cols$X4)

# simplify by assuming "n_eff" is same as "n_complete_samples" (not true for binary trait like t2d)
colnames(ss)[2] = "n_eff"

# Note that we needed to give specific names to columns in the sumstats file to
# avoid the following error when using `snp_match` later:
# "Error: Please use proper names for variables in 'sumstats'. Expected 'chr, pos, a0, a1, beta'."

# Restrict to 1000 variants to speed up test
# sumstats_head = sumstats[c(1:1000),] # TODO: Remove when done with testing

# Match variants
# New code: Match all variants on chr22
df_beta <- snp_match(ss, map)
# 179,766 variants to be matched.
# 23,302 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 156,131 variants have been matched; 0 were flipped and 156,131 were reversed.

# Subset df_beta to 10,000 variants for simplicity 
df_beta = df_beta = df_beta[c(1:10000),]

# Old code that tried subsetting to 1000 variants
# df_beta <- snp_match(sumstats_head, map)
# 1,000 variants to be matched.
# 3 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 993 variants have been matched; 0 were flipped and 993 were reversed.

# Old code that did not subset to 1000 variants
# df_beta <- snp_match(sumstats, map)
# > df_beta <- snp_match(sumstats, map)
# 180,078 variants to be matched.
# 23,363 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 156,036 variants have been matched; 0 were flipped and 156,036 were reversed.

# Get genomic position
# In a login node:
# $ cd /gpfs3/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/tmp
# $ wget https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_OMNI/chr22.OMNI.interpolated_genetic_map.gz
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/well/lindgren/UKBIOBANK/nbaya/resources/prs/LDpred2/gmap", ncores = NCORES)

######################
# Calculate correlation
######################

tmp <- tempfile(tmpdir = "/gpfs3/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/tmp")

# Only use chr22 for this test
chr=22

ind.chr <- which(df_beta$chr == chr)
## indices in 'G'
ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]

corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
		 infos.pos = POS2[ind.chr2], ncores = NCORES)

ld <- Matrix::colSums(corr0^2)
corr <- as_SFBM(corr0, tmp, compact = TRUE)


cat(paste("File size (GB):", file.size(corr$sbk) / 1024^3))  # file size in GB




######################
# Calculate heritability
######################

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
				sample_size = n_eff, blocks = NULL)))

h2_est = ldsc[["h2"]]


######################
# Run LDpred2-auto
######################

multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
			       vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
			       ncores = NCORES)
str(multi_auto, max.level = 1)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

saveRDS(beta_auto, file = "/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data/tmp/ldpred2/beta_auto.rds")
