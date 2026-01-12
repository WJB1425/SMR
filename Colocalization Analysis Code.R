library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubr)
library(data.table)
library(vcfR)


# Colocalization Analysis Example: eQTLs Data - [QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]----
dat1 <- fread('[RESULTS_DIR]/smr/smr_eQTLs_[EXPOSURE_ID]_[OUTCOME_ID].msmr') # Filtered SMR analysis results

## If the GWAS data has duplicate SNP entries, remove the duplicates
GWAS <- fread('[RAW_DATA_DIR]/[GWAS_INPUT_FILE].tsv',sep='\t',header=TRUE,na.strings = "",fill=TRUE )
# The column names here must match the actual input file. The user's code renames them to the structure below:
# colnames(GWAS) = c('CHR','pos','SNP','A2','A1','beta','se','freq','n','p') 
colnames(GWAS) = c('CHR','pos','SNP','A1','A2','freq','beta','se','p')
GWAS <- GWAS[which(GWAS$SNP!=''),]
GWAS <- GWAS[which(GWAS$SNP !='NA'),]
GWAS <- unique(GWAS)
GWAS <- GWAS[!duplicated(GWAS$SNP),]
GWAS <- GWAS[which(!is.na(GWAS$SNP)),]
GWAS <- separate_rows(GWAS, SNP, sep = ',')

gwas1 = GWAS[, c('pos','CHR','SNP','beta','se','A1','A2','freq','p')]
colnames(gwas1) = c('base_pair_location','chromosome','ID','beta','se','effect_allele','other_allele','effect_allele_frequency','p_value')
gwas1$MAF = ifelse(gwas1$effect_allele_frequency<0.5, gwas1$effect_allele_frequency, 1 - gwas1$effect_allele_frequency)
gwas1$beta = as.numeric(gwas1$beta)
gwas1$se = as.numeric(gwas1$se)

head(gwas1)


gwas2 <- fread('[QTL_COLOC_DATA].txt') # The complete QTL data from which key results were extracted
## Add column names
colnames(gwas2) <- c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","SE","p")
head(gwas2)
gwas2$MAF1 = ifelse(gwas2$Freq<0.5, gwas2$Freq, 1 - gwas2$Freq)
gwas2 = gwas2[order(gwas2$p),]
# gwas2 = gwas2[!duplicated(gwas2$SNP),]  # A single SNP may correspond to multiple genes/sites/proteins; do not remove duplicates in QTL data

output_table <- data.frame()
for (i in unique(dat1$probeID)) { # Loop through each filtered result for colocalization
  SNP <- dat1$topSNP[dat1$probeID == i]
  CHR <- dat1$ProbeChr[dat1$probeID == i]
  BP <- dat1$Probe_bp[dat1$probeID == i]
  gwas1_c <- gwas1[gwas1$chromosome == CHR,]
  head(gwas1_c)
  
  gwas2_c <- gwas2[gwas2$Chr == CHR,]
  # Filter SNPs within 1000 kb upstream and downstream of the key probeID in the QTL data
  gwas2_c <- gwas2_c[gwas2_c$BP >  BP - 1000000 &
                       gwas2_c$BP < BP + 1000000,]
  head(gwas2_c)
  gwas2_c <- gwas2_c[which(gwas2_c$Probe == i),] ## Ensure that the SNPs in the QTL data used for colocalization affect the target probe
  dat_merge <- merge(gwas1_c,gwas2_c,
                     by.x = 'variant_id',by.y = 'SNP', # Note: 'variant_id' might be an alias for the previously defined 'ID' column
                     suffixes = c("_gwas1","_gwas2"))
  head(dat_merge)
  
  # Handling GWAS P-values of 0
  table(is.na(dat_merge$p_value))
  min(dat_merge$p_value)
  dat_merge$p_value[dat_merge$p_value==0] <- min(dat_merge$p_value[dat_merge$p_value!=0])
  
  # Handling QTL P-values of 0
  table(is.na(dat_merge$p))
  min(dat_merge$p)
  dat_merge$p[dat_merge$p==0] <- min(dat_merge$p[dat_merge$p!=0])
  
  
  # Colocalization Analysis
  library("coloc")
  
  ngwas = [GWAS_N_CASES] + [GWAS_N_CONTROLS]  # GWAS sample size (modify according to the GWAS cohort)
  case = [GWAS_N_CASES] # Number of cases (modify according to the GWAS cohort)
  neqtl = [QTL_N_SAMPLES] # QTL sample size (e.g., eQTL N=31684, mQTL N=1980, pQTL N=10,708). Must be confirmed based on the data used.
  
  # Cohort information example (for reference/modification only):
  # [QTL_STUDY_NAME] | [PMID] | [AUTHOR] et al. | [POPULATION] | [N_SAMPLES]
  
  result <- coloc.abf(
    dataset1 = list(pvalues=dat_merge$p_value, # GWAS pvalue
                    snp=dat_merge$variant_id, # Note: 'variant_id' might be an alias for the previously defined 'ID' column
                    type="cc", N=ngwas, s=case/ngwas,
                    MAF = dat_merge$MAF),
    dataset2 = list(pvalues=dat_merge$p, # QTL pvalue
                    snp=dat_merge$variant_id, # Note: 'variant_id' might be an alias for the previously defined 'ID' column
                    type="quant", N=neqtl,
                    MAF = dat_merge$MAF1),
    p12 = 5e-05) # p12 default 1e-05. Setting to 5e-05 may yield better PPH4 results than the default.
  
  tmp <- data.frame(ProbID = i,
                    PPH4 = result$summary[6], # This value is used for the PPH4 in the forest plot
                    SNP.PP.H4 = max(result$results$SNP.PP.H4))
  
  
  output_table <- rbind(tmp, output_table) # Save the PPH4 value for each result
  
  if (!file.exists('[RESULTS_DIR]/coloc_result')) {
    dir.create('[RESULTS_DIR]/coloc_result')
  }
  
  if(result$summary[6] > 0.5){
    
    need_result=result$results # Output all SNPs involved in the colocalization analysis
    write.table(need_result,
                paste0("[RESULTS_DIR]/coloc_result/coloc_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",SNP,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    
    
    gwas = cbind(dat_merge$variant_id,dat_merge$p_value) # Note: 'variant_id' might be an alias for the previously defined 'ID' column
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(dat_merge$variant_id,dat_merge$p) # Note: 'variant_id' might be an alias for the previously defined 'ID' column
    colnames(eqtl) = c("rsid","pval")
    if (!file.exists('[DATA_DIR]/locusdata')) {
      dir.create('[DATA_DIR]/locusdata')
    }
    write.table(gwas,
                paste0("[DATA_DIR]/locusdata/gwas_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",i,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,
                paste0("[DATA_DIR]/locusdata/qtl_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",i,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    
    # Note: locuscompare function needs to be loaded or part of a package
    locuscompare(paste0("[DATA_DIR]/locusdata/gwas_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",i,".tsv"),
                 paste0("[DATA_DIR]/locusdata/qtl_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",i,".tsv"),
                 title1 = 'GWAS',title2 = 'eQTL',# title1 is the x-axis label. title2 should be modified based on the actual QTL type (e.g., 'eQTL', 'pQTL')
                 legend =T)
    ggsave(paste0("[FIGURE_DIR]/locus_cis-[QTL_TYPE]_[PROJECT_ID]_[OUTCOME_ID]-",i,".jpeg"),
           width=10,height = 6)
  }
  
}