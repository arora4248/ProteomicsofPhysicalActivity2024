
### MR of traits (PA) on proteins

library(data.table)
library(R.utils)
library(TwoSampleMR)
library(ieugwasr)

protein_files<-list.files("filepath_proteomics_gwas_summary_stats")
gwas_files<-list.files("filepath_PAtraits_gwas_summary_stats")

# j is a vector for number of PA tarits in order
j=9

#outmarix for results output 
Out = matrix(nrow=(2940), ncol=28, NA)
Out[,1] = protein_files
Out[,2] = sub("\\_.*", "", protein_files)
Out[,3] = sub("\\..*", "", gwas_files[j])
colnames(Out)= c("file_name","exposure","outcome","method1","nsnp1",
                 "b1","se1","pval1", "method2","nsnp2",
                 "b2","se2","pval2", "method3","nsnp3",
                 "b3","se3","pval3", "method4","nsnp4",
                 "b4","se4","pval4",
                 "method5","nsnp5","b5","se5","pval5")

### loop for MR
for (i in 1:length(protein_files)){
  protein<-fread(paste0("filepath_proteomics_gwas_summary_stats/",protein_files[i]),header=T)
  gwas<-fread(paste0(file="filepath_traits_gwas_summary_stats/",gwas_files[j]),header=T)
  
  ### select instrument variables
  gwas<-gwas[-which(gwas$P>5e-8),]
  
  if (nrow(gwas>0)){
    
    ## outcome
    out<-format_data(
      protein,
      type = "outcome",
      header = TRUE,
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P",
      samplesize_col = "N")
    
    ##exposure
    exp<-format_data(
      gwas,
      type = "exposure",
      header = TRUE,
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "P",
      samplesize_col = "N")
    
    
    ### clump locally
    
    possibleError<-tryCatch(tmp<- ld_clump(
      dplyr::tibble(rsid=exp$SNP, pval=exp$pval.exposure, id=exp$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "filepath_to_local_library/TwoSampleMR_local_ld/EUR"
    ), error=function(e) e)
    
    if(inherits(possibleError, "error")) next
    
    exp2<-exp[(exp$SNP%in%tmp$rsid),]
    if(nrow(exp2>0)) {
      
      dat <- harmonise_data(exp2,out,action = 1)
      
      if (nrow(dat)>0) {
        
        mr_results <- mr(dat, method_list = c("mr_egger_regression", "mr_wald_ratio"))
        mr_results
        
        mr_results2 <- mr(dat, method_list = c("mr_weighted_median", "mr_wald_ratio"))
        mr_results2
        
        mr_results3 <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
        mr_results3
        
        mr_results4 <- mr(dat, method_list = c("mr_simple_mode", "mr_wald_ratio"))
        mr_results4
        
        mr_results5 <- mr(dat, method_list = c("mr_weighted_mode", "mr_wald_ratio"))
        mr_results5
        
        if (nrow(mr_results>0)){
          
          Out[i,4] <- mr_results[1,5]
          Out[i,5] <- mr_results[1,6]
          Out[i,6] <- mr_results[1,7]
          Out[i,7] <- mr_results[1,8]
          Out[i,8] <- mr_results[1,9]
          Out[i,9] <- mr_results2[1,5]
          Out[i,10] <- mr_results2[1,6]
          Out[i,11] <- mr_results2[1,7]
          Out[i,12] <- mr_results2[1,8]
          Out[i,13] <- mr_results2[1,9]
          Out[i,14] <- mr_results3[1,5]
          Out[i,15] <- mr_results3[1,6]
          Out[i,16] <- mr_results3[1,7]
          Out[i,17] <- mr_results3[1,8]
          Out[i,18] <- mr_results3[1,9]
          Out[i,19] <- mr_results4[1,5]
          Out[i,20] <- mr_results4[1,6]
          Out[i,21] <- mr_results4[1,7]
          Out[i,22] <- mr_results4[1,8]
          Out[i,23] <- mr_results4[1,9]
          Out[i,24] <- mr_results5[1,5]
          Out[i,25] <- mr_results5[1,6]
          Out[i,26] <- mr_results5[1,7]
          Out[i,27] <- mr_results5[1,8]
          Out[i,28] <- mr_results5[1,9]
          
          
        }
        
      }
      
    }
    
  }
  print(i)
}

Outord = Out[order(Out[,7]),]

write.csv(Outord,paste0(file="/filepath_output_directory/MR_",gwas_files[j],"_on_2940proteins.csv"),row.names=F,quote=F)

#####
#####
