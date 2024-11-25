
### MR of protein cis-variants on PA traits

### ALLELE 1 is effect allele in both proteins and PA traits

library(data.table)
library(R.utils)
library(TwoSampleMR)
library(ieugwasr)

protein_files<-list.files("filepath_to_proteomics_gwas_sumstats")
gwas_files<-list.files("filepath_to_traits_gwas_sumstats")

# j is a vector for number of traits
j=10

##  Out is a matrix for results output
Out = matrix(nrow=(2940), ncol=23, NA)
Out[,1] = protein_files
Out[,2] = sub("\\_.*", "", protein_files)
Out[,3] = sub("\\..*", "", gwas_files[j])
colnames(Out)= c("file_name","exposure","outcome","method1","nsnp1",
                 "b1","se1","pval1", "method2","nsnp2",
                 "b2","se2","pval2", "method3","nsnp3",
                 "b3","se3","pval3", "method4","nsnp4",
                 "b4","se4","pval4")


## loop for MR's

for (i in 1:length(protein_files)){
  
  protein<-fread(paste0("filepath_to_proteomics_gwas_sumstats/",protein_files[i]),header=T)
  protein<-protein[-which(protein$P>5e-8),]
  
  gwas<-fread(paste0(file="filepath_to_traits_gwas_sumstats/",gwas_files[j]),header=T)

  
  
  ## pick up cis-variants +/-500kb from transcription start site (TSS)
  ## use TSS file
  TSS<-fread(file="filepath/Gene_transcription-start_position.txt.gz",header=TRUE)
  tmp<-sub("\\_.*", "", protein_files[i])
  TSS$gene_name=TSS$`Gene name`
  TSS$bp_range<-(TSS$`Transcript end (bp)` - TSS$`Transcript start (bp)`)
  gene<-TSS[(TSS$gene_name%in%tmp[1]),]
  gene2<-gene[which.min(gene$`Transcript start (bp)`),]
  myGWAS<-protein[protein$CHR==gene2$`Chromosome/scaffold name`,]
  downstream=myGWAS[(myGWAS$BP_grch38>=gene2$`Transcript start (bp)`) & (myGWAS$BP_grch38<(gene2$`Transcript start (bp)`+1000000)),]
  upstream=myGWAS[(myGWAS$BP_grch38>(gene2$`Transcript start (bp)`-1000000)) & (myGWAS$BP_grch38<gene2$`Transcript start (bp)`),]
  protein2=rbind(downstream,upstream)
  
  if (nrow(protein2>0)){
    
    ## outcome
    out<-format_data(
      gwas,
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
      protein2,
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
    
    possible_error<-tryCatch(tmp<- ld_clump(
      dplyr::tibble(rsid=exp$SNP, pval=exp$pval.exposure, id=exp$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "filepath_for_local_Rpackage/TwoSampleMR_local_ld/EUR"
    ), error=function(e) e)
    
    if(inherits(possible_error, "error")) next
    
    
    if(nrow(tmp>0)){
      
      exp2<-exp[(exp$SNP%in%tmp$rsid),]
      if(nrow(exp2>0)) {
        
        dat <- harmonise_data(exp2,out,action = 1)
        
        if (nrow(dat)>0) {
          
          mr_results <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
          mr_results
          
          mr_results2 <- mr(dat, method_list = c("mr_weighted_median", "mr_wald_ratio"))
          mr_results2
          
          mr_results3 <- mr(dat, method_list = c("mr_weighted_mode", "mr_wald_ratio"))
          mr_results3
          
          mr_results4 <- mr(dat, method_list = c("mr_egger_regression", "mr_wald_ratio"))
          mr_results4
          
          
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
            
          }
          
        }
        
      }
      
    }
    
    
  }
  print(i)
  
}

Outord = Out[order(Out[,8]),]

write.csv(Outord,paste0(file="filepath_to_output_directory/MR_TSS1mb_2940proteins_on_",gwas_files[j],".csv"),row.names=F,quote=F)

#####
#####