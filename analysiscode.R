
rm(list = ls(all=TRUE))

library(devtools)
library(TwoSampleMR)
library(ggplot2)

ao <- available_outcomes()


##analysis of set of exposures on skintone and hair colour
#set up set of exposures  - selected 50 exposures filtered by;
# - NOT biobank
# - European/Mixed population
# - Disease or risk factor (NOT Immune system or Metabolites)
# - Males and Females only
# - extra filtering to remove duplicates of the same exposures 
#filtering done in .csv file of all outcome and then resulting exposure ids saved in csv file. 

ex_list<-read.csv("MR-negativecontrols/exposurelist.csv", sep = ",")
exposures <- as.character(ex_list[c(1:dim(ex_list)[1]),])


#MR analysis of all the exposures on tanning
exposure_dat <- extract_instruments(exposures)
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('UKB-b:533'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
#write.csv(res, file = "IVWresults.csv")

##obtain plot of all exposures with more than 5 snps
ex_dat_5 <- extract_instruments(c("1096","975","966","85","44","1007","26","91","996","1058","93","90","298","798","72","7","86","60","1025","833","1084","302","31","1001","22","2","300","970","301","299","12","89"))
ex_dat_5 <- clump_data(ex_dat_5)
out_dat_5 <- extract_outcome_data(ex_dat_5$SNP, c('UKB-b:533'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat_5 <- harmonise_data(ex_dat_5, out_dat_5, action = 2)
res_5 <- mr(dat_5, method_list=c("mr_ivw"))
res_5<-subset_on_method(res_5) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res_5<-sort_1_to_many(res_5,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
res_5<-split_exposure(res_5)

plot1 <- forest_plot_1_to_many(res_5,b="b",se="se",
                      exponentiate=F,ao_slc=F,lo=-0.1,up=0.1,
                      TraitM="exposure",col1_width=2.2,by=NULL,
                      xlab="Decrease in tanning per SD/Log OR increase in risk factor")
ggsave(plot1, file="forest_skintone.pdf", width=7, height=7)



##rerun to get full set of results including sensitivity analyses for significant exposures
ex_dat_sig <- extract_instruments(c("1037","1058","302","31","1001","93","86","299","91","60","90","1096", "2"))
ex_dat_sig <- clump_data(ex_dat_sig)
out_dat_sig <- extract_outcome_data(ex_dat_sig$SNP, c('UKB-b:533'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token=NULL)
dat_sig <- harmonise_data(ex_dat_sig, out_dat_sig, action = 2)
#mr_results <- mr(dat)
res_sigall <- mr(dat_sig)
write.csv(res_sigall, file = "sigresults_skintone.csv")



##analysis with hair colour as the outcome

outcome_dat_hair <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "HairGWAS.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_LINREG"
)

dat_hair <- harmonise_data(exposure_dat, outcome_dat_hair, action = 2)
res_hair <- mr(dat_hair, method_list=c("mr_ivw", "mr_wald_ratio"))
write.csv(res_hair, file = "IVWresults_hair.csv")


##obtain plot of all exposures with more than 5 snps

out_dat_5_hair <- read_outcome_data(
  snps = ex_dat_5$SNP,
  filename = "HairGWAS.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_LINREG"
)

dat_5 <- harmonise_data(ex_dat_5, out_dat_5_hair, action = 2)
res_5 <- mr(dat_5, method_list=c("mr_ivw"))
res_5<-subset_on_method(res_5) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res_5<-sort_1_to_many(res_5,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
res_5<-split_exposure(res_5)

plot2 <- forest_plot_1_to_many(res_5,b="b",se="se",
                      exponentiate=F,ao_slc=F,
                      TraitM="exposure",col1_width=2.2,by=NULL,
                      xlab="Increase in darkness of hair per SD/Log OR increase in risk factor")



##rerun to get full set of results with significant exposures
ex_dat_sig <- extract_instruments(c("1058","302","1001","85","300","91","301"))
ex_dat_sig <- clump_data(ex_dat_sig)

outcome_dat_hair_sig <- read_outcome_data(
  snps = ex_dat_sig$SNP,
  filename = "HairGWAS.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_LINREG"
)

dat_sig <- harmonise_data(ex_dat_sig, outcome_dat_hair_sig, action = 2)
res_sigall <- mr(dat_sig)
write.csv(res_sigall, file = "sigresults_hair.csv")

