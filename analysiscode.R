

library(devtools)
library(TwoSampleMR)


ao <- available_outcomes()

##MR analysis of education on bmi
exposure_dat <- extract_instruments(c('1001'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('2'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token=NULL)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results_edubmiheight <- mr(dat)

##analysis of set of exposures on skintone and hair colour
#set up set of exposures  - selected 50 exposures filtered by;
# - NOT biobank
# - European/Mixed population
# - Disease or risk factor (NOT Immune system or Metabolites)
# - Males and Females only
# - extra filtering to remove duplicates of the same exposures 
#filtering done in .csv file of all outcome and then resulting exposure ids saved in csv file. 

ex_list<-read.csv("outcomes_formrbase.csv", sep = ",")
exposures <- as.character(ex_list[c(1:dim(ex_list)[1]),])


#MR analysis of all the exposures on tanning
exposure_dat <- extract_instruments(exposures)
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('UKB-b:533'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)
write.csv(mr_results, file = "allresults.csv")

res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
write.csv(res, file = "IVWresults.csv")



##rerun to get plot of significant exposures
ex_dat_sig <- extract_instruments(c("1037","1058","302","31","1001","93","86","299","91","60","90","1096", "2"))
ex_dat_sig <- clump_data(ex_dat_sig)
out_dat_sig <- extract_outcome_data(ex_dat_sig$SNP, c('UKB-b:533'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token=NULL)
dat_sig <- harmonise_data(ex_dat_sig, out_dat_sig, action = 2)
#mr_results <- mr(dat)
res_sigall <- mr(dat_sig)
write.csv(res_sigall, file = "MR base analysis/sigresults_allmethods.csv")
res_sig <- mr(dat_sig, method_list=c("mr_ivw", "mr_wald_ratio"))
res_sig<-subset_on_method(res_sig) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res_sig<-sort_1_to_many(res_sig,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
res_sig<-split_exposure(res_sig)

forest_plot_1_to_many(res_sig,b="b",se="se",
                      exponentiate=F,ao_slc=F,lo=-0.2,up=0.25,
                      TraitM="exposure",col1_width=2.2,by=NULL,
                      xlab="Decrease in tanning per SD increase in risk factor")
ggsave(p1[[1]], file="sigexposures_forest.pdf", width=7, height=7)

#plots of each significant exposure 
plot_allsig <- mr_scatter_plot(res_sigall, dat_sig)



