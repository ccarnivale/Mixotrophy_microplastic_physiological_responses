##Primary Production code for going from umol to ugrams per L per hour
##none of this data is normalized yet However, you *could* normalize the initial
#start values given you have the first counts. I have already calculated the DIC
#for you and the transformations are included. However, how you want to visualize
# or report the data will be up to you.

mixo_phys_pp <- read.csv("mixo_mp_pp.csv")

str(mixo_phys_pp)

mixo_phys_pp$Time_days <- factor(mixo_phys_pp$Time_days)

mixo_phys_pp$Treatment_plastic <- factor(mixo_phys_pp$Treatment_plastic,levels = c("control", "10","1000", "100000"))

mixo_phy_pp_calc <- mixo_phys_pp %>% mutate(pp_DPM = pp_CPM/C_14_efficiency, bkgr_DPM = background_CPM/C_14_efficiency/0.1, ugram_Cfixed_per_L_hr = pp_DPM*DIC.umol.kg./bkgr_DPM*12*1.029/2)

#12 is the molar mass of carbon o convert from umol to ugrams
#1.029 is to convert from kg to L of Seawater
#2 is to get an hourly uptake since it was exposed for 2 hours