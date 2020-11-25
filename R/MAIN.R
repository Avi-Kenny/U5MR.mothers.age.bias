# Title: "U5MR Time Bias"
# Author: Avi Kenny, Jon Wakefield



################.
##### SETUP #####
################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jon Wakefield/Project - U5MR Bias/z.u5mr.bias/R")
} else {
  setwd("z.u5mr.bias/R")
}

# Load packages
{
  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
  library(magrittr)
  library(dplyr)
  library(SUMMER)
  library(ggplot2)
  library(patchwork)
  # library(rdhs)
  # library(haven)
  library(rgdal)
  library(readstata13)
}

# Load functions
{
  source("smoothDirect_mod.R")
}

# Set code blocks to run
{
  run_main <- TRUE
}



###############################################.
##### MAIN: SUMMER smooth direct vignette #####
###############################################.

if (run_main) {
  
  # https://cran.r-project.org/web/packages/SUMMER/vignettes/u5mr-vignette.html
  # `dat` is a DF of birth data (one row per child-period)
  # `direct0` is a DF of Horvitz-Thompson direct estimates of U5MR
  
}



###################################################.
##### MAIN: Direct estimates, by mother's age #####
###################################################.

if (run_main) {
  
  load_data <- TRUE
  
  if (load_data) {
    
    births <- readRDS("../data/births.rds")
    geo <- readRDS("../data/geo.rds")
    Amat <- readRDS("../data/Amat.rds")
    dat <- readRDS("../data/dat.rds")
    direct0 <- readRDS("../data/direct0.rds")
    direct1 <- readRDS("../data/direct1.rds")
    direct2 <- readRDS("../data/direct2.rds")
    
  } else {
    
    # KEBR71FL.DTA    Kenya 2014
    # TDBR71FL.DTA    Chad 2014
    # GHBR72FL.DTA    Ghana 2014
    # RWBR70FL.DTA    Rwanda 2014
    # LSBR71FL.DTA    Lesotho 2014
    
    # Read in DHS births data
    filename <- "../data/DHS/KEBR71FL.DTA"
    births <- read.dta13(filename, generate.factors=TRUE)
    saveRDS(births, file="../data/births.rds")
    
    # Read in Kenya 2014 GIS data
    mapfilename <- "../data/shps/sdr_subnational_boundaries.shp"
    geo <- readOGR(mapfilename, verbose = FALSE)
    Amat <- getAmat(geo, geo$REGNAME)
    saveRDS(geo, file="../data/geo.rds")
    saveRDS(Amat, file="../data/Amat.rds")
    
    # Use getBirths() to extract birth data
    dat <- getBirths(
      data = births,
      strata = "v023",
      variables = c("caseid","v001","v002","v004","v005","v008","v011","v012",
                    "v021","v022","v023","v024","v025","v139","bidx","b3"),
      year.cut = seq(1978, 2013, by=5)
    )
    dat <- dat[, c("v001","v002","v008","v011","v012","v024","time","age",
                   "v005","strata","died","b3")]
    colnames(dat) <- c("clustid","id","survey_date_cmc","dob_cmc","woman_age","region","time","age","weights","strata",
                       "died","dob_child_cmc")
    dat %<>% mutate(
      # !!!!! Can get more precise than this by using woman's CMC DOB
      woman_age_ab = round(woman_age-(survey_date_cmc-dob_child_cmc)/12)
    )
    dat %<>% mutate(
      woman_age_bin = case_when(
        woman_age >= 15 & woman_age <= 19 ~ "15-19",
        woman_age >= 20 & woman_age <= 24 ~ "20-24",
        woman_age >= 25 & woman_age <= 29 ~ "25-29",
        woman_age >= 30 & woman_age <= 34 ~ "30-34",
        woman_age >= 35 & woman_age <= 39 ~ "35-39",
        woman_age >= 40 & woman_age <= 44 ~ "40-44",
        woman_age >= 45 & woman_age <= 49 ~ "45-49"
      ),
      woman_age_ab_bin = case_when(
        woman_age_ab >= 10 & woman_age_ab <= 14 ~ "10-14",
        woman_age_ab >= 15 & woman_age_ab <= 19 ~ "15-19",
        woman_age_ab >= 20 & woman_age_ab <= 24 ~ "20-24",
        woman_age_ab >= 25 & woman_age_ab <= 29 ~ "25-29",
        woman_age_ab >= 30 & woman_age_ab <= 34 ~ "30-34",
        woman_age_ab >= 35 & woman_age_ab <= 39 ~ "35-39",
        woman_age_ab >= 40 & woman_age_ab <= 44 ~ "40-44",
        woman_age_ab >= 45 & woman_age_ab <= 49 ~ "45-49",
        woman_age_ab < 10 | woman_age_ab > 49 ~ "other"
      )
    )
    
    # Remove outliers (age at birth equal to 6 or 9)
    dat %<>% filter(woman_age_ab_bin!="other")
    saveRDS(dat, file="../data/dat.rds")
    
    # Get direct estimates using getDirect()
    years <- levels(dat$time)
    direct0 <- getDirect(
      births=dat, years=years, regionVar="region", timeVar="time",
      clusterVar="~clustid + id", ageVar="age", weightsVar="weights",
      national.only = TRUE
    )
    direct1 <- getDirect(
      births=dat, years=years, regionVar="woman_age_bin", timeVar="time",
      clusterVar="~clustid + id", ageVar="age", weightsVar="weights"
    )
    direct2 <- getDirect(
      births=dat, years=years, regionVar="woman_age_ab_bin", timeVar="time",
      clusterVar="~clustid + id", ageVar="age", weightsVar="weights"
    )
    direct1 %<>% rename("woman_age_bin" = region)
    direct2 %<>% rename("woman_age_ab_bin" = region)
    saveRDS(direct0, file="../data/direct0.rds")
    saveRDS(direct1, file="../data/direct1.rds")
    saveRDS(direct2, file="../data/direct2.rds")
    
  }
  
  # !!!!! Skipped HIV adjustment for now; see vignette
  
  # U5MR national
  ggplot(
    data = direct0 %>%
      mutate(years = factor(years, levels=levels(dat$time))),
    aes(x=years, y=mean, group=1)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    labs(x="Years", y="U5MR")
  
  # U5MR national, by woman age group (at time of survey)
  # Filtered out one age group with huge CI
  direct1_mod <- direct1 %>% mutate(subset="Age-specific")
  direct1_all <- direct1_mod %>% filter(woman_age_bin=="All") %>%
    mutate(subset="All")
  direct1_mod %<>% filter(woman_age_bin!="All")
  for (bin in unique(direct1$woman_age_bin)) {
    direct1_all2 <- direct1_all %>% mutate(woman_age_bin=bin)
    direct1_mod <- rbind(direct1_mod, direct1_all2)
  }
  ggplot(
    data = direct1_mod %>%
      mutate(years = factor(years, levels=levels(dat$time))) %>%
      filter(!(woman_age_bin=="25-29" & years=="93-97")),
    aes(x=years, y=mean, group=subset, color=subset)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=subset), alpha=0.2) +
    facet_wrap(~woman_age_bin, ncol=4) +
    labs(x="Years", y="U5MR", fill="Subset",
         color="Subset", title="U5MR, by woman's age (at time of survey)")
  
  # U5MR national, by woman age group (at time of birth)
  direct2_mod <- direct2 %>% mutate(subset="Age-specific")
  direct2_all <- direct2_mod %>% filter(woman_age_ab_bin=="All") %>%
    mutate(subset="All")
  direct2_mod %<>% filter(woman_age_ab_bin!="All")
  for (bin in unique(direct2$woman_age_ab_bin)) {
    direct2_all2 <- direct2_all %>% mutate(woman_age_ab_bin=bin)
    direct2_mod <- rbind(direct2_mod, direct2_all2)
  }
  ggplot(
    data = direct2_mod %>%
      mutate(years = factor(years, levels=levels(dat$time))) %>%
      filter(years %in% c("98-02","03-07","08-12")) %>%
      filter(woman_age_ab_bin!="All"),
    aes(x=years, y=mean, group=subset, color=subset)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=subset), alpha=0.2) +
    facet_wrap(~woman_age_ab_bin, ncol=4) +
    labs(x="Years", y="U5MR", fill="Subset",
         color="Subset", title="U5MR, by woman's age (at birth)")
  
  # What percentage of deaths are in each age band (by time period)
  dat_deaths <- dat %>% filter(died==1)
  xtabs(~woman_age_ab_bin+time, data=dat_deaths, addNA=TRUE) %>% prop.table(2)
  #                  time
  # woman_age_ab_bin 78-82 83-87 88-92 93-97 98-02 03-07 08-12
  #            10-14  0.58  0.13  0.09  0.05  0.05  0.03  0.01
  #            15-19  0.42  0.71  0.38  0.28  0.23  0.21  0.19
  #            20-24  0.00  0.16  0.46  0.40  0.32  0.27  0.27
  #            25-29  0.00  0.00  0.08  0.24  0.25  0.21  0.21
  #            30-34  0.00  0.00  0.00  0.04  0.14  0.18  0.17
  #            35-39  0.00  0.00  0.00  0.00  0.01  0.08  0.09
  #            40-44  0.00  0.00  0.00  0.00  0.00  0.01  0.04
  #            45-49  0.00  0.00  0.00  0.00  0.00  0.00  0.00
  
  # !!!!! Also do % of deaths

}



##################################################.
##### MAIN: Modified smooth direct estimates #####
##################################################.

if (run_main) {
  
  # RW2 smoothed direct model (national)
  fit0 <- smoothDirect(
    data=direct0, Amat=NULL, year_label=years, year_range=c(1978,2012),
    time.model="rw2", is.yearly=FALSE
    # time.model="rw2", is.yearly=TRUE, m=5
  )
  out0 <- getSmoothed(fit0)
  
  # Combine smoooth and direct estimates into one dataframe
  c1 <- out0 %>% subset(select=c(years,upper,lower,median)) %>%
    mutate(type="Smoothed") %>% rename('estimate'=median)
  c2 <- direct0 %>% subset(select=c(years,upper,lower,mean)) %>%
    mutate(type="Direct") %>% rename('estimate'=mean)
  combine <- rbind(c1,c2)
  
  # Plot estimates
  ggplot(
    data = combine %>% mutate(years = factor(years, levels=levels(out0$years))),
    aes(x=years, y=estimate, group=type, color=type, fill=type)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    labs(x="Years", y="U5MR", title="U5MR over time", color="Type", fill="Type")

  # !!!!! Continue: implement smoothDirect_mod()
  
}
