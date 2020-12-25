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
  library(rgdal)
  library(readstata13)
  library(INLA)
  library(tidyr)
}

# Load functions
{
  source("helpers.R")
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
    # MWBR7HFL.DTA    Malawi 2015-16
    
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
  
  # !!!!! Currently ignoring the survey design
  
  # Select dataset and set time range
  {
    country <- "Malawi"
    # MWBR7HFL.DTA    Malawi 2015-16
    # births <- read.dta13("../data/DHS/MWBR7HFL.DTA", generate.factors=TRUE)
    births <- readRDS("../data/births_Malawi.rds")
    
    time_range <- c(1986:2015)
    time_blocks <- c("1986-1990", "1991-1995", "1996-2000", "2001-2005",
                     "2006-2010", "2011-2015")
    
  }
  
  # Import and process data
  {
    vars_to_keep = c("caseid","v001","v002","v005","v008","v011","v012",
                     "v024","v025","bidx","b2","b3","b5","b7")
    
    births %<>% subset(select=vars_to_keep)
    
    births %<>% rename(
      "clustid" = v001,
      "hhid" = v002,
      "weights" = v005,
      "survey_date_cmc" = v008,
      "dob_cmc" = v011,
      "age" = v012,
      "region" = v024,
      "ur" = v025,
      "birth_col_num" = bidx,
      "child_birth_year" = b2,
      "child_dob_cmc" = b3,
      "child_alive" = b5,
      "child_age_at_death" = b7
      # Defaults omitted: "ultimate_area_unit" = v004, "psu" = v021, "sample_strata" = v022, "stratification_design" = v023, "region2" = v139
    )
    
    births %<>% mutate(
      neo_death = ifelse(child_age_at_death==0, 1, 0),
      time_5yr = case_when(
        child_birth_year %in% time_range[1:5] ~ time_blocks[1],
        child_birth_year %in% time_range[6:10] ~ time_blocks[2],
        child_birth_year %in% time_range[11:15] ~ time_blocks[3],
        child_birth_year %in% time_range[16:20] ~ time_blocks[4],
        child_birth_year %in% time_range[21:25] ~ time_blocks[5],
        child_birth_year %in% time_range[26:30] ~ time_blocks[6],
        child_birth_year %in% c((time_range[1]-100):(time_range[1]-1)) ~ "other"
      )
    )
    births$neo_death <- replace_na(births$neo_death, 0)
    births %<>% filter(time_5yr!="other")
    
    births %<>% mutate(
      # !!!!! Can get more precise than this by using woman's CMC DOB
      woman_age_ab = round(age-(survey_date_cmc-child_dob_cmc)/12)
    )
    age_bins <- c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49")
    births %<>% mutate(
      woman_age_ab_bin = case_when(
        woman_age_ab >= 10 & woman_age_ab <= 14 ~ age_bins[1],
        woman_age_ab >= 15 & woman_age_ab <= 19 ~ age_bins[2],
        woman_age_ab >= 20 & woman_age_ab <= 24 ~ age_bins[3],
        woman_age_ab >= 25 & woman_age_ab <= 29 ~ age_bins[4],
        woman_age_ab >= 30 & woman_age_ab <= 34 ~ age_bins[5],
        woman_age_ab >= 35 & woman_age_ab <= 39 ~ age_bins[6],
        woman_age_ab >= 40 & woman_age_ab <= 44 ~ age_bins[7],
        woman_age_ab >= 45 & woman_age_ab <= 49 ~ age_bins[8],
        woman_age_ab < 10 | woman_age_ab > 49 ~ "other"
      )
    )
    births %<>% filter(woman_age_ab_bin!="other")
    
    # Convert characters to factors
    births %<>% mutate(
      time_5yr = factor(time_5yr, levels=time_blocks),
      woman_age_ab_bin = factor(woman_age_ab_bin, levels=age_bins)
    )
    
  }
  
  # Calculate summary dataframes
  {
    summ_1 <- births %>% group_by(child_birth_year) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_2 <- births %>% group_by(time_5yr, woman_age_ab_bin) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_3 <- births %>% group_by(time_5yr) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_4 <- births %>% group_by(child_birth_year, woman_age_ab_bin) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
  }
  
  # NMR: yearly vs. smoothed
  {
    model_1 <- inla(
      deaths ~ 1 + f(child_birth_year, model="rw2"),
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor = list(compute=TRUE, link=1),
      data = summ_1
    )
    
    summ_1$nmr_smoothed <- model_1$summary.fitted.values[,1]
    
    # Export: 600x400
    ggplot(
      data.frame(
        x = rep(summ_1$child_birth_year, 2),
        y = c(summ_1$nmr_direct, summ_1$nmr_smoothed),
        type = rep(c("direct", "smoothed"), each=nrow(summ_1))
      ),
      aes(x=x, y=y, group=type, color=type)
    ) +
      geom_line() +
      labs(title=paste0("Neonatal mortality rates (",country,")"), x="Year",
           y="NMR", color="Type", group="Type")
    
  }
  
  # NMR, by woman's age at birth
  {
    # !!!!! TO DO; copy code from above
  }
  
  # New model
  {
    # Calculate crosstab to show selection bias issue and fill missing values
    crosstab <- xtabs(~time_5yr+woman_age_ab_bin, data=summ_2, addNA=T)
    
    # Augment births_summ_2 with "missing values" for INLA prediction
    summ_2_aug <- summ_2
    for (time in unique(summ_2$time_5yr)) {
      for (age in unique(summ_2$woman_age_ab_bin)) {
        if (crosstab[time,age]==0) {
          summ_2_aug <- rbind(
            as.data.frame(summ_2_aug),
            data.frame(time_5yr=time, woman_age_ab_bin=age, births=NA,
                       deaths=NA, nmr_direct=NA, stringsAsFactors=FALSE)
          )
        }
      }
    }
    
    # Run INLA model to get smoothed estimates over woman's age and time
    model_2 <- inla(
      deaths ~ 1 + f(time_5yr, model="rw2") + f(woman_age_ab_bin, model="rw2"),
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor=list(compute=TRUE, link=1),
      data = summ_2_aug
    )
    summ_2_aug$nmr_smoothed <- model_2$summary.fitted.values[,1]
    
    # Calculate "gamma" values (method 1)
    # !!!!! Note: "2006-2010" is survey/country-specific
    birth_probs1 <- summ_2_aug %>% filter(time_5yr=="2006-2010") %>%
      subset(select=c(woman_age_ab_bin,births)) %>%
      mutate(births=ifelse(is.na(births),0,births))
    birth_probs1 %<>% mutate(prob_m1=births/sum(birth_probs1$births, na.rm=T)) %>%
      subset(select=-births)
    
    # Calculate "gamma" values (method 2; INLA model + average over time)
    model_2 <- inla(
      births ~ 1 + f(time_5yr, model="rw2") + f(woman_age_ab_bin, model="rw2"),
      family = "poisson",
      control.family = list(link="log"),
      control.predictor=list(compute=TRUE, link=1),
      data = summ_2_aug
    )
    summ_2_aug$births_smoothed <- round(model_2$summary.fitted.values[,1])
    
    # Smoothed; constant over time
    birth_probs2 <- summ_2_aug %>% group_by(woman_age_ab_bin) %>%
      summarize(births=round(sum(births_smoothed)))
    birth_probs2 %<>% mutate(prob_m2=births/sum(birth_probs2$births)) %>%
      subset(select=-births)
    
    # Smoothed; changes over time
    yearly_total_births <- summ_2_aug %>% group_by(time_5yr) %>%
      summarize(births_total=sum(births_smoothed))
    summ_2_aug <- inner_join(summ_2_aug, yearly_total_births, by="time_5yr")
    summ_2_aug %<>% mutate(prob_m3=births_smoothed/births_total)
    summ_2_aug %<>% inner_join(birth_probs2, by="woman_age_ab_bin")
    summ_2_aug %<>% inner_join(birth_probs1, by="woman_age_ab_bin")

    # Export: 600x400
    ggplot(
      data.frame(
        x = rep(summ_2_aug$time_5yr, 3),
        y = c(summ_2_aug$prob_m1, summ_2_aug$prob_m2, summ_2_aug$prob_m3),
        age = rep(summ_2_aug$woman_age_ab_bin, 3),
        method = rep(c("Method 1", "Method 2", "Method 3"), each=nrow(summ_2_aug))
      ),
      aes(x=x, y=y, group=age, color=age)
    ) +
      geom_line() +
      facet_wrap(~method) +
      labs(title=paste0("Probability of birth over mother's age (",country,")"),
           x="Year", y="Probability", color="Type", group="Type") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
    # Calculate "corrected" NMR estimates
    summ_2_aug %<>% mutate(
      nmr_t_m1 = nmr_smoothed * prob_m1,
      nmr_t_m3 = nmr_smoothed * prob_m3
    )
    
    df_corrected <- summ_2_aug %>% group_by(time_5yr) %>% summarize(
      nmr_m1 = sum(nmr_t_m1),
      nmr_m3 = sum(nmr_t_m3)
    )
    
    summ_3$nmr_corrected_m1 <- df_corrected$nmr_m1
    summ_3$nmr_corrected_m3 <- df_corrected$nmr_m3
    
    # Add smoothed estimates (over time only) to births_summ_3
    model_3 <- inla(
      deaths ~ 1 + f(time_5yr, model="rw2"),
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor = list(compute=TRUE, link=1),
      data = summ_3
    )
    summ_3$nmr_smoothed <- model_3$summary.fitted.values[,1]
    
    # Export: 600x400
    ggplot(
      data.frame(
        x = rep(summ_3$time_5yr, 3),
        y = c(summ_3$nmr_direct, summ_3$nmr_smoothed, summ_3$nmr_corrected_m3),
        type = rep(c("direct", "smoothed", "corrected"), each=6)
      ),
      aes(x=x, y=y, group=type, color=type)
      ) +
      geom_line() +
      labs(title=paste0("Neonatal mortality rates (",country,")"), x="Year",
           y="NMR", color="Type", group="Type") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
  }
  
  # Annual proportion of births in each mother's age bin
  {
    summ_4 <- inner_join(
      summ_4,
      subset(summ_1, select=c(child_birth_year,births)),
      by = "child_birth_year"
    ) %>% rename(
      "births" = births.x,
      "births_total" = births.y,
    ) %>% mutate(prop_of_year_births = births/births_total)
    
    # Export: 600x400
    ggplot(
      summ_4,
      aes(x=child_birth_year, y=prop_of_year_births, group=woman_age_ab_bin,
          color=woman_age_ab_bin)
    ) +
      geom_line() +
      scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2", "green3", 
                                  "darkorchid2", "orange", "red", "black")) +
      labs(
        title=paste0("Annual proportion of births in each mother's age bin (",
          country,")"), x="Year", y="Proportion", color="Mother's age at birth",
        group="Mother's age at birth"
      )
      
  }

}



##################################################.
##### MAIN: Modified smooth direct estimates #####
##################################################.

if (FALSE) {
  
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
