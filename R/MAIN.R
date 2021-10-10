# Title: "U5MR Mother's Age Bias"
# Author: Avi Kenny, Jon Wakefield



################.
##### SETUP #####
################.

# Load packages
{
  pkgs <- c("magrittr", "dplyr", "SUMMER", "ggplot2", "patchwork", "rgdal",
            "readstata13", "INLA", "tidyr")
  for (pkg in pkgs) { do.call("library", list(pkg)) }
}

# Helper functions
{
  #' Convert dates (year+month) to CMC
  #' @param year Calendar year
  #' @param month Calendar month
  #' @return Date in CMC format
  dates_to_cmc <- function(year, month) {
    12*(year-1900)+month
  }
  
  #' Convert CMC to dates (year+month)
  #' @param cmc CMC-formatted date
  #' @return Date in CMC format
  #' @return A list containing:
  #'     * `year`: calendar year
  #'     * `month`: calendar month
  cmc_to_dates <- function(cmc) {
    month <- ifelse(mod(cmc, 12)!=0, mod(cmc, 12), 12)
    return(list(
      year = 1900 + (cmc-month)/12,
      month = month
    ))
  }
}



###################################################.
##### MAIN: Direct estimates, by mother's age #####
###################################################.

if (run_main) {
  
  # On first run, set this to FALSE. Later, set this to TRUE to speed things up
  #     by reloading files so that they don't have to be recreated each time.
  load_data <- FALSE
  
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
    filename <- "../data/DHS/MWBR7HFL.DTA"
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
      year.cut = seq(1980, 2015, by=5)
    )
    dat <- dat[, c("v001","v002","v008","v011","v012","v024","time","age",
                   "v005","strata","died","b3")]
    colnames(dat) <- c("clustid","id","survey_date_cmc","dob_cmc","woman_age",
                       "region","time","age","weights","strata","died",
                       "dob_child_cmc")
    dat %<>% mutate(
      # Note: we can get more precise than this by using woman's CMC DOB
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
    
    # Get direct estimates of U5MR
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
  
  # 1. U5MR five-year direct estimates, Malawi
  # Export: 4" x 6"
  ggplot(
    data = direct0 %>%
      mutate(years = factor(years, levels=levels(dat$time))),
    aes(x=years, y=mean, group=1)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    labs(x="Years", y="U5MR",
         title="1. U5MR five-year direct estimates, Malawi")
  
  # 2. U5MR, by woman's age (at time of survey), Malawi
  # Filtered out one age group with huge CI
  # Export: 6" x 9"
  plotdata_2 <- direct1 %>% mutate(subset="Age-specific")
  plotdata_2b <- plotdata_2 %>% filter(woman_age_bin=="All") %>%
    mutate(subset="All")
  plotdata_2 %<>% filter(woman_age_bin!="All")
  for (bin in unique(direct1$woman_age_bin)) {
    plotdata_2c <- plotdata_2b %>% mutate(woman_age_bin=bin)
    plotdata_2 <- rbind(plotdata_2, plotdata_2c)
  }
  plotdata_2 %<>% mutate(years = factor(years, levels=levels(dat$time))) %>%
    filter(!(woman_age_bin=="25-29" & years=="95-99"))
  ggplot(
    data = plotdata_2,
    aes(x=years, y=mean, group=subset, color=subset)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=subset), alpha=0.2) +
    facet_wrap(~woman_age_bin, ncol=4) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    labs(x="Years", y="U5MR", fill="Subset", color="Subset",
         title="U5MR, by woman's age (at time of survey), Malawi")
  
  # 3. U5MR, by woman's age (at birth), Malawi
  # Export: 6" x 9"
  plotdata_3 <- direct2 %>% mutate(subset="Age-specific")
  plotdata_3b <- plotdata_3 %>% filter(woman_age_ab_bin=="All") %>%
    mutate(subset="All")
  plotdata_3 %<>% filter(woman_age_ab_bin!="All")
  for (bin in unique(direct2$woman_age_ab_bin)) {
    plotdata_3c <- plotdata_3b %>% mutate(woman_age_ab_bin=bin)
    plotdata_3 <- rbind(plotdata_3, plotdata_3c)
  }
  ggplot(
    data = plotdata_3 %>%
      mutate(years = factor(years, levels=levels(dat$time))) %>%
      filter(years %in% c("00-04","05-09","10-14")) %>%
      filter(woman_age_ab_bin!="All"),
    aes(x=years, y=mean, group=subset, color=subset)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=subset), alpha=0.2) +
    facet_wrap(~woman_age_ab_bin, ncol=4) +
    labs(x="Years", y="U5MR", fill="Subset",
         color="Subset", title="U5MR, by woman's age (at birth), Malawi")
  
  # What percentage of deaths are in each age band (by time period)
  # !!!!! Also do % of deaths
  xtabs(~woman_age_ab_bin+time, data=filter(dat, died==1), addNA=TRUE) %>%
    prop.table(2)

}



##################################################.
##### MAIN: Modified smooth direct estimates #####
##################################################.

if (run_main) {
  
  # !!!!! Currently ignoring the survey design
  
  # Select dataset and set time range
  {
    country <- "Malawi"
    # # MWBR7HFL.DTA    Malawi 2015-16
    # # births <- read.dta13("../data/DHS/MWBR7HFL.DTA", generate.factors=TRUE)
    # births <- readRDS("../data/births_Malawi.rds")
    # 
    births2 <- births
    time_range <- c(1986:2015)
    time_blocks <- c("1986-1990", "1991-1995", "1996-2000", "2001-2005",
                     "2006-2010", "2011-2015")
    
  }
  
  # Import and process data
  {
    vars_to_keep = c("caseid","v001","v002","v005","v008","v011","v012",
                     "v024","v025","bidx","b2","b3","b5","b7")
    
    births2 %<>% subset(select=vars_to_keep)
    
    births2 %<>% rename(
      "clustid" = v001,
      "hhid" = v002,
      "weights" = v005,
      "survey_date_cmc" = v008,
      "dob_cmc" = v011,
      "age" = v012,
      "region" = v024,
      "rural" = v025,
      "birth_col_num" = bidx,
      "child_birth_year" = b2,
      "child_dob_cmc" = b3,
      "child_alive" = b5,
      "child_age_at_death" = b7
      # Defaults omitted: "ultimate_area_unit" = v004, "psu" = v021, "sample_strata" = v022, "stratification_design" = v023, "region2" = v139
    )
    
    births2 %<>% mutate(
      rural = as.numeric(rural=="rural"),
      neo_death = as.numeric(child_age_at_death==0, 1, 0),
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
    births2$neo_death <- replace_na(births2$neo_death, 0)
    births2 %<>% filter(time_5yr!="other")
    
    births2 %<>% mutate(
      # !!!!! Can get more precise than this by using woman's CMC DOB
      woman_age_ab = round(age-(survey_date_cmc-child_dob_cmc)/12)
    )
    age_bins <- c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49")
    births2 %<>% mutate(
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
    
    # Check: should be 0
    print(nrow(births2 %>% filter(woman_age_ab_bin=="other")))
    
    # Convert characters to factors
    births2 %<>% mutate(
      time_5yr = factor(time_5yr, levels=time_blocks),
      woman_age_ab_bin = factor(woman_age_ab_bin, levels=age_bins)
    )
    
  }
  
  # Calculate summary dataframes
  {
    summ_1 <- births2 %>% group_by(child_birth_year) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_2 <- births2 %>% group_by(time_5yr, woman_age_ab_bin, rural) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_3 <- births2 %>% group_by(time_5yr, rural) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
    
    summ_4 <- births2 %>%
      group_by(child_birth_year, woman_age_ab_bin, rural) %>%
      summarize(births=n(), deaths=sum(neo_death)) %>%
      mutate(nmr_direct=deaths/births)
  }
  
  # NMR: yearly vs. smoothed
  {
    # Created smoothed direct annual estimates
    model_1 <- inla(
      deaths ~ 1 + f(child_birth_year, model="rw2"),
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor = list(compute=TRUE, link=1),
      data = summ_1
    )
    
    summ_1$nmr_smoothed <- model_1$summary.fitted.values[,1]
    
    # 4. Neonatal mortality rates, Malawi
    # Export: 4" x 6"
    ggplot(
      data.frame(
        x = rep(summ_1$child_birth_year, 2),
        y = c(summ_1$nmr_direct, summ_1$nmr_smoothed),
        type = rep(c("direct", "smoothed"), each=nrow(summ_1))
      ),
      aes(x=x, y=y, group=type, color=type)
    ) +
      geom_line() +
      labs(title=paste0("Neonatal mortality rates, ",country), x="Year",
           y="NMR", color="Type", group="Type")
    
  }
  
  # NMR, by woman's age at birth
  {
    # !!!!! TO DO; copy code from above
  }
  
  # New model
  {
    # Calculate crosstab to fill missing values
    crosstab <- xtabs(~time_5yr+woman_age_ab_bin, data=summ_2, addNA=T)
    
    # Augment births_summ_2 with "missing values" for INLA prediction
    # Assumes crosstab is the same for urban/rural
    summ_2_aug <- summ_2
    for (time in unique(summ_2$time_5yr)) {
      for (age in unique(summ_2$woman_age_ab_bin)) {
        if (crosstab[time,age]==0) {
          summ_2_aug <- rbind(
            as.data.frame(summ_2_aug),
            data.frame(time_5yr=time, woman_age_ab_bin=age, births=NA,
                       deaths=NA, rural=0, nmr_direct=NA,
                       stringsAsFactors=FALSE),
            data.frame(time_5yr=time, woman_age_ab_bin=age, births=NA,
                       deaths=NA, rural=1, nmr_direct=NA,
                       stringsAsFactors=FALSE)
          )
        }
      }
    }
    
    # Run INLA model to get smoothed estimates over woman's age and time
    model_2 <- inla(
      deaths ~ 1 + f(time_5yr, model="rw2") + f(woman_age_ab_bin, model="rw2") +
                   rural,
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor = list(compute=TRUE, link=1),
      data = summ_2_aug
    )
    summ_2_aug$nmr_smoothed <- model_2$summary.fitted.values[,1]
    
    # !!!!! Plot smoothed vs actual
    
    # Calculate "gamma" values (method 1)
    # !!!!! Note: "2006-2010" is survey/country-specific
    birth_probs1 <- summ_2_aug %>% filter(time_5yr=="2006-2010") %>%
      subset(select=c(woman_age_ab_bin,rural,births)) %>%
      mutate(births=ifelse(is.na(births),0,births))
    births_rural <- sum(birth_probs1$births[birth_probs1$rural==1])
    births_urban <- sum(birth_probs1$births[birth_probs1$rural==0])
    birth_probs1 %<>% mutate(
      gamma_m1=births/(ifelse(rural==1,births_rural,births_urban))
    ) %>%
      subset(select=-births)
    
    # A few transformations of summ_2_aug
    summ_2_aug$id <- c(1:nrow(summ_2_aug))
    summ_2_aug_urban <- filter(summ_2_aug, rural==0)
    summ_2_aug_rural <- filter(summ_2_aug, rural==1)
    
    # Calculate "gamma" values (method 2; INLA model)
    {
      model_2bi <- inla(
        births ~ 1 + f(time_5yr, model="rw1") + f(woman_age_ab_bin, model="rw1") +
                     rural,
        family = "poisson",
        control.family = list(link="log"),
        control.predictor = list(compute=TRUE, link=1),
        data = summ_2_aug
      )
      summ_2_aug$births_smoothed_i <- round(model_2bi$summary.fitted.values[,1])
      
      model_2bii_urban <- inla(
        births ~ 1 + f(time_5yr, model="rw1") + f(woman_age_ab_bin, model="rw1"),
        family = "poisson",
        control.family = list(link="log"),
        control.predictor = list(compute=TRUE, link=1),
        data = summ_2_aug_urban
      )
      model_2bii_rural <- inla(
        births ~ 1 + f(time_5yr, model="rw1") + f(woman_age_ab_bin, model="rw1"),
        family = "poisson",
        control.family = list(link="log"),
        control.predictor = list(compute=TRUE, link=1),
        data = summ_2_aug_rural
      )
      df_temp <- data.frame(
        id = c(summ_2_aug_urban$id, summ_2_aug_rural$id),
        v = c(round(model_2bii_urban$summary.fitted.values[,1]),
              round(model_2bii_rural$summary.fitted.values[,1]))
      )
      summ_2_aug$births_smoothed_ii <- (arrange(df_temp, id))$v
      
    }
    
    # Smoothed; changes over time
    if (FALSE) { # !!!!! temporarily skip this
      yearly_total_births <- summ_2_aug %>% group_by(time_5yr, rural) %>%
        summarize(births_total=sum(births_smoothed))
      summ_2_aug <- inner_join(summ_2_aug, yearly_total_births,
                               by=c("time_5yr", "rural"))
      summ_2_aug %<>% mutate(gamma_m2=births_smoothed/births_total)
      summ_2_aug %<>% inner_join(birth_probs1, by=c("woman_age_ab_bin","rural"))
    }
    summ_2_aug %<>% mutate(rural_ch=ifelse(rural==1, "rural", "urban"))
    
    # !!!!! START TEMP PLOT
    ggplot(
      data.frame(
        x = rep(summ_2_aug$time_5yr, 3),
        y = c(summ_2_aug$births, summ_2_aug$births_smoothed_i, summ_2_aug$births_smoothed_ii),
        age = rep(summ_2_aug$woman_age_ab_bin, 3),
        rural_ch = rep(summ_2_aug$rural_ch, 3),
        method = rep(c("Observed", "Smoothed (rw1)", "Smoothed (rw1; separate urban/rural)"),
                     each=nrow(summ_2_aug))
      ),
      aes(x=x, y=y, group=age, color=age)
    ) +
      geom_line() +
      facet_grid(rows=vars(rural_ch), cols=vars(method), scales="free_y") +
      labs(title=paste0("Number of births, by mother's age (smoothed vs. obser",
                        "ved), ", country),
           x="Year", y="Number of births", color="Group", group="Age bin") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      scale_color_manual(values=c("turquoise","salmon","orange","darkblue",
                                  "green3","darkorchid2","black","red"))
    # !!!!! END TEMP PLOT
    
    
    # 5b. Number of births, by mother's age (smoothed vs. observed), Malawi
    ggplot(
      data.frame(
        x = rep(summ_2_aug$time_5yr, 2),
        y = c(summ_2_aug$births, summ_2_aug$births_smoothed),
        age = rep(summ_2_aug$woman_age_ab_bin, 2),
        rural_ch = rep(summ_2_aug$rural_ch, 2),
        method = rep(c("Observed", "Smoothed"),
                     each=nrow(summ_2_aug))
      ),
      aes(x=x, y=y, group=age, color=age)
    ) +
      geom_line() +
      facet_grid(rows=vars(rural_ch), cols=vars(method), scales="free_y") +
      labs(title=paste0("Number of births, by mother's age (smoothed vs. obser",
                        "ved), ", country),
           x="Year", y="Number of births", color="Group", group="Age bin") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      scale_color_manual(values=c("turquoise","salmon","orange","darkblue",
                                  "green3","darkorchid2","black","red"))
    
    # 5. Probability of birth, by mother's age (gamma), Malawi
    # Export: 6" x 9"
    ggplot(
      data.frame(
        x = rep(summ_2_aug$time_5yr, 2),
        y = c(summ_2_aug$gamma_m1, summ_2_aug$gamma_m2),
        age = rep(summ_2_aug$woman_age_ab_bin, 2),
        rural_ch = rep(summ_2_aug$rural_ch, 2),
        method = rep(c("Method 1", "Method 2"),
                     each=nrow(summ_2_aug))
      ),
      aes(x=x, y=y, group=age, color=age)
    ) +
      geom_line() +
      facet_grid(rows=vars(method), cols=vars(rural_ch)) +
      labs(title=paste0("Probability of birth, by mother's age (gamma), ",
           country), x="Year", y="Probability", color="Group", group="Age bin") +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      scale_color_manual(values=c("turquoise","salmon","orange","darkblue",
                                  "green3","darkorchid2","black","red"))
    
    # Calculate "corrected" NMR estimates
    summ_2_aug %<>% mutate(
      nmr_t_m1 = nmr_smoothed * gamma_m1,
      nmr_t_m2 = nmr_smoothed * gamma_m2
    )
    
    df_corrected <- summ_2_aug %>% group_by(time_5yr) %>% summarize(
      nmr_m1 = sum(nmr_t_m1),
      nmr_m2 = sum(nmr_t_m2)
    )
    
    summ_3$nmr_corrected_m1 <- df_corrected$nmr_m1
    summ_3$nmr_corrected_m2 <- df_corrected$nmr_m2
    
    # Add smoothed estimates (over time only) to births_summ_3
    model_3 <- inla(
      deaths ~ 1 + f(time_5yr, model="rw2") + rural,
      family = "binomial",
      Ntrials = births,
      control.family = list(link="logit"),
      control.predictor = list(compute=TRUE, link=1),
      data = summ_3
    )
    summ_3$nmr_smoothed <- model_3$summary.fitted.values[,1]
    
    # 6. Neonatal mortality rates (three methods), Malawi
    # Export: 4" x 6"
    ggplot(
      data.frame(
        x = rep(summ_3$time_5yr, 3),
        y = c(summ_3$nmr_direct, summ_3$nmr_smoothed, summ_3$nmr_corrected_m2),
        type = rep(c("direct", "smoothed", "corrected"), each=6)
      ),
      aes(x=x, y=y, group=type, color=type)
      ) +
      geom_line() +
      labs(title=paste0("Neonatal mortality rates (three methods), ",country),
           x="Year", y="NMR", color="Type", group="Type") +
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
    
    # 7. Annual proportion of births in each mother's age bin, Malawi
    # Export: 4" x 6"
    ggplot(
      summ_4,
      aes(x=child_birth_year, y=prop_of_year_births, group=woman_age_ab_bin,
          color=woman_age_ab_bin)
    ) +
      geom_line() +
      scale_color_manual(
        values=c("turquoise", "salmon", "dodgerblue2", "green3", "darkorchid2",
                 "orange", "red", "black")) +
      labs(
        title=paste0("Annual proportion of births in each mother's age bin, ",
          country), x="Year", y="Proportion", color="Mother's age at birth",
        group="Mother's age at birth"
      )
    
  }

}



##################################################.
##### MAIN: Modified smooth direct estimates #####
##################################################.

if (FALSE) {
  
  # ?????
  
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
