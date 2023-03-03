# load libraries
rm(list=ls())
library(tidyverse)
library(rms)

# output directories
output_dir <- "C:/Users/jmd237/OneDrive - University of Exeter/John/Projects/2023_tsmodel/results/" 
data_dir <- "C:/Users/jmd237/OneDrive - University of Exeter/John/CPRD/mastermind22/"

# load data
load("C:/Users/jmd237/OneDrive - University of Exeter/John/CPRD/mastermind22/20230213_t2d_1stinstance.Rda")
load("C:/Users/jmd237/OneDrive - University of Exeter/John/CPRD/mastermind22/20230213_t2d_all_drug_periods.Rda")

#### Set up data function #####
set_up_data <- function(dataset.type, drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"), earliestdrugstartdate) {
  ##### Input variables
  # dataset.type: a character string mentioning the type of dataset required
  
  # initial checks
  if (missing(dataset.type)) {stop("'dataset.type' needs to be supplied")}
  if (!is.character(dataset.type)) {stop("'dataset.type' must be a character string")}
  if (!(dataset.type %in% c("diagnostics", "synthetic", "full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset", "ckd.dataset" , "cvd.dataset", "hf.dataset", "no_co.dataset", "semaglutide.dataset"))) {
    stop("'dataset.type' must be one of: diagnostics / synthetic / full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset / ckd.dataset / cvd.dataset / hf.dataset / no_co.dataset / semaglutide.dataset")
  }
  if (missing(drugs)) {stop("'drugs' needs to be supplied")}
  if (!is.character(drugs)) {stop("'drugs' must be a character string")}
  for (i in 1:length(drugs)) {
    if (!(drugs[i] %in% c("DPP4", "GLP1", "INS", "MFN", "SGLT2", "SU", "TZD"))) {
      stop("'drugs' must be one of: DPP4 / GLP1 / INS / MFN / SGLT2 / SU / TZD")
    }
  }
  
  if (dataset.type == "ckd.dataset") {stop("outcome variables for 'ckd.dataset' hasn't been coded")}
  
    
  # load original dataset # name - t2d_1stinstance
  cprd <- t2d_1stinstance
  
  
  ################################################
  ##### Select only input drugs
  ################################################
  
  cprd <- cprd %>% 
    filter(drugclass %in% drugs)
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Select only input drugs")
    print("################################################")
    print(nrow(cprd))
    print(table(cprd$drugclass))
    
  }
  
  
  ################################################
  ##### Drop patients initiating before 1/1/2013
  ################################################
  
  #######################
  # Explore adjusted HbA1c repsonse by calendar year
  
  cprd  <- cprd %>%
    mutate(yrdrugstart = format(dstartdate, format = "%Y")) %>%
    mutate(yrdrugstart = as.numeric(yrdrugstart))


  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < earliestdrugstartdate, 1 , NA_real_))

  # printing inclusion patients
  if (dataset.type == "diagnostics") {

    print("################################################")
    print("##### Drop patients initiating before 1/1/2013")
    print("################################################")
    print(table(cprd$dstartdate_cutoff))
    print(table(cprd$dstartdate_cutoff, cprd$drugclass))

  }

  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))
  
  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if treated with insulin when starting new drug")
    print("################################################")
    print(table(cprd$INS))
    print(table(cprd$INS, cprd$drugclass))
    
  }
  
  cprd <- cprd %>% 
    filter(INS == 0)      
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients with ESRD")
    print("################################################")
    print(table(cprd$preckdstage))
    print(table(cprd$preckdstage, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if first-line treatment")
    print("################################################")
    print(table(cprd$drugline_all))
    print(table(cprd$drugline_all, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(drugline_all != 1)
  
  
  ################################################
  ##### Drop if semaglutide
  ################################################
  
  cprd <- cprd %>%
    mutate(semaglutide_drug = ifelse(str_detect(drugsubstances, "Semaglutide"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if semaglutide")
    print("################################################")
    print(table(cprd$semaglutide_drug))
    
  }
  
  if (dataset.type == "semaglutide.dataset") {
    cprd <- cprd %>%
      filter(!is.na(semaglutide_drug))
  } else {
    cprd <- cprd %>%
      filter(is.na(semaglutide_drug))
  }
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  ### Add variable that identifies an individual entry in the data
  
  cprd <- cprd %>%
    mutate(pated = paste(patid, drugclass, dstartdate, sep = ".")) %>%
    
    ################################################
  ##### Drug of interest
  ################################################
  
  mutate(drugclass = factor(drugclass, levels = drugs)) 
  
  ################################################
  ##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
  #####   - posthba1c12m but if missing
  #####     - posthba1c6m
  ################################################
  
  cprd <- cprd %>%
    mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%
    mutate(posthba1cfinal = as.numeric(posthba1cfinal))
  
  
  ################################################
  ##### Sociodemographic variables
  ################################################
  
  cprd <- cprd %>%
    #####   - Age: agetx (new var)
    mutate(agetx = as.numeric(dstartdate_age)) %>%
    #####   - Sex: sex
    mutate(sex = factor(ifelse(gender == 1, "Male", "Female"))) %>%
    
    #####   - Duration of diabetes: t2dmduration
    mutate(t2dmduration = as.numeric(dstartdate_dm_dur_all)) %>%
    #####   - Ethnicity: ethnicity
    mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed"))) %>%
    
    #####   - Deprivation: deprivation
    mutate(deprivation = factor(imd2015_10)) %>%
    
    #####   - Smoking Status: smoke
    mutate(smoke = factor(smoking_cat)) %>%
    
    #####   - Line Therapy: drugline: turn all > 4 to 5+
    mutate(drugline = ifelse(drugline_all > 4, 5, drugline_all)) %>%
    mutate(drugline = factor(drugline, levels = c(2, 3, 4, 5), labels = c("2", "3", "4", "5+"))) %>%
    
    #####   - Hospitalisations in previous year
    mutate(prehospitalisation = factor(hosp_admission_prev_year, levels = c(0, 1), labels = c("No", "Yes")))
  
  ################################################
  ##### Diabetes treatment
  ################################################
  
  cprd <- cprd %>%
    #####   - Drugs taken alongside treatment
    # #####     - SU
    # #####     - MFN
    # #####     - DPP4
    # #####     - TZD
    # mutate(SU = factor(SU, levels = c(0, 1), labels = c("No", "Yes")),
    #        MFN = factor(MFN, levels = c(0, 1), labels = c("No", "Yes")),
    #        DPP4 = factor(DPP4, levels = c(0, 1), labels = c("No", "Yes")),
    #        TZD = factor(TZD, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    mutate(ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN) %>%
    mutate(ncurrtx = ifelse(ncurrtx > 4, 5, ncurrtx)) %>%
    mutate(ncurrtx = factor(ncurrtx, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+"))) %>%
    
    #####   - Outcome month: hba1cmonth
    mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
    mutate(hba1cmonth = as.numeric(hba1cmonth))
  
  
  ################################################
  ##### Biomarkers
  ################################################
  
  #####   - hba1c: prehba1c (Nothing to do)
  #####   - BMI: prebmi (Nothing to do)
  #####   - eGFR: preegfr (Nothing to do)
  #####   - Albumin:Creatine ratio: preacr (Nothing to do)
  #####   - Serum albumin: prealbumin_blood (remove the _ from the name)
  
  cprd <- cprd %>%
    rename("prealbuminblood" = "prealbumin_blood",
           "prealbuminblooddate" = "prealbumin_blooddate",
           "prealbuminblooddrugdiff" = "prealbumin_blooddrugdiff")
  
  #####   - Alanine aminotransferase: prealt (Nothing to do)
  #####   - Aspartate aminotransferase: preast (Nothing to do)
  #####   - Bilirubin: prebilirubin (Nothing to do)
  #####   - Fasting glucose: prefastingglucose (Nothing to do)
  #####   - Fasting haematocrit: prehaematocrit (Nothing to do)
  #####   - Fasting haemoglobin: prehaemoglobin (Nothing to do)
  #####   - High-density lipoprotein (HDL): prehdl (Nothing to do)
  
  #####   - Mean arterial BP: premap
  
  cprd <- cprd %>%
    mutate(premap = predbp + ((presbp - predbp) / 3))
  
  #####   - Total cholesterol: pretotalcholesterol (Nothing to do)
  #####   - Triglycerides: pretriglyceride (Nothing to do)
  
  
  
  ################################################
  ##### Comorbidities
  ################################################
  
  cprd <- cprd %>%
    #####   - Angina: predrug_earliest_angina
    mutate(preangina = factor(predrug_angina, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = factor(predrug_cld, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = factor(predrug_diabeticnephropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = factor(predrug_hypertension, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = factor(predrug_myocardialinfarction, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = factor(predrug_revasc, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = factor(predrug_stroke, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = factor(predrug_tia, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = factor(predrug_af, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - CKD stages
    mutate(preckd = factor(preckdstage)) %>%
    #####   - Pre-existing CVD
    mutate(predrug_cvd = ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1, 1, 0)) %>%
    mutate(predrug_cvd = factor(predrug_cvd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Date of CV / HF death for later outcomes
    mutate(cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA))) %>%
    #####   - Define 5 years post-drug start date for later censoring
    mutate(five_years_post_dstart=dstartdate+(365.25*5))
  
  
  
  ################################################
  ##### Add in later GLP1/SGLT2/TZD drug starts needed for censoring
  ################################################
  
  
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort") {
    
    # Add in later GLP1/SGLT2/TZD drug starts needed for censoring
    later_sglt2 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="SGLT2") %>%
                    select(patid, next_sglt2=dstartdate)), by="patid") %>%
      filter(next_sglt2>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
      ungroup()
    
    
    later_glp1 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="GLP1") %>%
                    select(patid, next_glp1=dstartdate)), by="patid") %>%
      filter(next_glp1>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
      ungroup()
    
    
    later_tzd <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="TZD") %>%
                    select(patid, next_tzd=dstartdate)), by="patid") %>%
      filter(next_tzd>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_tzd_start=min(next_tzd, na.rm=TRUE)) %>%
      ungroup()
    
    
    cprd <- cprd %>%
      left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
      left_join(later_glp1, by=c("patid", "dstartdate")) %>%
      left_join(later_tzd, by=c("patid", "dstartdate"))
    
  }
  
  
  
  
  ###############################################################################
  ###############################################################################
  #################### Final dataset - all patients #############################
  ###############################################################################
  ###############################################################################
  #
  # Add all variables necessary for ALL analysis in the paper.
  #
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort") {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # CKD
        preckd, predrug_cvd,
        # CVD
        predrug_cvd, postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, 
        five_years_post_dstart, death_date, next_sglt2_start, next_tzd_start, gp_record_end, next_glp1_start,
        # HF
        postdrug_first_primary_hhf, hf_death_date_primary_cause,
        # No comorbidities
        postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, postdrug_first_heartfailure, 
        hf_death_date_any_cause, qrisk2_10yr_score
      ) %>%
      as.data.frame()
    
  } else {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # No comorbidities
        qrisk2_10yr_score
      ) %>%
      as.data.frame()
    
  }
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Final dataset - all patients")
    print("################################################")
    print(nrow(final.dataset))
    print(table(final.dataset$drugclass))
    
  }
  
  # if full cohort was requested
  if (dataset.type == "full.cohort" | dataset.type == "semaglutide.dataset") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Synthetic dataset ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # Create synthetic dataset
  if (dataset.type == "synthetic") {
    # load package
    require(synthpop)
    
    set.seed(123)
    syn.dataset <- synthpop::syn(final.dataset %>%
                                   select(
                                     # response hba1c
                                     posthba1cfinal,
                                     # therapies of interest
                                     drugclass,
                                     # Sociodemographic features
                                     agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
                                     # Diabetes treatment 
                                     drugline, ncurrtx, hba1cmonth, yrdrugstart,
                                     # Biomarkers
                                     prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
                                     prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
                                     # Comorbidities
                                     preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
                                     preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
                                   ) %>%
                                   sample_n(520),
                                 print.flag = FALSE
    )
    
    return(syn.dataset$syn)
    
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Propensity score model ###############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  ps.model.dataset <- final.dataset
  
  #:----------------------------------------------------
  # Select variables needed
  
  ps.model.dataset <- ps.model.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
      # Diabetes treatment 
      drugline, ncurrtx, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, 
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # Training dataset
  set.seed(123)
  ps.model.dataset.train <- ps.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6) %>%
    ungroup() %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score training cohort")
    print("################################################")
    print(nrow(ps.model.dataset.train))
    print(table(ps.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "ps.model.train") {
    return(ps.model.dataset.train)
  }
  
  
  # Testing dataset
  ps.model.dataset.test <- subset(ps.model.dataset, !(pated %in% ps.model.dataset.train$pated)) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score testing cohort")
    print("################################################")
    print(nrow(ps.model.dataset.test))
    print(table(ps.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "ps.model.test") {
    return(ps.model.dataset.test)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model")
    print("################################################")
    
  }
  
  hba1c.model.dataset.train <- ps.model.dataset.train %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  hba1c.model.dataset.test <- ps.model.dataset.test %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$multi_drug_start))
    print(table(hba1c.model.dataset.train$multi_drug_start, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$multi_drug_start))
    print(table(hba1c.model.dataset.test$multi_drug_start, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(multi_drug_start == 0)
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$timeprevcombo_less61))
    print(table(hba1c.model.dataset.train$timeprevcombo_less61, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$timeprevcombo_less61))
    print(table(hba1c.model.dataset.test$timeprevcombo_less61, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(timeprevcombo_less61))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(timeprevcombo_less61))
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$hb_extreme_53))
    print(table(hba1c.model.dataset.train$hb_extreme_53, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$hb_extreme_53))
    print(table(hba1c.model.dataset.test$hb_extreme_53, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(hb_extreme_53))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$prehba1c)))
    print(table(is.na(hba1c.model.dataset.train$prehba1c), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$prehba1c)))
    print(table(is.na(hba1c.model.dataset.test$prehba1c), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(prehba1c))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post HbA1c missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(posthba1cfinal))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(posthba1cfinal))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  # Training dataset
  final.hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Training cohort")
    print("################################################")
    print("Training Cohort")
    print(nrow(final.hba1c.model.dataset.train))
    print(table(final.hba1c.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "hba1c.train") {
    return(final.hba1c.model.dataset.train)
  }
  
  # Testing dataset
  final.hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Testing cohort")
    print("################################################")
    print("Testing Cohort")
    print(nrow(final.hba1c.model.dataset.test))
    print(table(final.hba1c.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "hba1c.test") {
    return(final.hba1c.model.dataset.test)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Weight population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model")
    print("################################################")
    
  }
  
  weight.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(weight.dataset$multi_drug_start))
    print(table(weight.dataset$multi_drug_start, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(weight.dataset$timeprevcombo_less61))
    print(table(weight.dataset$timeprevcombo_less61, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(weight.dataset$hb_extreme_53))
    print(table(weight.dataset$hb_extreme_53, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(weight.dataset$prehba1c)))
    print(table(is.na(weight.dataset$prehba1c), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$preweight)))
    print(table(is.na(weight.dataset$preweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$postweight)))
    print(table(is.na(weight.dataset$postweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(postweight))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.weight.dataset <- weight.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model - final")
    print("################################################")
    print(nrow(final.weight.dataset))
    print(table(final.weight.dataset$drugclass))
    
  }
  
  if (dataset.type == "weight.dataset") {
    return(final.weight.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Discontinuation population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model")
    print("################################################")
    
  }
  
  discontinuation.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(discontinuation.dataset$multi_drug_start))
    print(table(discontinuation.dataset$multi_drug_start, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(discontinuation.dataset$timeprevcombo_less61))
    print(table(discontinuation.dataset$timeprevcombo_less61, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(discontinuation.dataset$hb_extreme_53))
    print(table(discontinuation.dataset$hb_extreme_53, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(discontinuation.dataset$prehba1c)))
    print(table(is.na(discontinuation.dataset$prehba1c), discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.discontinuation.dataset <- discontinuation.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # discontinuation
      stopdrug_6m_3mFU
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model - final")
    print("################################################")
    print(nrow(final.discontinuation.dataset))
    print(table(final.discontinuation.dataset$drugclass))
    
  }
  
  if (dataset.type == "discontinuation.dataset") {
    return(final.discontinuation.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################## eGFR population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model")
    print("################################################")
    
  }
  
  egfr.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(egfr.dataset$multi_drug_start))
    print(table(egfr.dataset$multi_drug_start, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(egfr.dataset$timeprevcombo_less61))
    print(table(egfr.dataset$timeprevcombo_less61, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(egfr.dataset$hb_extreme_53))
    print(table(egfr.dataset$hb_extreme_53, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$prehba1c)))
    print(table(is.na(egfr.dataset$prehba1c), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if baseline eGRF is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if baseline eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$preegfr)))
    print(table(is.na(egfr.dataset$preegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(preegfr))
  
  
  ################################################
  ##### Drop if post eGRF is missing
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(postegfr = ifelse(is.na(postegfr12m), postegfr6m, postegfr12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$postegfr)))
    print(table(is.na(egfr.dataset$postegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(postegfr))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.egfr.dataset <- egfr.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # eGFR analysis
      postegfr
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model - final")
    print("################################################")
    print(nrow(final.egfr.dataset))
    print(table(final.egfr.dataset$drugclass))
    
  }
  
  if (dataset.type == "egfr.dataset") {
    return(final.egfr.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CKD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD outcome model")
    print("################################################")
    
  }
  
  ckd.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(ckd.dataset$multi_drug_start))
    print(table(ckd.dataset$multi_drug_start, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(ckd.dataset$timeprevcombo_less61))
    print(table(ckd.dataset$timeprevcombo_less61, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(ckd.dataset$hb_extreme_53))
    print(table(ckd.dataset$hb_extreme_53, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(ckd.dataset$prehba1c)))
    print(table(is.na(ckd.dataset$prehba1c), ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if CKD
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(no_ckd = ifelse(!is.na(preckd) & (preckd=="stage_3a" | preckd=="stage_3b" | preckd=="stage_4"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CKD")
    print("################################################")
    print(table(ckd.dataset$no_ckd))
    print(table(ckd.dataset$no_ckd, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(no_ckd))
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(ckd.dataset$TZD, ckd.dataset$drugclass))
    print("GLP1 treated")
    print(table(ckd.dataset$GLP1, ckd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(ckd.dataset$SGLT2, ckd.dataset$drugclass))
    
  }
  
  
  ckd.dataset <- ckd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  # ckd.dataset <- ckd.dataset %>%
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.ckd.dataset <- ckd.dataset %>%
    # select(
    #   # information regarding patient
    #   patid, pated,
    #   # response hba1c
    #   posthba1cfinal,
    #   # therapies of interest
    #   drugclass,
    #   # Sociodemographic features
    #   agetx, sex, t2dmduration, 
    #   # Diabetes treatment 
    #   drugline, ncurrtx, hba1cmonth,
  #   # Biomarkers
  #   prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
  #   prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
  #   # Comorbidities
  #   preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
  #   preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
  #   # eGFR analysis
  #   postegfr
  # ) %>%
  as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD model - final")
    print("################################################")
    print(nrow(final.ckd.dataset))
    print(table(final.ckd.dataset$drugclass))
    
  }
  
  if (dataset.type == "ckd.dataset") {
    return(final.ckd.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CVD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD outcome model")
    print("################################################")
    
  }
  
  cvd.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(cvd.dataset$multi_drug_start))
    print(table(cvd.dataset$multi_drug_start, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(cvd.dataset$timeprevcombo_less61))
    print(table(cvd.dataset$timeprevcombo_less61, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(cvd.dataset$hb_extreme_53))
    print(table(cvd.dataset$hb_extreme_53, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(cvd.dataset$prehba1c)))
    print(table(is.na(cvd.dataset$prehba1c), cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if CVD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CVD")
    print("################################################")
    print(table(cvd.dataset$predrug_cvd))
    print(table(cvd.dataset$predrug_cvd, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(predrug_cvd == "No")
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(cvd.dataset$TZD, cvd.dataset$drugclass))
    print("GLP1 treated")
    print(table(cvd.dataset$GLP1, cvd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(cvd.dataset$SGLT2, cvd.dataset$drugclass))
    
  }
  
  
  cvd.dataset <- cvd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  ## MACE: narrow definition (hospitalisation/death): narrow ('incident') MI/stroke HES codes (primary cause only) + CV death in ONS death (primary cause only)
  
  cvd.dataset <- cvd.dataset %>%
    
    mutate(postdrug_mace=pmin(postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.cvd.dataset <- cvd.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.cvd.dataset))
    print(table(final.cvd.dataset$drugclass))
    
  }
  
  if (dataset.type == "cvd.dataset") {
    return(final.cvd.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ##################### Heart Failure outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Heart Failure outcome model")
    print("################################################")
    
  }
  
  hf.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(hf.dataset$multi_drug_start))
    print(table(hf.dataset$multi_drug_start, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(hf.dataset$timeprevcombo_less61))
    print(table(hf.dataset$timeprevcombo_less61, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(hf.dataset$hb_extreme_53))
    print(table(hf.dataset$hb_extreme_53, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(hf.dataset$prehba1c)))
    print(table(is.na(hf.dataset$prehba1c), hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if heart failure
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if heart failure")
    print("################################################")
    print(table(hf.dataset$preheartfailure))
    print(table(hf.dataset$preheartfailure, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(preheartfailure == "No")
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(hf.dataset$TZD, hf.dataset$drugclass))
    print("GLP1 treated")
    print(table(hf.dataset$GLP1, hf.dataset$drugclass))
    print("SGLT2 treated")
    print(table(hf.dataset$SGLT2, hf.dataset$drugclass))
    
  }
  
  
  hf.dataset <- hf.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  hf.dataset <- hf.dataset %>%
    
    ## HF: narrow definition (hospitalisation/death): HF HES codes (primary cause only) + HF death in ONS death (primary cause only)
    
    
    mutate(postdrug_hf=pmin(postdrug_first_primary_hhf, hf_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.hf.dataset <- hf.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.hf.dataset))
    print(table(final.hf.dataset$drugclass))
    
  }
  
  if (dataset.type == "hf.dataset") {
    return(final.hf.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ################# No comorbidities outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### No comorbidities outcome model")
    print("################################################")
    
  }
  
  no_co.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(no_co.dataset$multi_drug_start))
    print(table(no_co.dataset$multi_drug_start, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(no_co.dataset$timeprevcombo_less61))
    print(table(no_co.dataset$timeprevcombo_less61, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(no_co.dataset$hb_extreme_53))
    print(table(no_co.dataset$hb_extreme_53, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(no_co.dataset$prehba1c)))
    print(table(is.na(no_co.dataset$prehba1c), no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if comorbidity
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(comorbidities = ifelse(predrug_cvd == "No" & preheartfailure == "No" & (is.na(preckd) | (preckd!="stage_3a" & preckd!="stage_3b" & preckd!="stage_4")), NA_real_, 1))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if comorbidities")
    print("################################################")
    print(table(no_co.dataset$comorbidities, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(comorbidities))
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(no_co.dataset$TZD, no_co.dataset$drugclass))
    print("GLP1 treated")
    print(table(no_co.dataset$GLP1, no_co.dataset$drugclass))
    print("SGLT2 treated")
    print(table(no_co.dataset$SGLT2, no_co.dataset$drugclass))
    
  }
  
  
  no_co.dataset <- no_co.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  
  no_co.dataset <- no_co.dataset %>%
    
    
    mutate(postdrug_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
           postdrug_hf=pmin(postdrug_first_heartfailure, hf_death_date_any_cause, na.rm=TRUE)) %>%
    
    
    # Mace: broad definition: MI/stroke GP codes + broad MI/stroke HES codes (any cause) + CV death in ONS death (any cause)
    
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25,
           
           
           ## HF: broad definition: HF GP codes + HF HES codes (any cause) + HF death in ONS death (any cause)
           
           
           postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25)
  
  ## CKD: TBD
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.no_co.dataset <- no_co.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.no_co.dataset))
    print(table(final.no_co.dataset$drugclass))
    
  }
  
  if (dataset.type == "no_co.dataset") {
    return(final.no_co.dataset)
  }
  
  
  
  
}

#### define study cohorts ####

#Training and test datasets
  #2013 onwards
  earliestdrugstartdate <- "2013-01-01"
  md.train.2013 <- set_up_data("hba1c.train",drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"),earliestdrugstartdate)
  md.test.2013 <- set_up_data("hba1c.test",drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"),earliestdrugstartdate)
  
  save(md.train.2013,file=paste0(data_dir,"md.train.2013.Rda"))
  save(md.test.2013,file=paste0(data_dir,"md.test.2013.Rda"))
  
  #2004 onwards
  earliestdrugstartdate <- "2004-01-01"
  md.train.2004 <- set_up_data("hba1c.train",drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"),earliestdrugstartdate)
  md.test.2004 <- set_up_data("hba1c.test",drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"),earliestdrugstartdate)
  
  save(md.train.2004,file=paste0(data_dir,"md.train.2004.Rda"))
  save(md.test.2004,file=paste0(data_dir,"md.test.2004.Rda"))
  


#### quick load cohorts ####
  
  #2013 onwards
  load(paste0(data_dir,"md.train.2013.Rda"))
  load(paste0(data_dir,"md.test.2013.Rda"))
  
  #2004 onwards
  load(paste0(data_dir,"md.train.2004.Rda"))
  load(paste0(data_dir,"md.test.2004.Rda"))
  
#### Global settings ####
#Set global factors for prediction
hba1cmonth <- 12 
ncurrtx <- 2
drugline <- 1

#### Model development ####

md.train <- md.train.2013

describe(md.train)

#Candidate features
  # drugline
  # ncurrtx
  # hba1cmonth
  # agetx
  # sex
  # t2dmduration
  # pretotalcholesterol
  # prehdl
  # prealt
  # preegfr
  # prebmi
  # prehba1c

md.train.cc <- md.train %>% filter(complete.cases(pretotalcholesterol,
                                                  prehdl,
                                                  prealt,
                                                  preegfr,
                                                  prebmi))
# md.test.cc <- md.test %>% filter(complete.cases(pretotalcholesterol,
#                                                   prehdl,
#                                                   prealt,
#                                                   preegfr,
#                                                   prebmi))

formula1 <- "posthba1cfinal ~ drugclass +
  drugline +
  ncurrtx +
  rcs(hba1cmonth,3)*drugclass +
  rcs(agetx,3)*drugclass +
  sex*drugclass +
  rcs(t2dmduration,3)*drugclass +
  rcs(pretotalcholesterol,3)*drugclass +
  rcs(prehdl,3)*drugclass +
  rcs(prealt,3)*drugclass +
  rcs(preegfr,3)*drugclass +
  rcs(prebmi,3)*drugclass +
  rcs(prehba1c,3)*drugclass"

formula2 <- "posthba1cfinal ~ drugclass +
  drugline*drugclass +
  ncurrtx*drugclass +
  rcs(hba1cmonth,3)*drugclass +
  rcs(agetx,3)*drugclass +
  sex*drugclass +
  rcs(t2dmduration,3)*drugclass +
  rcs(pretotalcholesterol,3)*drugclass +
  rcs(prehdl,3)*drugclass +
  rcs(prealt,3)*drugclass +
  rcs(preegfr,3)*drugclass +
  rcs(prebmi,3)*drugclass +
  rcs(prehba1c,3)*drugclass"

formula3 <- 
  "posthba1cfinal ~ concordant +
  drugline +
  ncurrtx +
  rcs(hba1cmonth,3) +
  rcs(agetx,3) +
  sex +
  rcs(t2dmduration,3) +
  rcs(pretotalcholesterol,3) +
  rcs(prehdl,3) +
  rcs(prealt,3) +
  rcs(preegfr,3) +
  rcs(prebmi,3) +
  rcs(prehba1c,3)"

formula3 <- 
  "posthba1cfinal ~ concordant +
  drugline +
  ncurrtx +
  rcs(hba1cmonth,3) +
  rcs(agetx,3) +
  sex +
  rcs(t2dmduration,3) +
  rcs(pretotalcholesterol,3) +
  rcs(prehdl,3) +
  rcs(prealt,3) +
  rcs(preegfr,3) +
  rcs(prebmi,3) +
  rcs(prehba1c,3)"

formula4 <- 
  "posthba1cfinal ~ drugclass +
  drugline +
  ncurrtx +
  rcs(hba1cmonth,3) +
  rcs(agetx,3) +
  sex +
  rcs(t2dmduration,3) +
  rcs(pretotalcholesterol,3) +
  rcs(prehdl,3) +
  rcs(prealt,3) +
  rcs(preegfr,3) +
  rcs(prebmi,3) +
  rcs(prehba1c,3)"

#Set reference category
md.train.cc$drugclass <- relevel(md.train.cc$drugclass,ref="DPP4") 

#Set data dist for rms
ddist <- datadist(md.train.cc); options(datadist='ddist') 
 

# Variable importance formula 1
m1 <- ols(as.formula(formula1),data=md.train.cc,x=TRUE,y=TRUE)
m1
nobs(m1)
anova(m1,indnl=FALSE)
plot(anova(m1), margin=c('chisq', 'proportion chisq'))
plot(anova(m1), what='proportion R2')
plot(anova(m1), what='partial')

# Variable importance formula 2
m2 <- ols(as.formula(formula2),data=md.train.cc,x=TRUE,y=TRUE)
m2
nobs(m2)
anova(m2,indnl=FALSE)
plot(anova(m2), margin=c('chisq', 'proportion chisq'))
plot(anova(m2), what='proportion R2')
plot(anova(m2), what='partial')

#use formula1 as the model
pen<- pentrace(m2,
               list(simple=10*c(0.05,0.1,0.2,0.3,0.4,0.5,1,5,10,100,1000,10000),#, 15, 20, 30, 40, 100, 200, 500, 1000),
                    nonlinear=10*c(0.05,0.1,0.2,0.3,0.4,0.5,1,5,10,100,1000,10000),#, 15, 20, 30, 40, 100, 200, 500, 1000),
                    interaction=10*c(0.05,0.1,0.2,0.3,0.4,0.5,1,5,10,100,1000,10000)))#, 15, 20, 30, 40, 100, 200, 500, 1000)))
pen

m1 <- update(m1, penalty=list(simple=0.5, nonlinear=50, interaction=50))
effective.df(m1)
plot(anova(m1))

m2 <- update(m2, penalty=list(simple=0.5, nonlinear=50, interaction=50))
effective.df(m2)
plot(anova(m2))

Function(m1)

##Plot cont. variables and save model summary
#Baseline Hba1c  
#Predict 1-99th centile of POPULATION
quantile(md.train.cc$prehba1c, c(.01, .99), na.rm=TRUE)
c1 <- quantile(md.train.cc$prehba1c, .01, na.rm=TRUE)
c99 <- quantile(md.train.cc$prehba1c, .99, na.rm=TRUE)

w <- Predict(m2, drugclass=levels(md.train.cc$drugclass),prehba1c=seq(c1,c99,by=1),drugline=2,ncurrtx=1, hba1cmonth=12)
w <- data.frame(w)
w$yhat<- w$yhat - w$prehba1c
w$upper<- w$upper - w$prehba1c
w$lower<- w$lower - w$prehba1c

rplot.prehba1c <- ggplot(data = w, aes(x = prehba1c, y = yhat, group=drugclass)) + geom_line(aes(colour=drugclass), size = 1.5) + theme_bw() +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.2) + #scale_color_manual(values=c("royalblue", "red3")) +
  ylab("HbA1c response (mmol/mol)") + xlab("Baseline Hba1c (mmol/mol)") + theme(axis.text=element_text(size=rel(1.5)))+ theme(axis.title=element_text(size=rel(1.5)))+
  ggtitle("Baseline Hba1c") + theme(legend.text = element_text(colour="black", size=rel(1.5))) + 
  theme(legend.title=element_blank()) + theme(plot.margin = margin()) + geom_hline(yintercept = 0) +
  #scale_y_continuous(breaks=c(seq(12.5,0,by=2.5)), limits=c(12.5,0)) + coord_cartesian(ylim=c(12.5,0)) + 
  theme(legend.position = c(0.2, 0.2)) + theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) 
rplot.prehba1c

#### Predict outcomes on each drug

  md.train.cc <- md.train.cc %>% 
    mutate(drug=drugclass,
           drugclass="DPP4")
  md.train.cc$DPP4 <- predict(m1,md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(drugclass="SGLT2")
  md.train.cc$SGLT2 <- predict(m1,md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(drugclass="SU")
  md.train.cc$SU <- predict(m1,md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(drugclass="TZD")
  md.train.cc$TZD <- predict(m1,md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(drugclass="GLP1")
  md.train.cc$GLP1 <- predict(m1,md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(drugclass=drug) %>% 
    select(-drug)
  
  #Find best drug
  #https://stackoverflow.com/questions/37195322/create-a-new-variable-from-the-minimum-in-r
  library(data.table)
  setDT(md.train.cc)[, lowest.hba1c := apply(.SD, 1, min), .SDcols=c("DPP4", "SGLT2", "SU", "TZD", "GLP1")]
  md.train.cc[, bestdrug := apply(.SD, 1, function(x) names(x)[which.min(x)]), .SDcols = c("DPP4", "SGLT2", "SU", "TZD", "GLP1")]
  
  #Define concordant and discordant
  md.train.cc <- data.frame(md.train.cc)
  md.train.cc <- md.train.cc %>% 
    mutate(concordant = if_else(drugclass==bestdrug,1,0))

  head(md.train.cc)
  table(md.train.cc$concordant)
  table(md.train.cc$bestdrug)
  
  #Define drug pair subsets - not right - want to compare all predictions not just best drug
  md.train.cc.dpp4glp1 <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "GLP1") & (drugclass == "DPP4" | drugclass == "GLP1") )
  md.train.cc.dpp4sglt2 <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "SGLT2") & (drugclass == "DPP4" | drugclass == "SGLT2") )
  md.train.cc.dpp4su <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "SU") & (drugclass == "DPP4" | drugclass == "SU") )
  md.train.cc.dpp4tzd <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "TZD") & (drugclass == "DPP4" | drugclass == "TZD") )
  md.train.cc.glp1sglt2 <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "SGLT2") & (drugclass == "GLP1" | drugclass == "SGLT2") )
  md.train.cc.glp1su <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "SU") & (drugclass == "GLP1" | drugclass == "SU") )
  md.train.cc.glp1tzd <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "TZD") & (drugclass == "GLP1" | drugclass == "TZD") )
  md.train.cc.sglt2su <- md.train.cc %>% 
    filter((bestdrug == "SGLT2" | bestdrug == "SU") & (drugclass == "SGLT2" | drugclass == "SU") )
  md.train.cc.sglt2tzd <- md.train.cc %>% 
    filter((bestdrug == "SGLT2" | bestdrug == "TZD") & (drugclass == "SGLT2" | drugclass == "TZD") )
  md.train.cc.sutzd <- md.train.cc %>% 
    filter((bestdrug == "SU" | bestdrug == "TZD") & (drugclass == "SU" | drugclass == "TZD") )
  
  
  #DPP4GLP1
  md.train.cc.dpp4glp1 <- md.train.cc %>% 
    filter(drugclass == "DPP4" | drugclass == "GLP1") %>% 
    mutate(hba1c_diff = GLP1-DPP4,
           bestdrug=ifelse(hba1c_diff<=0,"GLP1","DPP4"),
           hba1c_diff.q = ntile(hba1c_diff, 10))
  
  #define dataset with predicted values
  t1 <- md.train.cc.dpp4glp1 %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarize(N=length(hba1c_diff),
              hba1c_diff.pred = mean(hba1c_diff))
    
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula4),data=md.train.cc.dpp4glp1,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.adj,lower.adj,upper.adj))
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  
  hte_plot <- function(data,pred,obs,obslowerci,obsupperci) {
    
    #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
    ymin  <- -16;  ymax <- 6
    
    ggplot(data=data,aes_string(x=pred,y=obs)) +
      geom_point(alpha=1) + theme_bw() +
      geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
      ylab("Observed HbA1c difference (mmol/mol)") + xlab("Predicted HbA1c difference (mmol/mol)") +
      scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
      scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
      # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
      # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
      theme_base() + geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
      geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") 
  }
  library(ggthemes)
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  
  
  #GLP1 SGLT2
  c1 <- md.train.cc %>% 
    filter(drugclass == "SGLT2" | drugclass == "GLP1") %>% 
    mutate(hba1c_diff = SGLT2-GLP1,
           bestdrug=ifelse(hba1c_diff<=0,"SGLT2","GLP1"),
           hba1c_diff.q = ntile(hba1c_diff, 10))
  
  #define dataset with predicted values
  t1 <- c1 %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarize(N=length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula4),data=c1,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.adj,lower.adj,upper.adj))
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  
  
  #DPP4 SGLT2
  c1 <- md.train.cc %>% 
    filter(drugclass == "SGLT2" | drugclass == "DPP4") %>% 
    mutate(hba1c_diff = SGLT2-DPP4,
           bestdrug=ifelse(hba1c_diff<=0,"SGLT2","DPP4"),
           hba1c_diff.q = ntile(hba1c_diff, 10))
  
  #define dataset with predicted values
  t1 <- c1 %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarize(N=length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula4),data=c1,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.adj,lower.adj,upper.adj))
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  
#SU TZD
  c1 <- md.train.cc %>% 
    filter(drugclass == "SU" | drugclass == "TZD") %>% 
    mutate(hba1c_diff = TZD-SU,
           bestdrug=ifelse(hba1c_diff<=0,"TZD","SU"),
           hba1c_diff.q = ntile(hba1c_diff, 10))
  
  #define dataset with predicted values
  t1 <- c1 %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarize(N=length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula4),data=c1,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.adj,lower.adj,upper.adj))
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  
#SU SGLT2
  c1 <- md.train.cc %>% 
    filter(drugclass == "SU" | drugclass == "SGLT2") %>% 
    mutate(hba1c_diff = SU-SGLT2,
           bestdrug=ifelse(hba1c_diff<=0,"SU","SGLT2"),
           hba1c_diff.q = ntile(hba1c_diff, 10))
  head(c1)
  #define dataset with predicted values
  t1 <- c1 %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarize(N=length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula4),data=c1,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.adj,lower.adj,upper.adj))
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  
  
  
  md.train.cc.dpp4sglt2 <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "SGLT2") & (drugclass == "DPP4" | drugclass == "SGLT2") )
  md.train.cc.dpp4su <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "SU") & (drugclass == "DPP4" | drugclass == "SU") )
  md.train.cc.dpp4tzd <- md.train.cc %>% 
    filter((bestdrug == "DPP4" | bestdrug == "TZD") & (drugclass == "DPP4" | drugclass == "TZD") )
  md.train.cc.glp1sglt2 <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "SGLT2") & (drugclass == "GLP1" | drugclass == "SGLT2") )
  md.train.cc.glp1su <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "SU") & (drugclass == "GLP1" | drugclass == "SU") )
  md.train.cc.glp1tzd <- md.train.cc %>% 
    filter((bestdrug == "GLP1" | bestdrug == "TZD") & (drugclass == "GLP1" | drugclass == "TZD") )
  md.train.cc.sglt2su <- md.train.cc %>% 
    filter((bestdrug == "SGLT2" | bestdrug == "SU") & (drugclass == "SGLT2" | drugclass == "SU") )
  md.train.cc.sglt2tzd <- md.train.cc %>% 
    filter((bestdrug == "SGLT2" | bestdrug == "TZD") & (drugclass == "SGLT2" | drugclass == "TZD") )
  md.train.cc.sutzd <- md.train.cc %>% 
    filter((bestdrug == "SU" | bestdrug == "TZD") & (drugclass == "SU" | drugclass == "TZD") )
  
    md.train.cc.dpp4su <- md.train.cc %>% 
    filter((bestdrug == "TZD" | bestdrug == "SU") & (drugclass == "TZD" | drugclass == "SU") )
  
#Function to fit a series of models and output the coefficient(s) of interest with CIs and p-value
    hte.model.coefs <- function(x,nmodels) {
      mnumber = c(1:nmodels)
      models <- as.list(1:nmodels)
      nobs <- vector()
      coef <- vector()
      lower <- vector()
      upper <- vector()
      pvalue <- vector()
      data <- x
      
      for(i in mnumber) {
        models[[i]] <- lm(as.formula(formula3),data=data)
        nobs <- append(nobs,nobs(models[[i]]))
        coef <- append(coef,models[[i]]$coefficients[2])
        confint_all <- confint(models[[i]], levels=0.95)
        lower <- append(lower,confint_all[2,1])
        upper <- append(upper,confint_all[2,2])
        pvalue <- append(pvalue,summary(models[[i]])$coefficients[2,4])
      }
      
      datasetname = c(deparse(substitute(x)),deparse(substitute(x)),deparse(substitute(x)))
      x <- data.frame(datasetname,modelname,cbind(nobs,coef,lower,upper,pvalue))
      rownames(x) <- c()
      return(x)
    }  
   
    hte.model.coefs(md.train.cc,1)
    
    md.train.cc.test <- md.train.cc %>% filter(drugclass!="DPP4")
    hte.model.coefs(md.train.cc.test,1)
    
    ols(as.formula(formula3),data=md.train.cc)
    
#### Model validation ####
md.test <- set_up_data("hba1c.test",drugs = c("GLP1", "SGLT2","DPP4","SU","TZD"))
