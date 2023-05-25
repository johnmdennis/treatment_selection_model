
#dataset
data <- md.train %>% filter(drugclass=="TZD")

#Alt 

#Add missing indicator
  describe(data$prealt)
  data <- data %>% mutate(prebmi.missing=if_else(is.na(data$prebmi), 1, 0),
                          preegfr.missing=if_else(is.na(data$preegfr), 1, 0),
                          preacr.missing=if_else(is.na(data$preacr), 1, 0),
                          prealt.missing=if_else(is.na(data$prealt), 1, 0),
                          prehdl.missing=if_else(is.na(data$prehdl), 1, 0),
                          premap.missing=if_else(is.na(data$premap), 1, 0),
                          pretotalcholesterol.missing=if_else(is.na(data$pretotalcholesterol), 1, 0),
                          pretriglyceride.missing=if_else(is.na(data$pretriglyceride), 1, 0))

  
#complete.case
  formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt"
  mi.cc <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE)
  mi.cc

#Regression imputation (following https://journals.sagepub.com/doi/10.1177/09622802231165001)
  
  #Predict ALT
    formula.imp <- 
      "prealt ~ drugline +
      ncurrtx +
      rcs(agetx,3) +
      sex +
      rcs(t2dmduration,3)+
      rcs(prehba1c,3)"
    
    m.alt <- ols(as.formula(formula.imp),data=data,x=TRUE,y=TRUE)
    m.alt
    
    data$prealt.i <- predict(m.alt,data)
    describe(data$prealt.i)
  
    #Define final predictor variable
    data <- data %>% mutate(prealt.c = round(if_else(alt.missing==0,prealt,prealt.i),0))
    describe(data$prealt.c)
     
    #RI: predictor only              
    formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt.c "
    mi.ri <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE)
    mi.ri
    
    #RI: predictor + missing indicator             
    formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt.c + alt.missing"
    mi.ri.mi <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE)
    mi.ri.mi
    
    #RI: predictor + missing indicator+ interaction prealt.c and missing indicator           
    formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt.c*alt.missing"
    mi.ri.mi.int <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE)
    mi.ri.mi.int
    
    lrtest(mi.ri.mi,mi.ri.mi.int)
    
#Regression imputation (following https://researchonline.lshtm.ac.uk/id/eprint/4655332/1/Estimating-treatment-effects-with-partially-observed-covariates-using-outcome-regression-with-missing-indicators.pdf)
    # For a continuous partially observed covariate, the missing indicator approach in outcome regression
    # replaces missing covariate values with some fixed value: the same value (for example, 0) is used for all
    # participants with that covariate missing. Both the modified covariate and the missing indicator R are then
    # included in the analysis model. 
    
    data <- data %>% mutate(prealt.mi = if_else(is.na(prealt),0,prealt))
    describe(data$prealt.mi)    
    
    formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt.mi + alt.missing"
    mi.ri.set0 <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE)
    mi.ri.set0
    
    #formula1 <- "posthba1cfinal ~ rcs(prehba1c,3) + prealt.mi*alt.missing"
    #mi.ri.set0.int <- ols(as.formula(formula1),data=data,x=TRUE,y=TRUE) #DNC
    #mi.ri.set0.int 
    
#Compare predictions
    data$hb.cc <- predict(mi.cc,data)
    data$mi.ri <- predict(mi.ri,data)
    data$mi.ri.mi <- predict(mi.ri.mi,data)
    data$mi.ri.mi.int <- predict(mi.ri.mi.int,data)
    data$mi.ri.set0 <- predict(mi.ri.set0,data)
    
    ggplot(data, aes(x = hb.cc, y = mi.ri)) + geom_point()
    ggplot(data, aes(x = mi.ri, y = mi.ri.mi)) + geom_point()
    ggplot(data, aes(x = mi.ri.mi, y = mi.ri.mi.int)) + geom_point()
    
    ggplot(data, aes(x = mi.ri.mi, y = mi.ri.set0)) + geom_point()
    ggplot(data, aes(x = mi.ri.mi.int, y = mi.ri.set0)) + geom_point()
 
    
    ggplot(subset(data,alt.missing==1), aes(x = mi.ri.mi, y = mi.ri.mi.int)) + geom_point()
    ggplot(subset(data,alt.missing==1), aes(x = mi.ri.mi.int, y = mi.ri.set0)) + geom_point()

#Conclusion
    # Based on the Sisk et al paper and the above, use:
    #     1) RI to define predicted values of biomarker X if missing
    #     2) Include missing indicator with interaction term in the outcome model (mi.ri.mi.int)
    
    #Steps for UK model
    
    #1 Define missing indicators for: prehba1c,prebmi,preegfr,preacr,prealt,prehdl,premap,pretotalcholesterol,pretriglyceride
    
    #2 Define a single RI model formula based on: agetx,sex,t2dmduration,
    #drugline,ncurrtx,prehba1c,prebmi,preegfr,preacr,prealt,prehdl,premap,pretotalcholesterol,pretriglyceride,
    #CVD, HF, Complication count,ethnicity, smok, deprivation quintile
    
    #3 Predict missing values for biomarkers in step 1
    
    #4 Fit a model with tx*feature*missingindicator(if required) using:
    #predictor set: agetx,sex,t2dmduration,prehba1c,prebmi,preegfr,preacr,prealt,prehdl,premap,pretotalcholesterol,pretriglyceride,
      #CVD, HF, Complication count
    #confounder set: ethnicity, smok, deprivation quintile, drugline,ncurrtx (as we assume these are not effect modifiers)
    
    #This will allow deployment for missing predictors
    #Sens analysis: redo with drugline+ncurrtx in predictor set
    
    #Cohort: geographical region 60:40
    
    #Dev: East Mid, East of Eng, London, South Central, SE Coast, SW
    #Val: NE, NW, WM, YH
    
    table(cprd$prac_region)
    cprd <- cprd %>% mutate(region=
                              ifelse(prac_region==1,"NorthEast",
                                     ifelse(prac_region==2,"NorthWest",
                                            ifelse(prac_region==3,"YorkshireHumber",
                                                   ifelse(prac_region==4,"EastMidlands",
                                                          ifelse(prac_region==5,"WestMidlands",
                                                                 ifelse(prac_region==6,"EastofEngland",
                                                                        ifelse(prac_region==7,"SouthWest",
                                                                               ifelse(prac_region==8,"SouthCentral",
                                                                                      ifelse(prac_region==9,"London",
                                                                                             ifelse(prac_region==10,"SouthEastCoast",NA
                                                                                             )))))))))))
    table(cprd$prac_region,cprd$region)
    describe(cprd$region)
    
    
    
#RI for multiple predictors
    
    
    #Predict ALT
    formula.imp <-"prealt ~ drugline + ncurrtx + rcs(agetx,3) + sex + rcs(t2dmduration,3)+ rcs(prehba1c,3)"
    m.alt <- ols(as.formula(formula.imp),data=data,x=TRUE,y=TRUE)
    data$prealt.i <- predict(m.alt,data)
    data <- data %>% mutate(prealt.c = round(if_else(alt.missing==0,prealt,prealt.i),0)) #Define final predictor variable

    describe(data$prealt.c)
    