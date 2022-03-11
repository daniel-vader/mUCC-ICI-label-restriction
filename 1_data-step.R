# OUR Data formatting
# Date created: 2021-07-14
# Author: Daniel Vader
################################################################################

# Load packages ################################################################
library(dplyr)
library(lubridate)
library(readr)

# Set base path for data files.
# Not setting WD to avoid accidentally outputting
# files to the wrong place.
path_base <- "/data/original/edm_bladder_062021/"

# Set date of source data in MMYYYY format. May append this to output file name.
source_date <- "062021"

# Path to output processed data set
path_out <- "/data/analyticdat/"


################################################################################
### Import data
################################################################################

### Main data ###
therapy <- readr::read_csv(
  paste(path_base, "lineoftherapy.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

biomark <- readr::read_csv(
  paste(path_base, "enhanced_advurothelialbiomarkers.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

demo <- readr::read_csv(
  paste(path_base, "demographics.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

death <- readr::read_csv(
  paste(path_base, "enhanced_mortality_v2.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

lab <- readr::read_csv(
  paste(path_base, "Lab.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

vital <- readr::read_csv(
  paste(path_base, "vitals.csv", sep=""),
             col_types = list(col_character(), col_date(format = ""),
                              col_character(),
                              col_character(), col_character(),
                              col_character(), col_character(),
                              col_character(), col_date(format = ""),
                              col_character(), col_double(),
                              col_double(), col_double(),
                              col_double(), col_double())) %>%
  dplyr::rename_all(.funs = tolower)

visit <- readr::read_csv(
  paste(path_base, "visit.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

canc <- readr::read_csv(
  paste(path_base, "enhanced_advurothelial.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

ecogd <- readr::read_csv(
  paste(path_base, "ecog.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

diagnosis <- readr::read_csv(
  paste(path_base, "diagnosis.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

insur <- readr::read_csv(
  paste(path_base, "insurance.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

meds <- readr::read_csv(
  paste(path_base, "medicationadministration.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

################################################################################
# Therapy (line of therapy, startdate, enddate), 
# demographics (age, gender, race/ethnicity, practicetype)
################################################################################

# Base N = 7078
length(unique(therapy$patientid))

# Grab first line therapy and limit therapy start to between 
   # April 1, 2017 and March 1, 2020
wdat <- therapy %>% 
  filter(linenumber == 1, startdate >= "2017-04-01", startdate < "2020-03-01") %>%
  mutate(tstartd = startdate, tendd = enddate) %>%
  select(patientid, linename, tstartd, tendd)

# N = 2581
nrow(wdat)

# Grab demographics. Create age & race/ethnicity variables. 
   # Limit cohort to age 18 +
wdat2 <- wdat %>% 
  left_join(demo, by="patientid") %>%
  mutate(birthdate = ymd(paste0(birthyear, "-01-01")),
         age = interval(birthdate, tstartd) %/% years(1),
         reth = ifelse(!is.na(ethnicity), "latinx",
                       ifelse(!is.na(race),
                              ifelse(race == "Hispanic or Latino", "latinx",
                              ifelse(race == "Black or African American", "black",
                              ifelse(race == "White", "white", "other"))), 
                              NA)
                       )
         ) %>%
  filter(age >= 18) %>%
  select(-primaryphysicianid, -birthdate)

# N = 2581
nrow(wdat2)  

# output therapy names for classification
#write.csv(table(demo.a$linename), "/data/classifiers/1st-line-therapy-table_raw.csv")

# Classify therapies and apply therapy exclusion
classify_therapy <- read_csv("/data/classifiers/1st-line-therapy-classified_2021-08-13.csv") %>%
  rename(linename = Var1,
         therapy_drop = drop) %>%
  select(-Freq)

wdat2t <- wdat2 %>% left_join(classify_therapy, by="linename") %>%
  filter(is.na(therapy_drop)) %>%
  mutate(therapyf = ifelse(immunotherapy == 1, "immuno",
                           ifelse(cisplatin_chemo == 1, "cisplatin",
                                  ifelse(carboplatin_chemo == 1, "carboplatin", NA)
                                  )
                           )
         ) %>%
  select(-therapy_drop)

# 2328
nrow(wdat2t)

################################################################################
# Add advanced diagnosis date, primary site, smoking status
################################################################################

wdat3 <- wdat2t %>% left_join(canc, by="patientid") %>%
  select(-diagnosisdate, -diseasegrade, -groupstage, -tstage, -nstage, -mstage, 
         -surgery, -surgerydate, -surgerytype) %>%
  mutate(psitebladder = ifelse(primarysite == "Bladder", 1, 0))


################################################################################
# Add last visit date, filter out subjects with < 2 encounters or no encounter
#  within 90 days of diagnosis
################################################################################

visit.a <- wdat3 %>% 
  select(patientid, advanceddiagnosisdate) %>% 
  full_join(visit, by="patientid") %>%
  group_by(patientid) %>%
  filter(visitdate > advanceddiagnosisdate) %>%
  summarize(lastvisit = max(visitdate),
            tfirstvisit = interval(advanceddiagnosisdate[1], min(visitdate)) %/% days(1),
            nvisit = sum(interval(advanceddiagnosisdate[1], visitdate) %/% days(1) > 0))

wdat4 <- wdat3 %>% 
  full_join(visit.a, by="patientid") %>%
  filter(tfirstvisit <= 90, nvisit >= 2)

# N = 2071
nrow(wdat4)

################################################################################
# Opioid and steroid status
################################################################################
steroid <- meds %>% 
  filter(drugcategory == "steroid", route == "Oral") %>%
  mutate(oralsteroid = 1) %>%
  rename(steroiddate = administereddate) %>%
  select(patientid, oralsteroid, steroiddate) 

painclassifier <- read.csv("/data/classifiers/opioid-classified_2021-08-05.csv") %>%
  rename(commondrugname = Var1) %>%
  select(-Freq)

pain <- meds %>% 
  filter(drugcategory == "pain agent") %>%
  right_join(painclassifier, by="commondrugname") %>%
  filter(opioid == 1) %>%
  rename(opioiddate = administereddate) %>%
  select(patientid, opioid, opioiddate) 

steroidsel <- wdat4 %>%
  left_join(steroid, by="patientid") %>%
  group_by(patientid) %>%
  mutate(sterdist = interval(tstartd, steroiddate) %/% days(1)) %>%
  filter(sterdist >= -60 & sterdist <= 30) %>%
  slice_min(abs(sterdist)) %>%
  summarize(steroid = max(oralsteroid),
            sterdist = sterdist[1]) %>%
  ungroup() 

painsel <- wdat4 %>%
  left_join(pain, by="patientid") %>%
  group_by(patientid) %>%
  mutate(opioiddist = interval(tstartd, opioiddate) %/% days(1)) %>%
  filter(opioiddist >= -60 & opioiddist <= 30) %>%
  slice_min(abs(opioiddist)) %>%
  summarize(opioid = opioid[1],
            opioiddist = opioiddist[1]) %>%
  ungroup() 

wdat5 <- wdat4 %>% 
  left_join(steroidsel, by="patientid") %>%
  mutate(steroid = ifelse(!is.na(steroid), steroid, 0)) %>%
  left_join(painsel, by="patientid") %>%
  mutate(opioid = ifelse(!is.na(opioid), opioid, 0)) %>%
  # High volume disease
  mutate(hvoldis = ifelse(steroid + opioid > 0, 1, 0))

# # output steroid names for classification
# pain <- meds %>% 
#   filter(drugcategory == "pain agent" | drugcategory == "steroid") %>%
#   select(patientid, commondrugname, drugcategory, administereddate, route) %>%
#   filter(patientid %in% wdat4$patientid)

# write.csv(table(painsteroid$commondrugname[painsteroid$drugcategory == "steroid"],
#                 painsteroid$route[painsteroid$drugcategory == "steroid"]),
#           "/data/classifiers/steroid-table_raw.csv")
# 
# # output steroid names for classification
# write.csv(table(painsteroid$commondrugname[painsteroid$drugcategory == "pain agent"]),
#           "/data/classifiers/opioid-table_raw.csv")



################################################################################
# ECOG PS
################################################################################
ecog6030 <- wdat5 %>% select(patientid, tstartd) %>% 
  full_join(ecogd, by="patientid") %>%
  group_by(patientid) %>%
  mutate(ecogdist = interval(tstartd, ecogdate) %/% days(1)) %>%
  filter(ecogdist >= -60 & ecogdist <= 30) %>%
  slice_min(abs(ecogdist)) %>%
  summarize(ecogvalue = mean(ecogvalue),
            ecogdist = ecogdist[1]) %>%
  ungroup() %>%
  mutate(ecogvalue = ceiling(ecogvalue))
  
wdat6 <- wdat5 %>% left_join(ecog6030, by="patientid")

################################################################################
# BMI (vitals)
################################################################################
bmigrab <- wdat6 %>% select(patientid, tstartd) %>%
  left_join(vital, by="patientid") %>%
  filter(testbasename %in% c("body height", "body weight")) %>%
  group_by(patientid) %>%
  mutate(weightkg = ifelse(testbasename == "body weight", testresultcleaned, NA),
         heightcm = ifelse(testbasename == "body height", testresultcleaned, NA),
         mdist = interval(tstartd, testdate) %/% days(1))  %>%
  filter(mdist >= -60 & mdist <= 30)

weightd <- bmigrab %>% filter(!is.na(weightkg)) %>%
  group_by(patientid) %>%
  slice_min(abs(mdist)) %>%
  summarize(weightkg = mean(weightkg),
            weightdist = mdist[1]) %>%
  ungroup() 

heightd <- bmigrab %>% filter(!is.na(heightcm)) %>%
  group_by(patientid) %>%
  slice_min(abs(mdist)) %>%
  summarize(heightcm = mean(heightcm),
            heightdist = mdist[1]) %>%
  ungroup()

weightheight <- full_join(weightd, heightd, by="patientid") %>%
  mutate(bmi = weightkg/((heightcm/100)^2)) 

wdat7 <- wdat6 %>% left_join(weightheight, by="patientid")

################################################################################
# Creatinine (lab), "creatinine", 
#                   "glomerular filtration rate/1.73 sq m.predicted",
#                   "glomerular filtration rate/1.73 sq m.predicted.black", <- I don't think we should use these
#                   "glomerular filtration rate/1.73 sq m.predicted.non black",
# CCr={((140-age) x weight)/(72 x creatinine)}x 0.85 (if female)
  # CCr (creatinine clearance) = mL/minute
  # Age = years
  # Weight = kg
  # Creatinine = mg/dL
################################################################################
creatgrab <- wdat7 %>% select(patientid, tstartd) %>%
  left_join(lab, by="patientid") %>%
  filter(testbasename %in% c("creatinine"#, 
                             # "glomerular filtration rate/1.73 sq m.predicted",
                             # "glomerular filtration rate/1.73 sq m.predicted.black",
                             # "glomerular filtration rate/1.73 sq m.predicted.non black"
                             )) %>%
  group_by(patientid) %>%
  mutate(mdist = interval(tstartd, testdate) %/% days(1)) %>%
  filter(mdist >= -60 & mdist <= 30) %>%
  slice_min(abs(mdist)) %>%
  summarize(creat = mean(testresultcleaned),
            creatunits = testunitscleaned[1],
            creatdist = mdist[1]
  ) %>%
  ungroup() 

# All units checkout (mg/dL))
#table(creatgrab$creatunits)

# Add creatinine info to overall data. Calculate clearance (mL/minute)
wdat8 <- wdat7 %>%
  left_join(creatgrab, by="patientid") %>%
  mutate(creat_clear = (((140-age)*weightkg)/(72*creat))*(1-.15*I(gender=="F")),
         creat_clear = ifelse(creat_clear > 150, 150, creat_clear))

# creatd <- egfrgrab %>% 
#   filter(testbasename == "creatinine") %>%
#   group_by(patientid) %>%
#   slice_min(abs(mdist)) %>%
#   summarize(creat = mean(testresultcleaned),
#             creatunits = testunitscleaned[1],
#             creatdist = mdist[1]
#   ) %>%
#   ungroup()

# gfrp <- egfrgrab %>% 
#   filter(testbasename == "glomerular filtration rate/1.73 sq m.predicted") %>%
#   group_by(patientid) %>%
#   slice_min(abs(mdist)) %>%
#   summarize(gfrp = mean(testresultcleaned),
#             gfrp.dist = mdist[1]
#             ) %>%
#   ungroup()
# 
# gfrp.b <- egfrgrab %>% 
#   filter(testbasename == "glomerular filtration rate/1.73 sq m.predicted.black") %>%
#   group_by(patientid) %>%
#   slice_min(abs(mdist)) %>%
#   summarize(gfrp.b = mean(testresultcleaned),
#             gfrp.b.dist = mdist[1]
#   ) %>%
#   ungroup()
# 
# gfrp.nb <- egfrgrab %>% 
#   filter(testbasename == "glomerular filtration rate/1.73 sq m.predicted.non black") %>%
#   group_by(patientid) %>%
#   slice_min(abs(mdist)) %>%
#   summarize(gfrp.nb = mean(testresultcleaned),
#             gfrp.nb.dist = mdist[1]
#   ) %>%
#   ungroup()
# 
# wdat8 <- wdat7 %>% 
#   left_join(creatd, by="patientid") %>%
#   left_join(gfrp, by="patientid") %>%
#   left_join(gfrp.b, by="patientid") %>%
#   left_join(gfrp.nb, by="patientid")


################################################################################
# PD-L1 (biomark) 
  # Dako combined + score >= 10%
  # Ventana tumor infiltrating immune cells >= 5%
################################################################################
pdlgrab <- biomark %>%
  filter(biomarkername == "PDL1", 
         assay %in% c("Dako PD-L1 IHC 22C3 pharmDx (Keytruda companion diagnostic)",
                      "Dako PD-L1 IHC 28-8 pharmDx (Opdivo complementary diagnostic)",
                      "Ventana PD-L1 (SP142) Assay (Tecentriq complementary diagnostic)",
                      "Ventana PD-L1 (SP263) Assay (Imfinzi complementary diagnostic)")) %>%
  mutate(assayns = ifelse(grepl('Dako', assay), "Dako", "Ventana"),
         pdlcps = ifelse(is.na(combinedpositivescore) | 
                                 combinedpositivescore == "Unknown/not documented", NA, 
                         ifelse(!(combinedpositivescore %in% c("0",
                                                               "<1",
                                                               "1",
                                                               "2-4",
                                                               "5-9")), 1, 0)),
         pdlps = ifelse(is.na(percentstaining), NA, 
                        ifelse(!(percentstaining %in% c("< 1%",
                                                 "0%",
                                                 "1%",
                                                 "2% - 4%")), 1, 0)),
         pdlstat = ifelse(assayns == "Dako", pdlcps,
                                 ifelse(assayns == "Ventana", pdlps, NA)),
         pdltest = 1,
         pdltextbin = ifelse(biomarkerstatus == "PD-L1 positive", 1, 
                             ifelse(biomarkerstatus == "PD-L1 negative/not detected", 0, NA)),
         pdlanyp = ifelse(!is.na(pdltextbin) & pdltextbin == 1, 1,
                          ifelse(!is.na(pdlcps) & pdlcps == 1, 1,
                                 ifelse(!is.na(pdlps) & pdlps == 1, 1,
                                        ifelse(is.na(pdltextbin) & 
                                                 is.na(pdlcps) &
                                                 is.na(pdlps), NA, 0))))
         ) %>% 
  left_join(wdat8, by="patientid") %>%
  group_by(patientid) %>%
  mutate(pdldist = interval(tstartd, specimencollecteddate) %/% days(1)) %>%
  ungroup()

pdlgrab3060 <- pdlgrab %>%
  group_by(patientid) %>%
  filter(pdldist >= -60 & pdldist <= 30) %>%
  #slice_min(abs(pdldist)) %>%
  summarize(pdlstat3060 = if(all(is.na(pdlstat))){NA}else{max(pdlstat, na.rm=T)},
            #pdldist3060 = pdldist[1],
            pdlcps3060 =  if(all(is.na(pdlcps))){NA}else{max(pdlcps, na.rm=T)},
            pdlps3060 = if(all(is.na(pdlps))){NA}else{max(pdlps, na.rm=T)},
            pdltest3060 = pdltest[1],
            pdlanyp3060 = if(all(is.na(pdlanyp))){NA}else{max(pdlanyp, na.rm=T)},
            #pdlbmstat3060 = biomarkerstatus
            
  ) %>%
  ungroup()

pdlgrabfull <- pdlgrab %>% filter(pdldist <= 30) %>%
  group_by(patientid) %>%
  #slice_min(abs(pdldist)) %>%
  summarize(pdlstat = if(all(is.na(pdlstat))){NA}else{max(pdlstat, na.rm=T)},
            #pdldist = pdldist[1],
            pdlcps =  if(all(is.na(pdlcps))){NA}else{max(pdlcps, na.rm=T)},
            pdlps = if(all(is.na(pdlps))){NA}else{max(pdlps, na.rm=T)},
            pdltest = pdltest[1],
            pdlanyp = if(all(is.na(pdlanyp))){NA}else{max(pdlanyp, na.rm=T)},
            #pdlbmstat = biomarkerstatus
            
  ) %>%
  ungroup()
  
wdat9 <- wdat8 %>% 
  left_join(pdlgrabfull, by="patientid") %>%
  left_join(pdlgrab3060, by="patientid")


################################################################################
# Insurance 
################################################################################

insgrab <- insur %>%
  rename(pc = payercategory) %>%
  mutate(
    # if start date but not end date is missing, set start date to 1 year before
    # end date
    isd = if_else(is.na(startdate) & !is.na(enddate), enddate - years(1), startdate),
    
    # if end date but not start date is missing, set end date to 2022.
    ied = if_else(is.na(enddate) & !is.na(startdate), ymd("2022/01/01"), enddate),
    coverp = interval(isd, ied),
    icat = ifelse(pc == "", "blank",
                  ifelse(pc == "Self Pay", "selfpay",
                         ifelse(pc %in% c("Medicaid", "Other Government Program"), "othergov",
                                ifelse(pc %in% c("Medicare"), "medicare",
                                  ifelse(pc %in% c("Other Payer - Type Unknown", 
                                                   "Patient Assistance Program",
                                                   "Workers Compensation"), "other",
                                         ifelse(pc == "Commercial Health Plan", 
                                                "commercial", "drop catch"))))))) %>%
  right_join(wdat9, by="patientid") %>%
  filter(tstartd %within% coverp) %>%
  mutate(insurey = ifelse(icat %in% c("selfpay","blank"), 0, 1)) %>%
  group_by(patientid) %>%
  
  # Collapse patients with multiple insurance types using an a priori hierarchy
  summarize(insure_y = max(insurey),
            insurecat = 
              if(all(is.na(icat) | icat == "blank")){"uncategorized"
                }else if(any(icat == "commercial")){"commercial" 
                  }else if(any(icat == "medicare")){"medicare"
                    }else if(any(icat == "othergov")){"othergov"
                    }else if(any(icat == "other")){"other"
                    }else if(any(icat == "selfpay")){"other"
                        }else{"drop catch"}
            )

wdat10 <- wdat9 %>% 
  left_join(insgrab, by="patientid") %>%
  mutate(insurecat =  ifelse(is.na(insurecat), "other", insurecat),
         tperiod = ifelse(tstartd %within% interval(ymd("2017-04-01"), ymd("2018-05-17")), 0,
                          ifelse(tstartd %within% interval(ymd("2018-06-20"), ymd("2020-03-01")), 1,
                                 9)),
         tperiodf = factor(tperiod, levels=c(0,1,9), 
                          labels=c("2017-04-01 to 2018-05-17", 
                                   "2018-06-20 to 2020-03-01", 
                                   "Washout")))

# Add survival variables
wdat11 <- wdat10 %>%
  left_join(death, by="patientid") %>%
  # Death only specific to month, convert to date
  mutate(dmonth = ifelse(!is.na(dateofdeath), paste0(dateofdeath, "-15"), NA),
         dmonth = ymd(dmonth),
         # If for some reason there is a visit that occurs after death "day", 
            # set death day to day of last visit. Use if_else to preseve date class
         dmonth = if_else(lastvisit > dmonth, lastvisit, dmonth), 
         # Set death status (event status)
         dead = ifelse(is.na(dmonth), 0, 1),
         # Summarize survival time in months
         smonths = if_else(dead == 0, interval(tstartd, lastvisit) %/% months(1),
                          interval(tstartd, dmonth) %/% months(1)),
         # treatment month
         calday = yday(tstartd),
         # clean smoking status
         smokey = ifelse(smokingstatus == "Unknown/not documented", NA,
                         ifelse(smokingstatus == "History of smoking", 1, 0)),
         # clean practice type
         acadprac = ifelse(practicetype == "ACADEMIC", 1, 0),
         # clean gender
         gendf = ifelse(gender == "F", 1, 0),
         # clean pdltest
         pdltest = ifelse(is.na(pdltest), 0, 1),
         # clean race/ethnicity
         reth = ifelse(reth == "white", 1, 
                       ifelse(reth == "black", 2,
                              ifelse(reth == "latinx", 3, 4))),
         reth = factor(reth, levels = c(1,2,3,4), 
                          labels=c("white", "black", "latinx", "other"))
         )


saveRDS(wdat11, paste0("/data/analytic/analyticdat_", Sys.Date(), ".rds"))
