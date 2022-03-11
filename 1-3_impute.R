# Purpose: Impute missing values using MICE
# Date created: 2021-08-31

library(dplyr)
library(mice)
library(survival)

# Load data and setup to work with MICE
ucc <- readRDS("/data/analytic/analyticdat_2021-09-14.rds")  %>%
  filter(tperiodf != "Washout")  %>%
  mutate(nahaz = nelsonaalen(., smonths, dead),
         treat = ifelse(cisplatin_chemo == 1, 0, 1),
         smokey = as.factor(smokey),
         gendf = as.factor(gendf),
         ecogvalue = factor(ecogvalue, levels=c(0,1,2,3,4), 
                            labels=c("0","1","2","3","4")),
         insurecatf = factor(ifelse(insurecat=="commercial", 0,
                                   ifelse(insurecat=="medicare", 1,
                                          ifelse(insurecat=="othergov", 2, 3))),
                             levels=c(0,1,2,3),
                             labels=c("Commercial", "Medicare", "Other gov", 
                                      "Other"))
         ) %>%
  select(reth, smokey, bmi, creat_clear, pdltest, acadprac,
         psitebladder, hvoldis, tperiod, treat, insurecatf,
         ecogvalue, gendf, age, smonths, dead, nahaz, calday,
         practiceid)

# Check missingness patterns
md.pattern(ucc)

# Set seed
#sample(1:100000, 1) # 40794

#Check missing data and set imputation methods
imp0 <- mice(ucc, m=1, maxit=0)
meth <- imp0$method
imp0$nmis

meth["reth"] <- "polyreg"
meth["smokey"] <- "logreg"
meth["bmi"] <- "pmm"
meth["creat_clear"] <- "pmm"
meth["ecogvalue"] <- "polr"
meth["gendf"] <- "logreg"
meth["insurecatf"] <- "polyreg"

#Setup predictor matrix
pm <- imp0$predictorMatrix
pm[,"smonths"] <- 0
pm[,"practiceid"] <- 0

# impute
imp10 <- mice(ucc, m=10, maxit=10, method=meth, predictorMatrix = pm, seed=40794)

# assess results
#mice::bwplot(imp10, ecogvalue + reth + smokey + gendf ~ .imp)
densityplot(imp10, .imp~ecogvalue)
densityplot(imp10, .imp~reth)
densityplot(imp10, .imp~bmi+creat_clear)

# Save imputed data prior to new variable creation
saveRDS(imp10, paste0("/data/analytic/imputedat-raw_", Sys.Date(), ".rds"))


# Dichotomize ECOG and save
imp10.l <- mice::complete(imp10, "long", include=T) %>%
  mutate(ecog24 = ifelse(ecogvalue %in% c("2", "3", "4"), 1, 0),
         creat_clear = ifelse(creat_clear > 150, 150, creat_clear),
         treatperiod = factor(ifelse(treat == 0 & tperiod == 0, 0, # cisp pre
                                        ifelse(treat == 0 & tperiod == 1, 1, # cisp post
                                               ifelse(treat == 1 & tperiod == 0, 2, # imm/carb pre
                                                      ifelse(treat == 1 & tperiod == 1, 3, NA)))), #imm/carb post
                                 levels=c(0,1,2,3),
                                 labels=c("Cisp, pre", "Cisp, post", 
                                          "Imm/Carbo, pre", "Imm/Carbo, post")),
         practiceid = as.numeric(factor(practiceid)))

imp10_2 <- as.mids(imp10.l)

saveRDS(imp10_2, paste0("/data/analytic/imputedat_", Sys.Date(), ".rds"))
