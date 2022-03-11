################################################################################
# Purpose: Generate tables of baseline descriptive statistics
################################################################################

library(dplyr)

ucc <- readRDS("/data/analytic/analyticdat_2021-08-31.rds")  %>%
  mutate(ecog24 = ifelse(ecogvalue > 1, 1, 0)) %>%
  filter(tperiodf != "Washout")


# Generate baseline tables #####################################################

# Set categorical and numerical variables to tabulate
cvars <- c("gender", "reth", "therapyf", "primarysite", "smokingstatus", "steroid", "opioid",
             "ecog24", "insurecat", 
             "pdltest", "pdlanyp")
             #"pdlcps", "pdlps", "pdlanyp",
             #"pdltest3060", "pdlcps3060", "pdlps3060", "pdlanyp3060")
nvars <- c("age", "bmi", "creat_clear")

# Initiate table generating function
btab.gen <- function(dat, catvars, numvars, una="ifany"){
  tabdat <- tibble(var = character(),
                   cat = character(),
                   npre = numeric(), 
                   ppre = numeric(),
                   qpre = numeric(),
                   npost = numeric(), 
                   ppost = numeric(), 
                   qpost = numeric(),
                   ntot = numeric(), 
                   ptot = numeric(),
                   qtot = numeric())
  
  for(i in 1:length(catvars)){
    v <- dat %>% pull(catvars[i])
    t <- dat %>% pull(tperiod)
    tn <- table(v, useNA = una)
    tp <- prop.table(tn)
    
    pn <- table(v, t, useNA = una)
    pp <- prop.table(pn, margin=2)
    
    for(j in 1:length(tn)){
      tabdat <- tabdat %>% 
        add_row(var=catvars[i], cat=as.character(names(tn)[j]), 
                npre=pn[j,1], ppre=pp[j,1], qpre=NA,
                npost=pn[j,2], ppost=pp[j,2], qpost=NA,
                ntot=tn[j], ptot=tp[j], qtot=NA)
    }
    if(una == "no"){
      ttotna <- table(is.na(v))
      pttot <- prop.table(ttotna)
      tna <- table(is.na(v), t)
      pna <- prop.table(tna)
      
      if(sum(is.na(v)) > 0){
        tabdat <- tabdat %>% 
          add_row(var=catvars[i], cat="Missing", 
                  npre=tna[2,1], ppre=NA,#pna0, 
                  qpre=NA,
                  npost=tna[2,2], ppost=NA,#pna1, 
                  qpost=NA,
                  ntot=ttotna[2], ptot=NA,#pnat, 
                  qtot=NA)
      }
    }
  }
  
  for(i in 1:length(numvars)){
    v <- dat %>% pull(numvars[i])
    v0 <- dat %>% filter(tperiod == 0) %>% pull(numvars[i])
    v1 <- dat %>% filter(tperiod == 1) %>% pull(numvars[i])
    
    t <- dat %>% pull(tperiod)
    
    mtot <- mean(v, na.rm=T)
    lqtot <- quantile(v, 0.25, na.rm=T)
    uqtot <- quantile(v, 0.75, na.rm=T)
    
    mv0 <- mean(v0, na.rm=T)
    lqv0 <- quantile(v0, 0.25, na.rm=T)
    uqv0 <- quantile(v0, 0.75, na.rm=T)
    
    mv1 <- mean(v1, na.rm=T)
    lqv1 <- quantile(v1, 0.25, na.rm=T)
    uqv1 <- quantile(v1, 0.75, na.rm=T)
    
    tabdat <- tabdat %>% 
      add_row(var=numvars[i], cat="mean(IQR)", 
              npre=mv0, ppre=lqv0, qpre=uqv0,
              npost=mv1, ppost=lqv1, qpost=uqv1,
              ntot=mtot, ptot=lqtot, qtot=uqtot)
    
    ttot <- table(is.na(v))
    pttot <- prop.table(ttot)
    t0 <- table(is.na(v0))
    pt0 <- prop.table(t0)
    t1 <- table(is.na(v1))
    pt1 <- prop.table(t1)
    
    if(length(ttot) > 1){
      nat <- ttot[2]
      pnat <- pttot[2]
      na0 <- t0[2]
      pna0 <- pt0[2]
      na1 <- t1[2]
      pna1 <- pt1[2]
      
      tabdat <- tabdat %>% 
        add_row(var=numvars[i], cat="Missing", 
                npre=na0, ppre=NA,#pna0, 
                qpre=NA,
                npost=na1, ppost=NA,#pna1, 
                qpost=NA,
                ntot=nat, ptot=NA,#pnat, 
                qtot=NA)
    }
  }
  return(tabdat)
}

# Main baseline table, all data
maintab <- btab.gen(ucc, catvars = cvars, numvars = nvars, una="no")
write.csv(maintab, file="/results/aim1/T1raw.csv", na = "")

# Stratified baseline table, only treatment == carbo
cvars2 <- cvars[!cvars == "therapyf"]
carbtab <- btab.gen(ucc[ucc$therapyf == "carboplatin",], catvars=cvars2, numvars=nvars, una="no")
write.csv(carbtab, file="/results/aim1/T1carbo.csv", na = "")

# Stratified baseline table, only treatment == cisplatin
cistab <- btab.gen(ucc[ucc$therapyf == "cisplatin",], catvars=cvars2, numvars=nvars, una="no")
write.csv(cistab, file="/results/aim1/T1cis.csv", na = "")

# Stratified baseline table, only treatment == immuno
immtab <- btab.gen(ucc[ucc$therapyf == "immuno",], catvars=cvars2, numvars=nvars, una="no")
write.csv(immtab, file="/results/aim1/T1immuno.csv", na = "")

# Generate risk table for a priori determination of probability cut ############
library(survminer)
library(survival)
library(ggplot2)
ggrisktable(fit=survfit(Surv(smonths, dead) ~ tperiodf + cispf, data=ucc.ni), 
            risk.table = T, data=ucc,
            break.time.by=6)
