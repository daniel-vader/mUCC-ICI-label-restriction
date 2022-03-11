# Functions
# Created 2021/09/14

# Load required packages
library(miceadds)
library(stdReg)
library(lmtest)

# Get predicted probs for each imputation (function)
# m = class mira (list of fitted models)
# d = imputed data sets
# cp = logical, calculate probs conditional on cancer if TRUE
# sts = logical, calculate cumulative probs instead of survival if TRUE
# cvar = conditional variable for setting up conditional by cancer, race / ethnicity,
# insurance, denovo
survpp <- function(m, d, cp = F, sts = T){
  m <- m$analyses
  d <- mids2datlist(d)
  d1 <- d[[1]]
  nimp <- m %>% length()
  
  # Create a data matrix with potential combination of values we want to explore.
  # Also set up a tibble to accommodate those combinations where we can store
  # probabilities.
  sdat <- tibble(est.pre.cisp = numeric(), 
                 lci.pre.cisp = numeric(), 
                 uci.pre.cisp = numeric(),
                 
                 est.post.cisp = numeric(), 
                 lci.post.cisp = numeric(), 
                 uci.post.cisp = numeric(), 
                 
                 est.pre.immcarb = numeric(), 
                 lci.pre.immcarb = numeric(), 
                 uci.pre.immcarb = numeric(),
                 
                 est.post.immcarb = numeric(), 
                 lci.post.immcarb = numeric(), 
                 uci.post.immcarb = numeric(),
                 
                 est.diff.cisp = numeric(), 
                 lci.diff.cisp = numeric(), 
                 uci.diff.cisp = numeric(),
                 p.diff.cisp = numeric(),
                 
                 est.diff.immcarb = numeric(), 
                 lci.diff.immcarb = numeric(),
                 uci.diff.immcarb = numeric(),
                 p.diff.immcarb = numeric(),
                 
                 est.dd = numeric(), 
                 lci.dd = numeric(),
                 uci.dd = numeric(),
                 p.dd = numeric())
  
  if(!cp){
    # Set # possible combinations (marginal, 1)
    combos <- 1
  } else { 
    #   yp <- d1 %>% tidyr::expand(ccrc, cnsc, cpan, cpro, crcc, cucc) %>%
    #     filter(ccrc + cnsc + cpan + cpro + crcc + cucc < 2)
    #   cvar.t <- tibble(ccrc = numeric(), cnsc = numeric(), cpan = numeric(),
    #                    cpro = numeric(), crcc = numeric(), cucc = numeric())
    # # Set # possible combinations
    # combos <- nrow(yp)
  }
  
  
  for(k in 1:combos){ # loop through potential combinations
    sp.store <- tibble(
      est.pre.cisp = numeric(), 
      se.pre.cisp = numeric(), 
      
      est.post.cisp = numeric(), 
      se.post.cisp = numeric(), 
      
      est.pre.immcarb = numeric(), 
      se.pre.immcarb = numeric(), 
      
      est.post.immcarb = numeric(), 
      se.post.immcarb = numeric(), 
      
      est.diff.cisp = numeric(), 
      se.diff.cisp = numeric(), 
      
      est.diff.immcarb = numeric(), 
      se.diff.immcarb = numeric(),
      
      est.dd = numeric(), 
      se.dd = numeric(),
    )
    
    sp.pooled <- sp.store
    
    
    
    for(i in 1:length(m)){ # loop through each imputation
      # Grab all model values in imputed data i where our combination set is met.
      tt <- d[[i]]
      
      # Create TF vector for subset and get n in subset
      # Since no years/cancers/period obs are missing this will be the same
      # in every imputation
      if(!cp){
        subtf <- T
      } else {
        # if(cvar == "cancer"){
        #   subtf <- ifelse(tt$ccrc == yp$ccrc[k] &
        #                     tt$cnsc == yp$cnsc[k] &
        #                     tt$cpan == yp$cpan[k] &
        #                     tt$cpro == yp$cpro[k] &
        #                     tt$crcc == yp$crcc[k] &
        #                     tt$cucc == yp$cucc[k], T, F)
        # } else if(cvar == "reth"){
        #   subtf <- ifelse(tt$rblack == yp$rblack[k] &
        #                     tt$rhisp == yp$rhisp[k] &
        #                     tt$rother == yp$rother[k], T, F)
        # } else if(cvar == "denov"){
        #   subtf <- ifelse(tt$denovo_met == yp$denovo_met[k], T, F)
        # } else if(cvar == "insur"){
        #   subtf <- ifelse(tt$igov == yp$igov[k] &
        #                     tt$iother == yp$iother[k], T, F)
        # } else if(cvar == "age"){
        #   subtf <- ifelse(tt$agegt75 == yp$agegt75[k], T, F)
        # }
        
      }
      
      tt$sub <- subtf
      
      sf <- stdCoxph(m[[i]], data=tt, X="treatperiod", clusterid = "practiceid", 
                     subsetnew = sub, t=24)
      
      covm <- sf$vcov[[1]]
      estv <- as.vector(sf$est)
      sev <- diag(covm)
      
      ## Get estimates and standard errors for diffs and diff in diffs
      # Diff 1: -cisp, pre + immcarbo, pre
      d1v <- c(-1, 1, 0, 0)
      est.dcisp <- estv[2] - estv[1]
      se.dcisp <- d1v %*% covm %*% d1v
      
      # Diff 2: -cisp, post + immcarbo, post
      d2v <- c(0, 0, -1, 1)
      est.dimmcarb <- estv[4] - estv[3]
      se.dimmcarb <- d2v %*% covm %*% d2v
      
      # Diff 3: -(-janmar2019 + aprjul2019) + (-janmar2020 + aprjul2020) 
      # = janmar2019 - aprjul2019 - janmar2020 + aprjul2020
      d3v <- c(1, -1, -1, 1)
      est.dd <- est.dimmcarb - est.dcisp
      se.dd <- d3v %*% covm %*% d3v
      
      # Store estimates and variances 
      sp.store <- sp.store %>% 
        add_row(
          est.pre.cisp = estv[1], 
          se.pre.cisp = sev[1], 
          
          est.post.cisp = estv[2], 
          se.post.cisp = sev[2], 
          
          est.pre.immcarb = estv[3], 
          se.pre.immcarb = sev[3], 
          
          est.post.immcarb = estv[4], 
          se.post.immcarb= sev[4], 
          
          est.diff.cisp = est.dcisp, 
          se.diff.cisp = se.dcisp, 
          
          est.diff.immcarb = est.dimmcarb, 
          se.diff.immcarb = se.dimmcarb,
          
          est.dd = est.dd, 
          se.dd = se.dd 
        )
      print(paste("Completed iteration", (k-1)*nimp + i, "of", nimp*combos))
    } # end imputation loop
    
    # pool results
    p <- NULL
    pval <- NULL
    subn <- sum(tt$sub)
    sp.store <- sp.store %>% as.data.frame()
    for(i in 1:(ncol(sp.store)/2)){
      ind <- i*2
      pf <- pool.scalar(sp.store[,ind-1], sp.store[,ind])
      #, n = subn, 
      #k=length(m[[1]]$coefficients))
      if(!sts){
        e <- pf$qbar
      } else if(i > 4){
        e <- -pf$qbar
      } else {
        e <- 1 - pf$qbar
      }
      
      lci <- e - 1.96*sqrt(pf$t)
      uci <- e + 1.96*sqrt(pf$t)
      
      if(i>4){
        pval <- c(pval, 2*pt(-abs(e/sqrt(pf$t)), nrow(tt)))
      }
      
      if(lci < 0 & i <= 4){lci <- 0}
      
      p <- c(p, e, lci, uci)
    }
    
    sdat <- sdat %>% add_row( 
      est.pre.cisp = p[1], 
      lci.pre.cisp = p[2], 
      uci.pre.cisp = p[3],
      
      est.post.cisp = p[4], 
      lci.post.cisp = p[5], 
      uci.post.cisp = p[6], 
      
      est.pre.immcarb = p[7], 
      lci.pre.immcarb = p[8], 
      uci.pre.immcarb = p[9],
      
      est.post.immcarb = p[10], 
      lci.post.immcarb = p[11], 
      uci.post.immcarb = p[12],
      
      est.diff.cisp = p[13], 
      lci.diff.cisp = p[14], 
      uci.diff.cisp = p[15],
      p.diff.cisp = pval[1],
      
      est.diff.immcarb = p[16], 
      lci.diff.immcarb = p[17],
      uci.diff.immcarb = p[18],
      p.diff.immcarb = pval[2],
      
      est.dd = p[19], 
      lci.dd = p[20],
      uci.dd = p[21],
      p.dd = pval[3]
      )
    
    # Store probabilities in tibble
    if(cp){
      # if(cvar == "cancer"){
      #   cvar.t <- cvar.t %>% add_row(ccrc = yp$ccrc[k], cnsc = yp$cnsc[k], cpan = yp$cpan[k],
      #                                cpro = yp$cpro[k], crcc = yp$crcc[k], cucc = yp$cucc[k])
      # } else if(cvar == "reth"){
      #   cvar.t <- cvar.t %>% add_row(rblack = yp$rblack[k], rhisp = yp$rhisp[k], 
      #                                rother = yp$rother[k])
      # } else if(cvar == "denov"){
      #   cvar.t <- cvar.t %>% add_row(denovo_met = yp$denovo_met[k])
      # } else if(cvar == "insur"){
      #   cvar.t <- cvar.t %>% add_row(igov = yp$igov[k], iother = yp$iother[k])
      # } else if(cvar == "age"){
      #   cvar.t <- cvar.t %>% add_row(agegt75 = yp$agegt75[k])
      # }
    } 
  }
  
  if(cp){
    sdat.f <- bind_cols(sdat, cvar.t)
  } else{
    sdat.f <- sdat
  }
  
  sdat.f <- sdat.f %>%
    relocate(est.pre.cisp, lci.pre.cisp, uci.pre.cisp, 
             est.post.cisp, lci.post.cisp, uci.post.cisp,
             est.diff.cisp, lci.diff.cisp, uci.diff.cisp,
             est.pre.immcarb, lci.pre.immcarb, uci.pre.immcarb,
             est.post.immcarb, lci.post.immcarb, uci.post.immcarb,
             est.diff.immcarb, lci.diff.immcarb, uci.diff.immcarb,
             est.dd, lci.dd, uci.dd, 
             p.diff.cisp, p.diff.immcarb, p.dd)
  
  
  return(sdat.f)
}