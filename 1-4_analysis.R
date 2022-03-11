# Purpose: Conduct primary analyses: 
  # 1) HR for relationship between policy time period and survival
  # 2) Predicted probabilities (24mos)
# Date Created: 2021/09/14

library(dplyr)
library(survival)
library(splines)
library(ggplot2)
library(survminer)
library(mice)

# Load non-imputed data
ucc <- readRDS("/data/analytic/imputedat_2021-09-16.rds")

# Model with interpretable interaction term
m1 <- with(ucc, 
               survival::coxph(Surv(smonths, dead) ~ treat + tperiod + treat*tperiod +
                                 gendf + smokey + reth + acadprac + 
                                 psitebladder + ecog24 + hvoldis + 
                                 insurecatf + 
                                 ns(age, df=4) +
                                 ns(bmi, df=4) + 
                                 ns(creat_clear, df=4) +
                                 ns(calday, df=3),
                               robust=T, 
                               cluster=practiceid,
                               ties = "breslow",
                               x = T
               )
)

# Same model but with interaction as factor to feed to predicted probabilities 
  # function
m2 <- with(ucc, 
           survival::coxph(Surv(smonths, dead) ~ treatperiod +
                             gendf + smokey + reth + acadprac + 
                             psitebladder + ecog24 + hvoldis + 
                             insurecatf + 
                             ns(age, df=4) +
                             ns(bmi, df=4) + 
                             ns(creat_clear, df=4) +
                             ns(calday, df=3),
                           robust=T, 
                           cluster=practiceid,
                           ties = "breslow",
                           x = T
           )
)


# Plot martingale residuals for continuous data
ra <- c(sapply(m1[[4]], FUN=residuals, type="martingale", simplify = T))
ucc.l <- mice::complete(ucc, "long", include=F)
ucc.l$m_resid <- ra

resplot.age <- ggplot(data=ucc.l, aes(x=age, y=m_resid)) +
  geom_point() +
  geom_smooth(method="loess") + 
  xlab("Age") +
  ylab("Martingale Residuals") +
  theme_bw()

resplot.bmi <- ggplot(data=ucc.l, aes(x=bmi, y=m_resid)) +
  geom_point() +
  geom_smooth(method="loess") + 
  xlab("BMI") +
  ylab("Martingale Residuals") + 
  theme_bw()

resplot.creat <- ggplot(data=ucc.l, aes(x=creat_clear, y=m_resid)) +
  geom_point() +
  geom_smooth(method="loess") + 
  xlab("Creatinine clearance") +
  ylab("Martingale Residuals") + 
  theme_bw()

resplot.calday <- ggplot(data=ucc.l, aes(x=calday, y=m_resid)) +
  geom_point() +
  geom_smooth(method="loess") + 
  xlab("Calendar day") +
  ylab("Martingale Residuals") + 
  theme_bw()

tiff("/results/aim1/figures/mart-resid-plots.tiff", res=300, width = 8, height = 8, units = 'in')
ggarrange(resplot.age, resplot.bmi, resplot.creat, resplot.calday, nrow=3, ncol=2)
dev.off()

# Pool results
library(parameters)
m1.p <- pool(m1)

m1.p.parms <- model_parameters(m1.p, ci=0.95, exponentiate=T)
write.csv(m1.p.parms, "/results/aim1/model-params.csv")

m2.p.parms <- model_parameters(pool(m2), ci=0.95, exponentiate=T)
write.csv(m2.p.parms, "/results/aim1/model-params2.csv")

# Contrast treatment pre vs treatment post
cperiod <- function(dat){
  pindex <- which(names(dat$analyses[[1]]$coefficients) == "tperiod")
  iindex <- which(names(dat$analyses[[1]]$coefficients) == "treat:tperiod")
  n <- dat$analyses[[1]]$n
  eff <- NULL
  pivar <- NULL
  
  for(i in 1:length(m1$analyses)){
    an <- dat$analyses[i][[1]]
    eff <- c(eff, an$coefficients[pindex] + an$coefficients[iindex])
    pvar <- an$var[pindex, pindex]
    ivar <- an$var[iindex, iindex]
    picov <- an$var[pindex, iindex]
    pivar <- c(pivar, pvar + ivar + 2*picov)
  }
  pf <- pool.scalar(eff, pivar)
  est <- pf$qbar
  tvar <- pf$t
  hrfinal <- tibble(HR=exp(est), 
                    lci95=exp(est - 1.96*sqrt(tvar)),
                    uci95=exp(est + 1.96*sqrt(tvar)),
                    p=2*pt(-abs(est/sqrt(pf$t)), n)
                    )
  return(hrfinal)
}
cperiod(m1)


# Generate 24 month survival probabilities
source("1-6_functions.R") 
survprobs24 <- survpp(m2, ucc, sts=F)
write.csv(survprobs24, "/results/aim1/survprobs.csv")
