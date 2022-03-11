# Purpose: Conduct descriptive and analyses
  # 1) Proportion patients receiving each treatment type over time,
  #              stratified by PD-L1 status
  # 2) KM Survival curves stratified by therapy, time period, and PD-L1 status
# Date Created: 2021/09/07

library(dplyr)
library(survival)
library(ggplot2)
library(viridis)
library(ggquickeda)

# Load non-imputed data for descriptive analysis and primary (1)
ucc.ni <- readRDS("/data/analytic/analyticdat_2021-09-14.rds") %>%
  mutate(cispf = factor(cisplatin_chemo, levels=c(0,1), 
                        labels=c("ICI or carboplatin-based\nchemotherapy",
                                 "Cisplatin-based\nchemotherapy")),
         pdl1f = factor(ifelse(is.na(pdlanyp), 9, pdlanyp), levels=c(0,1,9), 
                        labels=c("Negative", "Positive", "Unknown")),
         tperiodf2 = factor(tperiod, levels=c(0,1), 
                            labels=c("Pre-label change","Post-label change")),
         efflab = rep_len(1:3, nrow(.)),
         efflab = factor(efflab, 
                         levels=c(1,2,3),
                         labels=c("Cisplatin: 1.03 (0.83, 1.29)",
                                  "Carboplatin/IO: 1.02 (0.88, 1.18)",
                                  "Diff-in-diffs: 0.99 (0.77, 1.27)")),
         fx = -1,
         fy = -1) %>%
  filter(tperiodf != "Washout")

# Survival by treatment and time period (all therapies)
surv.alltreat <- ggplot(ucc.ni, aes(time=smonths, 
                                    status = dead, 
                                    fill=therapyf, 
                                    color=therapyf,
                                    linetype=tperiodf)) +
  geom_km() + 
  #geom_kmband() +
  #geom_kmticks() +
  #facet_wrap(~tperiodf, ncol=2) + 
  ylim(0,1) +
  xlab("Time (months)") +
  ylab("Cumulative survival probability") +
  scale_color_viridis(discrete = T, end=0.8, option="C") +
  scale_fill_viridis(discrete = T, option="C") +
  theme_bw()

surv.alltreat$labels$group <- "Treatment"
surv.alltreat$labels$colour <- "Treatment"
surv.alltreat$labels$fill <- "Treatment"
surv.alltreat$labels$linetype <- "Time Period"

tiff("/results/aim1/figures/surv_alltreat.tiff", res=300, width = 8, height = 5, units = 'in')
surv.alltreat
dev.off()

# Survival by treatment and time period (Cisplatin vs immuno/carbo)
surv.2treat <- ggplot(ucc.ni, aes(time=smonths, 
                                    status = dead, 
                                    #fill=cispf, 
                                    color=cispf,
                                  
                                    #group=cispf,
                                  linetype=tperiodf2)) +
  geom_km() + 
  #geom_kmband() +
  #geom_kmticks() +
  #facet_wrap(~tperiodf, ncol=2) + 
  ylim(0,1) +
  xlab("Time (months)") +
  ylab("Cumulative survival probability") +
  labs(group = "Treatment",
       colour = "Treatment",
       fill = "Treatment",
       linetype = "Time Period",
       shape = "Adjusted mortality HR (95% CI)") +
  scale_color_viridis(discrete = T, end=0.8, option="C") +
  scale_fill_viridis(discrete = T, option="C")  +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

tiff("/results/aim1/figures/surv_2treat.tiff", res=300, width = 8, height = 5, units = 'in')
surv.2treat 
dev.off()

png("/results/aim1/figures/surv_2treat.png", res=300, width = 5, height = 4, units = 'in')
surv.2treat 
dev.off()

# Survival by PD-L1 status
surv.pdl <- ggplot(ucc.ni, aes(time=smonths, 
                                  status = dead, 
                                  fill=cispf, 
                                  color=cispf,
                                  #group=cispf,
                                  linetype=tperiodf)) +
  geom_km() + 
  #geom_kmband() +
  #geom_kmticks() +
  facet_wrap(~pdl1f, ncol=2) + 
  ylim(0,1) +
  xlab("Time (months)") +
  ylab("Cumulative survival probability") +
  ggtitle("Survival by time period, therapy, and PD-L1 status") +
  scale_color_viridis(discrete = T, end=0.8, option="C") +
  scale_fill_viridis(discrete = T, option="C") +
  theme_bw() 

surv.pdl$labels$group <- "Therapy"
surv.pdl$labels$colour <- "Therapy"
surv.pdl$labels$fill <- "Therapy"
surv.pdl$labels$linetype <- "Time Period"

tiff("/results/aim1/figures/surv_pdl1.tiff", res=300, width = 8, height = 5, units = 'in')
surv.pdl
dev.off()

# Treatment type over time ##################################################
library(lubridate)

# Bar plots of treatment type by time period. 
ucc <- readRDS("/data/analytic/analyticdat_2021-09-14.rds")  %>%
  filter(tperiodf != "Washout") %>%
  mutate(cispf = factor(cisplatin_chemo, levels=c(0,1), 
                        labels=c("ICI or \ncarboplatin-based chemotherapy",
                                 "Cisplatin-based chemotherapy")),
         pdlf = factor(ifelse(is.na(pdlanyp), 9, pdlanyp), levels=c(0,1,9),
                       labels=c("PDL\1 Negative", "PDL\1 Positive", "Unknown")),
         therapyf = factor(ifelse(therapyf == "carboplatin", 0,
                                  ifelse(therapyf == "cisplatin", 1, 2)),
                           levels=c(0,1,2),
                           labels=c("Carboplatin-based\nchemotherapy", 
                                    "Cisplatin-based\nchemotherapy", 
                                    "ICI"))) %>%
  group_by(tperiodf, pdlf, therapyf) %>%
  summarize(n=length(tperiodf)) %>%
    group_by(tperiodf, pdlf) %>%
    mutate(pdl.prop = n/sum(n)) %>%
    ungroup() 
levels(ucc$tperiodf) <- c("Pre-label\nchange", "Post-label\nchange", "Washout")

ucc2 <- ucc %>% 
  group_by(tperiodf, therapyf) %>%
  summarize(n = sum(n)) %>%
  group_by(tperiodf) %>%
  mutate(tp.prop = n/sum(n)) %>%
  ungroup()

ucc.t <- readRDS("/data/analytic/analyticdat_2021-09-14.rds")  %>%
  filter(tperiodf != "Washout") %>%
  group_by(tperiodf) %>%
  summarize(n = sum(pdltest),
            tp.prop=sum(pdltest)/length(tperiodf)) %>%
  ungroup()
levels(ucc.t$tperiodf) <- c("Pre-label\nchange", "Post-label\nchange", "Washout")
ucc.t$therapyf <- "PDL/1 Tested"

ucc.t2 <- rbind(ucc2, ucc.t)

dens.tp <- ggplot(ucc.t2, aes(x=tperiodf, y=tp.prop, fill=therapyf)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=paste0(sprintf("%.1f", round(tp.prop*100,1)), "%")), 
            position=position_dodge(width=0.9), vjust=-.25) +
  xlab("Time Period") +
  ylab("Proportion") +
  #ggtitle("Therapy distribution by time period") +
  scale_fill_viridis(discrete=T, end=0.85, option="C") +
  scale_y_continuous(labels=scales::percent) +
  labs(fill="Therapy")  +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())

tiff("/results/aim1/figures/therapy-dist-by-time-period.tiff", res=300, 
     width = 8, height = 5, units = 'in')
dens.tp
dev.off()

png("/results/aim1/figures/therapy-dist-by-time-period.png", res=300, 
     width = 5, height = 4, units = 'in')
dens.tp
dev.off()


dens.pdl <- ggplot(ucc, aes(x=tperiodf, y=pdl.prop, fill=therapyf)) +
  geom_bar(position=position_dodge(), stat="identity") +
  xlab("Time Period") +
  ylab("Proportion") +
  ggtitle("Therapy distribution by time period and PD-L1 status") +
  facet_wrap(~pdlf) +
  scale_fill_viridis(discrete=T, end=0.85, option="C") +
  labs(fill="Therapy") +
  theme_minimal()
tiff("/results/aim1/figures/therapy-dist-by-time-period-and-pdl1.tiff", res=300, 
     width = 8, height = 5, units = 'in')
dens.pdl
dev.off()

# Combine plots
library(ggpubr)
dens.tp.g <- dens.tp + theme(plot.margin=unit(c(.5, 0, 0, 1.1), "cm")) +
  scale_y_continuous(labels=scales::percent, limits=c(0,0.6))
surv.2treat.g <- surv.2treat + theme(plot.margin=unit(c(.5, 0, 0, 1.1), "cm"))

tiff("/results/aim1/figures/combfig.tiff", res=300, 
     width = 6, height = 10, units = 'in')
ggarrange(dens.tp.g, surv.2treat.g, nrow=2, vjust=2, 
          labels=c("(A)", "(B)"), heights=c(1.7,2))
dev.off()
