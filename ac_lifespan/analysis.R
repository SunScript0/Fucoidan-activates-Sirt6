library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

setwd("/scratch/fmorandi/internal/Ali/merging/ac_lifespan_new")

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

ccodes = c("#cccccc", "#666666", "#ccff66", "#99cc00")

#### LOAD DATA ####

data = read.table("./data/survival_data.tsv", header = T, fill=T) %>%
  mutate(Birth = as.Date(Birth, format="%d-%b-%Y"))%>%
  mutate(Death = as.Date(Death, format="%d-%b-%Y"))%>%
  mutate(AgeAtDeath = as.numeric(Death-Birth))%>%
  mutate(AgeAtDeathMonths = AgeAtDeath / 30.4) %>%
  mutate(Status = 1)

data_s6ko = read.table("./data/survival_data_s6ko.tsv", header = T, sep="\t") %>%
  mutate(Birth = as.Date(Birth, format="%Y-%m-%d"))%>%
  mutate(Death = as.Date(Death, format="%Y-%m-%d"))%>%
  mutate(AgeAtDeath = as.numeric(Death-Birth))%>%
  mutate(AgeAtDeathMonths = AgeAtDeath / 30.4) %>%
  mutate(Status = 1)

table(data$Sex, data$Diet)
any(duplicated(data$MouseID))

dataM = subset(data, Sex == "M")
dataF = subset(data, Sex == "F")

dataK_ko = subset(data_s6ko, Sex != "Unknown")
dataM_ko = subset(data_s6ko, Sex == "M")
dataF_ko = subset(data_s6ko, Sex == "F")

data = list(
  "all" = data,
  "F" = dataF,
  "M" = dataM,
  "s6ko_all" = data_s6ko,
  "s6ko_known" = dataK_ko,
  "s6ko_F" = dataF_ko,
  "s6ko_M" = dataM_ko
)
rm(data_s6ko, dataF, dataF_ko, dataM, dataM_ko, dataK_ko)

#### PLOT ####

fit = list()

##### Wild-type #####

# Males
fit$M = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$M)
p=ggsurvplot(fit$M, data = data$M)
p1=p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Males")

# Females
fit$F = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$F)
p=ggsurvplot(fit$F, data = data$F)
p2=p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none",
        axis.title.y=element_blank())+
  labs(x="Age (Months)")+
  ggtitle("Females")
p1+p2
ggsave("./results/survival.pdf", width = 1*w, height = 0.3*h, units="mm")

##### Sirt6 KO #####

# Combined all
fit$s6ko_all = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$s6ko_all)
p=ggsurvplot(fit$s6ko_all, data = data$s6ko_all)
p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Sirt6 KO (Combined)")
ggsave("./results/survival_s6ko_combined_all.pdf", width = 0.6*w, height = 0.3*h, units="mm")

# Combined, known sex only
fit$s6ko_known = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$s6ko_known)
p=ggsurvplot(fit$s6ko_known, data = data$s6ko_known)
p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Sirt6 KO (Combined, known)")
ggsave("./results/survival_s6ko_combined_known.pdf", width = 0.6*w, height = 0.3*h, units="mm")

# Males
fit$s6ko_M = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$s6ko_M)
p=ggsurvplot(fit$s6ko_M, data = data$s6ko_M)
p1=p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Males")

# Females
fit$s6ko_F = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data$s6ko_F)
p=ggsurvplot(fit$s6ko_F, data = data$s6ko_F)
p2=p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none",
        axis.title.y=element_blank())+
  labs(x="Age (Months)")+
  ggtitle("Females")
p1+p2
ggsave("./results/survival_s6ko.pdf", width = 1*w, height = 0.3*h, units="mm")

#### STATS ####

sink("./results/survival_stats.txt")
cat("======= Wild-types =======\n")
cat("Males:\n")
surv_pvalue(fit$M, method="n", data = data$M)
cat("Females:\n")
surv_pvalue(fit$F, method="n", data = data$F)
cat("======= Sirt6 KOs =======\n")
cat("Combined:\n")
surv_pvalue(fit$s6ko_all, method="n", data = data$s6ko_all)
cat("Combined, known sex:\n")
surv_pvalue(fit$s6ko_known, method="n", data = data$s6ko_known)
cat("Males:\n")
surv_pvalue(fit$s6ko_M, method="n", data = data$s6ko_M)
cat("Females:\n")
surv_pvalue(fit$s6ko_F, method="n", data = data$s6ko_F)
closeAllConnections()

#### EXTRA CHECKS ####

ggplot(data$M, aes(Birth, fill=Diet))+
  geom_bar(position="fill")
ggplot(data$F, aes(Birth, fill=Diet))+
  geom_bar(position="fill")

data$all %>%
  group_by(Birth, Sex) %>%
  summarize(n = n())

summary(coxph(Surv(AgeAtDeath, Status)~as.factor(Birth)+Diet, data = data$M))
summary(coxph(Surv(AgeAtDeath, Status)~as.factor(Birth)+Diet, data = data$F))
summary(coxph(Surv(AgeAtDeath, Status)~Diet, data = data$M))
summary(coxph(Surv(AgeAtDeath, Status)~Diet, data = data$F))
# Including DoB does not affect results -> no cohort effect

median(data$M[data$M$Diet == "Ctrl", "AgeAtDeath"]) / 365
median(data$M[data$M$Diet == "Fuco", "AgeAtDeath"]) / 365
median(data$M[data$M$Diet == "Fuco", "AgeAtDeath"]) / median(data$M[data$M$Diet == "Ctrl", "AgeAtDeath"])
median(data$F[data$F$Diet == "Fuco", "AgeAtDeath"]) / median(data$F[data$F$Diet == "Ctrl", "AgeAtDeath"])
