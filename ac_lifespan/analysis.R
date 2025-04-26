library(survival)
library(survminer)

setwd("/scratch/fmorandi/internal/Ali/github/ac_lifespan")

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

data_s6ko = read.table("./data/survival_data_s6ko.tsv", header = T) %>%
  mutate(Birth = as.Date(Birth, format="%d-%b-%Y"))%>%
  mutate(Death = as.Date(Death, format="%d-%b-%Y"))%>%
  mutate(AgeAtDeath = as.numeric(Death-Birth))%>%
  mutate(AgeAtDeathMonths = AgeAtDeath / 30.4) %>%
  mutate(Status = 1)

table(data$Sex, data$Diet)
any(duplicated(data$MouseID))

dataM = subset(data, Sex == "M")
dataF = subset(data, Sex == "F")

#### PLOT ####

# Males
fitM = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = dataM)
p=ggsurvplot(fitM, data = dataM)
p1=p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Males")

# Females
fitF = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = dataF)
p=ggsurvplot(fitF, data = dataF)
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

# S6KO
fit_s6ko = survfit(Surv(AgeAtDeathMonths, Status) ~ Diet, data = data_s6ko)
p=ggsurvplot(fit_s6ko, data = data_s6ko)
p$plot+
  scale_x_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_color_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")+
  labs(x="Age (Months)", y="Survival (%)")+
  ggtitle("Sirt6 KO (Males and Females combined)")
ggsave("./results/survival_s6ko.pdf", width = 0.6*w, height = 0.3*h, units="mm")

#### STATS ####

sink("./results/survival_stats.txt")
cat("Males:\n")
surv_pvalue(fitM, method="n")
cat("Females:\n")
surv_pvalue(fitF, method="n")
cat("S6KO:\n")
surv_pvalue(fit_s6ko, method="n")
closeAllConnections()

#### EXTRA CHECKS ####

# Are there any cohort effects?
ggplot(dataM, aes(Birth, fill=Diet))+
  geom_bar(position="fill")
ggplot(dataF, aes(Birth, fill=Diet))+
  geom_bar(position="fill")
data %>%
  group_by(Birth, Sex) %>%
  summarize(n = n())
# Treatments are fairly balanced over the 2 birth cohorts but best to check

summary(coxph(Surv(AgeAtDeath, Status)~as.factor(Birth)+Diet, data = dataM))
summary(coxph(Surv(AgeAtDeath, Status)~as.factor(Birth)+Diet, data = dataF))
summary(coxph(Surv(AgeAtDeath, Status)~Diet, data = dataM))
summary(coxph(Surv(AgeAtDeath, Status)~Diet, data = dataF))
# Including DoB does not affect results and DoB is never sig

median(dataM[dataM$Diet == "Ctrl", "AgeAtDeath"]) / 365
median(dataM[dataM$Diet == "Fuco", "AgeAtDeath"]) / 365
median(dataM[dataM$Diet == "Fuco", "AgeAtDeath"]) / median(dataM[dataM$Diet == "Ctrl", "AgeAtDeath"])
median(dataF[dataF$Diet == "Fuco", "AgeAtDeath"]) / median(dataF[dataF$Diet == "Ctrl", "AgeAtDeath"])
max(data$AgeAtDeath)
