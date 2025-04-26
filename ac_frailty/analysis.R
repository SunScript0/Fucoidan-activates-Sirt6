library(lmerTest)

setwd("/scratch/fmorandi/internal/Ali/github/ac_frailty")

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701
ccodes=c("#666666", "#99CC00")

#### LOAD DATA ####

n_per_cat = data.frame(
  Category = c("WeightScore", "Integument", "PhysicalMusculoskeletal",
               "VestibulocochlearAuditory", "OcularNasal", "DigestiveUrogenital",
               "Respiratory", 'Discomfort'),
  N = c(1, 5, 8, 2, 7, 4, 1, 2))

# Normalize categories by number of subscores
# Note: it would also be fine to calc overall score by summing sub-scores and 
# dividing by total, but this way sub-scores are brought to the same scale

# Excluding weightscore because unsure of validity, also sus cages
data = read.csv("./data/frailty_raw.csv") %>%
  mutate(Cage = str_extract(AnimalID, "^(\\d+)", group=1), .after=AnimalID) %>%
  mutate(WeightScore = WeightScore / 1) %>%
  mutate(Integument = Integument / 5) %>%
  mutate(PhysicalMusculoskeletal = PhysicalMusculoskeletal / 8) %>%
  mutate(VestibulocochlearAuditory = VestibulocochlearAuditory / 2) %>%
  mutate(OcularNasal = OcularNasal / 7) %>%
  mutate(DigestiveUrogenital = DigestiveUrogenital / 4) %>%
  mutate(Respiratory = Respiratory / 1) %>%
  mutate(Discomfort = Discomfort / 2) %>%
  mutate(Total = 
           Integument*5/29+
           PhysicalMusculoskeletal*8/29+
           VestibulocochlearAuditory*2/29+
           OcularNasal*7/29+
           DigestiveUrogenital*4/29+
           Respiratory*1/29+
           Discomfort*2/29)
write.csv(data, "./data/frailty_norm.csv", quote = F)

#### MAIN PLOT ####

ggplot(data, aes(as.factor(Month), Total, fill=Diet))+
  geom_boxplot()+
  scale_fill_manual(values=ccodes)+
  facet_grid(cols=vars(Sex))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))+
  labs(x = "Month", y="Frailty Score")
ggsave("./results/frailty_main.pdf", width = w, height = 0.3*h, units="mm")

#### WEIGHT PLOT ####

ggplot(data, aes(as.factor(Month), Weight, fill=Diet))+
  geom_boxplot()+
  scale_fill_manual(values=ccodes)+
  facet_grid(cols=vars(Sex))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))+
  labs(x = "Month", y="Body Mass (g)")
ggsave("./results/weight_plot.pdf", width = w, height = 0.3*h, units="mm")

#### CATEGORY BREAKDOWN ####

data_long = data %>%
  mutate(Month = as.factor(paste(Month, "mo"))) %>%
  pivot_longer(cols=-c(Month:Cage)) %>%
  dplyr::filter(!name %in% c("Weight", "WeightScore", "Temperature"))
data_summary = data_long %>%
  group_by(Month, Diet, Sex, name) %>%
  summarize(mean = mean(value), sd=sd(value), n=n()) %>%
  mutate(sem = sd / sqrt(n))%>%
  mutate(name = factor(name, c("WeightScore", "Integument", "PhysicalMusculoskeletal", "VestibulocochlearAuditory", "OcularNasal", "DigestiveUrogenital", "Respiratory", "Discomfort", "Total"))) %>%
  ungroup()%>%
  mutate(name = fct_recode(name,
                           "Physical\nMusculoskeletal" = "PhysicalMusculoskeletal", 
                           "Vestibulocochlear\nAuditory" = "VestibulocochlearAuditory", 
                           "Ocular\nNasal" = "OcularNasal", 
                           "Digestive\nUrogenital" = "DigestiveUrogenital"))
dodge = position_dodge(width=0.7)
p1 = data_summary %>%
  dplyr::filter(Sex == "M") %>%
  ggplot(., aes(name, mean, fill=Diet, ymin=mean-sem, ymax=mean+sem))+
  geom_bar(stat="identity",position=dodge, width = 0.6, color="black")+
  geom_errorbar(position=dodge, width=0.2)+
  coord_flip()+
  facet_grid(cols = vars(Month))+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)), breaks=c(0,0.3,0.6))+
  scale_fill_manual(values=c("#666666", "#99CC00"))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))+
  labs(y="Normalized Score")+
  ggtitle("Males")
p2 = data_summary %>%
  dplyr::filter(Sex == "F") %>%
  ggplot(., aes(name, mean, fill=Diet, ymin=mean-sem, ymax=mean+sem))+
  geom_bar(stat="identity",position=dodge, width = 0.6, color="black")+
  geom_errorbar(position=dodge, width=0.2)+
  coord_flip()+
  facet_grid(cols = vars(Month))+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)), breaks=c(0,0.3,0.6))+
  scale_fill_manual(values=c("#666666", "#99CC00"))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))+
  labs(y="Normalized Score")+
  ggtitle("Females")

data_summary2 = data_summary %>%
  dplyr::filter(name == "Total")
p3 = data_summary2 %>%
  dplyr::filter(Sex == "M") %>%
  mutate(Diet = paste(Diet, "(n)")) %>%
  ggplot(., aes(1, Diet, label=n, color=Diet))+
  geom_text()+
  facet_grid(cols = vars(Month))+
  scale_color_manual(values=c("#666666", "#99CC00"))+
  theme_void()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_text())
p4 = data_summary2 %>%
  dplyr::filter(Sex == "F") %>%
  mutate(Diet = paste(Diet, "(n)")) %>%
  ggplot(., aes(1, Diet, label=n, color=Diet))+
  geom_text()+
  facet_grid(cols = vars(Month))+
  scale_color_manual(values=c("#666666", "#99CC00"))+
  theme_void()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_text())

(p1+p3+p2+p4)+plot_layout(guides="collect",ncol = 1, heights = c(9,1,9,1))
ggsave("./results/frailty_categories.pdf", width = w, height = 1*h, units="mm")

#### STATS ####

summary(lmer(Total~Month+Diet+Sex+(1|AnimalID), data=data))
summary(lmer(Total~Month+Diet+Sex+(1+Month|AnimalID), data=data))
summary(lmer(Total~Month*Diet+Sex+(1|AnimalID), data=data))
summary(lmer(Total~Month*Diet+Sex+(1+Month|AnimalID), data=data))

ggplot(data, aes(as.factor(Month), Total, fill=Diet))+
  geom_boxplot()

# When doing Month*Diet, the effect of DietFuc. is estimated at Month 0
# But because the slope is more shallow for Fuco group, at M0 Fuc. appears
# more frail. I subtract 15 to month so intercept refers to start of treatment
data = data %>%
  mutate(TreatmentMonth = Month - 15, .after=Month)
summary(lmer(Total~TreatmentMonth*Diet+Sex+(1|AnimalID), data=data)) # Oh very cool! No sig difference in intercept at M15!
summary(lm(Total~TreatmentMonth*Diet+Sex, data=data)) # Simple lm gives similar results, turns out

# Should we include random slope or not?
summary(lmer(Total~TreatmentMonth*Diet+Sex+(1|AnimalID), data=data))
summary(lmer(Total~TreatmentMonth*Diet+Sex+(1+TreatmentMonth|AnimalID), data=data))
summary(lmer(Total~TreatmentMonth*Diet+Sex+(0+TreatmentMonth|AnimalID), data=data))
# All are sig, random slope has higher pval
# Im not convinced by how well this behaves when considering that the trend
# over time is not linear. Lets go for intercept only or maybe report both.

# What about Sex interactions?
summary(lmer(Total~TreatmentMonth*Diet*Sex+(1|AnimalID), data=data))
summary(lmer(Total~TreatmentMonth*Diet*Sex+(1+TreatmentMonth|AnimalID), data=data))
# None are significant, and visually it seems correct
# Fucoidan effect on slope is always sig no matter what
# Sex interactions make for a very complicated model so let's run 2 tests instead

dataM = dplyr::filter(data, Sex == "M")
dataF = dplyr::filter(data, Sex == "F")

summary(lmer(Total~TreatmentMonth*Diet+(1|AnimalID), data=dataM))
summary(lmer(Total~TreatmentMonth*Diet+(1|AnimalID), data=dataF))
summary(lmer(Total~TreatmentMonth*Diet+Sex+(1|AnimalID), data=data))

summary(lmer(Total~TreatmentMonth*Diet+(1+TreatmentMonth|AnimalID), data=dataM))
summary(lmer(Total~TreatmentMonth*Diet+(1+TreatmentMonth|AnimalID), data=dataF))
summary(lmer(Total~TreatmentMonth*Diet+Sex+(1+TreatmentMonth|AnimalID), data=data))

# In conclusion, lets test separately for male and female, with random intercept
# Reporting the effect of fucoidan on the slope.

# Ah! I understand why females are more significant despite the differences 
# being clearly larger for males in the plot: the female trajectory is much more 
# linear and simply fits the model better.
# Therefore we will report both pairwise comparisons and overall with mixed model

sink("./results/frailty_stats.txt")
cat("===================== Test frailty month by month =====================\n")
compare_means(Total~Diet, group.by = c("Month", "Sex"), data=data)
cat("\n===================== Test weight month by month =====================\n")
compare_means(Weight~Diet, group.by = c("Month", "Sex"), data=data)
cat("\n===================== Mixed model: males =====================\n")
summary(lmer(Total~TreatmentMonth*Diet+(1|AnimalID), data=dataM))
cat("\n===================== Mixed model: females =====================\n")
summary(lmer(Total~TreatmentMonth*Diet+(1|AnimalID), data=dataF))
closeAllConnections()

#### ZERO CROSSING ####

# We did not start measuring frailty right away so it would be interesting to
# check if the point the trajectories of Ctrl and Fuco meet is the start of treatment
summary(lmer(Total~Month*Diet+Sex+(1|AnimalID), data=data))

ggplot(data, aes(Month, Total, color=Diet))+
  geom_jitter(width = 0.2, size=0.5)+
  geom_abline(intercept = -0.261615, slope = 0.020173, color="#666666")+
  geom_abline(intercept = -0.261615+0.106097, slope = 0.020173-0.008444, color="#99CC00")+
  geom_hline(yintercept=0, linetype="dashed")+
  coord_cartesian(xlim=c(10,26))
ggsave("./results/zero_crossing.pdf", width = 0.8*w, height = 0.4*h, units="mm")

# Quite close actually!

  

