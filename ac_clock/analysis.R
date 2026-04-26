library(tidyverse)

setwd("/scratch/fmorandi/internal/Ali/github/ac_clock")

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### LOAD DATA ####

data_lifespan = read.table("../ac_lifespan/data/survival_data.tsv", header = T, sep="\t", fill=T) %>%
  mutate(MouseID = paste0("c", str_replace(MouseID, "-", "m")))

data_clock = read.table("./data/mammalian_array_predictions.tsv", header = T, sep = "\t")

data_clock = merge(data_lifespan, data_clock, by=c("MouseID", "Sex", "Diet")) %>%
  mutate(PredictedAge = UniversalClock3Blood*12)

#### PAPER FIGURES ####

data_clock_F = data_clock %>%
  dplyr::filter(Sex=="F")
data_clock_M = data_clock %>%
  dplyr::filter(Sex=="M")

summary(lm(PredictedAge~TreatmentMonths*Diet, data = data_clock))
summary(lm(PredictedAge~TreatmentMonths*Diet, data = data_clock_F))
summary(lm(PredictedAge~TreatmentMonths*Diet, data = data_clock_M))

ccodes = c("#cccccc", "#666666", "#ccff66", "#99cc00")

ggplot(data_clock, aes(Timepoint, PredictedAge, fill=Diet))+
  geom_boxplot()+
  scale_y_continuous(breaks=c(14, 18, 22, 26, 30))+
  stat_compare_means(method = "wilcox", label="p.format")+
  scale_y_continuous(expand = expansion(mult=c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")
ggsave("./results/preds_final.pdf", width = 0.4*w, height=0.4*h, units="mm")

ggplot(data_clock, aes(Timepoint, PredictedAge, fill=Diet))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(breaks=c(14, 18, 22, 26, 30))+
  stat_compare_means(method = "wilcox", label="p.format")+
  scale_y_continuous(expand = expansion(mult=c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")
ggsave("./results/preds_final_wpoints.pdf", width = 0.4*w, height=0.4*h, units="mm")

ggplot(data_clock, aes(Timepoint, PredictedAge, fill=Diet))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = 0.1))+
  facet_wrap(~Sex)+
  scale_y_continuous(breaks=c(14, 18, 22, 26, 30))+
  stat_compare_means(method = "wilcox", label="p.format")+
  scale_y_continuous(expand = expansion(mult=c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes[c(2,4)])+
  theme(legend.position = "none")
ggsave("./results/preds_final_wpoints_separate.pdf", width = 0.6*w, height=0.4*h, units="mm")
