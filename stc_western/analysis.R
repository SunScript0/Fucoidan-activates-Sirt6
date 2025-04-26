library(tidyverse)
library(readxl)
library(ggpubr)

setwd("/scratch/fmorandi/internal/Ali/github/stc_western")

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### PLOT ####

data = read_excel("./data/Western-Blot-Quantification.xlsx")%>%
  mutate(Treatment = ifelse(grepl("Fuco", Sample), "Fuco", "Ctrl")) %>%
  pivot_wider(names_from=Target, values_from=Intensity) %>%
  mutate(H3K9ac_norm = H3K9ac/ H3) %>%
  dplyr::select(-H3, -H3K9ac) %>%
  pivot_longer(cols=c(H3K9ac_norm), names_to = "Target", values_to = "Value")

data_summary = data %>%
  group_by(Target, Sex, Tissue, Treatment) %>%
  summarize(mean = mean(Value, na.rm=T), sd = sd(Value, na.rm=T), n=n()) %>%
  mutate(sem = sd / sqrt(n))

ggplot(data, aes(Treatment, Value, fill=Treatment))+
  geom_boxplot()+
  facet_grid(rows=vars(Target), cols = vars(Tissue, Sex))+
  stat_compare_means(method="t.test")

dodge = position_dodge(width = 0.8)
jitterdodge = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
ggplot(data_summary)+
  geom_bar(aes(Sex, mean, fill=Treatment), stat="identity", color="black", position=dodge, width=0.7)+
  geom_errorbar(aes(x=Sex, y=mean, ymin=mean-sem, ymax=mean+sem, group=Treatment), position=dodge, width=0.35)+
  geom_jitter(data=data, aes(Sex, Value, color=Treatment), position=jitterdodge)+
  facet_grid(rows=vars(Target), cols = vars(Tissue))+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)))+
  scale_fill_manual(values=c("#666666", "#99CC00"))+
  scale_color_manual(values=c("#000000", "#000000"))+
  theme(axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "bottom")
ggsave("./data/western_quant.pdf", width = 1*w, height=0.35*h, units="mm") 

sink("./data/western_stats.txt")
compare_means(Value ~ Treatment, data, group.by = c("Target","Tissue", "Sex"), method="t.test")
closeAllConnections()
