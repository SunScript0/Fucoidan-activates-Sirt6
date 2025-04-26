library(tidyverse)
library(ggpubr)

setwd("/scratch/fmorandi/internal/Ali/github/biochem")

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### PANEL A ####

data = read_excel("./biochem_data.xlsx", sheet = "Fig1A") %>%
  mutate(Norm1 = Intensity / Coomassie)
data2 = data %>%
  dplyr::select(Replicate, Condition, Intensity, Norm1) %>%
  dplyr::filter(Condition %in% c("H2O", "DMSO")) %>%
  pivot_wider(names_from = Condition, values_from = c(Intensity, Norm1))
data = merge(data, data2, by="Replicate") %>%
  mutate(Norm2 = ifelse(Condition %in% c("Fucoidan", "H2O"), Norm1 / Norm1_H2O, Norm1 / Norm1_DMSO)) %>%
  mutate(Norm3 = ifelse(Condition %in% c("Fucoidan", "H2O"), Intensity / Intensity_H2O, Intensity / Intensity_DMSO))
cor(data[, c("Intensity", "Norm1", "Norm2", "Norm3")])
# Norm1 is Intensity normalized by coomassie
# Norm2 is Norm1 norm by respective control
# Norm3 is Intensity norm by respective control
  
ggbarplot(data, x="Condition", y="Norm3", fill="Condition", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
ggsave("./results/fig1a/plot.pdf", width = 0.5*w_in, height = 0.3*h_in)

ggbarplot(data, x="Condition", y="Norm2", fill="Condition", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
ggsave("./results/fig1a/plot_by_coomassie.pdf", width = 0.5*w_in, height = 0.3*h_in)
# Very similar to other normalization

# For stats I do not normalize by respective control
# That would artificially zero the variance of the control
sink("./results/fig1a/stats.txt")
cat("=== Comparing intensity directly ===\n")
compare_means(Intensity~Condition, data, method = "t.test")%>%
  dplyr::filter(group1 == "H2O" & group2 == "Fucoidan" | group1 == "DMSO" & group2 != "Fucoidan")%>%
  dplyr::select(-p.adj)
cat("\n=== Comparing intensity normalized by coomassie ===\n")
compare_means(Norm1~Condition, data, method = "t.test")%>%
  dplyr::filter(group1 == "H2O" & group2 == "Fucoidan" | group1 == "DMSO" & group2 != "Fucoidan")%>%
  dplyr::select(-p.adj)
closeAllConnections()

#### PANEL B ####
                     
data = read_excel("./biochem_data.xlsx", sheet = "Fig1B") %>%
  mutate(Norm1 = Intensity / Coomassie)
data2 = data %>%
  dplyr::select(Replicate, Time, Condition, Intensity, Norm1) %>%
  dplyr::filter(Time == 0) %>%
  pivot_wider(names_from = Time, values_from = c(Intensity, Norm1))
data = merge(data, data2, by=c("Condition", "Replicate")) %>%
  mutate(Norm2 = Norm1 / Norm1_0)%>%
  mutate(Norm3 = Intensity / Intensity_0)

data2 = data %>%
  group_by(Condition, Time) %>%
  summarize(Norm2mean = mean(Norm2), Norm3mean = mean(Norm3),
            Norm2sd = sd(Norm2), Norm3sd = sd(Norm3), n = n()) %>%
  mutate(Norm2sem = Norm2sd / sqrt(n), Norm3sem = Norm3sd / sqrt(n))

# Here normalizing by coomassie looks weird, so I'll stick to Norm3 for everything
ggplot(data2, aes(Time, Norm2mean, ymin=Norm2mean-Norm2sem, ymax=Norm2mean+Norm2sem, color=Condition))+
  geom_point()+
  geom_line()+
  geom_errorbar(width=2)
ggplot(data2, aes(Time, Norm3mean, ymin=Norm3mean-Norm3sem, ymax=Norm3mean+Norm3sem, color=Condition))+
  geom_point()+
  geom_line()+
  geom_errorbar(width=2)
ggsave("./results/fig1b/plot.pdf", width = 0.6*w_in, height = 0.3*h_in)

ggplot(data2, aes(Time, Norm3mean, ymin=Norm3mean-Norm3sem, ymax=Norm3mean+Norm3sem, color=Condition))+
  geom_point()+
  geom_line()+
  geom_errorbar(width=2)+
  theme(legend.position="none")
ggsave("./results/fig1b/plot_no_legend.pdf", width = 0.5*w_in, height = 0.3*h_in)

#### PANEL C ####

data = read_excel("./biochem_data.xlsx", sheet = "Fig1C") %>%
  mutate(Norm1 = Intensity / Coomassie)
data2 = data %>%
  dplyr::select(Replicate, Condition, Intensity, Norm1) %>%
  dplyr::filter(Condition == "Control") %>%
  pivot_wider(names_from = Condition, values_from = c(Intensity, Norm1))
data = merge(data, data2, by="Replicate") %>%
  mutate(Norm2 = Norm1 / Norm1_Control) %>%
  mutate(Norm3 = Intensity / Intensity_Control)
cor(data[, c("Intensity", "Norm1", "Norm2", "Norm3")])

ggbarplot(data, x="Condition", y="Norm3", fill="Condition", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.2)))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
ggsave("./results/fig1c/plot.pdf", width = 0.5*w_in, height = 0.3*h_in)

sink("./results/fig1c/stats.txt")
cat("=== Comparing intensity directly ===\n")
compare_means(Intensity~Condition, data, method = "t.test")%>%
  dplyr::select(-p.adj) %>%
  print(n=100)
cat("\n=== Comparing intensity normalized by coomassie ===\n")
compare_means(Norm1~Condition, data, method = "t.test")%>%
  dplyr::select(-p.adj)%>%
  print(n=100)
closeAllConnections()

#### PANEL D ####

data = read_excel("./biochem_data.xlsx", sheet = "Fig1D") %>%
  mutate(Norm1 = Intensity / Coomassie)
data2 = data %>%
  dplyr::select(Replicate, Concentration, Intensity, Norm1) %>%
  dplyr::filter(Concentration == 0) %>%
  pivot_wider(names_from = Concentration, values_from = c(Intensity, Norm1))
data = merge(data, data2, by="Replicate") %>%
  mutate(Norm2 = Norm1 / Norm1_0) %>%
  mutate(Norm3 = Intensity / Intensity_0)
cor(data[, c("Intensity", "Norm1", "Norm2", "Norm3")])

data2 = data %>%
  group_by(Concentration) %>%
  summarize(Norm2mean = mean(Norm2), Norm3mean = mean(Norm3),
            Norm2sd = sd(Norm2), Norm3sd = sd(Norm3), n = n()) %>%
  mutate(Norm2sem = Norm2sd / sqrt(n), Norm3sem = Norm3sd / sqrt(n))
# ggplot() + 
#   geom_point(data=data2, aes(Concentration, Norm3mean)) + 
#   geom_errorbar(data=data2, aes(Concentration, Norm3mean, ymin=Norm3mean-Norm3sem,ymax=Norm3mean+Norm3sem)) + 
#   geom_smooth(data=data, aes(Concentration, Norm3), method = "nls", se = FALSE,
#               method.args = list(formula = y ~ Vmax * x / (Km + x),
#                                  start = list(Km = 1, Vmax = 2)))

ggplot() +
  geom_point(data=data2, aes(Concentration, Norm3mean)) +
  geom_errorbar(data=data2, aes(Concentration, Norm3mean, ymin=Norm3mean-Norm3sem,ymax=Norm3mean+Norm3sem)) +
  geom_smooth(data=data, aes(Concentration, Norm3), method = "nls", se = FALSE,
              method.args = list(formula = y ~ (Vmax*x^h)/(Khalf^h + x^h),
                                 start = list(h = 1, Vmax = 4, Khalf=0.1)))
ggsave("./results/fig1d/plot.pdf", width = 0.5*w_in, height = 0.3*h_in)
