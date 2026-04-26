library(tidyverse)
library(patchwork)
library(contrast)

setwd("/scratch/fmorandi/internal/Ali/merging/ac_cytokines_new")

paths = list()
paths$data = "./data"
paths$results = "./results"
dir.create(paste0(paths$results), showWarnings = F)
dir.create(paste0(paths$results, "/volcanos"), showWarnings = F)
dir.create(paste0(paths$results, "/individual_cytokines"), showWarnings = F)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701
ccodes = c("#cccccc", "#666666", "#ccff66", "#99cc00")

#### LOAD DATA ####

data = read.table(paste0(paths$data, "/data_cytokines.tsv"), header=T, sep="\t") %>%
  mutate(MouseID = str_extract(Sample, "(\\d+-\\d+) 2", group=1), .after=Sample)
meta = read.table("../../github/ac_lifespan/data/survival_data.tsv", header=T, sep="\t", fill=T)
cyto_classes = read.table("./data/cytokine_classes.tsv", header = T, sep="\t")

setdiff(data$MouseID, meta$MouseID)
data = merge(meta, data, by="MouseID")
meta_cols = 1:7
total_event_col = 44

tmp = data
# Some values are below or above dynamic range, change values to ceiling/roof
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    tmp[i, j] = str_replace_all(tmp[i, j], "[<>]", "")
  }
}
# Convert to numeric
tmp =  mutate_if(tmp, is.character, as.numeric)
tmp[, meta_cols] = data[, meta_cols]
data = tmp

# Prettify cytokine names
colnames(data) = str_replace_all(colnames(data), "\\.", "_")
colnames(data) = str_replace_all(colnames(data), "__", "_")
colnames(data) = str_remove(colnames(data), "_$")

#### FILTER LOW ####

hist(colSums(data[, -c(meta_cols, total_event_col)]), breaks = 100) # Range of values
plot(rowSums(data[, -c(meta_cols, total_event_col)]), data[,total_event_col]) # Does the row sum to total_events? No

# Do any columns have zero variance? No
sort(sapply(data[, -c(meta_cols, total_event_col)], sd)) # Least to most variable
data = data[, -which(colnames(data) == "Total_Events")] # Won't use this

#### OUTLIER REMOVAL #####

# Below I look for outliers based on PCA. I won't actually remove them, however:
# The data has very outliery distributions in linear scale and removing them until
# the cloud becomes ~ gaussian would remove too much. Instead, I will work in log scale

pdf(paste0(paths$results, "/outlier_annotation.pdf"), width = w_in, height = 0.5*h_in)
pca = prcomp(data[, -c(meta_cols)], scale = T)
pca = cbind(data[, meta_cols], pca$x)
ggplot(pca, aes(x=PC1, y=PC2, label=MouseID))+
  geom_text()
# outliers = c("8-1", "9-2", "10-3") # Based on unscaled PCA
outliers = c("24-2")
data_noout = subset(data, !data$MouseID %in% outliers)

# Any more outliers? Yes
pca = prcomp(data_noout[, -c(meta_cols)], scale = T)
pca = cbind(data_noout[, c(meta_cols)], pca$x)
ggplot(pca, aes(x=PC1, y=PC2, label=MouseID))+
  geom_text()
outliers = c(outliers, "19-3")
data_noout = subset(data_noout, !data_noout$MouseID %in% outliers)

# Any more outliers? Yes
pca = prcomp(data_noout[, -c(meta_cols)], scale = T)
pca = cbind(data_noout[, c(meta_cols)], pca$x)
ggplot(pca, aes(x=PC1, y=PC2, label=MouseID))+
  geom_text()
outliers = c(outliers, "6-1", "15-4")
data_noout = subset(data_noout, !data_noout$MouseID %in% outliers)

# Any more outliers? Yes, it just keeps going
pca = prcomp(data_noout[, -c(meta_cols)], scale = T)
pca = cbind(data_noout[, c(meta_cols)], pca$x)
ggplot(pca, aes(x=PC1, y=PC2, label=MouseID))+
  geom_text()
dev.off()

#### LOG SCALE ####

data[, -c(meta_cols)] = log1p(data[, -c(meta_cols)])

#### EXPLORATION ####

# PCA
pca = prcomp(data[, -c(meta_cols)], scale = T)
pca = cbind(data[, c(meta_cols)], pca$x)
pdf(paste0(paths$results, "/pca.pdf"), width = w_in, height = 0.5*h_in)
ggplot(pca, aes(x=PC1, y=PC2, color=Sex))+
  geom_point()
ggplot(pca, aes(x=PC1, y=PC2, color=Diet))+
  geom_point()
dev.off()

#### TESTING #####

data_M = subset(data, Sex == "M")
data_F = subset(data, Sex == "F")

vars = colnames(data)[-meta_cols]

##### Linear model #####

resM = data.frame()
for (v in vars) {
  frml = paste(v, "~ Diet")
  frml = as.formula(frml)
  model = lm(frml, data=data_M)
  cont = contrast(model, list(Diet = "Fuco"), list(Diet = "Ctrl"), fcType="log")
  resM[v, c("logFC_FucoVsCtrl", "pval_FucoVsCtrl")] = c(cont$foldChange, cont$Pvalue)
}

resF = data.frame()
for (v in vars) {
  frml = paste(v, "~ Diet")
  frml = as.formula(frml)
  model = lm(frml, data=data_F)
  cont = contrast(model, list(Diet = "Fuco"), list(Diet = "Ctrl"), fcType="log")
  resF[v, c("logFC_FucoVsCtrl", "pval_FucoVsCtrl")] = c(cont$foldChange, cont$Pvalue)
}

# Multiple testing correction
resM = resM %>%
  mutate(padj_FucoVsCtrl = p.adjust(pval_FucoVsCtrl, method="BH"), .after=pval_FucoVsCtrl)
resF = resF %>%
  mutate(padj_FucoVsCtrl = p.adjust(pval_FucoVsCtrl, method="BH"), .after=pval_FucoVsCtrl)

##### U test #####

for (v in vars) {
  frml = paste(v, "~ Diet")
  frml = as.formula(frml)
  test = wilcox.test(frml, data=data_M)
  resM[v, "pvalU_FucoVsCtrl"] = test$p.value
}

for (v in vars) {
  frml = paste(v, "~ Diet")
  frml = as.formula(frml)
  test = wilcox.test(frml, data=data_F)
  resF[v, "pvalU_FucoVsCtrl"] = test$p.value
}

# Multiple testing correction
resM = resM %>%
  mutate(padjU_FucoVsCtrl = p.adjust(pvalU_FucoVsCtrl, method="BH"), .after=pvalU_FucoVsCtrl)
resF = resF %>%
  mutate(padjU_FucoVsCtrl = p.adjust(pvalU_FucoVsCtrl, method="BH"), .after=pvalU_FucoVsCtrl)

##### Compare #####

ggplot(resM, aes(pval_FucoVsCtrl, pvalU_FucoVsCtrl))+
  geom_point()

#### VOLCANOS ####

pval_to_sig = function(logFCs, pvals, th_logc, th_pval) {
  out = rep("Not sig.", length(pvals))
  out[logFCs > th_logc & pvals < th_pval] = "Up"
  out[logFCs < -th_logc & pvals < th_pval] = "Down"
  return(out)
}

my_volcano = function(table, col, v1, v2, title=NULL) {
  maxx = max(abs(range(table[v1])))
  maxy = min(max(-log10(table[[v2]])), 50)
  colmap = c("Down" = "#5555cc", "Not sig." = "#999999", "Up" = "#cc5555")
  tmp2 = table[table[col] != "Not sig.", ]
  p = ggplot(table, aes(x=.data[[v1]], y=-log10(.data[[v2]]), color=.data[[col]]))+
    geom_point(size=0.5)+
    lims(x=c(-maxx, maxx), y=c(0,maxy))+
    ggtitle(title)+
    # geom_text(data=tmp, aes(x=x, y=0.9*maxy, label=n, vjust=0))+
    geom_text_repel(data=tmp2, aes(x=.data[[v1]], y=-log10(.data[[v2]]), label=Cytokine))+
    scale_color_manual(values=colmap)
  if (!is.null(title)) p = p + ggtitle(title)
  return(p)
}

resM %>%
  rownames_to_column("Cytokine")%>%
  mutate(sig_FucoVsCtrl = pval_to_sig(logFC_FucoVsCtrl, pval_FucoVsCtrl, 0, 0.05)) %>%
  my_volcano(., "sig_FucoVsCtrl", "logFC_FucoVsCtrl", "pval_FucoVsCtrl")
ggsave(sprintf("%s/volcanos/FucoVsCtrl_M.png", paths$results), width=0.8*w, height=0.5*h, units="mm")
resF %>%
  rownames_to_column("Cytokine")%>%
  mutate(sig_FucoVsCtrl = pval_to_sig(logFC_FucoVsCtrl, pval_FucoVsCtrl, 0, 0.05)) %>%
  my_volcano(., "sig_FucoVsCtrl", "logFC_FucoVsCtrl", "pval_FucoVsCtrl")
ggsave(sprintf("%s/volcanos/FucoVsCtrl_F.png", paths$results), width=0.8*w, height=0.5*h, units="mm")

#### GROUP COMPARISON HEATMAP ####

heat = data %>%
  group_by(Diet, Sex) %>%
  dplyr::select(-c(MouseID:Death), -c(Note:Sample)) %>%
  summarize_all(mean)
heat[, vars] = apply(heat[, vars], 2, scale)
colMeans(heat[, vars])
varclust = hclust(dist(t(heat[, vars])))
plot(varclust)
varord = vars[varclust$order]

heat %>%
  mutate(group = paste(Diet, Sex, sep = "\n")) %>%
  mutate(group = factor(group, levels=c("Ctrl\nM", "Fuco\nM", "Ctrl\nF", "Fuco\nF"))) %>%
  pivot_longer(cols=all_of(vars)) %>%
  mutate(name = factor(name, levels=varord)) %>%
  ggplot(., aes(x=group, y=name, fill=value)) +
  geom_tile()+
  theme(axis.title=element_blank())+
  scale_fill_gradient2(high="red", low="blue")+
  scale_x_discrete(expand = c(0, 0))+
  labs(fill="Z-score\n(Row-wise)")
ggsave(paste0(paths$results, "/heatmap_clustered.pdf"), width=0.8*w, height=0.5*h, units="mm")

#### MALE VS FEMALE ####

res_comp = merge(resM, resF, by=0, suffix=c("_M", "_F"))
ggplot(res_comp, aes(logFC_FucoVsCtrl_M, logFC_FucoVsCtrl_F))+
  geom_point()+
  stat_cor()

#### INDIVIDUAL FEATURES ####

th_pval = 0.1

resM_sig = resM %>%
  rownames_to_column("Cytokine") %>%
  dplyr::filter(pval_FucoVsCtrl < th_pval)%>%
  mutate(sig = stars.pval(pval_FucoVsCtrl))
for (ctk in unique(resM_sig$Cytokine)) {
  p=ggplot(data_M, aes(Diet, .data[[ctk]], fill=Diet))+
    geom_boxplot()+
    labs(x="Age [years]", color="Treatment")+
    scale_fill_manual(values=ccodes[c(2,4)])+
    stat_compare_means(method="t.test")
  ggsave(plot=p, sprintf("%s/individual_cytokines/%s_M.png", paths$results, ctk), width=0.8*w, height=0.5*h, units="mm")
}

resF_sig = resF %>%
  rownames_to_column("Cytokine") %>%
  dplyr::filter(pval_FucoVsCtrl < th_pval)%>%
  mutate(sig = stars.pval(pval_FucoVsCtrl))
for (ctk in unique(resF_sig$Cytokine)) {
  p=ggplot(data_F, aes(Diet, .data[[ctk]], fill=Diet))+
    geom_boxplot()+
    labs(x="Age [years]", color="Treatment")+
    scale_fill_manual(values=ccodes[c(2,4)])+
    stat_compare_means(method="t.test")
  ggsave(plot=p, sprintf("%s/individual_cytokines/%s_F.png", paths$results, ctk), width=0.8*w, height=0.5*h, units="mm")
}

#### PAPER FIGS ####

ccodes = c("#cccccc", "#666666", "#ccff66", "#99cc00")

resM2 = merge(cyto_classes, resM, by.x="ID", by.y=0)
resF2 = merge(cyto_classes, resF, by.x="ID", by.y=0)

resM2 %>%
  mutate(Name = fct_reorder(Name, logFC_FucoVsCtrl)) %>%
  ggplot(., aes(1, Name, fill=logFC_FucoVsCtrl))+
  geom_tile()+
  facet_grid(rows=vars(Classes), scale="free", space="free")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0(paths$results, "/heatmap_logFC_by_classM.pdf"), width=0.4*w, height=1*h, units="mm")

resF2 %>%
  mutate(Name = fct_reorder(Name, logFC_FucoVsCtrl)) %>%
  ggplot(., aes(1, Name, fill=logFC_FucoVsCtrl))+
  geom_tile()+
  facet_grid(rows=vars(Classes), scale="free", space="free")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0(paths$results, "/heatmap_logFC_by_classF.pdf"), width=0.4*w, height=1*h, units="mm")

#### WRITE TABLES ####

write.csv(resM, paste0(paths$results, "/resultsM.csv"))
write.csv(resF, paste0(paths$results, "/resultsF.csv"))
