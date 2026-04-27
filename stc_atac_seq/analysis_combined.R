library(tidyverse)
library(data.table)
library(edgeR)
library(rasterpdf)
library(patchwork)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(BiocParallel)

setwd("~/Desktop/fmorandi/internal/Ali/merging/stc_atac_seq_new")

paths = list()
paths$data = "./results"
paths$results = "./results/combined"
dir.create(paths$results, showWarnings = F)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### FUNCTIONS ####

my_volcano = function(table, logFC_col, pval_col, sig_col, clip_axes=T) {
  if (clip_axes) {
    xmax = quantile(abs(table[, logFC_col]), 0.999)
    ymax = quantile(-log10(table[, pval_col]), 0.999)
  } else {
    xmax = max(abs(table[, logFC_col]))
    ymax = max(-log10(table[, pval_col]))
  }
  tmp = table %>%
    group_by_at(sig_col, .drop=F) %>%
    dplyr::summarize(n=n()) %>%
    mutate(x = c(-xmax*0.9, 0, xmax*0.9)) %>%
    mutate(hjust = c(0, 0.5, 1))
  p = ggplot(table, aes(x=.data[[logFC_col]], y=-log10(.data[[pval_col]]), color=.data[[sig_col]]))+
    geom_point(size=0.1)+
    lims(x=c(-xmax, xmax), y=c(0,ymax))+
    geom_text(data=tmp, aes(x=x, y=ymax*0.95, label=n))+
    scale_color_manual(values=c("#5555cc", "#999999", "#cc5555"), drop=F)+
    labs(x="logFC", y="-log10(pval)", color="Significance")
  return(p)
}

#### LOAD DATA ####

load(paste0(paths$data, "/prepro_combined.Rdata"))
meta$Replicate = as.factor(meta$Replicate)

meta = read.table("./meta.txt", sep="\t", header = T)
counts_ocr = counts_ocr[, meta$SampleID]
counts_re = counts_re[, meta$SampleID]
norm_ocr = norm_ocr[, meta$SampleID]
norm_re = norm_re[, meta$SampleID]

#### PCA ####

##### Ocrs ####

pca = prcomp(t(norm_ocr), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Treatment", "Replicate", "Tissue", "Batch")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(paste0(paths$results, "/pca_ocrs.pdf"), width=1*w_in, height=0.5*h_in)
ps
dev.off()

##### Repeats ####

tmp = norm_re[apply(norm_re, 1, sd) > 0, ]
pca = prcomp(t(tmp), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Treatment", "Replicate", "Tissue", "Batch")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(paste0(paths$results, "/pca_res.pdf"), width=1*w_in, height=0.5*h_in)
ps
dev.off()

#### CHECK SOME MARKER GENES ####

all.equal(meta$SampleID, colnames(norm_ocr))

pinfo$Symbol = mapIds(org.Mm.eg.db, keys=pinfo$geneId, keytype="ENSEMBL", column = "SYMBOL")

ps = list()

# Liver
subset(pinfo, Symbol == "Alb")
Alb = "p215321" # p215321, p215322
ps[["liver"]] = meta %>%
  mutate(marker = unlist(t(norm_ocr[Alb, ]))) %>%
  ggplot(., aes(x=Tissue, y=marker, fill=Treatment))+
  geom_boxplot()+
  ggtitle("Liver marker: Alb", subtitle = "(Albumin)")+
  labs(y="Promoter accessibility")

# Lung
subset(pinfo, geneId == "Sftpd")
SftpdP = "p86428" # p86428, p86429
meta$marker = t(norm_ocr[SftpdP, ])
ps[["lung"]] = ggplot(meta, aes(x=Tissue, y=marker, fill=Treatment))+
  geom_boxplot()+
  ggtitle("Lung marker: Sftpd",
          subtitle="(Surfactant associated protein D)")+
  labs(y="Promoter accessibility")

wrap_plots(ps, nrow=1)+plot_layout(guides="collect")
ggsave(paste0(paths$results, "/markers_tissue.png"), width=1*w, height=0.4*h, units="mm")

#### CHECK SEXES ####

all.equal(colnames(norm_ocr), meta$SampleID)

peaksY = pinfo[pinfo$chr == "Y", ]
accY = colMeans(norm_ocr[rownames(peaksY), ])

meta %>%
  mutate(AccY = accY) %>%
  ggplot(., aes(Sex, AccY, fill=Sex))+
  geom_boxplot()+
  facet_wrap(~Tissue)

sex_markers = c("ENSMUSG00000086503", "ENSMUSG00000068457", "ENSMUSG00000056673") # Xist, Uty, Kdm5d
sex_genes = c("Xist", "Uty", "Kdm5d")
peaks_sex = pinfo[pinfo$geneId %in% sex_markers, ] %>%
  dplyr::filter(annotation == "Promoter") %>%
  dplyr::filter(distanceToTSS == 0)
tmp = data.frame(t(norm_ocr[rownames(peaks_sex), ]))
colnames(tmp) = sex_genes
rownames(meta) = meta$SampleID
tmp$Sex = meta[rownames(tmp), "Sex"]
tmp$Tissue = meta[rownames(tmp), "Tissue"]
tmp %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols=-c(Sample, Sex, Tissue), values_to = "Accessibility") %>%
  ggplot(., aes(x=Sex, y=Accessibility, fill=Sex))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1)+
  facet_grid(rows=vars(Tissue), cols=vars(name))
ggsave(paste0(paths$results, "/markers_sex.png"), width=1*w, height=0.4*h, units="mm")
  
#### REMOVE SAMPLES PROCESSED IN DIFFERENT BATCH ####

meta = subset(meta, Replicate == 1)
counts_ocr = counts_ocr[, meta$SampleID]
counts_re = counts_re[, meta$SampleID]
norm_re = norm_re[, meta$SampleID]

#### FILTER REPEATS ####

ggplot(reinfo, aes(class, meanAcc))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0(paths$results, "/TE_acc_by_class.png"), width=1*w, height=0.6*h, units="mm")

reinfo = reinfo %>%
  dplyr::filter(class %in% c("LINE", "SINE", "LTR", "DNA")) %>%
  dplyr::filter(count > 5, length > 200)
table(reinfo$class)
counts_re = counts_re[rownames(reinfo), ]
norm_re = norm_re[rownames(reinfo), ]

#### DIFFERENTIAL EXPRESSION OF REPEATS ####

##### Main #####

# Note: sexes are in different batches -> Sex effect cant be separated from batch effect

all.equal(colnames(counts_re), meta$SampleID)

meta$Group = paste(meta$Treatment, meta$Sex, sep="_")

# dge = DGEList(counts_re, samples = meta, lib.size = meta$cut_sites) # <!> Originally lib.sizes not specified manually , lib.size = meta$annotated
dge = DGEList(rbind(counts_ocr, counts_re), samples = meta)

all_results_re = list()
for (tissue in unique(meta$Tissue)) {
  dge_sset = dge[, dge$samples$Tissue == tissue]
  design = model.matrix(~0+Group, data=dge_sset$samples)
  contrasts = makeContrasts(
    FucovsCtrl_F=GroupFuco_F-GroupCtrl_F, 
    FucovsCtrl_M=GroupFuco_M-GroupCtrl_M,
    levels = design)
  # Filter further
  keep = filterByExpr(dge_sset, design=design)
  # perc_kept = 100 * colSums(counts_re[keep, colnames(dge_sset)]) / colSums(counts_re[, colnames(dge_sset)])
  # hist(perc_kept, breaks=100)
  dge_sset = dge_sset[keep, ]
  dim(dge_sset)
  
  # TMM normalization
  dge_sset = calcNormFactors(dge_sset)
  ggplot(dge_sset$samples, aes(x=norm.factors))+
    geom_histogram()
  
  dge_sset = estimateDisp(dge_sset, design)
  fit = glmFit(dge_sset, design)
  
  resultsDE_re = list()
  for (cont in colnames(contrasts)) {
    print(cont)
    resultsDE_re[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
      dplyr::filter(rownames(.) %in% rownames(reinfo)) %>%
      dplyr::select(-logCPM, -LR) %>%
      dplyr::rename(pval = PValue) %>%
      mutate(padj = p.adjust(pval, method="BH")) %>%
      rename_with(function(x) paste0(x, "_", cont)) %>%
      rownames_to_column("repID") 
      
  }
  resultsDE_re = Reduce(function(x, y) merge(x, y, by="repID"), resultsDE_re)
  resultsDE_re = merge(reinfo, resultsDE_re, by.x=0, by.y="repID") %>%
    dplyr::rename("repID"="Row.names")
  all_results_re[[tissue]] = resultsDE_re
}

##### Save #####

for (tiss in names(all_results_re)) {
  write.csv(all_results_re[[tiss]], file=paste0(paths$results, "/da_res_", tolower(tiss), ".csv"))
}

#### DEFINE SIGNIFICANCE ####

max_padj = 0.05
min_logFC = 0
for (tissue in names(all_results_re)) {
  for (cont in colnames(contrasts)) {
    logFC = all_results_re[[tissue]][paste0("logFC_", cont)]
    padj = all_results_re[[tissue]][paste0("padj_", cont)]
    sig = abs(logFC) > min_logFC & padj < max_padj
    all_results_re[[tissue]][paste0("sig_", cont)] = "Not Sig."
    all_results_re[[tissue]][sig & logFC > 0, paste0("sig_", cont)] = "Sig. Up"
    all_results_re[[tissue]][sig & logFC < 0, paste0("sig_", cont)] = "Sig. Down"
    all_results_re[[tissue]][paste0("sig_", cont)] = factor(
      all_results_re[[tissue]][, paste0("sig_", cont)],
      levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
    all_results_re[[tissue]] = relocate(all_results_re[[tissue]], paste0("sig_", cont), .after=paste0("padj_", cont))
  }
}

#### VOLCANOS ####

pdf(paste0(paths$results, "/TE_volcanos_joined.pdf"), width=2*w_in, height=0.5*h_in)
for (tissue in unique(meta$Tissue)) {
  for (cont in colnames(contrasts)) {
    ps = list()
    # bps = list()
    for (re_class in c("LINE", "SINE", "LTR", "DNA")) {
      tmp = all_results_re[[tissue]] %>%
        dplyr::filter(class == re_class)
      pval = round(wilcox.test(tmp[, paste0("logFC_", cont)])$p.value, 3)
      pval = ifelse(pval == 0, "<0.001", paste0("=", pval))
      median_logFC = median(tmp[, paste0("logFC_", cont)])
      ps[[re_class]] = my_volcano(
        tmp, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
        ggtitle(re_class, subtitle=sprintf("Median logFC=%f, p%s", median_logFC, pval))
    }
    p=(wrap_plots(ps, nrow=1)&labs(color="Significance\npadj < 0.05"))+ # , |logFC| > 1"
      plot_layout(guides="collect")+
      plot_annotation(title=paste(cont, tissue))
    print(p)
  }
}
dev.off()

#### LINE LENGTH EFFECT ####

all_resultsL1 = list()
for (tissue in names(all_results_re)) {
  all_resultsL1[[tissue]] = all_results_re[[tissue]] %>%
    dplyr::select(class:meanLen, starts_with("logFC")) %>%
    mutate(tissue = tissue) %>%
    dplyr::rename("F" = "logFC_FucovsCtrl_F", "M" = "logFC_FucovsCtrl_M")
}
all_resultsL1 = Reduce(function(x, y) rbind(x,y), all_resultsL1) %>%
  dplyr::filter(superf == "LINE/L1") %>%
  pivot_longer(cols = c("F", "M"), names_to = "sex", values_to = "logFC") %>%
  # merge(re_stats, ., by=c("class", "superf", "fam")) %>%
  mutate(young = grepl("L1Md", fam)) %>%
  mutate(fam=fct_reorder(fam, meanLen))

all_resultsL1 %>%
  dplyr::filter(count > 40) %>%
  ggplot(., aes(x=fam, y=logFC, fill=young))+
  geom_bar(stat="identity")+
  facet_grid(rows = vars(tissue), cols=vars(sex), scale="free_y")+
  theme(axis.text.x = element_text(angle=90, hjust=1, size=7))+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank())
ggsave(paste0(paths$results, "/LINE_length_effect_joined.png"), width=2*w, height=0.5*h, units="mm")
