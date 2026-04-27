library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(BiocParallel)
library(rasterpdf)
library(patchwork)
library(enrichplot)
library(GOSemSim)
library(org.Mm.eg.db)
library(data.table)

setwd("/scratch/fmorandi/internal/Ali/github/stc_medip")

paths = list()
paths$data = "./pipeline_out"
paths$prot_coding_ids = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.protein_coding.ids" 
paths$results = "./results"
paths$tables = paste0(paths$results, "/tables")

dir.create(paths$results, showWarnings = F)
dir.create(paths$tables, showWarnings = F)

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

load(paste0(paths$results, "/prepro.Rdata"))

protein_coding_ids = readLines(paths$prot_coding_ids)
ginfo = ginfo[ginfo$Geneid %in% protein_coding_ids, ]
counts_ge = counts_ge[rownames(ginfo), ]

meta$Treatment = factor(meta$Treatment, levels=c("Ctrl", "Fuco"))
meta$Group = paste(meta$Treatment, meta$Sex, sep="_")

table(meta$Group, meta$Tissue)

#### PCA ####

all.equal(colnames(counts_ge), meta$SampleID)
dge = DGEList(counts_ge, samples = meta)

##### Global #####

norm_ge = cpm(dge, log = T)
pca = prcomp(t(norm_ge), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Tissue", "Treatment", "Sex", "Batch")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(sprintf("%s/pca.pdf", paths$results), width=0.65*w_in, height=0.35*h_in)
print(ps)
dev.off()

##### Within tissue #####

for (tissue in c("Lung", "Liver")) {
  this_dge = dge[, dge$samples$Tissue == tissue]
  this_dge = calcNormFactors(this_dge)
  this_norm = cpm(this_dge, log = T)
  
  pca = prcomp(t(this_norm), center=T, scale.=T)
  pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
  ps = list()
  for (var in c("Treatment", "Sex", "Batch")) { 
    ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
      geom_point()
  }
  
  ps = align_patches(ps)
  pdf(sprintf("%s/pca_%s.pdf", paths$results, tissue), width=0.65*w_in, height=0.35*h_in)
  print(ps)
  dev.off()
}

#### DIFFERENTIAL EXPRESSION ####

##### Genes #####

all.equal(colnames(counts_ge), meta$SampleID)
dge = DGEList(counts_ge, samples = meta)

all_results = list()
for (tissue in unique(meta$Tissue)) {
  dge_sset = dge[, dge$samples$Tissue == tissue]
  design = model.matrix(~0+Group, data=dge_sset$samples)
  contrasts = makeContrasts(
    FucovsCtrl_F=GroupFuco_F-GroupCtrl_F, 
    FucovsCtrl_M=GroupFuco_M-GroupCtrl_M,
    levels = design)

  # TMM normalization
  dge_sset = calcNormFactors(dge_sset)
  ggplot(dge_sset$samples, aes(x=norm.factors))+
    geom_histogram()
  
  dge_sset = estimateDisp(dge_sset, design)
  fit = glmFit(dge_sset, design)
  
  resultsDE = list()
  for (cont in colnames(contrasts)) {
    print(cont)
    resultsDE[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
      dplyr::select(-logCPM, -LR) %>%
      dplyr::rename(pval = PValue) %>%
      mutate(padj = p.adjust(pval, method="BH")) %>%
      rename_with(function(x) paste0(x, "_", cont)) %>%
      rownames_to_column("geneID")
  }
  resultsDE = Reduce(function(x, y) merge(x, y, by="geneID"), resultsDE)
  resultsDE = merge(ginfo, resultsDE, by.x=0, by.y="geneID") %>%
    dplyr::rename("geneID"="Row.names")
  all_results[[tissue]] = resultsDE
}

View(all_results$Liver)
  
##### Repeats #####

counts = rbind(counts_ge, counts_re)
all.equal(colnames(counts), meta$SampleID)
dge = DGEList(counts, samples = meta)

all_results_re1 = list()
for (tissue in unique(meta$Tissue)) {
  dge_sset = dge[, dge$samples$Tissue == tissue]
  design = model.matrix(~0+Group, data=dge_sset$samples)
  contrasts = makeContrasts(
    FucovsCtrl_F=GroupFuco_F-GroupCtrl_F, 
    FucovsCtrl_M=GroupFuco_M-GroupCtrl_M,
    levels = design)
  
  # TMM normalization
  dge_sset = calcNormFactors(dge_sset)
  ggplot(dge_sset$samples, aes(x=norm.factors))+
    geom_histogram()
  
  dge_sset = estimateDisp(dge_sset, design)
  fit = glmFit(dge_sset, design)
  
  resultsDE = list()
  for (cont in colnames(contrasts)) {
    print(cont)
    resultsDE[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
      dplyr::select(-logCPM, -LR) %>%
      dplyr::rename(pval = PValue) %>%
      mutate(padj = p.adjust(pval, method="BH")) %>%
      rename_with(function(x) paste0(x, "_", cont)) %>%
      rownames_to_column("repID")
  }
  resultsDE = Reduce(function(x, y) merge(x, y, by="repID"), resultsDE)
  resultsDE = merge(rinfo, resultsDE, by.x=0, by.y="repID") %>%
    dplyr::rename("repID"="Row.names")
  all_results_re1[[tissue]] = resultsDE
}

#### DEFINE SIGNIFICANCE ####

max_padj = 0.05
min_logFC = 0.0
for (tissue in names(all_results)) {
  for (cont in colnames(contrasts)) {
    logFC = all_results[[tissue]][paste0("logFC_", cont)]
    padj = all_results[[tissue]][paste0("padj_", cont)]
    sig = abs(logFC) > min_logFC & padj < max_padj
    all_results[[tissue]][paste0("sig_", cont)] = "Not Sig."
    all_results[[tissue]][sig & logFC > 0, paste0("sig_", cont)] = "Sig. Up"
    all_results[[tissue]][sig & logFC < 0, paste0("sig_", cont)] = "Sig. Down"
    all_results[[tissue]][paste0("sig_", cont)] = factor(
      all_results[[tissue]][, paste0("sig_", cont)],
      levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
    all_results[[tissue]] = relocate(all_results[[tissue]], paste0("sig_", cont), .after=paste0("padj_", cont))
  }
}

for (tissue in names(all_results)) {
  for (cont in colnames(contrasts)) {
    logFC = all_results_re1[[tissue]][paste0("logFC_", cont)]
    padj = all_results_re1[[tissue]][paste0("padj_", cont)]
    sig = abs(logFC) > min_logFC & padj < max_padj
    all_results_re1[[tissue]][paste0("sig_", cont)] = "Not Sig."
    all_results_re1[[tissue]][sig & logFC > 0, paste0("sig_", cont)] = "Sig. Up"
    all_results_re1[[tissue]][sig & logFC < 0, paste0("sig_", cont)] = "Sig. Down"
    all_results_re1[[tissue]][paste0("sig_", cont)] = factor(
      all_results_re1[[tissue]][, paste0("sig_", cont)],
      levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
    all_results_re1[[tissue]] = relocate(all_results_re1[[tissue]], paste0("sig_", cont), .after=paste0("padj_", cont))
  }
}

#### VOLCANOS ####

##### Genes #####

ps = list()
for (cont in colnames(contrasts)) {
  ps[[cont]] = my_volcano(
    all_results$Lung, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
    ggtitle(cont, subtitle=sprintf("padj < %.2f, |logFC| > %.2f", max_padj, min_logFC))
}

ps = align_patches(ps)
raster_pdf(paste0(paths$results, "/volcanos_ge_Lung.pdf"), width=1*w_in, height=0.5*h_in, res=300)
ps
dev.off()

ps = list()
for (cont in colnames(contrasts)) {
  ps[[cont]] = my_volcano(
    all_results$Liver, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
    ggtitle(cont, subtitle=sprintf("padj < %.2f, |logFC| > %.2f", max_padj, min_logFC))
}

ps = align_patches(ps)
raster_pdf(paste0(paths$results, "/volcanos_ge_Liver.pdf"), width=1*w_in, height=0.5*h_in, res=300)
ps
dev.off()

##### Repeats #####

pdf(paste0(paths$results, "/volcanos_re1_Lung.pdf"), width=2*w_in, height=0.5*h_in)
for (cont in colnames(contrasts)) {
  ps = list()
  # bps = list()
  for (re_class in c("LINE", "SINE", "LTR", "DNA")) {
    tmp = all_results_re1$Lung %>%
      dplyr::filter(class == re_class)
    pval = round(wilcox.test(tmp[, paste0("logFC_", cont)])$p.value, 3)
    pval = ifelse(pval == 0, "<0.001", paste0("=", pval))
    median_logFC = median(tmp[, paste0("logFC_", cont)])
    ps[[re_class]] = my_volcano(
      tmp, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
      ggtitle(re_class, subtitle=sprintf("Median logFC=%f, p%s", median_logFC, pval))
  }
  p=(wrap_plots(ps, nrow=1)&labs(color=sprintf("Significance\npadj < %.2f, |logFC| > %.2f", max_padj, min_logFC)))+
    plot_layout(guides="collect")+
    plot_annotation(title=cont)
  print(p)
}
dev.off()

pdf(paste0(paths$results, "/volcanos_re1_Liver.pdf"), width=2*w_in, height=0.5*h_in)
for (cont in colnames(contrasts)) {
  ps = list()
  # bps = list()
  for (re_class in c("LINE", "SINE", "LTR", "DNA")) {
    tmp = all_results_re1$Liver %>%
      dplyr::filter(class == re_class)
    pval = round(wilcox.test(tmp[, paste0("logFC_", cont)])$p.value, 3)
    pval = ifelse(pval == 0, "<0.001", paste0("=", pval))
    median_logFC = median(tmp[, paste0("logFC_", cont)])
    ps[[re_class]] = my_volcano(
      tmp, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
      ggtitle(re_class, subtitle=sprintf("Median logFC=%f, p%s", median_logFC, pval))
  }
  p=(wrap_plots(ps, nrow=1)&labs(color=sprintf("Significance\npadj < %.2f, |logFC| > %.2f", max_padj, min_logFC)))+
    plot_layout(guides="collect")+
    plot_annotation(title=cont)
  print(p)
}
dev.off()

#### SAVE ####

for (tissue in c("Lung", "Liver")) {
  write.csv(all_results[[tissue]], paste0(paths$tables, "/resultsDE_ge_", tissue, ".csv"))
  write.csv(all_results_re1[[tissue]], paste0(paths$tables, "/resultsDE_re1_", tissue, ".csv"))
}

