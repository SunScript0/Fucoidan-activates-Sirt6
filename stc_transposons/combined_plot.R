library(tidyverse)
library(patchwork)

setwd("/scratch/fmorandi/internal/Ali/github/stc_transposons")

dir.create("./results", showWarnings = F)
dir.create("./results/source_data", showWarnings = F)

#### PLOTTING SETTINGS ####

# Letter format
w_in = 8.5
h_in = 11

#### COPY DATA ####

# file.copy(from="../stc_rna_seq/results_merged/tables/resultsDE_Lung_re_joined.csv", to="./results/source_data/rna_lung.csv")
# file.copy(from="../stc_rna_seq/results_merged/tables/resultsDE_Liver_re_joined.csv", to="./results/source_data/rna_liver.csv")
# 
# file.copy(from="../stc_atac_seq/results/combined/da_res_lung.csv", to="./results/source_data/atac_lung.csv")
# file.copy(from="../stc_atac_seq/results/combined/da_res_liver.csv", to="./results/source_data/atac_liver.csv")
# 
# file.copy(from="../stc_medip/results/tables/resultsDE_re1_Lung.csv", to="./results/source_data/medip_lung.csv")
# file.copy(from="../stc_medip/results/tables/resultsDE_re1_Liver.csv", to="./results/source_data/medip_liver.csv")

#### LUNG NO BH ####

res_rna = read.csv("./results/source_data/rna_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_atac = read.csv("./results/source_data/atac_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_medip = read.csv("./results/source_data/medip_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")

all_res = list(
  "rna" = res_rna,
  "atac" = res_atac,
  "medip" = res_medip
)

ps = list()
for (omic in names(all_res)) {
  for (cont in c("FucovsCtrl_M", "FucovsCtrl_F")) {
    tmp = all_res[[omic]]
    tmp$logFC = tmp[, paste0("logFC_", cont)]
    tmp$pval = tmp[, paste0("pval_", cont)]
    tmp$sig = "Not Sig."
    tmp[tmp$logFC > 0 & tmp$pval < 0.05, "sig"] = "Sig. Up"
    tmp[tmp$logFC < 0 & tmp$pval < 0.05, "sig"] = "Sig. Down"
    tmp$sig = factor(tmp$sig, levels=c("Sig. Down", "Not Sig.", "Sig. Up"))
    n_up = sum(tmp$logFC > 0)
    n_down = sum(tmp$logFC < 0)
    # xrange = ifelse(max(abs(tmp$logFC)) < 0.5, 0.5, 1)
    xrange = 1
    xrange = c(-xrange, xrange)
    test = binom.test(sum(tmp$logFC < 0), nrow(tmp))
    p=ggplot(tmp, aes(logFC, -log10(pval), fill=sig))+
      geom_point(pch=21, size=3)+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed", color="grey")+
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey")+
      scale_fill_manual(values=c("blue", "gray", "red"), drop =F)+
      lims(x=xrange, y=c(0, 3))+
      annotate("text", x=xrange, y=c(0,0), label=c(n_down, n_up), color=c("blue", "red"), hjust=c(0, 1), vjust=0)+
      ggtitle(paste(omic, cont), subtitle = sprintf("p-value = %.0e", test$p.value))+
      theme_classic()+
      theme(axis.line = element_blank(),
            panel.border = element_rect())
    ps[[paste(omic, cont, sep="_")]] = p
  }
}
ps$rna_FucovsCtrl_M + ps$atac_FucovsCtrl_M + ps$medip_FucovsCtrl_M + 
ps$rna_FucovsCtrl_F + ps$atac_FucovsCtrl_F + ps$medip_FucovsCtrl_F + plot_layout(nrow = 2, guides="collect")
ggsave("./results/lung.pdf", width = w_in, height = 0.5*h_in)

#### LIVER NO BH ####

res_rna = read.csv("./results/source_data/rna_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_atac = read.csv("./results/source_data/atac_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_medip = read.csv("./results/source_data/medip_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")

all_res = list(
  "rna" = res_rna,
  "atac" = res_atac,
  "medip" = res_medip
)

ps = list()
for (omic in names(all_res)) {
  for (cont in c("FucovsCtrl_M", "FucovsCtrl_F")) {
    tmp = all_res[[omic]]
    tmp$logFC = tmp[, paste0("logFC_", cont)]
    tmp$pval = tmp[, paste0("pval_", cont)]
    tmp$sig = "Not Sig."
    tmp[tmp$logFC > 0 & tmp$pval < 0.05, "sig"] = "Sig. Up"
    tmp[tmp$logFC < 0 & tmp$pval < 0.05, "sig"] = "Sig. Down"
    tmp$sig = factor(tmp$sig, levels=c("Sig. Down", "Not Sig.", "Sig. Up"))
    n_up = sum(tmp$logFC > 0)
    n_down = sum(tmp$logFC < 0)
    xrange = ifelse(max(abs(tmp$logFC)) < 0.5, 0.5, 1)
    xrange = c(-xrange, xrange)
    test = binom.test(sum(tmp$logFC < 0), nrow(tmp))
    p=ggplot(tmp, aes(logFC, -log10(pval), fill=sig))+
      geom_point(pch=21, size=3)+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed", color="grey")+
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey")+
      scale_fill_manual(values=c("blue", "gray", "red"), drop =F)+
      lims(x=xrange, y=c(0, 3))+
      annotate("text", x=xrange, y=c(0,0), label=c(n_down, n_up), color=c("blue", "red"), hjust=c(0, 1), vjust=0)+
      ggtitle(paste(omic, cont), subtitle = sprintf("p-value = %.0e", test$p.value))+
      theme_classic()+
      theme(axis.line = element_blank(),
            panel.border = element_rect())
    ps[[paste(omic, cont, sep="_")]] = p
  }
}
ps$rna_FucovsCtrl_M + ps$atac_FucovsCtrl_M + ps$medip_FucovsCtrl_M + 
  ps$rna_FucovsCtrl_F + ps$atac_FucovsCtrl_F + ps$medip_FucovsCtrl_F + plot_layout(nrow = 2, guides="collect")
ggsave("./results/liver.pdf", width = w_in, height = 0.5*h_in)

#### LUNG BH WITHIN TEs ####

res_rna = read.csv("./results/source_data/rna_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_atac = read.csv("./results/source_data/atac_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_medip = read.csv("./results/source_data/medip_lung.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")

all_res = list(
  "rna" = res_rna,
  "atac" = res_atac,
  "medip" = res_medip
)

ps = list()
for (omic in names(all_res)) {
  for (cont in c("FucovsCtrl_M", "FucovsCtrl_F")) {
    tmp = all_res[[omic]]
    tmp$logFC = tmp[, paste0("logFC_", cont)]
    tmp$pval = tmp[, paste0("pval_", cont)]
    tmp$padj = p.adjust(tmp$pval, method = "BH")
    tmp$sig = "Not Sig."
    tmp[tmp$logFC > 0 & tmp$padj < 0.05, "sig"] = "Sig. Up"
    tmp[tmp$logFC < 0 & tmp$padj < 0.05, "sig"] = "Sig. Down"
    tmp$sig = factor(tmp$sig, levels=c("Sig. Down", "Not Sig.", "Sig. Up"))
    n_up = sum(tmp$logFC > 0)
    n_down = sum(tmp$logFC < 0)
    # xrange = ifelse(max(abs(tmp$logFC)) < 0.5, 0.5, 1)
    xrange = 1
    xrange = c(-xrange, xrange)
    test = binom.test(sum(tmp$logFC < 0), nrow(tmp))
    p=ggplot(tmp, aes(logFC, -log10(pval), fill=sig))+
      geom_point(pch=21, size=3)+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed", color="grey")+
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey")+
      scale_fill_manual(values=c("blue", "gray", "red"), drop =F)+
      lims(x=xrange, y=c(0, 3))+
      annotate("text", x=xrange, y=c(0,0), label=c(n_down, n_up), color=c("blue", "red"), hjust=c(0, 1), vjust=0)+
      ggtitle(paste(omic, cont), subtitle = sprintf("p-value = %.0e", test$p.value))+
      theme_classic()+
      theme(axis.line = element_blank(),
            panel.border = element_rect())
    ps[[paste(omic, cont, sep="_")]] = p
  }
}
ps$rna_FucovsCtrl_M + ps$atac_FucovsCtrl_M + ps$medip_FucovsCtrl_M + 
  ps$rna_FucovsCtrl_F + ps$atac_FucovsCtrl_F + ps$medip_FucovsCtrl_F + plot_layout(nrow = 2, guides="collect")
ggsave("./results/lung_bh_tes.pdf", width = w_in, height = 0.5*h_in)

#### LIVER BH WITHIN TEs ####

res_rna = read.csv("./results/source_data/rna_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_atac = read.csv("./results/source_data/atac_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")
res_medip = read.csv("./results/source_data/medip_liver.csv", row.names = 1) %>%
  dplyr::filter(class == "LINE")

all_res = list(
  "rna" = res_rna,
  "atac" = res_atac,
  "medip" = res_medip
)

ps = list()
for (omic in names(all_res)) {
  for (cont in c("FucovsCtrl_M", "FucovsCtrl_F")) {
    tmp = all_res[[omic]]
    tmp$logFC = tmp[, paste0("logFC_", cont)]
    tmp$pval = tmp[, paste0("pval_", cont)]
    tmp$padj = p.adjust(tmp$pval, method = "BH")
    tmp$sig = "Not Sig."
    tmp[tmp$logFC > 0 & tmp$padj < 0.05, "sig"] = "Sig. Up"
    tmp[tmp$logFC < 0 & tmp$padj < 0.05, "sig"] = "Sig. Down"
    tmp$sig = factor(tmp$sig, levels=c("Sig. Down", "Not Sig.", "Sig. Up"))
    n_up = sum(tmp$logFC > 0)
    n_down = sum(tmp$logFC < 0)
    xrange = ifelse(max(abs(tmp$logFC)) < 0.5, 0.5, 1)
    xrange = c(-xrange, xrange)
    test = binom.test(sum(tmp$logFC < 0), nrow(tmp))
    p=ggplot(tmp, aes(logFC, -log10(pval), fill=sig))+
      geom_point(pch=21, size=3)+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed", color="grey")+
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey")+
      scale_fill_manual(values=c("blue", "gray", "red"), drop =F)+
      lims(x=xrange, y=c(0, 3))+
      annotate("text", x=xrange, y=c(0,0), label=c(n_down, n_up), color=c("blue", "red"), hjust=c(0, 1), vjust=0)+
      ggtitle(paste(omic, cont), subtitle = sprintf("p-value = %.0e", test$p.value))+
      theme_classic()+
      theme(axis.line = element_blank(),
            panel.border = element_rect())
    ps[[paste(omic, cont, sep="_")]] = p
  }
}
ps$rna_FucovsCtrl_M + ps$atac_FucovsCtrl_M + ps$medip_FucovsCtrl_M + 
  ps$rna_FucovsCtrl_F + ps$atac_FucovsCtrl_F + ps$medip_FucovsCtrl_F + plot_layout(nrow = 2, guides="collect")
ggsave("./results/liver_bh_tes.pdf", width = w_in, height = 0.5*h_in)
