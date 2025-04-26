library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(edgeR)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(BiocParallel)
library(patchwork)
library(enrichplot)
library(GOSemSim)
library(immunedeconv)
library(GenomicFeatures)

setwd("/scratch/fmorandi/internal/Ali/github/ac_rna_seq")

paths = list()
paths$data = "./pipeline_out"
paths$prepro = "./results/prepro.Rdata"
paths$results = "./results"
paths$lists = "./results/lists"
paths$tables = "./results/tables"
paths$paper_figs = "./results/paper_figures"
paths$gsea = "./results/gsea"
paths$prot_coding_ids = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.protein_coding.ids" 

dir.create(paths$lists)
dir.create(paths$tables)
dir.create(paths$paper_figs)
dir.create(paths$gsea)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### FUNCTIONS ####

my_volcano = function(table, logFC_col, pval_col, sig_col, xmax, ymax) {
  tmp = table %>%
    mutate_at(c(sig_col), factor, levels=c("Down", "Not sig.", "Up"))%>%
    group_by_at(sig_col, .drop=F) %>%
    dplyr::summarize(n=n()) %>%
    mutate(x = c(-xmax*0.85, 0, xmax*0.85)) %>%
    mutate(hjust = c(0, 0.5, 1))
  p = ggplot(table, aes(x=.data[[logFC_col]], y=-log10(.data[[pval_col]]), color=.data[[sig_col]]))+
    geom_point(size=0.1)+
    lims(x=c(-xmax, xmax), y=c(0,ymax))+
    geom_text(data=tmp, aes(x=x, y=ymax*0.95, label=n))+
    scale_color_manual(values=c("#5555cc", "#999999", "#cc5555"), drop=F)+
    labs(x="logFC", y="-log10(pval)", color="Significance")
  return(p)
}

my_gseGO = function(table, logFC_column, orgdb) {
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$Symbol
  gene_list = sort(gene_list, decreasing = TRUE)
  # Add a tiny amount to break ties
  # gene_list[duplicated(gene_list)] = gene_list[duplicated(gene_list)] - 1e-5
  # Run GSEA
  set.seed(1337)
  res = gseGO(
    geneList=gene_list,
    ont ="BP",
    keyType = "SYMBOL",
    verbose = TRUE,
    OrgDb = orgdb,
    pvalueCutoff = 1.1,
    BPPARAM = SerialParam())
  return(res)
}

my_gseGOcustom = function(table, logFC_column, gmt) {
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$Symbol
  gene_list = sort(gene_list, decreasing = TRUE)
  # Add a tiny amount to break ties
  # gene_list[duplicated(gene_list)] = gene_list[duplicated(gene_list)] - 1e-5
  # Run GSEA
  set.seed(1337)
  res = GSEA(
    geneList=gene_list,
    TERM2GENE = gmt,
    verbose = TRUE,
    pvalueCutoff = 1.1,
    maxGSSize = 1400,
    BPPARAM = SerialParam())
  return(res)
}

my_gseReactome = function(table, logFC_column, org, ginfo) {
  options(warn = 1)
  # Make gene list
  table = merge(ginfo[c("Symbol", "Entrez")], table, by="Symbol") %>%
    dplyr::filter(!is.na(Entrez))
  gene_list = table[[logFC_column]]
  names(gene_list) = table$Entrez
  gene_list = sort(gene_list, decreasing = TRUE)
  # Add a tiny amount to break ties
  # gene_list[duplicated(gene_list)] = gene_list[duplicated(gene_list)] - 1e-5
  # Run GSEA
  set.seed(1337)
  res = gsePathway(
    geneList=gene_list,
    verbose = TRUE,
    organism = org,
    pvalueCutoff = 1.1,
    BPPARAM = SerialParam())
  return(res)
}

my_GO = function(fg, bg, orgdb, ontology="BP") {
  options(warn = 1)
  # Run GSEA
  set.seed(1337)
  res = enrichGO(
    gene = fg,
    universe = bg,
    OrgDb = orgdb,
    keyType ="SYMBOL",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1.1,
    qvalueCutoff = 1.1,
    readable = TRUE)
  return(res)
}

find_by_keyword = function(res, max_p, keywords) {
  res = data.frame(res)
  if (is.null(res) | nrow(res) == 0) return(NULL)
  res = dplyr::filter(res, p.adjust < max_p)
  out = data.frame()
  for (kw in keywords) {
    tmp = res %>%
      dplyr::filter(grepl(kw, Description, ignore.case = T)) %>%
      mutate(keyword = kw) %>%
      dplyr::select(keyword, Description, p.adjust, any_of(c("Direction", "NES"))) #, "geneID", "core_enrichment"
    out = rbind(out, tmp)
  }
  return(out)
}

#### LOAD DATA ####

load(paths$prepro)

table(meta$Treatment, meta$Sex)

protein_coding_ids = readLines(paths$prot_coding_ids)

ginfo = ginfo[ginfo$Geneid %in% protein_coding_ids, ]
counts_ge = counts_ge[ginfo$Symbol, ]
norm = norm[, ginfo$Symbol]

#### EXPLORATION ####

all.equal(rownames(norm), meta$FileName)
pca = prcomp(norm, scale=T)
pca = cbind(meta, pca$x)

# PCA, all vars on one
p = pca %>%
  mutate(Treatment = ifelse(Treatment == "F", "Fuco", "Ctrl")) %>%
  mutate(Group = paste(Treatment, Sex, sep="_")) %>%
  ggplot(., aes(x=PC1, y=PC2, color=Group))+
  geom_point()+
  stat_ellipse(level = 0.67)+
  # facet_wrap(~Sex)+
  scale_color_manual(values=c("#fdbf6f", "#ff7f00", "#a6cee3", "#1f78b4"))
ggsave(paste0(paths$results, "/pca.pdf"), plot = p, width=1*w, height=0.5*h, units="mm")

table(meta$Treatment, meta$Sex)

#### DECONVOLUTION ####

##### Prep #####

all.equal(rownames(counts_ge), ginfo$Symbol)
all.equal(colnames(counts_ge), meta$FileName)
tpm = counts_ge / ginfo$Length
tpm = t(t(tpm) * 1e6 / colSums(tpm))

ortho = ginfo %>%
  dplyr::filter(!duplicated(Symbol)) %>%
  dplyr::filter(!is.na(SymbolHSA))
tpm_hsa = tpm[ortho$Symbol, ]
rownames(tpm_hsa) = ortho[rownames(tpm_hsa), "SymbolHSA"]

##### Run #####

deconv = list()
# ctypes_oi = c("lymphoid", "B cell", "T cell", "T cell CD8+", "T cell CD4+", "NK cell", 
#               "myeloid", "Basophil", "Neutrophil", "Eosinophil", "Mast cell", "Monocytic lineage", "Macrophage/Monocyte")
ctypes_oi = c("B cell", "T cell", "NK cell", "Neutrophil", "Monocytic lineage", "Macrophage/Monocyte")

# mMCP counter
tmp = immunedeconv::deconvolute_mouse(tpm, "mmcp_counter")%>%
  map_result_to_celltypes(ctypes_oi, "mmcp_counter") %>%
  drop_na()
deconv[["mmcp_counter"]] = merge(meta, t(tmp), by.x="FileName", by.y=0)

# DCQ
tmp = immunedeconv::deconvolute_mouse(tpm, "dcq") %>%
  map_result_to_celltypes(ctypes_oi, "dcq") %>%
  drop_na()
deconv[["dcq"]] = merge(meta, t(tmp), by.x="FileName", by.y=0)

# Quantiseq
tmp = immunedeconv::deconvolute(tpm_hsa, "quantiseq", tumor=F, array=F, scale_mrna=F) %>%
  map_result_to_celltypes(ctypes_oi, "quantiseq") %>%
  drop_na()
deconv[["quantiseq"]] = merge(meta, t(tmp), by.x="FileName", by.y=0)

# EPIC
tmp = immunedeconv::deconvolute(tpm_hsa, "epic", tumor=F, scale_mrna=F)%>%
  map_result_to_celltypes(ctypes_oi, "epic") %>%
  drop_na()
deconv[["epic"]] = merge(meta, t(tmp), by.x="FileName", by.y=0)

##### Males #####

# MMCP counter
cell_vars = setdiff(colnames(deconv$mmcp_counter), colnames(meta))
p1 = deconv$mmcp_counter %>%
  dplyr::filter(Sex == "M") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: mMCP counter")

# DCQ
cell_vars = setdiff(colnames(deconv$dcq), colnames(meta))
p2 = deconv$dcq %>%
  dplyr::filter(Sex == "M") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: DCQ")

# quanTIseq
cell_vars = setdiff(colnames(deconv$quantiseq), colnames(meta))
p3 = deconv$quantiseq %>%
  dplyr::filter(Sex == "M") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: quanTIseq")

# EPIC
cell_vars = setdiff(colnames(deconv$epic), colnames(meta))
p4 = deconv$epic %>%
  dplyr::filter(Sex == "M") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: EPIC")

pdf(paste0(paths$results, "/deconvM.pdf"), width=1*w_in, height=0.7*h_in)
list(p1, p2, p3, p4)
dev.off()

##### Females #####

# MMCP counter
cell_vars = setdiff(colnames(deconv$mmcp_counter), colnames(meta))
p1 = deconv$mmcp_counter %>%
  dplyr::filter(Sex == "F") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: mMCP counter")

# DCQ
cell_vars = setdiff(colnames(deconv$dcq), colnames(meta))
p2 = deconv$dcq %>%
  dplyr::filter(Sex == "F") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: DCQ")

# quanTIseq
cell_vars = setdiff(colnames(deconv$quantiseq), colnames(meta))
p3 = deconv$quantiseq %>%
  dplyr::filter(Sex == "F") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: quanTIseq")

# EPIC
cell_vars = setdiff(colnames(deconv$epic), colnames(meta))
p4 = deconv$epic %>%
  dplyr::filter(Sex == "F") %>%
  mutate(Age2 = ifelse(Age == "1.8YO", "21mo", "31mo")) %>%
  pivot_longer(cols=cell_vars)%>%
  ggplot(., aes(x=Age2, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", ncol=3)+
  geom_boxplot()+
  stat_compare_means()+
  ggtitle("Cell type deconvolution: EPIC")

pdf(paste0(paths$results, "/deconvF.pdf"), width=1*w_in, height=0.7*h_in)
list(p1, p2, p3, p4)
dev.off()

##### Compact for paper #####

ccodes = c("#cccccc", "#666666", "#ccff66", "#99cc00")

tmp = deconv$mmcp_counter %>%
  mutate(Sex = factor(Sex, levels= c("M", "F"))) %>%
  pivot_longer(cols=cell_vars) %>%
  dplyr::filter(name != "NK cell") %>%
  mutate(name = str_replace(name, "Macrophage/Monocyte", "Mϕ, Monocyte"))
tmp %>%
  dplyr::filter(Sex == "M") %>%
  ggplot(., aes(x=Treatment, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", nrow=1)+
  # facet_wrap(~name+Sex, scales = "free_y")+
  geom_boxplot(outlier.size = 0.5)+
  theme(axis.title.x = element_blank())+
  labs(y="Score")+
  stat_compare_means(aes(group = Treatment), label = "p.format", size=2.4)+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.2)))+
  scale_fill_manual(values=ccodes[c(2,4)])+
  guides(fill="none")
ggsave(paste0(paths$paper_figs, "/cell_comp_M.pdf"), width = 1*w, height = 0.25*h, units="mm")
tmp %>%
  dplyr::filter(Sex == "F") %>%
  ggplot(., aes(x=Treatment, y=value, fill=Treatment))+
  facet_wrap(~name, scales = "free_y", nrow=1)+
  # facet_wrap(~name+Sex, scales = "free_y")+
  geom_boxplot(outlier.size = 0.5)+
  theme(axis.title.x = element_blank())+
  labs(y="Score")+
  stat_compare_means(aes(group = Treatment), label = "p.format", size=2.4)+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.2)))+
  scale_fill_manual(values=ccodes[c(2,4)])+
  guides(fill="none")
ggsave(paste0(paths$paper_figs, "/cell_comp_F.pdf"), width = 1*w, height = 0.25*h, units="mm")

#### DIFFERENTIAL EXPRESSION ####

all.equal(meta$FileName, colnames(counts_ge))
frml = as.formula("~0+Treatment")

##### Females #####

all.equal(colnames(counts_ge), meta$FileName)
dge = DGEList(counts_ge, samples = meta)
dge = dge[, dge$samples$Sex == "F"]
design = model.matrix(frml, data=dge$samples)
colnames(design) = c("Ctrl", "Fuco")
contrasts = makeContrasts(
  FucoVsCtrl = Fuco - Ctrl,
  levels=design)

# Filter further
keep = filterByExpr(dge, design=design)
perc_kept = 100 * colSums(dge$counts[keep, ]) / colSums(dge$counts)
hist(perc_kept, breaks=100)
dge = dge[keep,,keep.lib.sizes=FALSE]

# TMM normalization
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=norm.factors))+
  geom_histogram()
# Very large deviations from 1

dge = estimateDisp(dge, design)
fit = glmFit(dge, design)

resultsF = list()
for (cont in colnames(contrasts)) {
  print(cont)
  resultsF[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
    dplyr::select(logFC, PValue) %>%
    rename(pval = PValue) %>%
    mutate(padj = p.adjust(pval, method="BH")) %>%
    rename_with(function(x) paste0(x, "_", cont)) %>%
    rownames_to_column("Symbol")
}
resultsF = Reduce(function(x, y) merge(x, y, by="Symbol"), resultsF)
rownames(resultsF) = resultsF$Symbol

##### Males #####

all.equal(colnames(counts_ge), meta$FileName)
dge = DGEList(counts_ge, samples = meta)
dge = dge[, dge$samples$Sex == "M"]
design = model.matrix(frml, data=dge$samples)
colnames(design) = c("Ctrl", "Fuco")
contrasts = makeContrasts(
  FucoVsCtrl = Fuco - Ctrl,
  levels=design)

# Filter further
keep = filterByExpr(dge, design=design)
perc_kept = 100 * colSums(dge$counts[keep, ]) / colSums(dge$counts)
hist(perc_kept, breaks=100)
dge = dge[keep,,keep.lib.sizes=FALSE]

# TMM normalization
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=norm.factors))+
  geom_histogram()
# Very large deviations from 1

dge = estimateDisp(dge, design)
fit = glmFit(dge, design)

resultsM = list()
for (cont in colnames(contrasts)) {
  print(cont)
  resultsM[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
    dplyr::select(logFC, PValue) %>%
    rename(pval = PValue) %>%
    mutate(padj = p.adjust(pval, method="BH")) %>%
    rename_with(function(x) paste0(x, "_", cont)) %>%
    rownames_to_column("Symbol")
}
resultsM = Reduce(function(x, y) merge(x, y, by="Symbol"), resultsM)
rownames(resultsM) = resultsM$Symbol

##### LogFC comparison #####

common_genes = intersect(resultsM$Symbol, resultsF$Symbol)
rownames(resultsM) = resultsM$Symbol
rownames(resultsF) = resultsF$Symbol
cor(resultsM[common_genes, "logFC_FucoVsCtrl"], resultsF[common_genes, "logFC_FucoVsCtrl"])

##### Define significance #####

th_logFC = 0
th_padj = 0.05

for (cont in colnames(contrasts)) {
  # Males
  resultsM[paste0("up_", cont)] = (resultsM[paste0("padj_", cont)] < th_padj) & (resultsM[paste0("logFC_", cont)] > th_logFC)
  resultsM[paste0("down_", cont)] = (resultsM[paste0("padj_", cont)] < th_padj) & (resultsM[paste0("logFC_", cont)] < -th_logFC)
  resultsM[paste0("sig_", cont)] = "Not sig."
  resultsM[resultsM[, paste0("up_", cont)], paste0("sig_", cont)] = "Up"
  resultsM[resultsM[, paste0("down_", cont)], paste0("sig_", cont)] = "Down"
  levels(resultsM[paste0("sig_", cont)]) = c("Down", "Not sig.", "Up")
  # Females
  resultsF[paste0("up_", cont)] = (resultsF[paste0("padj_", cont)] < th_padj) & (resultsF[paste0("logFC_", cont)] > th_logFC)
  resultsF[paste0("down_", cont)] = (resultsF[paste0("padj_", cont)] < th_padj) & (resultsF[paste0("logFC_", cont)] < -th_logFC)
  resultsF[paste0("sig_", cont)] = "Not sig."
  resultsF[resultsF[, paste0("up_", cont)], paste0("sig_", cont)] = "Up"
  resultsF[resultsF[, paste0("down_", cont)], paste0("sig_", cont)] = "Down"
  levels(resultsF[paste0("sig_", cont)]) = c("Down", "Not sig.", "Up")
}

##### Volcanos #####

# Many more sig for male than female, which makes sense
p1 = my_volcano(resultsM, "logFC_FucoVsCtrl", "pval_FucoVsCtrl", "sig_FucoVsCtrl", 4, 15)+
  ggtitle("22mo males")
p2 = my_volcano(resultsF, "logFC_FucoVsCtrl", "pval_FucoVsCtrl", "sig_FucoVsCtrl", 4, 15)+
  ggtitle("22mo females")

p=p1+p2 & theme(legend.position = "none") & plot_layout(nrow=1) #, design = layout
ggsave(plot = p, filename = paste0(paths$paper_figs, "/volcanos.pdf"), width=1*w, height=0.3*h, units="mm")

##### Gene lists #####

gene_lists = list()

for (cont in colnames(contrasts)) {
  gene_lists[[paste0("up_", cont, "_F")]] = resultsF %>%
    dplyr::filter(!!sym(paste0("padj_", cont)) < 0.05) %>%
    dplyr::filter(!!sym(paste0("logFC_", cont)) > 0) %>%
    pull(Symbol)
  gene_lists[[paste0("down_", cont, "_F")]] = resultsF %>%
    dplyr::filter(!!sym(paste0("padj_", cont)) < 0.05) %>%
    dplyr::filter(!!sym(paste0("logFC_", cont)) < 0) %>%
    pull(Symbol)
  gene_lists[[paste0("up_", cont, "_M")]] = resultsM %>%
    dplyr::filter(!!sym(paste0("padj_", cont)) < 0.05) %>%
    dplyr::filter(!!sym(paste0("logFC_", cont)) > 0) %>%
    pull(Symbol)
  gene_lists[[paste0("down_", cont, "_M")]] = resultsM %>%
    dplyr::filter(!!sym(paste0("padj_", cont)) < 0.05) %>%
    dplyr::filter(!!sym(paste0("logFC_", cont)) < 0) %>%
    pull(Symbol)
}
gene_lists$bgF = rownames(resultsF)
gene_lists$bgM = rownames(resultsM)

for (lst in names(gene_lists)) {
  writeLines(gene_lists[[lst]], paste0(paths$lists, "/", lst, ".txt"))
}

#### CHECKPOINT 1 ####

# save.image(paste0(paths$results, "/checkpoint1.Rdata"))
load(paste0(paths$results, "/checkpoint1.Rdata"))

#### GSEA ####

gmt = read.gmt("../custom_gsets/mouse.GO-Aging-MLS_fixed.gmt")
gmtTMS = gmt %>%
  dplyr::filter(grepl(".tissues.", term))

##### Run #####

# gsea = list()
# 
# # === GO ===
# 
# gsea$GO = list()
# gsea$GO$F = list()
# gsea$GO$M = list()
# 
# for (cont in colnames(contrasts)) {
#   gsea$GO$F[[cont]] = my_gseGO(resultsF, paste0("logFC_", cont), org.Mm.eg.db)
#   gsea$GO$M[[cont]] = my_gseGO(resultsM, paste0("logFC_", cont), org.Mm.eg.db)
# }
# 
# # === Custom list: genes associated with aging (TMS) ===
# 
# gsea$TMS = list()
# gsea$TMS$F = list()
# gsea$TMS$M = list()
# 
# for (cont in colnames(contrasts)) {
#   gsea$TMS$F[[cont]] = my_gseGOcustom(resultsF, paste0("logFC_", cont), gmtTMS)
#   gsea$TMS$M[[cont]] = my_gseGOcustom(resultsM, paste0("logFC_", cont), gmtTMS)
# }
# 
# # === Reactome ===
# 
# gsea$Reactome = list()
# gsea$Reactome$F = list()
# gsea$Reactome$M = list()
# 
# for (cont in colnames(contrasts)) {
#   gsea$Reactome$F[[cont]] = my_gseReactome(resultsF, paste0("logFC_", cont), "mouse", ginfo)
#   gsea$Reactome$M[[cont]] = my_gseReactome(resultsM, paste0("logFC_", cont), "mouse", ginfo)
# }
# 
# save(gsea, file=paste0(paths$results, "/gsea.Rdata"))
load(paste0(paths$results, "/gsea.Rdata"))

##### Simplify #####
# 
# gsea_simple = list()
# gsea_simple$GO = list()
# for (sex in names(gsea$GO)) {
#   gsea_simple$GO[[sex]] = list()
#   for (contr in names(gsea$GO[[sex]])) {
#     gsea_simple$GO[[sex]][[contr]] = gsea$GO[[sex]][[contr]] %>%
#       dplyr::filter(p.adjust < 0.05) %>%
#       simplify()
#   }
# }
# 
# save(gsea, gsea_simple, file=paste0(paths$results, "/gsea.Rdata"))
load(paste0(paths$results, "/gsea.Rdata"))

##### Graphs #####

gobp = GOSemSim::godata(org.Mm.eg.db, ont="BP")

for (sex in names(gsea_simple$GO)) {
  pdf(paste0(paths$gsea, "/GOsimplified_graphs", sex, ".pdf"), width=w_in*1, height=h_in*0.5)
  for (contr in names(gsea_simple$GO[[sex]])) {
    tmp = dplyr::filter(gsea_simple$GO[[sex]][[contr]], NES > 0)
    if (nrow(tmp) > 0) {
      tmp = pairwise_termsim(tmp, method="Wang", semData = gobp)
      p = emapplot(tmp, cex_label_category=0.5, cex_category=0.7, cex_line=0.5, max.overlaps=3)+
        ggtitle(sprintf("Activated pathways, %s", contr))+
        scale_x_continuous(expand = expansion(mult=0.2))
      print(p)
    }
    tmp = dplyr::filter(gsea_simple$GO[[sex]][[contr]], NES < 0)
    if (nrow(tmp) > 0) {
      tmp = pairwise_termsim(tmp, method="Wang", semData = gobp)
      p = emapplot(tmp, cex_label_category=0.5, cex_category=0.7, cex_line=0.5, max.overlaps=3)+
        ggtitle(sprintf("Repressed pathways, %s", contr))+
        scale_x_continuous(expand = expansion(mult=0.2))
      print(p)
    }
  }
  dev.off()
}

##### Dotplots #####

# Convert to data.frames and note the comparison as a column
for (ontology in names(gsea)) {
  for (cont in colnames(contrasts)) {
    gsea[[ontology]]$F[[cont]] = data.frame(gsea[[ontology]]$F[[cont]]) %>%
      mutate(Comp = cont)
    gsea[[ontology]]$M[[cont]] = data.frame(gsea[[ontology]]$M[[cont]]) %>%
      mutate(Comp = cont)
  }
}
for (ontology in names(gsea_simple)) {
  for (cont in colnames(contrasts)) {
    gsea_simple[[ontology]]$F[[cont]] = data.frame(gsea_simple[[ontology]]$F[[cont]])%>%
      mutate(Comp = cont)
    gsea_simple[[ontology]]$M[[cont]] = data.frame(gsea_simple[[ontology]]$M[[cont]])%>%
      mutate(Comp = cont)
  }
}

for (ontology in c("GO", "Reactome")) {
  for (sex in names(gsea[[ontology]])) {
    ps = list()
    for (cont in colnames(contrasts)) {
      tmp = gsea[[ontology]][[sex]][[cont]] %>%
        dplyr::filter(p.adjust < 0.05) %>%
        mutate(Description = fct_reorder(Description, abs(NES))) %>%
        mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
        group_by(Sign) %>%
        slice_max(abs(NES), n=20)
      ps[[cont]] = ggplot(tmp, aes(x=abs(NES), y=Description, color=p.adjust))+
        geom_point(size=5)+
        scale_color_gradient(low="red", high="blue")+
        ggtitle(paste(ontology, "GSEA on", cont, ",", sex, "only"))+
        theme(axis.text.y = element_text(size = 10))
      if (nrow(tmp) > 0) {
        ps[[cont]] =  ps[[cont]] + facet_wrap(~Sign)
      }
    }
    ps = align_patches(ps)
    pdf(paste0(paths$gsea, "/gsea_", ontology, sex, ".pdf"), width=w_in*2, height=h_in)
    print(ps)
    dev.off()
  }
}

for (ontology in names(gsea_simple)) {
  for (sex in names(gsea_simple[[ontology]])) {
    ps = list()
    for (cont in colnames(contrasts)) {
      tmp = gsea_simple[[ontology]][[sex]][[cont]] %>%
        dplyr::filter(p.adjust < 0.05) %>%
        mutate(Description = fct_reorder(Description, abs(NES))) %>%
        mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
        group_by(Sign) %>%
        slice_max(abs(NES), n=20)
      ps[[cont]] = ggplot(tmp, aes(x=abs(NES), y=Description, color=p.adjust))+
        geom_point(size=5)+
        scale_color_gradient(low="red", high="blue")+
        ggtitle(paste(ontology, "GSEA on", cont, ",", sex, "only"))+
        theme(axis.text.y = element_text(size = 10))
      if (nrow(tmp) > 0) {
        ps[[cont]] =  ps[[cont]] + facet_wrap(~Sign)
      }
    }
    ps = align_patches(ps)
    pdf(paste0(paths$gsea, "/gsea_GOsimplified", sex, ".pdf"), width=w_in*2, height=h_in)
    print(ps)
    dev.off()
  }
}

##### Custom terms (TMS) #####

gseaTMS = rbind(
  do.call(rbind, gsea$TMS$F) %>%
    pivot_wider(names_from = Comp, values_from=-c(ID, Description, setSize, Comp))%>%
    mutate(Description = substr(Description, 1, 60)) %>%
    mutate(Description = fct_reorder(Description, NES_FucoVsCtrl))%>%
    mutate(Sex = "F"),
  do.call(rbind, gsea$TMS$M) %>%
    pivot_wider(names_from = Comp, values_from=-c(ID, Description, setSize, Comp))%>%
    mutate(Description = substr(Description, 1, 60)) %>%
    mutate(Description = fct_reorder(Description, NES_FucoVsCtrl))%>%
    mutate(Sex = "M")
)

gseaTMS %>%
  dplyr::filter(grepl("WBC|Marrow", ID)) %>%
  arrange(grepl("young", ID), as.character(Description)) %>%
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  dplyr::select(Description, starts_with("NES"), starts_with("p.adj"), Sex) %>%
  pivot_longer(cols=-c(Description, Sex), names_to = c(".value", "Comp"), names_pattern = "^([^_]+)_(.*)") %>%
  mutate(sig = stars.pval(p.adjust)) %>%
  mutate(Comp = paste(Comp, Sex)) %>%
  ggplot(., aes(x=Comp, y=Description, fill=NES, label=sig))+
  geom_tile()+
  geom_text(color="white")+
  ggtitle("Age-associated genes (TMS)") +
  scale_fill_gradient2(low="blue", high="red")+
  scale_x_discrete(expand=expansion(mult=0))+
  scale_y_discrete(expand=expansion(mult=0))+
  theme(axis.title=element_blank())
ggsave(paste0(paths$gsea, "/gsea_TMS.pdf"), width=0.65*w, height=0.2*h, units="mm")

##### Compact dot matrix - Prep #####

gseaGO_comp = merge(
  do.call(rbind, gsea$GO$F) %>%
    pivot_wider(names_from = Comp, values_from=-c(ID, Description, setSize, Comp))%>%
    mutate(Description = substr(Description, 1, 60)) %>%
    mutate(Description = fct_reorder(Description, NES_FucoVsCtrl)),
  do.call(rbind, gsea$GO$M) %>%
    pivot_wider(names_from = Comp, values_from=-c(ID, Description, setSize, Comp))%>%
    mutate(Description = substr(Description, 1, 60)) %>%
    mutate(Description = fct_reorder(Description, NES_FucoVsCtrl)),
  by="ID", suffix=c("_F", "_M"))

names(gseaGO_comp)

# For GO all 
# Get number of sig
tmp = gseaGO_comp[, c("p.adjust_FucoVsCtrl_F", "p.adjust_FucoVsCtrl_M")]
gseaGO_comp = mutate(gseaGO_comp, nsig = rowSums(tmp < 0.05), .after=Description_F)
# Get number of concordant NES sign
tmp = gseaGO_comp[, c("NES_FucoVsCtrl_F", "NES_FucoVsCtrl_M")]
gseaGO_comp = mutate(gseaGO_comp, signsum = rowSums(sign(tmp)), .after=nsig)
# Get mean pvalue and mean of pvalue rank
tmp = gseaGO_comp[, c("p.adjust_FucoVsCtrl_F", "p.adjust_FucoVsCtrl_M")] 
gseaGO_comp = mutate(gseaGO_comp, meanp = rowMeans(-log10(tmp)), .after=signsum)
gseaGO_comp = mutate(gseaGO_comp, meanprank = rowMeans(sapply(tmp, rank)), .after=meanp)

write.csv(gseaGO_comp, paste0(paths$tables, "/gseaGO_comp.csv"))

##### Compact dot matrix - Systematic filtering ####

table(gseaGO_comp$nsig)
table(gseaGO_comp$signsum)
table(gseaGO_comp$nsig, gseaGO_comp$signsum)

# Get terms significant in at least 1 (2<!>) comparison, NES having same sign in 2 (3<!>) or more
# Actually not even because i do want to see things that are sig in female
gseaGO_comp_filt = gseaGO_comp %>%
  dplyr::filter(nsig >= 1)%>%
  dplyr::filter(abs(signsum) >= 0) # Doesnt filter anything, placeholder

# Sort by mean p value or mean pvalue rank?
# Actually i only care so much about pvalue, once it gets below a certain amount
plot(gseaGO_comp_filt$meanp, gseaGO_comp_filt$meanprank)

# terms_oi1 = gseaGO_comp_filt %>%
#   dplyr::filter(p.adjust_FucoVsCtrl_F < 0.05) # Quite many this time, lets just go overall
# Lots of redundancy in the next, so select more than wanted by top_n and then pick to remove redundancy
terms_oi = gseaGO_comp %>%
  dplyr::filter(abs(signsum) == 2) %>%
  group_by(overallsign = sign(signsum)) %>%
  top_n(meanp, n=10) %>%
  ungroup() %>%
  dplyr::select(-overallsign)
# terms_oi = rbind(terms_oi1, terms_oi2) %>%
#   distinct()

maxchar = 50
gseaGO_comp_filt %>%
  dplyr::filter(ID %in% terms_oi$ID) %>%
  mutate(Description_F = if_else(str_length(Description_F) > maxchar, str_c(str_sub(Description_F, 1, maxchar-3), "..."), as.character(Description_F))) %>%
  mutate(Description_F = fct_reorder(Description_F, -sign(signsum)*-log(meanp))) %>%
  dplyr::select(Description_F, nsig:meanprank, starts_with("NES"), starts_with("p.adjust")) %>%
  dplyr::select(Description_F, nsig:meanprank, contains("FucoVsCtrl")) %>%
  pivot_longer(cols=-c(Description_F:meanprank), names_pattern = "(NES|p.adjust)_FucoVsCtrl_(F|M)", names_to = c(".value", "Sex")) %>%
  mutate(Group = paste("22mo", Sex)) %>%
  mutate(Group = factor(Group, levels=c("22mo M", "22mo F"))) %>%
  mutate(sig = stars.pval(p.adjust)) %>%
  ggplot(., aes(Group, Description_F, fill=NES, size=-log10(p.adjust), label=sig))+
  geom_point(pch=21)+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  guides(alpha="none")+
  theme(axis.title = element_blank())
ggsave(filename = paste0(paths$paper_figs, "/gseaGO_summary_systematic.pdf"), width=0.8*w, height=0.5*h, units="mm")

##### Compact dot matrix - Manual filtering ####

# List of interesting terms, with some redundancies
terms_gsea = readLines(paste0(paths$lists, "/terms_oi_gseaGO.txt"))
# tmp =  dplyr::filter(gseaGO_comp, ID %in% terms_gsea)
terms_gsea = setdiff(terms_gsea, c("GO:0072539", "GO:0045064")) # Redundant T helper terms
terms_gsea = setdiff(terms_gsea, c("GO:0008033", "GO:0006400")) # Redundant RNA mod terms

maxchar = 50
gseaGO_comp %>%
  dplyr::filter(ID %in% terms_gsea) %>%
  mutate(Description_F = if_else(str_length(Description_F) > maxchar, str_c(str_sub(Description_F, 1, maxchar-3), "..."), as.character(Description_F))) %>%
  mutate(Description_F = fct_reorder(Description_F, NES_FucoVsCtrl_M)) %>%
  dplyr::select(Description_F, nsig:meanprank, starts_with("NES"), starts_with("p.adjust")) %>%
  dplyr::select(Description_F, nsig:meanprank, contains("FucoVsCtrl")) %>%
  pivot_longer(cols=-c(Description_F:meanprank), names_pattern = "(NES|p.adjust)_FucoVsCtrl_(F|M)", names_to = c(".value", "Sex")) %>%
  mutate(Group = paste("22mo", Sex)) %>%
  mutate(Group = factor(Group, levels=c("22mo M", "22mo F"))) %>%
  mutate(sig = stars.pval(p.adjust)) %>%
  ggplot(., aes(Group, Description_F, fill=NES, size=-log10(p.adjust), label=sig))+
  geom_point(pch=21)+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  guides(alpha="none")+
  theme(axis.title = element_blank())
ggsave(filename = paste0(paths$paper_figs, "/gseaGO_summary_manual.pdf"), width=0.8*w, height=0.4*h, units="mm")

##### Expanded visualization #####

maxchar = 50

# === Males ===
tmp = gsea$GO$M$FucoVsCtrl %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(ID %in% terms_gsea)%>%
  mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  mutate(Direction = factor(Direction, levels=c("Upregulated", "Downregulated"))) %>%
  mutate(Description = ifelse(nchar(Description) > maxchar, paste0(substr(Description, 1, maxchar), "..."), Description))%>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))
ggplot(tmp, aes(Description, -log10(p.adjust), fill=NES))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  coord_flip()+
  facet_grid(rows=vars(Direction), scales = "free_y", space="free_y")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.text = element_text(size=9),
        axis.title.y = element_blank())
ggsave(paste0(paths$paper_figs, "/gseaGO_bars_M.pdf"), width=1*w, height=0.3*h, units="mm")

# === Females ===
tmp = gsea$GO$F$FucoVsCtrl %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(ID %in% terms_gsea)%>%
  mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  mutate(Direction = factor(Direction, levels=c("Upregulated", "Downregulated"))) %>%
  mutate(Description = ifelse(nchar(Description) > maxchar, paste0(substr(Description, 1, maxchar), "..."), Description))%>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))
ggplot(tmp, aes(Description, -log10(p.adjust), fill=NES))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  coord_flip()+
  facet_grid(rows=vars(Direction), scales = "free_y", space="free_y")+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.text = element_text(size=9),
        axis.title.y = element_blank())
ggsave(paste0(paths$paper_figs, "/gseaGO_bars_F.pdf"), width=1*w, height=0.2*h, units="mm")

#### ORA ####

##### Run #####

go_ora = list()
go_ora$F = list()
go_ora$M = list()
for (cont in c("FucoVsCtrl")) {
  bg = resultsF$Symbol
  go_ora$F[[cont]] = rbind(
    resultsF %>%
      dplyr::filter(get(paste0("padj_", cont)) < 0.05) %>%
      dplyr::filter(get(paste0("logFC_", cont)) > 0) %>%
      pull(Symbol) %>%
      my_GO(., bg, org.Mm.eg.db) %>%
      as.data.frame() %>%
      mutate(Direction = "Up"),
    resultsF %>%
      dplyr::filter(get(paste0("padj_", cont)) < 0.05) %>%
      dplyr::filter(get(paste0("logFC_", cont)) < 0) %>%
      pull(Symbol) %>%
      my_GO(., bg, org.Mm.eg.db) %>%
      as.data.frame()%>%
      mutate(Direction = "Down"))
  bg = resultsM$Symbol
  go_ora$M[[cont]] = rbind(
    resultsM %>%
      dplyr::filter(get(paste0("padj_", cont)) < 0.05) %>%
      dplyr::filter(get(paste0("logFC_", cont)) > 0) %>%
      pull(Symbol) %>%
      my_GO(., bg, org.Mm.eg.db) %>%
      as.data.frame() %>%
      mutate(Direction = "Up"),
    resultsM %>%
      dplyr::filter(get(paste0("padj_", cont)) < 0.05) %>%
      dplyr::filter(get(paste0("logFC_", cont)) < 0) %>%
      pull(Symbol) %>%
      my_GO(., bg, org.Mm.eg.db) %>%
      as.data.frame()%>%
      mutate(Direction = "Down"))
}

go_ora$M$FucoVsCtrl %>%
  dplyr::filter(p.adjust < 0.05) %>%
  View()
# Nothing
# go_ora$logFC1$F$FucoVsCtrl %>% 
#   dplyr::filter(p.adjust < 0.05) %>%
#   View()

##### Save #####

# save(go_ora, file=paste0(paths$results, "/go_ora.Rdata"))
load(paste0(paths$results, "/go_ora.Rdata"))

for (sex in names(go_ora)) {
  tmp = go_ora[[sex]]$FucoVsCtrl
  if (nrow(tmp>0)) {
    tmp %>%
      dplyr::filter(p.adjust < 0.05) %>%
      write.csv(., file = sprintf("%s/oraGO_%s.csv", paths$tables, sex))
  }
}

##### Plot selected terms #####

terms_ora = readLines(paste0(paths$lists, "/terms_oi_oraGO.txt"))
terms_ora = setdiff(terms_ora, c("GO:0050658")) # Redundant RNA terms

tmp = go_ora$M$FucoVsCtrl %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(ID %in% terms_ora)%>%
  mutate(Direction = factor(Direction, levels=c("Up", "Down"))) %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))
ggplot(tmp, aes(Description, -log10(p.adjust), fill=Direction))+
  geom_bar(stat="identity")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  coord_flip()+
  facet_grid(rows=vars(Direction), scales = "free_y", space="free_y")+
  scale_fill_manual(values = c("#dd3311", "#1155bb"))+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  theme(axis.text = element_text(size=9),
        axis.title.y = element_blank())
ggsave(paste0(paths$paper_figs, "/go_bars_M.pdf"), width=1*w, height=0.4*h, units="mm")

#### CHECKPOINT 2 ####

# save.image(paste0(paths$results, "/checkpoint2.Rdata"))
load(paste0(paths$results, "/checkpoint2.Rdata"))

#### WRITE TABLES ####

write.csv(counts_ge, paste0(paths$tables, "/counts_ge.csv"))
write.csv(resultsF, paste0(paths$tables, "/de_resultsF.csv"))
write.csv(resultsM, paste0(paths$tables, "/de_resultsM.csv"))
write.csv(deconv$mmcp_counter, paste0(paths$tables, "/deconv_mmcp_counter.csv"))

for (ont in names(gsea)) {
  for (sex in names(gsea[[ont]])) {
    for (cont in names(gsea[[ont]][[sex]])) {
      tmp = gsea[[ont]][[sex]][[cont]] %>%
        dplyr::filter(p.adjust < 0.05) %>%
        as.data.frame() %>%
        mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed"))
      write.csv(tmp, paste0(paths$tables, "/gsea_", ont, "_", cont, "_", sex, ".csv"))
    }
  }
}

