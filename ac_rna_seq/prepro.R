library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(edgeR)
library(biomaRt)

setwd("/scratch/fmorandi/internal/Ali/github/ac_rna_seq")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
gtf_path = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.gtf"

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### LOAD DATA ####

# Read counts
files = dir(paste0(paths$data, "/05_counts"), pattern="*cntTable", full.names = T)
counts = list()
for (f in files) {
  sname = str_extract(f, "/([^/]+).cntTable", group=1)
  counts[[sname]] = fread(f, col.names = c("Feature", sname))
}
counts = Reduce(function(x,y) merge(x, y, by="Feature"), counts)
counts = column_to_rownames(counts, "Feature")

# Take out gene info
ginfo = data.frame(
  "Geneid" = rownames(counts), 
  "is_re" = grepl(":", rownames(counts)),
  row.names = rownames(counts))

# Read meta
meta = read.table("meta.txt", sep="\t", header=T) # Actual real metadata

# Keep genes exressed in at least 3 samples (More filtering at later stages)
counts = counts[, meta$FileName]
counts = counts[rowSums(counts > 0) > 3, ]
ginfo = ginfo[rownames(counts), ]

# Read core qc
core_qc = read.csv("core_qc.csv")
core_qc$Sample.Name = str_replace_all(core_qc$Sample.Name, "-", "_")
core_qc = core_qc[, c(1,4,6:9)]
colnames(core_qc) = c("SampleName", "Conc", "RQN", "A260_230", "A260_280", "ratio28s_18s")
meta = merge(meta, core_qc, by.x="FileName", by.y="SampleName")
setdiff(core_qc$SampleName, meta$FileName) # These were not sequenced

# Read qc info
qc = list()
files = dir(paste0(paths$data, "/04_mapped"), pattern="*Log.final.out", full.names = T)
for (f in files) {
  sname = str_extract(f, "/([^/]+)_Log.final.out", group=1)
  tmp = read_file(f)
  input_reads = str_extract(tmp, "Number of input reads \\|\\t([[:digit:]]+)", group=1)
  uniquely_mapped = str_extract(tmp, "Uniquely mapped reads number \\|\\t([[:digit:]]+)", group=1)
  multi_mapped = str_extract(tmp, "Number of reads mapped to multiple loci \\|\\t([[:digit:]]+)", group=1)
  unmapped = str_extract_all(tmp, "Number of reads unmapped.*([[:digit:]]+)", simplify = T)
  unmapped = str_extract_all(unmapped, "[[:digit:]]+", simplify = T)
  qc[[sname]] = list(
    input_reads = as.numeric(input_reads),
    uniquely_mapped = as.numeric(uniquely_mapped),
    multi_mapped = as.numeric(multi_mapped),
    unmapped = sum(as.numeric(unmapped))
  )
}
files = dir(paste0(paths$data, "/B_TEcounts_logs"), pattern="*txt", full.names = T)
for (f in files) {
  sname = str_extract(f, "/\\d+-\\d+-\\d+_([^/]+).txt", group=1)
  tmp = read_file(f)
  annotated = str_extract(tmp, "Total annotated reads = (\\d+) ", group=1)
  unannotated = str_extract(tmp, "Total unannotated reads = (\\d+) ", group=1)
  qc[[sname]][["annotated"]] = as.numeric(annotated)
  qc[[sname]][["unannotated"]] = as.numeric(unannotated)
}
rm(unmapped, annotated, input_reads, multi_mapped, unannotated, uniquely_mapped)
qc = as.data.frame(do.call(rbind, qc)) %>%
  mutate_all(as.numeric)
meta = merge(meta, qc, by.x="FileName", by.y=0)

meta = arrange(meta, Timepoint, Treatment, Sex)
counts = counts[, meta$FileName]
all.equal(meta$FileName, colnames(counts))

rm(core_qc, qc)

meta$Sex = as.factor(meta$Sex)
meta$Treatment = as.factor(meta$Treatment)
meta$Timepoint = as.factor(meta$Timepoint)

#### SPLIT ginfo and re info ####

counts_re = counts[ginfo$is_re, ]
counts_ge = counts[!ginfo$is_re, ]

rinfo = ginfo[ginfo$is_re, ]
ginfo = ginfo[!ginfo$is_re, ]

all.equal(ginfo$Geneid, rownames(counts_ge))
all.equal(rinfo$Geneid, rownames(counts_re))

rm(counts)

#### CONVERT GENE IDS ####

# Fetch gene symbols and entrez ids
mart = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror="www")
conv = getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
             filter = "ensembl_gene_id",
             values = ginfo$Geneid, mart = mart)

# Note down symbols
conv2 = conv %>%
  dplyr::filter(!duplicated(ensembl_gene_id)) %>%
  dplyr::filter(!duplicated(mgi_symbol))
rownames(conv2) = conv2$ensembl_gene_id
ginfo$Symbol = conv2[ginfo$Geneid, "mgi_symbol"]

# Note down entrez
conv2 = conv %>%
  dplyr::filter(!duplicated(ensembl_gene_id)) %>%
  dplyr::filter(!duplicated(entrezgene_id))
rownames(conv2) = conv2$ensembl_gene_id
ginfo$Entrez = conv2[ginfo$Geneid, "entrezgene_id"]

# Discard genes with no gene symbol
rownames(ginfo) = ginfo$Geneid
ginfo = ginfo[!is.na(ginfo$Symbol), ]
ginfo = ginfo[!ginfo$Symbol == "", ]
counts_ge = counts_ge[ginfo$Geneid, ]
rownames(counts_ge) = ginfo[rownames(counts_ge), "Symbol"]
rownames(ginfo) = ginfo$Symbol

#### GET HUMAN ORTHOLOGS ####

human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

ortho = getLDS(
  attributes = "mgi_symbol",
  filters = "mgi_symbol", values = ginfo$Symbol,
  mart = mouse, attributesL = "hgnc_symbol", martL = human, uniqueRows = T)
ortho = ortho %>%
  dplyr::filter(!duplicated(MGI.symbol)) %>%
  column_to_rownames("MGI.symbol")
ginfo$SymbolHSA = ortho[ginfo$Symbol, "HGNC.symbol"]

#### GET GENE LENGTHS ####

# I checked and this matched what featureCounts returns
gtf = makeTxDbFromGFF(gtf_path, format="gtf")
transcript_lengths2 = transcriptLengths(gtf)
exons_list = exonsBy(gtf, by="gene")
exonic_lengths = data.frame(
  Length = sum(width(reduce(exons_list))))
ginfo$Length = exonic_lengths[ginfo$Geneid, "Length"]

#### OUTLIER REMOVAL ####

dge = DGEList(counts_ge, samples = meta)
dge = calcNormFactors(dge)
norm = data.frame(t(cpm(dge))) # Not log scale to make outliers stick out

pca = prcomp(norm, scale=T)
pca = cbind(meta, pca$x[, c("PC1", "PC2")])
p1 = ggplot(pca, aes(x=PC1, y=PC2, label=FileName))+
  geom_text(size=2)

outliers = c() # Nothing
meta$PassesQC = !meta$FileName %in% outliers

write.table(meta, paste0(paths$results, "/qc_summary.tsv"), sep="\t")

meta = subset(meta, meta$PassesQC)
counts_ge = counts_ge[, meta$FileName]
counts_re = counts_re[, meta$FileName]
dge = dge[, meta$FileName]
dge = calcNormFactors(dge)
norm = data.frame(t(cpm(dge, log = T)), check.names = F)

pca = prcomp(norm, scale=T)
pca = cbind(meta, pca$x[, c("PC1", "PC2")])
p2 = ggplot(pca, aes(x=PC1, y=PC2, label=FileName))+
  geom_text(size=2)
p1 + p2
ggsave(paste0(paths$results, "/pca_outliers.png"), width=1*w, height=0.4*h, units="mm")

#### CHECKPOINT ####

save(counts_ge, counts_re, norm, ginfo, rinfo, meta, file = paste0(paths$results, "/prepro.Rdata"))
