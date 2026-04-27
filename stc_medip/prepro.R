library(tidyverse)
library(data.table)
library(biomaRt)
library(edgeR)

setwd("/gpfs/fs2/scratch/fmorandi/internal/Ali/github/stc_medip")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
dir.create("./results", showWarnings = F)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### OUR DATA ####

##### Load data #####

# Read counts
files = dir(paste0(paths$data, "/05_counts"), pattern="*cntTable", full.names = T)
counts = list()
for (f in c(files)) {
  sname = str_extract(f, "/([^/]+).bam.cntTable", group=1)
  counts[[sname]] = fread(f, col.names = c("Feature", sname))
}
counts = Reduce(function(x,y) merge(x, y, by="Feature"), counts)
counts = column_to_rownames(counts, "Feature")

# Take out gene info
ginfo = data.frame(
  "Geneid" = rownames(counts), 
  "is_re" = grepl(":", rownames(counts)),
  row.names = rownames(counts))

# Keep genes exressed in at least 3 samples (More filtering at later stages)
counts = counts[rowSums(counts > 0) > 3, ]
ginfo = ginfo[ginfo$Geneid %in% rownames(counts), ]

# Read meta
meta = read.table("./meta.txt", sep="\t", header=T)

# # Read mapping logs
# logs = dir(paste0(paths$data, "/A_mapping_logs"), full.names=T)
# qc = data.frame()
# for (log in logs) {
#   f = read_file(log)
#   sname = str_extract(log, "\\/[0-9]{4}-[0-9]{2}-[0-9]{2}_([^\\/]+)\\.txt", group=1)
#   # TrimGalore! prints the total number of pairs processed after validation
#   #   This is the number of raw pairs
#   #   Sequences that were trimmed too much are removed
#   s = str_extract(f, "Total number of sequences analysed: ([0-9]+)", group=1)
#   qc[sname, "raw"] = as.numeric(s)
#   # Bowtie2 summary includes the number of pairs processed
#   #   This corresponds to the number of pairs after trimming
#   s = str_extract(f, "([0-9]+) reads; of these:", group=1)
#   qc[sname, "trimmed"] = as.numeric(s)
# }

# Read qc info from counting
qc = list()
files = dir(paste0(paths$data, "/B_TEcounts_logs"), pattern="*txt", full.names = T)
for (f in files) {
  sname = str_extract(f, "/\\d+-\\d+-\\d+_([^/]+).bam.txt", group=1)
  tmp = read_file(f)
  annotated = str_extract(tmp, "Total annotated reads = (\\d+) ", group=1)
  unannotated = str_extract(tmp, "Total unannotated reads = (\\d+) ", group=1)
  qc[[sname]][["mapped"]] = as.numeric(annotated) + as.numeric(unannotated)
  qc[[sname]][["annotated"]] = as.numeric(annotated)
  qc[[sname]][["unannotated"]] = as.numeric(unannotated)
}

# Merge qc and meta
qc = as.data.frame(do.call(rbind, qc)) %>%
  mutate_all(as.numeric)
meta = merge(meta, qc, by.x="SampleID", by.y=0)
rownames(meta) = meta$SampleID
meta = meta[colnames(counts), ]
write.table(meta, paste0(paths$results, "/qc_summary.tsv"), sep="\t", quote=F, row.names=F)

##### Split ginfo and re info ####

all.equal(rownames(counts), ginfo$Geneid)
counts_re = counts[ginfo$is_re, ]
counts_ge = counts[!ginfo$is_re, ]

rinfo = ginfo[ginfo$is_re, ]
ginfo = ginfo[!ginfo$is_re, ]

rm(counts)

##### Convert gene ids #####

# Fetch gene symbols and entrez ids
mart = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://oct2024.archive.ensembl.org")
conv = getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
             filter = "ensembl_gene_id",
             values = ginfo$Geneid, mart = mart)

for (id_type in c("mgi_symbol", "entrezgene_id")) {
  conv2 = conv %>%
    dplyr::filter(!duplicated(ensembl_gene_id)) %>%
    dplyr::filter(!duplicated(id_type))
  rownames(conv2) = conv2$ensembl_gene_id
  ginfo[[id_type]] = conv2[ginfo$Geneid, id_type]
}

# Discard genes with no mgi_symbol
ginfo = drop_na(ginfo)
ginfo = ginfo[!ginfo$mgi_symbol == "", ]
ginfo = ginfo[!duplicated(ginfo$mgi_symbol), ]
counts_ge = counts_ge[ginfo$Geneid, ]
rownames(ginfo) = ginfo$Geneid
rownames(counts_ge) = ginfo[rownames(counts_ge), "mgi_symbol"]
rownames(ginfo) = ginfo$mgi_symbol

##### Split repeat names #####

tmp = str_split_fixed(rinfo$Geneid, ":", n = 3)
rinfo$class = tmp[, 3]
rinfo$superf = paste(tmp[, 3], tmp[, 2], sep="/")
rinfo$fam = tmp[, 1]

##### SAVE #####

save(counts_ge, counts_re, ginfo, rinfo, meta,
     file=paste0(paths$results, "/prepro.Rdata"))

