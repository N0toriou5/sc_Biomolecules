### Setup working environment
### set your home directory
### load libraries (if library is not available, it will be installed)

if (!requireNamespace("BiocManager",quietly=TRUE)){
  install.packages("BiocManager")
}
packages<-c("AnnotationDbi","biomaRt","corto","DESeq2",
            "dplyr","ggplot2","grid","gridExtra",
            "msigdbr","org.Hs.eg.db","org.Mm.eg.db",
            "pander","Rtsne","Seurat","slalom",
            "stringr","TeachingDemos",
            "tidyverse","xlsx"
)
for(p in packages){
  if (!p %in% rownames(installed.packages())){
    BiocManager::install(p)
  }
  library(p,character.only=TRUE)
}

## Load useful functions to perform the whole analysis
source("code/area.R") # Functions to perform advanced rank enrichment analysis (area) for single cell GSEA
source("code/geneids.R") # Functions to convert gene ids
source("code/heatmaps.R") # Functions to plot Heatmaps
source("code/expmerger.R") # Functions to sum multiple expression tracks
source("code/qol.R") # Various "quality of life" functions
source("code/vst.R") # Variance Stabilizing Transformation functions

# Prepare data (define a color for each cell line)
col_be2c <- "salmon"
col_kelly <- "cornflowerblue"

if(!file.exists("data/rawcounts_symbols.rda")){
  ### Load counts (from CellRanger) ----
  rawcounts_be2c <- read.csv(gzfile("data/BE2C.csv.gz"), as.is =  TRUE, row.names =  1)
  ncol(rawcounts_be2c) # 962 BE2C cells
  rawcounts_kelly <- read.csv(gzfile("data/Kelly.csv.gz"), as.is =  TRUE, row.names =  1)
  ncol(rawcounts_kelly) # 1105 Kelly cells
  # R loads "-" as ".", we must fix this
  colnames(rawcounts_be2c) <- gsub("\\.", "-", colnames(rawcounts_be2c))
  colnames(rawcounts_be2c) <- paste0("be2c_", colnames(rawcounts_be2c))
  colnames(rawcounts_kelly) <- gsub("\\.", "-", colnames(rawcounts_kelly))
  colnames(rawcounts_kelly) <- paste0("kelly_", colnames(rawcounts_kelly))
  # Merge the two datasets in a single raw counts matrix
  # rows are genes, columns are single cells
  rawcounts <- cbind(rawcounts_be2c, rawcounts_kelly)
  rawcounts <- as.matrix(rawcounts)
  dim(rawcounts) # 32738 ENSEMBL genes, 2067 samples
  save(rawcounts, file =  "data/rawcounts_ensg.rda")
  
  ## Convert to gene symbols
  # Ensembl ids mapping to the same gene will see their raw counts summed up
  ensgmat <- rawcounts
  tmp <- ens2eg(rownames(ensgmat))
  convlist <- eg2sym(tmp)
  names(convlist) <- names(tmp)
  rawcounts <- squish(ensgmat, convlist =  convlist, method =  "sum", verbose =  TRUE)
  dim(rawcounts) # 22135 genes, 2067 samples
  save(rawcounts, file =  "data/rawcounts_symbols.rda")
  
  ## Save as CSV
  ksym <- rawcounts[ , grep("kelly", colnames(rawcounts))]
  bsym <- rawcounts[ , grep("be2c", colnames(rawcounts))]
  write.csv(bsym,file =  "results/Supp_File_S1_BE2C_rawcounts.csv")
  write.csv(ksym,file =  "results/Supp_File_S2_Kelly_rawcounts.csv")
} else {load("data/rawcounts_symbols.rda")}

## Normalization
# Extract gene lengths (more precisely, transcript lengths)
fname <- "data/genelengths.rda"
if(!file.exists(fname)){
  library(GenomicFeatures)
  supportedUCSCtables(genome =  "hg38", url =  "http://genome.ucsc.edu/cgi-bin/")
  hg <- makeTxDbFromUCSC(genome =  "hg38", tablename =  "refGene")
  exonic <- exonsBy(hg,by =  "gene")
  redexonic <- reduce(exonic)
  genelengths <- sum(width(redexonic))
  names(genelengths) <- eg2sym(names(genelengths))
  genelengths <- genelengths[!duplicated(names(genelengths))]
  genelengths <- genelengths[genelengths>0]
  genelengths <- genelengths[!is.na(genelengths)]
  save(genelengths,file =  fname)
}else{load(fname)}

# Function to calculate FPKM (https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
fpkm <- function(counts,genelengths){
  common <- intersect(rownames(counts), names(genelengths))
  counts <- counts[common, ]
  lengths <- genelengths[common]
  fpms <- apply(counts,2,function(x){1E6*x/sum(x)})
  fpkms <- fpms
  for(i in 1:nrow(fpms)){
    fpkms[i,] <- 1E3*fpms[i,] / lengths[i]
  }
  return(fpkms)
}
fname <- "data/fpkms.rda"
if(!file.exists(fname)){
  fpkms <- fpkm(rawcounts, genelengths)
  save(fpkms,file =  fname)
}
# Function to calculate TPM (https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
tpm <- function(counts,genelengths){
  common <- intersect(rownames(counts), names(genelengths))
  counts <- counts[common,]
  lengths <- genelengths[common]
  intermediate <- counts
  for(i in 1:nrow(counts)){
    intermediate[i,] <- 1E3*counts[i,] / lengths[i]
  }
  tpms <- apply(intermediate, 2, function(x){1E6*x/sum(x)})
  return(tpms)
}
fname <- "data/tpms.rda"
if(!file.exists(fname)){
  tpms <- tpm(rawcounts, genelengths)
  save(tpms, file =  fname)
}

# Variance-Stabilizing Transformation (expmat)
fname <- "data/tpms.rda"
if(!file.exists(fname)){
  vstmat <- vst(rawcounts)
  save(vstmat, file =  fname)
}

## Log-Normalize with Seurat
# By default, Seurat employs a global-scaling normalization method LogNormalize that
# normalizes the gene expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result:
fname <- "data/expmat.rda"
if(!file.exists(fname)){
  seuset <- CreateSeuratObject(counts =  rawcounts,project =  "nblcells", min.cells =  3, min.features =  1000)
  seuset # 15782 genes, 2067 samples
  seuset <- NormalizeData(seuset, normalization.method =  "LogNormalize", scale.factor =  10000)
  save(seuset, file =  "data/seuset.rda")
  
  ## Mean variability plot showing most expressed and variable genes
  expmat <- as.matrix(seuset[["RNA"]]@data)
  dim(expmat) # 15782 genes, 2067 samples
  save(expmat, file =  fname)
}

### Descriptive plots using TPMs ----
load("data/tpms.rda")
fname <- "results/tpms.csv"
if(!file.exists("results/tpms.csv.gz")){
  if(!file.exists(fname)){
    write.csv(tpms, file =  fname)
  }}
bmat <- tpms[, grep("be2c", colnames(tpms))]
kmat <- tpms[, grep("kelly", colnames(tpms))]


### Variance vs. Expression ----
bexpmean <- apply(bmat, 1, mean)
kexpmean <- apply(kmat, 1, mean)
bexpvar <- apply(bmat, 1, var)
kexpvar <- apply(kmat, 1, var)


# Define pseudovalues for log10 transform
min(kexpmean[kexpmean!=  0]) # 0.0006598291
min(bexpmean[bexpmean!=  0]) # 0.0003172931
pseudo <- 0.0001
min(kexpvar[kexpvar!=  0]) # 0.0003762165
min(bexpvar[bexpvar!=  0]) # 0.00009684925
pseudov <- 0.0001


# Calculate fitted variance
bx <- log10(bexpmean+pseudo)
by <- log10(bexpvar+pseudov)
bres <- loess(formula =  by~bx)$resid
kx <- log10(kexpmean+pseudo)
ky <- log10(kexpvar+pseudov)
kres <- loess(formula =  ky~kx)$resid 



### Genes to show ----
## Housekeeping
toshow <- c("ACTB", "GAPDH", "B2M", "GUSB")
## MYCN network
# Moreover, as existing literature refers to PRDM8, MYBL8, HMGB2 and TEAD4
# as crucial factors in the regulatory network triggered by MYCN, we examined how
# their transcription level varied across the two divergent neuroblastoma cell groups
# as a function of MYCN expression.
toshow <- c(toshow,c("MYCN", "PRDM8", "MYBL2", "HMGB2", "TEAD4"))
## Highly variable
bhighvar <- names(sort(bres, dec =  TRUE))[1:4]
khighvar <- names(sort(kres, dec =  TRUE))[1:4]
toshow <- c(toshow, khighvar, bhighvar)
## Low variance, high expression (candidate housekeeping)
bhk <- intersect(names(sort(bres, dec =  FALSE))[1:500], names(sort(bx, dec =  TRUE))[1:1000])
khk <- intersect(names(sort(kres, dec =  FALSE))[1:200], names(sort(kx, dec =  TRUE))[1:1000])
toshow <- c(toshow, bhk, khk)
## Other genes
toshow <- c(toshow, "MYC", "MYCL", "ALK", "LMO1")
toshow <- unique(toshow)

# Variance vs. Expression plots
png("plots/003_var_vs_exp.png", w =  3000, h =  3000, res =  300)
set.seed(3)
par(mfrow =  c(2, 2))
scatter2(bx, by, main =  "BE2C", xlab =  "Log10 Average Expression (TPM)", ylab =  "Log10 Variance",col =  col_be2c)
scatter2(bx, bres, main =  "BE2C", xlab =  "Log10 Average Expression (TPM)", ylab =  "Log10 Variance (Loess Residuals)",
         col =  col_be2c, showLine =  FALSE, extendXlim =  TRUE)
textplot3(bx[toshow], bres[toshow], words =  toshow,font =  2)

scatter2(kx, ky, main =  "Kelly", xlab =  "Log10 Average Expression (TPM)", ylab =  "Log10 Variance", col =  col_kelly)
scatter2(kx, kres, main =  "Kelly", xlab =  "Log10 Average Expression (TPM)", ylab =  "Log10 Variance (Loess Residuals)",
         col =  col_kelly, showLine =  FALSE, extendXlim =  TRUE)
textplot3(kx[toshow], kres[toshow], words =  toshow, font =  2)
dev.off()

### TPM vs. TPM ----
# Plot and compare genes in common (they should be the same already, but better to be sure)
png("plots/003_tpm_vs_tpm_exp.png", w =  4000,  =  3000, res =  600)
set.seed(6)
common <- intersect(names(bexpmean), names(kexpmean))
corto::scatter(bx, kx, xlab =  "BE2C Log10 Average Expression (TPM)", ylab =  "Kelly Log10 Average Expression (TPM)",
               main =  "Average Expression", col =  "gainsboro", extendXlim =  TRUE)
textplot3(bx[toshow], kx[toshow], words =  toshow, font =  2)
dev.off()

png("plots/003_tpm_vs_tpm_var.png", w =  4000, h =  3000, res =  600)
set.seed(6)
corto::scatter(bres, kres, xlab =  "BE2C Log10 Variance", ylab =  "Kelly Log10 Variance", main =  "Variance (Loess Residuals)",
               col =  "gainsboro", extendXlim =  TRUE)
textplot3(bres[toshow], kres[toshow], words =  toshow, font =  2)
dev.off()

### Housekeeping and MYCN expression ----
genes <- c("ACTB", "GAPDH", "B2M", "GUSB")
png("plots/003_exp_vs_mycn.png", w =  4000, h =  2500, res =  300)
par(mfrow =  c(2, 4))
for(gene in genes){
  corto::scatter(log10(bmat["MYCN",]+pseudo), log10(bmat[gene,]+pseudo), main =  "BE2C",
                 col =  col_be2c, xlab =  "MYCN", ylab =  gene)
}
for(gene in genes){
  corto::scatter(log10(kmat["MYCN",]+pseudo), log10(kmat[gene,]+pseudo), main =  "Kelly",
                 col =  col_kelly, xlab =  "MYCN", ylab =  gene)
}
dev.off()

### Expression vs. Cells with Gene ----
# Nr. cells with gene > 0
bcells <- apply(bmat, 1, function(x){sum(x>0)})
kcells <- apply(kmat, 1, function(x){sum(x>0)})

png("plots/003_nrcells_be2c.png", w =  4000, h =  3000, res =  600)
set.seed(2)
plot(bcells, bx, pch =  20, col =  col_be2c, xlim =  c(-200,1400), xlab =  "Nr. Cells with Gene > 0 TPM",
     ylab =  "Log10 Average Expression (TPM)", main =  "BE2C")
mtext(paste0("Total: ", ncol(bmat), " cells"), cex =  1, font =  2)
textplot3(bcells[toshow], bx[toshow], words =  toshow, font =  2)
dev.off()

png("plots/003_nrcells_kelly.png", w =  4000, h =  3000, res =  600)
set.seed(2)
plot(kcells, kx, pch =  20, col =  col_kelly, xlim =  c(-200,1400), xlab =  "Nr. Cells with Gene > 0 TPM",
     ylab =  "Log10 Average Expression (TPM)", main =  "Kelly")
textplot3(kcells[toshow], kx[toshow], words =  toshow, font =  2)
mtext(paste0("Total: ", ncol(kmat)," cells"), cex =  1, font =  2)
dev.off()

