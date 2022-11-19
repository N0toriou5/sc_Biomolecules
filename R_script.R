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
min(kexpmean[kexpmean!=0]) # 0.0006598291
min(bexpmean[bexpmean!=0]) # 0.0003172931
pseudo <- 0.0001
min(kexpvar[kexpvar!=0]) # 0.0003762165
min(bexpvar[bexpvar!=0]) # 0.00009684925
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

### Selected genes with error bars
error.bar  <-  function(x, y, upper, length =  0.1,...){
  arrows(x, y+abs(upper), x, y, angle =  90, code =  3, length =  length, ...)
}
# BE2C
png("plots/003_somegenes_be2c.png", w =  4000, h =  3000, res =  600)
par(las =  2)
value <- sort(bx[toshow], dec =  TRUE)+4
error <- bres[toshow]
bp <- barplot(value, main =  "BE2C",col =  col_be2c, ylab =  c("Log10 Average Expression (TPM)"), ylim =  c(0, 9), yaxt =  "n")
axis(2, at =  c(0:9), labels =  c(-4:5))
error.bar(bp[,1], value, error)
dev.off()

# Kelly
png("plots/003_somegenes_kelly.png", w =  4000, h =  3000, res =  600)
par(las =  2)
value <- sort(kx[toshow], dec =  TRUE)+4
error <- kres[toshow]
bp <- barplot(value, main =  "Kelly", col =  col_kelly, ylab =  c("Log10 Average Expression (TPM)"),
              ylim =  c(0,9), yaxt =  "n")
axis(2, at =  c(0:9), labels =  c(-4:5))
error.bar(bp[, 1], value, error)
dev.off()

# Expression by chromosome band
### Chromosome bands ----
fname <- "data/mlist.rda"
if(!file.exists(fname)){
  mdf <- msigdbr(species =  "Homo sapiens") # Retrieve all human gene sets
  mlist <- mdf %>% split(x =  .$gene_symbol, f =  .$gs_name)
  save(mlist, file =  fname)
}else{load(fname)}
chrom_bands  <-  mlist[grep("chr", names(mlist))] 
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
        "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# Find gene
gene <- "MYCN"
for(i in 1:length(chrom_bands)){
  band <- names(chrom_bands)[i]
  genes_here <- chrom_bands[[band]]
  if(gene%in%genes_here){
    message(gene, " is in ",b and)
  }
}
# MYCN is in chr2p24


# Kelly
cols <- c()
pickcols <- c("royalblue4", "skyblue")
means <- c()
coords <- c()
for(i in 1:length(chrs)){
  chr <- chrs[i]
  herecol <- pickcols[(i%%2)+1]
  sub <- chrom_bands[grep(paste0(chr, "(p|q)"), names(chrom_bands))]
  herecols <- rep(herecol,length(sub))
  cols <- c(cols,herecols)
  coords <- c(coords,rep(chr,length(sub)))
  for(band in names(sub)){
    genes_here <- sub[[band]]
    genes_here <- intersect(rownames(kmat),genes_here)
    allexp <- kmat[genes_here, ]
    mean <- mean(allexp)
    means <- c(means, mean)
    names(means)[length(means)] <- band
  }
}

png("plots/003_chrom_kelly.png", w =  11000, h =  1000, res =  300)
par(las =  2)
bp <- barplot(means, col =  cols, ylab =  "Mean expression in band (TPM)", xaxt =  "n", 
              main =  "Kelly gene expression by chromosome band", ylim =  c(0, 650))
for(i in 1:length(chrs)){
  chr <- chrs[i]
  xwhere <- mean(bp[which(coords == chr)])
  text(xwhere, 600, labels =  chr, font =  2)
}
axis(1, at =  bp, labels =  names(means), cex.axis =  0.7)
dev.off()

# BE2C
cols <- c()
pickcols <- c("red3", "salmon")
means <- c()
coords <- c()
for(i in 1:length(chrs)){
  chr <- chrs[i]
  herecol <- pickcols[(i%%2)+1]
  sub <- chrom_bands[grep(paste0(chr, "(p|q)"), names(chrom_bands))]
  herecols <- rep(herecol, length(sub))
  cols <- c(cols, herecols)
  coords <- c(coords, rep(chr, length(sub)))
  for(band in names(sub)){
    genes_here <- sub[[band]]
    genes_here <- intersect(rownames(bmat), genes_here)
    allexp <- bmat[genes_here, ]
    mean <- mean(allexp)
    means <- c(means, mean)
    names(means)[length(means)] <- band
  }
}

png("plots/003_chrom_be2c.png", w =  11000, h =  1000, res =  300)
par(las =  2)
bp <- barplot(means, col =  cols, ylab =  "Mean expression in band (TPM)", xaxt =  "n", main =  "BE2C gene expression by chromosome band", ylim =  c(0, 650))
for(i in 1:length(chrs)){
  chr <- chrs[i]
  xwhere <- mean(bp[which(coords == chr)])
  text(xwhere, 600, labels =  chr, font =  2)
}
axis(1, at =  bp, labels =  names(means), cex.axis =  0.7)
dev.off()

# Comparison with bulk RNA-Seq data
# Loading the Harenza dataset
fname <- "data/harenza/rawcounts_symbols.rda"
if(!file.exists(fname)){
  ### Load counts (from CellRanger) ----
  rawcounts <- read.delim("data/harenza/harenza.counts.txt.gz", as.is =  TRUE, skip =  1, row.names =  1)
  rawcounts <- as.matrix(rawcounts[, 6:ncol(rawcounts)])
  colnames(rawcounts) <- gsub("\\.sorted\\.bam", "", colnames(rawcounts))
  colnames(rawcounts) <- gsub("harenza_", "", colnames(rawcounts))
  colnames(rawcounts) <- gsub("_.+", "", colnames(rawcounts))
  colnames(rawcounts) <- gsub("\\.", "-", colnames(rawcounts))
  save(rawcounts, file =  "data/harenza/rawcounts_ensg.rda")
  dim(rawcounts) # 58721 genes,  40 samples
  ## Convert to gene symbols
  ensgmat <- rawcounts
  rownames(ensgmat) <- gsub("\\..+", "", rownames(ensgmat))
  tmp <- ens2eg(rownames(ensgmat))
  convlist <- eg2sym(tmp)
  names(convlist) <- names(tmp)
  rawcounts <- squish(ensgmat, convlist =  convlist, method =  "sum", verbose =  TRUE)
  dim(rawcounts) # 26131 genes,  40 samples
  save(rawcounts, file =  fname)
} else {load(fname)}


### Normalization ----
# Extract gene lengths (more precisely,  transcript lengths)
fname <- "data/genelengths.rda"
if(!file.exists(fname)){
  library(GenomicFeatures)
  supportedUCSCtables(genome =  "hg38",  url =  "http://genome.ucsc.edu/cgi-bin/")
  hg <- makeTxDbFromUCSC(genome =  "hg38", tablename =  "refGene")
  exonic <- exonsBy(hg, by =  "gene")
  redexonic <- reduce(exonic)
  genelengths <- sum(width(redexonic))
  names(genelengths) <- eg2sym(names(genelengths))
  genelengths <- genelengths[!duplicated(names(genelengths))]
  genelengths <- genelengths[genelengths>0]
  genelengths <- genelengths[!is.na(genelengths)]
  save(genelengths, file =  fname)
}else{load(fname)}



# Function to calculate FPKM (https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
fpkm <- function(counts, genelengths){
  common <- intersect(rownames(counts), names(genelengths))
  counts <- counts[common, ]
  lengths <- genelengths[common]
  fpms <- apply(counts, 2, function(x){1E6*x/sum(x)})
  fpkms <- fpms
  for(i in 1:nrow(fpms)){
    fpkms[i, ] <- 1E3*fpms[i, ]/lengths[i]
  }
  return(fpkms)
}
fname <- "data/harenza/fpkms.rda"
if(!file.exists(fname)){
  fpkms <- fpkm(rawcounts, genelengths)
  save(fpkms, file =  fname)
}
# Function to calculate TPM (https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
tpm <- function(counts, genelengths){
  common <- intersect(rownames(counts), names(genelengths))
  counts <- counts[common, ]
  lengths <- genelengths[common]
  intermediate <- counts
  for(i in 1:nrow(counts)){
    intermediate[i, ] <- 1E3*counts[i, ]/lengths[i]
  }
  tpms <- apply(intermediate, 2, function(x){1E6*x/sum(x)})
  return(tpms)
}
fname <- "data/harenza/tpms.rda"
if(!file.exists(fname)){
  tpms <- tpm(rawcounts, genelengths)
  save(tpms, file =  fname)
}

# Variance-Stabilizing Transformation (expmat)
fname <- "data/harenza/vstmat.rda"
if(!file.exists(fname)){
  vstmat <- vst(rawcounts)
  save(vstmat, file =  fname)
}

# comparison between sc and bulk
### Load data ----
# Harenza
load("data/harenza/tpms.rda")
harenza <- tpms
# Our own
load("data/tpms.rda")
data <- tpms
bmat <- data[, grep("be2c", colnames(data))]
kmat <- data[, grep("kelly", colnames(data))]
b <- apply(bmat, 1, sum)
k <- apply(kmat, 1, sum)
# Correct names
colnames(harenza)[colnames(harenza) == "SK-N-BE-2--C"] <- "BE2C"
colnames(harenza)[colnames(harenza) == "SK-N-BE-2-"] <- "BE2"

### Correlation matrices ----
output <- matrix(NA, nrow =  ncol(harenza), ncol =  2)
colnames(output) <- c("SCC", "SCC p-value")
rownames(output) <- colnames(harenza)

# Kelly
grpm <- k
for (i in 1:ncol(harenza)){
    cellname <- colnames(harenza)[i]
    message("Doing ", cellname)
    hrpm <- harenza[, cellname]
    common <- intersect(names(grpm), names(hrpm))
    grpm <- grpm[common]
    hrpm <- hrpm[common]
    x <- grpm
    y <- hrpm
    
    ii <- 1
    for(cortype in c("s")){
        cortest <- cor.test(x, y, method =  cortype)
        coeff <- signif(cortest$estimate, 4)
        p <- signif(cortest$p.value, 3)
        output[i, ii] <- coeff
        ii <- ii+1
        output[i, ii] <- p
        ii <- ii+1
    }
}
output <- output[order(-output[, 1]), ]

pander(output, style =  "rmarkdown") # Correlation between Kelly single cell dataset and other datasets
png("plots/005_cortests_kelly_vs_harenza.png", w =  2000, h =  5000, res =  350)
grid.newpage()
grid.table(output, theme =  ttheme_default(base_colour =  "navy"))
dev.off()
write.xlsx2(output, file =  "results/Supp_Table_S2_Kelly_vs_Bulk.xlsx")

# BE2C
output <- matrix(NA, nrow =  ncol(harenza), ncol =  2)
colnames(output) <- c("SCC", "SCC p-value")
rownames(output) <- colnames(harenza)
grpm <- b
for (i in 1:ncol(harenza)){
    cellname <- colnames(harenza)[i]
    message("Doing ", cellname)
    hrpm <- harenza[, cellname]
    common <- intersect(names(grpm), names(hrpm))
    grpm <- grpm[common]
    hrpm <- hrpm[common]
    x <- grpm
    y <- hrpm
    
    ii <- 1
    for(cortype in c("s")){
        cortest <- cor.test(x, y, method =  cortype)
        coeff <- signif(cortest$estimate, 4)
        p <- signif(cortest$p.value, 3)
        output[i, ii] <- coeff
        ii <- ii+1
        output[i, ii] <- p
        ii <- ii+1
    }
}
output <- output[order(-output[, 1]), ]

pander(output, style =  "rmarkdown") # Correlation between BE2C single cell dataset and other datasets
png("plots/005_cortests_be2c_vs_harenza.png", w =  2000, h =  5000, res =  350)
grid.table(output, theme =  ttheme_default(base_colour =  "red3"))
dev.off()
write.xlsx2(output, file =  "results/Supp_Table_S1_BE2C_vs_Bulk.xlsx")

### Scatter plots bulk vs sc ----
pseudo <- 0.0001
toshow <- c("MYCN", "GAPDH", "GUSB", "ACTB", "B2M")
# Kelly
x <- log10(k+pseudo)
y <- log10(harenza[, "KELLY"]+pseudo)
png("plots/005_scatter_kelly.png", w =  4000, h =  3000, res =  600)
scatter(x, y, xlab =  "Single Cell RNA-Seq (Sum of TPMs)", ylab =  "Bulk RNA-Seq (TPM)", 
        main =  "Kelly cells", col =  col_kelly, method =  "spearman")
textplot3(x[toshow], y[toshow], words =  toshow, font =  2)
dev.off()

# BE2C
x <- log10(b+pseudo)
y <- log10(harenza[, "BE2C"]+pseudo)

png("plots/005_scatter_be2c.png", w =  4000, h =  3000, res =  600)
scatter(x, y, xlab =  "Single Cell RNA-Seq (Sum of TPMs)", ylab =  "Bulk RNA-Seq (TPM)", 
        main =  "BE2C cells", col =  col_be2c, method =  "spearman")
textplot3(x[toshow], y[toshow], words =  toshow, font =  2)
dev.off()

### Rtsne ----
# Load cell line annotation
cl <- read.csv("data/NBLcellLines.csv", header =  TRUE)
cl <- setNames(cl[, 4], cl[, 1])

# Everything together
lh <- log10(harenza+pseudo)
lk <- log10(k+pseudo)
lb <- log10(b+pseudo)
common <- intersect(rownames(lh), intersect(names(lk), names(lb)))
tsnemat <- cbind(lh[common, ], lk[common], lb[common])
colnames(tsnemat)[(ncol(tsnemat)-1):ncol(tsnemat)] <- c("scKELLY", "scBE2C")
# Prepare the matrix
topvars <- names(sort(apply(tsnemat, 1, var), decreasing =  TRUE))[1:5000]
tsnemat <- tsnemat[topvars, ]
# Seed for TSNE and calculate TSNE
set.seed(4)
ttt <- Rtsne(t(tsnemat), perplexity =  10)


# Start preparing the plot
x <- setNames(ttt$Y[, 1], colnames(tsnemat))
y <- setNames(ttt$Y[, 2], colnames(tsnemat))
# Shapes
shapes <- rep(15, length(x))
names(shapes) <- names(x)
shapes[names(cl[cl == "notAmplified"])] <- 16
png("plots/005_tsne_with_bulk.png", w =  3000, h =  3000, res =  500)
plot(x, y, pch =  shapes, xlab =  "TSNE1", ylab =  "TSNE2", main =  "TSNE Representation of NBL Cell Lines", 
     xlim =  1.1*c(min(x), max(x)), ylim =  1.1*c(min(y), max(y)), 
     col =  c(rep("#00000099", length(x)-2), col_kelly, col_be2c)
)
set.seed(4)
textplot3(x, y, words =  names(x), cex =  0.8, 
          font =  c(rep(1, length(x)-2), 2, 2), 
          col =  c(rep("black", length(x)-2), col_kelly, col_be2c), padding =  "  ", pos =  3, offset =  0.1
)
legend("bottomleft", pch =  c(15, 16), c("Amplified", "Not Amplified"), title =  "MYCN status")
dev.off()

# Again without BE2C
x <- setNames(ttt$Y[, 1], colnames(tsnemat))
y <- setNames(ttt$Y[, 2], colnames(tsnemat))
x <- x[1:(length(x)-1)]
y <- y[1:(length(y)-1)]
# Shapes
shapes <- rep(15, length(x))
names(shapes) <- names(x)
shapes[names(cl[cl == "notAmplified"])] <- 16

png("plots/005_tsne_with_bulk_Kelly.png", w =  6000, h =  3000, res =  600)
plot(x, y, pch =  shapes, xlab =  "TSNE1", ylab =  "TSNE2", main =  "Clustering of NBL Cell Lines", 
     xlim =  1.1*c(min(x), max(x)), ylim =  1.1*c(min(y), max(y)), 
     col =  c(rep("#00000099", length(x)-1), col_kelly)
)
set.seed(4)
textplot3(x, y, words =  names(x), cex =  0.8, 
          font =  c(rep(1, length(x)-1), 2), 
          col =  c(rep("black", length(x)-1), col_kelly), padding =  "  ", pos =  3, offset =  0.1
)
legend("bottomleft", pch =  c(15, 16), c("Amplified", "Not Amplified"), title =  "MYCN status")
dev.off()

# Clustering single cell data
### Load Seurat object (already LogNormalized)
load("data/seuset.rda")
load("data/rawcounts_symbols.rda")

### Further processing ----
seuset <- FindVariableFeatures(seuset, selection.method =  "vst", nfeatures =  Inf)
all.genes <- rownames(seuset)
seuset <- ScaleData(seuset, features =  all.genes)
expmat <- as.matrix(seuset[["RNA"]]@scale.data)
dim(expmat) # 15782 genes,  2067 cells

### Seurat-based clustering ----
PCA <- RunPCA(seuset, features =  VariableFeatures(seuset))

# TSNE
set.seed(1)
TSNE <- RunTSNE(PCA)
png("plots/006a_tsne_seurat_clustering.png", width =  3000, height =  3000,  res =  600)
DimPlot(TSNE, reduction =  "tsne", cols =  c(col_be2c, col_kelly))+
ggtitle('Seurat TSNE clustering of KELLY and BE2C')+
theme(plot.title =  element_text(hjust =  0.5))
dev.off()

# UMAP
set.seed(1)
UMAP <- RunUMAP(PCA,  dims  =   1:10)
png("plots/006a_umap_seurat_clustering.png", width =  3000, height =  3000,  res =  600)
DimPlot(UMAP,  reduction  =   'umap', cols =  c(col_be2c, col_kelly))+
ggtitle('Seurat UMAP clustering of KELLY and BE2C')+
theme(plot.title  =   element_text(hjust  =   0.5))
dev.off()

## Cell Cycle clustering
### Effects of cell cycle and read numbers on clustering ----
# Cell Cycle Markers,  from Tirosh et al,  2015
ccgenes <- readLines("data/regev_lab_cell_cycle_genes.txt")
ccgenes <- eg2sym(sym2eg(ccgenes))
s.genes  <-  ccgenes[1:43]
g2m.genes  <-  ccgenes[44:97]
s.genes <- intersect(s.genes, rownames(seuset))
g2m.genes <- intersect(g2m.genes, rownames(seuset))
# Apply scoring
seuset <- CellCycleScoring(seuset, s.features =  s.genes, g2m.features =  g2m.genes, set.ident =  TRUE)

# Rtsne,  top var genes
allvars <- apply(expmat, 1, var)
topvars <- names(sort(allvars, dec =  TRUE))[1:1000]
topvarmat <- expmat[topvars, ]
fname <- "data/rtsne.rda"
if(!file.exists(fname)){
    set.seed(1)
    rtsne <- Rtsne(t(topvarmat))
    save(rtsne, file =  fname)
}else{load(fname)}

# Coloring cell lines
x <- setNames(rtsne$Y[, 1], colnames(topvarmat))
y <- setNames(rtsne$Y[, 2], colnames(topvarmat))
mycols <- rep("black", ncol(topvarmat))
mycols[grep("be2c", colnames(topvarmat))] <- col_be2c
mycols[grep("kelly", colnames(topvarmat))] <- col_kelly
png("plots/006b_tsne_seurat_celllines.png", w =  3000, h =  3000, res =  600)
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "TSNE clustering of KELLY and BE2C", 
     col =  mycols, xlim =  c(min(x), max(x)*1.5))
mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
grid()
legend("bottomright", pch =  19, col =  c(col_be2c, col_kelly), legend =  c("BE2C", "Kelly"))
dev.off()

# Color cell cycle phases
phases <- setNames(as.character(seuset@meta.data$Phase), rownames(seuset@meta.data))
mycols <- phases
mycols[phases == "G1"] <- "salmon"
mycols[phases == "S"] <- "cornflowerblue"
mycols[phases == "G2M"] <- "seagreen"
mycols <- mycols[colnames(topvarmat)]
png("plots/006b_tsne_seurat_celllcycle.png", w =  3000, h =  3000, res =  600)
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "TSNE clustering of KELLY and BE2C", 
     col =  mycols, xlim =  c(min(x), max(x)*1.5))
mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
grid()
legend("bottomright", pch =  19, col =  c("salmon", "cornflowerblue", "seagreen"), 
       legend =  c("G1", "S", "G2M"))
dev.off()

# Barplot cell cycle distribution
kphases <- phases[grep("kelly", names(phases))]
bphases <- phases[grep("be2c", names(phases))]
ktab <- table(kphases)[c("G1", "S", "G2M")]
btab <- table(bphases)[c("G1", "S", "G2M")]

png("plots/006b_phases_barplots.png", w =  3000, h =  3000, res =  500)
par(mfrow =  c(1, 2))
max <- max(c(ktab, btab))*1.2
bp <- barplot(ktab, col =  c("salmon", "cornflowerblue", "seagreen"), main =  "Kelly", 
              ylim =  c(0, max), ylab =  "nr. cells")
mtext(paste0(length(kphases), " cells"), cex =  0.8)
perc <- round(100*ktab/length(kphases), 2)
text(bp, ktab, labels =  paste0(perc, "%"), pos =  3, font =  3)
bp <- barplot(btab, col =  c("salmon", "cornflowerblue", "seagreen"), main =  "BE2C", 
              ylim =  c(0, max))
mtext(paste0(length(bphases), " cells"), cex =  0.8)
perc <- round(100*btab/length(bphases), 2)
text(bp, btab, labels =  paste0(perc, "%"), pos =  3, font =  3)
dev.off()
par(mfrow =  c(1, 1))

## Read number clustering
# Color nr. reads coloring
nreads <- apply(rawcounts, 2, sum)
colfunc  <-  colorRampPalette(c("cornflowerblue", "red3", "orange"))
mycols <- colfunc(100)[as.numeric(cut(nreads, breaks =  100))]
png("plots/006b_tsne_seurat_nreads.png", w =  3000, h =  3000, res =  600)
layout(matrix(1:2, ncol =  2),  width  =   c(3, 1), height  =   c(1, 1))
par(mar =  c(5.1, 4.1, 4.1, 0.1))
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "TSNE clustering", col =  mycols, 
     xlim =  c(min(x), max(x)*1.5))
mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
grid()
par(mar =  c(5.1, 0.1, 4.1, 0.1))
plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', main =  'Nr. Reads (k)')
legend_image <- as.raster(rev(matrix(colfunc(100),  ncol =  1)))
rasterImage(legend_image, 0, 0, 1, 1)
text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(quantile(nreads/1E3)))
dev.off()

# Read counts distributions
png("plots/006b_nreads_lineplots.png", w =  3000, h =  3000, res =  600)
knreads <- nreads[grep("kelly", names(nreads))]
bnreads <- nreads[grep("be2c", names(nreads))]
plot(density(knreads), col =  col_kelly, lwd =  3, ylim =  c(0, 4e-5), main =  "Nr. reads/cell distribution")
lines(density(bnreads), col =  col_be2c, lwd =  3)
legend("topright", legend =  c(
    paste0("Kelly (mean =  ", round(mean(knreads), 2), ")"), 
    paste0("BE2C (mean =  ", round(mean(bnreads), 2), ")")
), lwd =  3, col =  c(col_kelly, col_be2c))
dev.off()

## Regress out cell cycle and n. UMI
