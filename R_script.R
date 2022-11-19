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
fname <- "data/seuset_regressed.rda"
if(!file.exists(fname)){
    seuset <- ScaleData(object =  seuset, vars.to.regress =  c("nUMI", "S.Score", "G2M.Score"))
    save(seuset, file =  fname)
}else{load(fname)}
regmat <- as.matrix(seuset[["RNA"]]@scale.data)
dim(regmat) # 15782 genes 2067 cells
save(regmat, file =  "data/regmat.rda")
topvarmat <- regmat

# Rtsne,  top var genes
allvars <- apply(topvarmat, 1, var)
topvars <- names(sort(allvars, dec =  TRUE))[1:1000]
topvarmat <- topvarmat[topvars, ]
set.seed(1)
rtsne <- Rtsne(t(topvarmat))
save(rtsne, file =  fname)

# Coloring cell lines
x <- setNames(rtsne$Y[, 1], colnames(topvarmat))
y <- setNames(rtsne$Y[, 2], colnames(topvarmat))
mycols <- rep("black", ncol(topvarmat))
mycols[grep("be2c", colnames(topvarmat))] <- col_be2c
mycols[grep("kelly", colnames(topvarmat))] <- col_kelly
png("plots/006c_tsne_seurat_celllines_postregr.png", w =  3000, h =  3000, res =  600)
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "Post-regression TSNE", col =  mycols, xlim =  c(min(x), max(x)*1.5))
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
png("plots/006c_tsne_seurat_celllcycle_postregr.png", w =  3000, h =  3000, res =  600)
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "Post-regression TSNE", 
     col =  mycols, xlim =  c(min(x), max(x)*1.5))
mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
grid()
legend("bottomright", pch =  19, col =  c("salmon", "cornflowerblue", "seagreen"), 
       legend =  c("G1", "S", "G2M"))
dev.off()

# Color nr. reads coloring
nreads <- apply(rawcounts, 2, sum)
colfunc  <-  colorRampPalette(c("cornflowerblue", "red3", "orange"))
mycols <- colfunc(100)[as.numeric(cut(nreads, breaks =  100))]
png("plots/006c_tsne_seurat_nreads_postregr.png", w =  3000, h =  3000, res =  600)
layout(matrix(1:2, ncol =  2),  width  =   c(2, 1), height  =   c(1, 1))
plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  "Post-regression TSNE", 
     col =  mycols, xlim =  c(min(x), max(x)*1.5))
mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
grid()
plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', main =  'Nr. Reads (k)')
legend_image <- as.raster(rev(matrix(colfunc(100),  ncol =  1)))
rasterImage(legend_image, 0, 0, 1, 1)
text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(quantile(nreads/1E3)))
dev.off()

### Perform differential expression analysis
load("data/seuset_regressed.rda")
load("data/expmat.rda") # LogNormalized by Seurat
load("data/regmat.rda") # LogNormalized by Seurat + Regressed out cycle and nUMI
regmat <- regmat+abs(min(regmat))
# Function to calculate Differential expression with Wilcoxon tests
wexp <- function(matx, maty){
    columns <- c("log2fc", "wstat", "p", "fdr")
    output <- matrix(NA, nrow =  nrow(matx), ncol =  length(columns))
    colnames(output) <- columns
    rownames(output) <- rownames(matx)
    pb <- txtProgressBar(0, nrow(matx), style =  3)
    for(i in 1:nrow(matx)){
        x <- matx[i, ]
        y <- maty[i, ]
        l2fc <- log2(mean(x)/mean(y))
        wt <- wilcox.test(x, y)
        p <- wt$p.value
        if(p == 0){p <- 1e-301}
        stat <- -log10(p)*sign(l2fc)
        output[i, ] <- c(l2fc, stat, p, NA)
        setTxtProgressBar(pb, i)
    }
    output[, "fdr"] <- p.adjust(output[, "p"], method =  "BH")
    return(as.data.frame(output, stringsAsFactors =  FALSE))
}

# be2c vs. kelly (regressed)
fname <- "results/007_res_be2c_vs_kelly.rda"
if(!file.exists(fname)){
    g1 <- regmat[, grep("kelly", colnames(regmat))]
    g2 <- regmat[, grep("be2c", colnames(regmat))]
    res <- wexp(g2, g1)
    res <- res[order(res[, "p"]), ]
    save(res, file =  fname)
    write.csv(res, file =  "results/007_res_be2c_vs_kelly.csv")
}else{load(fname)}


# be2c vs. kelly (not regressed)
fname <- "results/007_res_be2c_vs_kelly_notregressed.rda"
if(!file.exists(fname)){
    g1 <- expmat[, grep("kelly", colnames(expmat))]
    g2 <- expmat[, grep("be2c", colnames(expmat))]
    res <- wexp(g2, g1)
    res <- res[order(res[, "p"]), ]
    save(res, file =  fname)
    write.csv(res, file =  "results/007_res_be2c_vs_kelly_notregressed.csv")
}else{load(fname)}

### Compare signatures (regressed and not regressed) ----
load("results/007_res_be2c_vs_kelly.rda")
sig_reg <- setNames(res$log2fc, rownames(res))
load("results/007_res_be2c_vs_kelly_notregressed.rda")
sig_notreg <- setNames(res$log2fc, rownames(res))
sig_reg <- sig_reg[!is.infinite(sig_reg)]
sig_notreg <- sig_notreg[!is.infinite(sig_notreg)]

png("plots/007_sigreg_vs_signotreg.png", w =  3000, h =  3000, res =  600)
scatter(sig_reg, sig_notreg, xlab =  "Signature,  regressed", 
        ylab =  "Signature,  not regressed", main =  "BE2C vs. Kelly Differential Expression")
dev.off()

### Bulk signature (from Harenza) ----
# We must compare log2FC,  since the bulk RNA-Seq has only a single sample per condition
load("data/harenza/vstmat.rda")
hsig <- log2(vstmat[, "SK-N-BE-2--C"]/vstmat[, "KELLY"])

load("data/expmat.rda")
min(expmat[expmat!= 0]) # 0.05
expmat <- expmat+0.01
scb <- apply(expmat[, grep("be2c", colnames(expmat))], 1, mean)
sck <- apply(expmat[, grep("kelly", colnames(expmat))], 1, mean)
scsig <- log2(scb/sck)

png("plots/007_sc_vs_bulk.png", w =  3000, h =  3000, res =  600)
scatter(hsig, scsig, xlab =  "log2FC,  Bulk RNA-Seq", ylab =  "log2FC,  Single Cell RNA-Seq",
        main =  "BE2C vs. Kelly Differential Expression")
dev.off()

### Volcano plot ----
load("results/007_res_be2c_vs_kelly.rda")
x <- setNames(res$log2fc, rownames(res))
y <- setNames(-log10(res$fdr), rownames(res))


png("plots/007_volcanoes.png", w =  6000, h =  3000, res =  600)
par(mfrow =  c(3, 2))
pthrs <- c(2, 5, 10, 20, 50, 100)
toshow <- c("RPS25", "RPL27", "MYCN")
for(pthr in pthrs){
    plot(x, y, pch =  20, xlab =  "log2 Fold Change", ylab =  "-log10(FDR)", 
         col =  "#00000011", main =  "BE2C vs. Kelly")
    up <- which(x>0&y>pthr)
    points(x[up], y[up], col =  "#FF000011", pch =  20)
    dn <- which(x<0&y>pthr)
    points(x[dn], y[dn], col =  "#0000FF11", pch =  20)
    mtext(paste0("Significant at FDR =  10^-", pthr, ": ", length(up), "(up),  ", length(dn), "(down)"), cex =  0.6)
    if(pthr == 2){
        text(x[toshow], y[toshow], labels =  toshow, font =  2)
    }
}
dev.off()

### Gene set enrichment analysis
fname <- "data/mlist.rda"
if(!file.exists(fname)){
  mdf <- msigdbr(species =  "Homo sapiens") # Retrieve all human gene sets
  mlist <- mdf %>% split(x =  .$gene_symbol, f =  .$gs_name)
  save(mlist, file =  fname)
}else{load(fname)}


### Get signature BE2C vs. Kelly ----
load("results/007_res_be2c_vs_kelly.rda")
signature <- setNames(-log10(res$p)*sign(res$log2fc), rownames(res))

### Check contrast orientation ----
signature[1:10]
load("data/tpms.rda")
mean(tpms["RPS25", grep("be2c", colnames(tpms))])
mean(tpms["RPS25", grep("kelly", colnames(tpms))])
mean(tpms["RPL27", grep("be2c", colnames(tpms))])
mean(tpms["RPL27", grep("kelly", colnames(tpms))])
# Positive is,  as planned,  higher in BE2C


### Run GSEA (fgsea function) ----
fname <- "results/008_gsea.rda"
if(!file.exists(fname)){
    library(fgsea)
    gsea <- fgsea(pathways =  mlist, stats =  signature, nperm =  1E6, minSize =  4, maxSize =  Inf, nproc =  7)
    gsea <- gsea[order(gsea$pval), ]
    gsea <- as.data.frame(gsea)
    save(gsea, file =  fname)
    # Save as supp table
    write.csv(gsea, file =  "results/008_gsea_BE2C_vs_Kelly.csv")
}else{load(fname)}

### Run GSEA (corto function) ----
fname <- "results/008_gsea_corto.rda"
if(!file.exists(fname)){
    gsea <- matrix(nrow =  length(mlist), ncol =  3)
    rownames(gsea) <- names(mlist)
    colnames(gsea) <- c("ES", "NES", "p")
    set.seed(1)
    pb <- txtProgressBar(0, length(mlist), style =  3)
    for(i in 1:length(mlist)){
        pname <- names(mlist)[i]
        pathway <- mlist[[pname]]
        if(length(intersect(names(signature), pathway))>1){
            obj <- gsea(signature, pathway, method =  "pareto")
            gsea[pname, ] <- c(obj$es, obj$nes, obj$p)
        }else{
            gsea[pname, ] <- c(0, 0, 1)
        }
        setTxtProgressBar(pb, i)
    }
    fdr <- p.adjust(gsea[, "p"], method =  "BH")
    gsea <- cbind(gsea, fdr)
    colnames(gsea)[ncol(gsea)] <- "fdr"
    save(gsea, file =  fname)
    # Save as supp table
    write.csv(gsea, file =  "results/008_gsea_BE2C_vs_Kelly_corto.csv")
}else{load(fname)}

# All pathways
gsea <- as.data.frame(gsea)
gsea <- gsea[order(-abs(gsea$NES)), ]


### Table of top pathways ----
top <- gsea[gsea$NES<0, ][1:15, ]
top <- rbind(top, gsea[gsea$NES>0, ][15:1, ])
top <- top[order(top$NES), ]
toplot <- setNames(top$NES, rownames(top))
# Format
#names(toplot) <- gsub("GO_", "", names(toplot))
names(toplot) <- gsub("_", " ", names(toplot))
names(toplot) <- str_to_title(names(toplot))
png("plots/008_gsea_be2c_vs_kelly.png", w =  6000, h =  3000, res =  500)
par(mar =  c(4, 1, 3, 1))
bp <- barplot(toplot, horiz =  TRUE, xlab =  "Normalized Enrichment Score", 
            xlim =  1.3*c(-max(abs(toplot)), max(abs(toplot))), 
            main =  "BE2C vs. Kelly,  top Pathways", 
            col =  rep(c("cornflowerblue", "salmon"), each =  15), 
            yaxt =  "n", cex.main =  2
)
text(0, bp[1:15, 1], names(toplot)[1:15], pos =  4)
text(0, bp[16:30, 1], names(toplot)[16:30], pos =  2)
abline(v =  c(-p2z(0.05), p2z(0.05)), lty =  2)
dev.off()

### Check for Neuroblastoma ----
gsea[grep("NEUROBLASTOMA|NBL", rownames(gsea)), ]

## Plot top GSEAs ----
pathways <- c("REACTOME_METABOLISM_OF_RNA", "CHICAS_RB1_TARGETS_CONFLUENT", "GO_RIBOSOME_BIOGENESIS", "LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_UP", "LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN")
for(pathway in pathways){
    set <- mlist[[pathway]]
    obj <- gsea(signature, set, method =  "pareto")
    png(paste0("plots/008_gsea_", pathway, ".png"), w =  3000, h =  3000, res =  600)
    plot_gsea(obj, ext_nes =  gsea[rownames(gsea)==pathway, "NES"], title =  pathway, colBarcode = "#00000033")
    dev.off()
}

## Double gsea (TEST)
path1 <- "LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_UP"
path2 <- "LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN"
set1 <- mlist[[path1]]
set2 <- mlist[[path2]]
obj <- gsea2(signature, set1, set2, method =  "pareto")
plot_gsea2(obj)

### Pathway analysis at single-cell level
load("data/seuset_regressed.rda")

### Calculate Rtsne ----
load("data/rtsne.rda")
load("data/regmat.rda")
x <- setNames(rtsne$Y[, 1], colnames(regmat))
y <- setNames(rtsne$Y[, 2], colnames(regmat))

### Overlay single genes over Rtsne ----
# Get TPM per cell
load("data/tpms.rda")
genes <- c("MYCN", "GAPDH", "ACTB", "GUSB", "B2M", "ALK", "LMO1")
for(gene in genes){
    vector <- rank(tpms[gene, ])
    
    ## Color by TPM rank
    colfunc <- colorRampPalette(c("navy", "grey", "red3"))
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
    # Plot
    png(paste0("plots/009_tsne_rank_", gene, ".png"), w =  4000, h =  2000, res =  450)
    layout(matrix(1:2, ncol =  2),  width  =   c(2, 1), height  =   c(1, 1))
    plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  paste0(gene), 
         col =  mycols, xlim =  c(min(x), max(x)*1.5))
    mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
    grid()
    plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', 
         main =  paste0(" expression (rank)"))
    legend_image <- as.raster(rev(matrix(colfunc(1000),  ncol =  1)))
    rasterImage(legend_image, 0, 0, 1, 1)
    text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(quantile(vector)))
    dev.off()
}

for(gene in genes){
    vector <- rank(tpms[gene, ])
    
    ## Color by TPM rank
    colfunc <- colorRampPalette(c("navy", "grey", "red3"))
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]

    ## Color by TPM
    colfunc <- colorRampPalette(c("navy", "navy", "navy", "grey", "salmon", "red", "red3"))
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
    vector <- tpms[gene, ]
    pseudo <- min(vector[vector!= 0]) # Pseudocount to ignore dropout effects
    vector <- log10(vector+pseudo)
    
    # Plot
    png(paste0("plots/009_tsne_log10_", gene, ".png"), w =  4500, h =  2000, res =  500)
    lmatrix <- t(c(1, 2, 2, 2, 3))
    layout(lmatrix)
    par(mar =  c(5.1, 4.1, 4.1, 2.1))
    hist(vector, main =  paste0(gene, ""), xlab =  "Log10 TPM")
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
    plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  paste0(gene), col =  mycols, xlim =  c(min(x), max(x)*1.5))
    mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
    par(mar =  c(5.1, 1.1, 4.1, 2.1))
    plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', main =  paste0("Log10 TPM"))
    legend_image <- as.raster(rev(matrix(colfunc(1000),  ncol =  1)))
    rasterImage(legend_image, 0, 0, 1, 1)
    text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(seq(min(vector), max(vector), length.out =  5), 2))
    dev.off()
}

### Load the msigdb database 7.1.1 ----
fname <- "data/mlist.rda"
if(!file.exists(fname)){
  mdf <- msigdbr(species =  "Homo sapiens") # Retrieve all human gene sets
  mlist <- mdf %>% split(x =  .$gene_symbol, f =  .$gs_name)
  save(mlist, file =  fname)
}else{load(fname)}

### Calculate single-cell gsea (better to use area) ----
fname <- "results/009_ssgsea_regexpmat.rda"
load("data/regmat.rda") # LogNormalized by Seurat + Regressed out cycle and nUMI
# regmat <- regmat+abs(min(regmat))
if(!file.exists(fname)){
    scalemat <- t(scale(t(regmat)))
    ssgsea <- area(signatures =  scalemat, groups =  mlist, minsize =  5)
    save(ssgsea, file =  fname)    
}else{load(fname)}

### Extract highest significance + highest variance pathways ----
dim(ssgsea) # 24472 pathways,  2067 cells
sumval <- apply(ssgsea, 1, sum)
sumabs <- apply(ssgsea, 1, function(x){sum(abs(x))})
var <- apply(ssgsea, 1, var)

png("plots/009b_ss_pathways.png", w =  6000, h =  2000, res =  600)
par(mfrow =  c(1, 3))
plot(sumval, var, xlab =  "NES sum", ylab =  "NES variance")
plot(sumabs, var, xlab =  "NES sum of absolute", ylab =  "NES variance", main =  "Pathways in Single Cell Dataset")
toshow <- names(sort(var, dec =  TRUE))[1:4]
text(sumabs[toshow], var[toshow], labels =  toshow, pos =  2, cex =  0.4)
plot(sumabs, var/sumabs, xlab =  "NES sum of absolute", ylab =  "NES variance / sum of abs")
toshow <- names(sort(var/sumabs, dec =  TRUE))[1:4]
text(sumabs[toshow], (var/sumabs)[toshow], labels =  toshow, pos =  2, cex =  0.4)
dev.off()

### Highest variable within Kelly and BE2C ----
ssk <- ssgsea[, grep("kelly", colnames(ssgsea))]
ssb <- ssgsea[, grep("be2c", colnames(ssgsea))]
# Kelly
sumval <- apply(ssk, 1, sum)
sumabs <- apply(ssk, 1, function(x){sum(abs(x))})
var <- apply(ssk, 1, var)
png("plots/009b_ss_pathways_kelly.png", w =  6000, h =  2000, res =  600)
par(mfrow =  c(1, 3))
plot(sumval, var, xlab =  "NES sum", ylab =  "NES variance")
plot(sumabs, var, xlab =  "NES sum of absolute", ylab =  "NES variance", main =  "Pathways in Kelly Cells")
toshow <- names(sort(var, dec =  TRUE))[1:6]
text(sumabs[toshow], var[toshow], labels =  toshow, pos =  2, cex =  0.4)
plot(sumabs, var/sumabs, xlab =  "NES sum of absolute", ylab =  "NES variance / sum of abs")
toshow <- names(sort(var/sumabs, dec =  TRUE))[1:4]
text(sumabs[toshow], (var/sumabs)[toshow], labels =  toshow, pos =  2, cex =  0.4)
dev.off()

# BE2C
sumval <- apply(ssb, 1, sum)
sumabs <- apply(ssb, 1, function(x){sum(abs(x))})
var <- apply(ssb, 1, var)
png("plots/009b_ss_pathways_be2c.png", w =  6000, h =  2000, res =  600)
par(mfrow =  c(1, 3))
plot(sumval, var, xlab =  "NES sum", ylab =  "NES variance")
plot(sumabs, var, xlab =  "NES sum of absolute", ylab =  "NES variance", main =  "Pathways in BE2C Cells")
toshow <- names(sort(var, dec =  TRUE))[1:6]
text(sumabs[toshow], var[toshow], labels =  toshow, pos =  2, cex =  0.4)
plot(sumabs, var/sumabs, xlab =  "NES sum of absolute", ylab =  "NES variance / sum of abs")
toshow <- names(sort(var/sumabs, dec =  TRUE))[c(1, 3, 4)]
text(sumabs[toshow], (var/sumabs)[toshow], labels =  toshow, pos =  4, cex =  0.4)
dev.off()

### Overlay selected pathways over TSNE ----
colfunc <- colorRampPalette(c("navy", "navy", "navy", "grey", "salmon", "red", "red3"))
mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
pathways <- c("KEGG_RIBOSOME", "MORF_SOD1", "HALLMARK_MYC_TARGETS_V1", "GO_OXIDATIVE_PHOSPHORYLATION")
for(pathway in pathways){
    vector <- ssgsea[pathway, ]
    # Plot
    png(paste0("plots/009b_ss_pathway_", pathway, ".png"), w =  4000, h =  2000, res =  450)
    layout(matrix(1:2, ncol =  2),  width  =   c(2, 1), height  =   c(1, 1))
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
    plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  paste0(pathway), col =  mycols, xlim =  c(min(x), max(x)*1.5))
    mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
    plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', main =  paste0("Normalized Enrichment Score"), cex.main =  0.8)
    legend_image <- as.raster(rev(matrix(colfunc(1000),  ncol =  1)))
    rasterImage(legend_image, 0, 0, 1, 1)
    text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(seq(min(vector), max(vector), length.out =  5), 2))
    dev.off()
}

### correlation MYCN vs. MYC pathway ----
path <- ssgsea["HALLMARK_MYC_TARGETS_V1", ]
gene <- tpms["MYCN", ]
cols <- rep("white", ncol(tpms))
cols[grep("kelly", colnames(tpms))] <- col_kelly
cols[grep("be2c", colnames(tpms))] <- col_be2c
png("plots/009c_MYCN_vs_MYCpathway.png", w =  3000, h =  3000, res =  600)
scatter(path, gene, col =  cols, xlab =  "MSIGDB Hallmark MYC V1 Targets (NES)", ylab =  "MYCN Expression (TPM)", 
        main =  "MYCN vs. MYC pathway"
)
legend("topright", pch =  19, col =  c(col_be2c, col_kelly), legend =  c("BE2C", "Kelly"))
dev.off()

### Master Regulator Analysis
load("data/regmat.rda")
be2c <- regmat[, grep("be2c", colnames(regmat))]
kelly <- regmat[, grep("kelly", colnames(regmat))]

### Get regulons from independent datasets ----
tfs <- read.delim("code/tfgenes_2020_09_11.txt", header =  FALSE)[, 2]

# Generate networks for all NBL datasets
ds <- c("kocak_NBL", "nrc_NBL", "target_NBL")
for(d in ds){
    message("Doing ", d)
    fname <- paste0("networks/", d, "-regulon.rda")
    if(!file.exists(fname)){
        load(paste0("../masterset/data/", d, "-expmat.rda")) # This step requires the original data,  not included here
        regulon <- corto(expmat, tfs, nbootstraps =  100, p =  1e-8, nthreads =  7, verbose =  TRUE)
        save(regulon, file =  fname)
    }
}

### Run BE2C vs. Kelly MRA ----
fname <- "results/010_mra.rda"
if(!file.exists(fname)){
    # Run master regulator analysis (target dataset)
    load("networks/target_NBL-regulon.rda")
    corto_tmra <- mra(be2c, kelly, regulon =  regulon, minsize =  10, nthreads =  8)
    # Run master regulator analysis (kocak dataset)
    load("networks/kocak_NBL-regulon.rda")
    corto_kmra <- mra(be2c, kelly, regulon =  regulon, minsize =  10, nthreads =  8)
    # Run master regulator analysis (nrc dataset)
    load("networks/nrc_NBL-regulon.rda")
    corto_nmra <- mra(be2c, kelly, regulon =  regulon, minsize =  10, nthreads =  8)
    #
    save(corto_tmra, corto_kmra, corto_nmra, file =  fname)
}else{load(fname)}


### Agreement between MRAs ----
target <- corto_tmra$nes
kocak <- corto_kmra$nes
nrc <- corto_nmra$nes

# Significant MRs
sig_target_up <- names(which(corto_tmra$pvalue<= 1e-180&corto_tmra$nes>0, useNames =  TRUE))
sig_kocak_up <- names(which(corto_kmra$pvalue<= 1e-180&corto_kmra$nes>0, useNames =  TRUE))
sig_nrc_up <- names(which(corto_nmra$pvalue<= 1e-180&corto_nmra$nes>0, useNames =  TRUE))
length(sig_target_up) #
length(sig_kocak_up) #
length(sig_nrc_up) #
sig_target_dn <- names(which(corto_tmra$pvalue<= 1e-110&corto_tmra$nes<0, useNames =  TRUE))
sig_kocak_dn <- names(which(corto_kmra$pvalue<= 1e-110&corto_kmra$nes<0, useNames =  TRUE))
sig_nrc_dn <- names(which(corto_nmra$pvalue<= 1e-110&corto_nmra$nes<0, useNames =  TRUE))
length(sig_target_dn) #
length(sig_kocak_dn) #
length(sig_nrc_dn) #

# Intersection
int_up <- intersect(sig_target_up, intersect(sig_kocak_up, sig_nrc_up))
int_dn <- intersect(sig_target_dn, intersect(sig_kocak_dn, sig_nrc_dn))
length(int_up) #
length(int_dn) #
genes <- unique(c(int_up, int_dn))
length(genes) # 
# NOTCH1 is not a TF
genes <- setdiff(genes, "NOTCH1")
png("plots/010_compareMRA.png", w =  4000, h =  2000, res =  450)
set.seed(1) # for reproducible label placement
par(mfrow =  c(1, 2))
y <- spread.labs(kocak[genes], mindiff =  4)
scatter(target, kocak, main =  "BE2C vs. Kelly Master Regulator Analysis", xlab =  "TARGET (NES)", 
        ylab =  "Kocak (NES)", col =  "#00000099", xlim =  c(-50, 50), ylim =  c(-50, 60))
shadowtext(target[genes], y, labels =  genes, col =  "white")
scatter(nrc, kocak, xlab =  "NRC (NES)", ylab =  "Kocak (NES)", col =  "#00000099", 
        xlim =  c(-50, 50), ylim =  c(-50, 60))
shadowtext(nrc[genes], y, labels =  genes, col =  "white")
dev.off()

png("plots/010_mraplot_target.png", w =  3000, h =  8000, res =  300)
mraplot(corto_tmra, mrs =  genes)
dev.off()

png("plots/010_mraplot_kocak.png", w =  3000, h =  8000, res =  300)
mraplot(corto_kmra, mrs =  genes)
dev.off()

png("plots/010_mraplot_nrc.png", w =  3000, h =  8000, res =  300)
mraplot(corto_nmra, mrs =  genes)
dev.off()

### MRA at a single-cell level
load("data/regmat.rda")
load("results/010_mra.rda")
load(paste0("networks/target_NBL-regulon.rda"))
load("data/rtsne.rda")
x <- setNames(rtsne$Y[, 1], colnames(regmat))
y <- setNames(rtsne$Y[, 2], colnames(regmat))

### Single-cell master regulator analysis ----
fname <- "results/011_scmra.rda"
if(!file.exists(fname)){
    scmra <- mra(regmat, regulon =  regulon)
    save(scmra, file =  fname)
}else{load(fname)}

### Show top MRs ----
genes <- c("MYCN", "DNAJC1", "TWIST1", "NOTCH1", "E2F3", "TEAD4")
for(gene in genes){
    vector <- scmra[gene, ]
    
    ## Color by TPM
    colfunc <- colorRampPalette(c("navy", "navy", "navy", "grey", "salmon", "red", "red3"))
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
 
    # Plot
    png(paste0("plots/011_tsne_mra_", gene, ".png"), w =  4500, h =  2000, res =  500)
    lmatrix <- t(c(1, 2, 2, 2, 3))
    layout(lmatrix)
    par(mar =  c(5.1, 4.1, 4.1, 2.1))
    hist(vector, main =  paste0(gene, ""), xlab =  "NES")
    mycols <- colfunc(1000)[as.numeric(cut(vector, breaks =  1000))]
    plot(x, y, pch =  20, xlab =  "TSNE1", ylab =  "TSNE2", main =  paste0(gene), col =  mycols, xlim =  c(min(x), max(x)*1.5))
    mtext(paste0(length(x), " cells"), cex =  0.8, font =  2)
    par(mar =  c(5.1, 1.1, 4.1, 2.1))
    plot(c(0, 2), c(0, 1), type =  'n', axes =  F, xlab =  '', ylab =  '', main =  paste0("NES"))
    legend_image <- as.raster(rev(matrix(colfunc(1000),  ncol =  1)))
    rasterImage(legend_image, 0, 0, 1, 1)
    text(x =  1.5, y =  seq(0, 1, l =  5), labels =  round(seq(min(vector), max(vector), length.out =  5), 2))
    dev.off()
}

### Heatmap of single cell MR activity ----
source("../shared/functions/heatmaps.R")

## Select genes
# Manually selected
genes <- c("MYCN", "E2F3", "TEAD4", "E2F1", "MYC", "E2F7")
# Most differentially active
newgenes <- c(
    "DNAJC1", "ETV4", "HEYL", "HINFP", "MBD3", "NFRKB", "NPAT", "SCYL1", "TAF10", "TAF6", "TWIST1", "ZCCHC24", "ZNF25", "ZNHIT1", 
    "SESN2", "TRIM28", "UXT", "ZNF581"
)
genes <- c(genes, newgenes)
# Highest sum of absolute NES
sums <- sort(apply(scmra, 1, function(x){sum(abs(x))}), dec =  TRUE)
newgenes <- names(sums)[1:5] # "ZNF429" "RNF10"  "ZNF264" "ZBTB43" "MAX"
genes <- c(genes, newgenes)
#
genes <- unique(genes)

## Prepare table
toshow <- scmra[genes, ]
toshow[toshow>5] <- 5
toshow[toshow<(-5)] <- -5
colside <- rep(col_be2c, ncol(scmra))
names(colside) <- colnames(scmra)
colside[grep("kelly", names(colside))] <- col_kelly

# Color function
colfun <- colorRampPalette(c("navy", "navy", "blue", "blue", "white", "red", "red", "red3", "red3"))


png("plots/011_heatmap.png", w =  4000, h =  3000, res =  300)
heatmap.3(toshow, mar =  c(0, 10), ColSideColors =  t(t(colside)), KeyValueName =  "NES", col =  colfun)
dev.off()

# Load Seurat object (already LogNormalized)
load("data/seuset.rda")
## Further processing
seuset <- FindVariableFeatures(seuset, selection.method =  "vst", nfeatures =  Inf)
all.genes <- rownames(seuset)
seuset <- ScaleData(seuset, features =  all.genes)
expmat <- as.matrix(seuset[["RNA"]]@scale.data)
dim(expmat) # 15782 genes,  2067 cells

## PCA is performed considering most variable features,  not the entire rows set
seuset <- RunPCA(seuset,  features  =   VariableFeatures(object  =   seuset))
ElbowPlot(seuset) # choose the elbow,  try +1 and -1

## Determine the optimal number of significant PCs

## Determine the optimal number of significant PCs for subsequent clustering
## Determine percent of variation associated with each PC
pct <- seuset@reductions$pca@stdev / sum(seuset@reductions$pca@stdev) * 100
## Calculate cumulative percents for each PC
cum <- cumsum(pct)
## Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum>90&pct<5)[1]
## Determine the difference between variation of PC and subsequent PC
co2  <-  sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),   decreasing  =   T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2 # 13
## Minimum of the two calculation
pcs  <-  min(co1,  co2) # change to any other number
pcs # 13
seuset <-  RunUMAP(seuset,  dims  =   1:pcs)
seuset <- FindNeighbors(object =  seuset, dims =  1:pcs, reduction =  'pca')

## Assign cells to clusters
seuset <- FindClusters(object =  seuset, resolution =  0.01)
table(seuset@meta.data$seurat_clusters)
png("plots/013b_umap_res001.png", w =  1700, h =  1500, res =  300)
DimPlot(seuset, cols =  c("cornflowerblue", "salmon"))+ggtitle("Resolution =  0.01")
dev.off()

fname <- "results/013_markers_BE2C_vs_Kelly.rda"
if(!file.exists(fname)){
  markers <- FindMarkers(seuset, 1, 0) # min.pct =  0.1,  the default,  will test only genes detected in at least 10% of either cell population
  save(markers, file =  fname)
}else{load(fname)}
pander(markers[1:30, ], style =  "rmarkdown")

### Try a resolution a little big higher
seuset  <-  FindClusters(object  =   seuset,  resolution  =   0.2)
table(seuset@meta.data$seurat_clusters)
png("plots/013b_umap_res02.png", w =  1700, h =  1500, res =  300)
DimPlot(seuset, cols =  c("cornflowerblue", "red", "red3"))+ggtitle("Resolution =  0.2")
dev.off()

fname <- "results/013_markers_BE2C_2_vs_BE2C_1.rda"
if(!file.exists(fname)){
  markers <- FindMarkers(seuset, 2, 1)
  save(markers, file =  fname)
}else{load(fname)}

fname <- "results/014_fclsvm_model.rda"
if(!file.exists(fname)){
  # Construct a SingleCellExperiment object
  load("data/rawcounts_symbols.rda")
  dim(rawcounts) # 22135 2067
  log2counts <- log2(rawcounts+0.001)
  vars <- sort(apply(log2counts, 1, var), dec =  TRUE)[1:5000]
  log2counts <- log2counts[names(vars), ]
  dim(log2counts) # 13265 2067
  sce <- SingleCellExperiment::SingleCellExperiment(assays =  list(logcounts =  log2counts))
  
  # We will supply f-scLVM with genesets in a GeneSetCollection object.
  # We will use a curated annotation as suggested in the f-scLCM vignette
  # The following file was downloaded from MSigDB on Jan 20,  2021
  gmtfile <- "data/c2.cp.wikipathways.v7.2.symbols.gmt"
  genesets <- GSEABase::getGmt(gmtfile)
  
  # Generate a f-scLVM model
  model <- newSlalomModel(sce, genesets)
  # 194 annotated factors retained;  393 annotated factors dropped.
  # 1616  genes retained for analysis.
  # Initialize it
  model <- initSlalom(model)
  # Train it
  model <- trainSlalom(model, nIterations =  10000) # train the model until it converges
  save(model, file =  fname)
} else {load(fname)}
