sample="Kelly"
cellranger mat2csv ${sample}/filtered_gene_bc_matrices_h5.h5 ${sample}.csv
sample="BE2C"
cellranger mat2csv ${sample}/filtered_gene_bc_matrices_h5.h5 ${sample}.csv
gzip *csv
