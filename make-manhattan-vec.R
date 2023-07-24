## Make example Manhattan plot and save as vector graphics

## Load data.table to read summary stats file
library(data.table)

## Load plotting function
source("fn-ggmanh_vec.R")



## Download UK Biobank height GWAS in European population (Watanabe et al. 2019)
system("wget https://atlas.ctglab.nl/ukb2_sumstats/f.50.0.0_res.EUR.sumstats.MACfilt.txt.gz")

## Read summary stats
gwas <- fread("f.50.0.0_res.EUR.sumstats.MACfilt.txt.gz", data.table = FALSE)

## -log10 transform p-values
gwas <- gwas[gwas$P > 0,]  # Remove p-values with value 0 before log-transformation
gwas$LOG10_P <- -log10(gwas$P)



## Make vector Manhattan plot and save as pdf
manh_vec <- ggmanh_vec(gwas,
                       chr = "CHR", pos = "BP", log_pv = "LOG10_P",
                       title = "GWAS of standing height in the UK Biobank (EUR population)",
                       ncores = 4)
ggsave(manh_vec, filename = paste0("ukb-height-gwas.pdf"),
       device = cairo_pdf, width = 14, height = 7, units = "cm")


## Make quick version (only plot SNPs with p-value <= 0.01)
manh_vec_pv2 <- ggmanh_vec(gwas[gwas$LOG10_P >= 2,],
                           chr = "CHR", pos = "BP", log_pv = "LOG10_P",
                           title = "GWAS of standing height in the UK Biobank (EUR population)",
                           ncores = 4)
ggsave(manh_vec_pv2, filename = paste0("ukb-height-gwas-pv2.pdf"),
       device = cairo_pdf, width = 14, height = 7, units = "cm")
