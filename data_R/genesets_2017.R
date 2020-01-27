# BiocManager::install("GSVA")
# BiocManager::install("GSVAdata")
# BiocManager::install("limma")
# install.packages("readxl") # CRAN version

# https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.R

library(readxl)
library(dplyr)
library(tibble)


#library(limma)

#library(ComplexHeatmap)


###

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'

###

genesets_2016 <- read.delim(paste0(DATA_DIR, 'genesets_2016.GSVA_scores.txt'), as.is = TRUE)$geneset

fpkm_df <- read_excel(paste0(DATA_DIR, '2017 GSE96619_PatientFPKM.xlsx'),
                      sheet = "PD1cohort.OT.fpkm.Feb17.2016.tx") %>%
  select('Gene', ends_with('-baseline')) %>%
  # Remove suffix: https://stackoverflow.com/a/45960434/310453
  rename_at(.vars = vars(ends_with("-baseline")),
            .funs = funs(sub("-baseline$", "", .)))

library(GSVA)
library(GSVAdata)

data(c2BroadSets)
c2BroadSets
length(names(c2BroadSets))


# https://www.biostars.org/p/242157/
detach("package:dplyr", unload=TRUE)
library(org.Hs.eg.db) # I used this database because I assumed we are dealing with Homo sapiens
anno_df <- select(org.Hs.eg.db, fpkm_df$Gene, "ENTREZID", "SYMBOL")
head(anno_df)

library(dplyr)
melanoma_es <- inner_join(anno_df, fpkm_df, by = c('SYMBOL' = 'Gene')) %>%
  filter(!is.na(ENTREZID)) %>%
  select(-SYMBOL) %>%
  column_to_rownames('ENTREZID') %>%
  as.matrix() %>%
  gsva(c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)
dim(melanoma_es)

melanoma_es[genesets_2016,] %>%
  as.data.frame() %>%
  rownames_to_column('geneset') %>%
  write.table(file = paste0(DATA_DIR, "genesets_2017.GSVA_scores.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")

