# BiocManager::install("GSVA")
# BiocManager::install("GSVAdata")
# BiocManager::install("limma")
# install.packages("readxl") # CRAN version

# https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.R

library(readxl)
library(dplyr)
library(tibble)

library(GSVA)
library(GSVAdata)

library(limma)

library(ComplexHeatmap)


###

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'

###

# https://stackoverflow.com/a/29138966/310453
pt_info <- read_excel(paste0(DATA_DIR, '2016 Table S1.xls'), sheet = "S1B", skip = 2) %>%
  # 'Patient ID'  =>  'Patient_ID'
  rename_all(function(x) gsub(" ", "_", x)) %>%
  filter(!is.na(Patient_ID)) %>%
  mutate(ID = Patient_ID) %>%
  column_to_rownames('ID')

# Add 2 more rows
tmp_df <- pt_info[c("Pt27", "Pt27"), ]
rownames(tmp_df)  <-  c("Pt27A", "Pt27B")
pt_info <- rbind(pt_info, tmp_df)
pt_info$ID <- rownames(pt_info)
tail(as.data.frame(pt_info))


fpkm_df <- read_excel(paste0(DATA_DIR, '2016 GSE78220_PatientFPKM.xlsx'), sheet = "FPKM")
# "Pt1.baseline"  =>  "Pt1"
colnames(fpkm_df) <- gsub('\\..+', '', colnames(fpkm_df))
head(as.data.frame(fpkm_df))


genesets_2016 <- read_excel(paste0(DATA_DIR, '2016 Table S2.xlsx'), sheet = "S2B", skip = 2)


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

# check
all(colnames(melanoma_es) %in% rownames(pt_info))

pt_subset <- pt_info[colnames(melanoma_es), ]

selected_2016 <- genesets_2016[genesets_2016$Geneset %in% rownames(melanoma_es), ]

melanoma_es[selected_2016$Geneset,] %>%
  as.data.frame() %>%
  rownames_to_column('geneset') %>%
  write.table(file = paste0(DATA_DIR, "genesets_2016.GSVA_scores.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")

