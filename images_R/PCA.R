
# https://www.r-bloggers.com/principal-component-analysis-in-r/
# wine <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data", sep=",")
# head(wine)
# 
# # Name the variables
# colnames(wine) <- c("Cvs","Alcohol","Malic acid","Ash","Alcalinity of ash", "Magnesium", "Total phenols", "Flavanoids", "Nonflavanoid phenols", "Proanthocyanins", "Color intensity", "Hue", "OD280/OD315 of diluted wines", "Proline")
# 
# wineClasses <- factor(wine$Cvs)
# winePCA <- prcomp(scale(wine[,-1]))
# plot(winePCA$x[,1:2], col = wineClasses)

library(readxl)
library(dplyr)
library(tibble)

library(ggplot2)

###

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'

###
# 2016

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
# head(as.data.frame(fpkm_df))

# head(t(fpkm_df[,-1]))
mtx <- t(fpkm_df[,-1])
# remove rows with identical values https://stackoverflow.com/a/40317343/310453
mtx <- mtx[ , apply(mtx, 2, var) != 0]

ptPCA <- prcomp(scale(mtx))

gg_df_2016 <- as.data.frame(ptPCA$x[,1:2]) %>%
  rownames_to_column(var = 'ID') %>%
  left_join(pt_info, by = 'ID') %>%
  mutate(ID = paste0('2016-', ID)) %>%
  select(ID, PC1, PC2, Response)


ggplot(gg_df_2016) +
  aes(x = PC1, y = PC2, col = Response) +
  geom_point() +
  # https://stackoverflow.com/a/15625149/310453
  geom_text(aes(label=ifelse(PC2>100, as.character(ID), '')), hjust=-0.2, vjust=0) +
  ggtitle('Hugo et al, 2016') +
  theme_bw() +
  theme(legend.position="bottom")
ggsave('PCA.2016.pdf', path = OUT_DIR, width = 4, height = 5)


###
# 2017

fpkm_df_2017 <- read_excel(paste0(DATA_DIR, '2017 GSE96619_PatientFPKM.xlsx'),
                           sheet = "PD1cohort.OT.fpkm.Feb17.2016.tx") %>%
  select('Gene', ends_with('-baseline')) %>%
  # Remove suffix: https://stackoverflow.com/a/45960434/310453
  rename_at(.vars = vars(ends_with("-baseline")),
            .funs = funs(sub("-baseline$", "", .)))

pt_info_2017 <- data.frame(
  ID = paste0("Pt", 1:5),
  Response = c(rep("R", 2), rep("NR", 3)),
  stringsAsFactors = FALSE)
rownames(pt_info) <- pt_info$ID

# head(t(fpkm_df[,-1]))
mtx_2017 <- t(fpkm_df_2017[,-1])
# remove rows with identical values https://stackoverflow.com/a/40317343/310453
mtx_2017 <- mtx_2017[ , apply(mtx_2017, 2, var) != 0]

ptPCA_2017 <- prcomp(scale(mtx_2017))

gg_df_2017 <- as.data.frame(ptPCA_2017$x[,1:2]) %>%
  rownames_to_column(var = 'ID') %>%
  left_join(pt_info_2017, by = 'ID') %>%
  mutate(ID = paste0('2017-', ID)) %>%
  select(ID, PC1, PC2, Response)

ggplot(gg_df_2017) +
  aes(x = PC1, y = PC2, col = Response) +
  geom_point() +
  # https://stackoverflow.com/a/15625149/310453
  geom_text(aes(label=ifelse(PC2>100, as.character(ID), '')), hjust=-0.2, vjust=0) +
  ggtitle('Garcia-Diaz et al, 2017') +
  theme_bw() +
  theme(legend.position="bottom")
ggsave('PCA.2017.pdf', path = OUT_DIR, width = 4, height = 5)




###
# 2016 and 2017

# Rename columns: Pt1 => 2016_Pt1
df_16 <- rename_all(fpkm_df, function(x) gsub("Pt", "y2016_Pt", x))
df_17 <- rename_all(fpkm_df_2017, function(x) gsub("Pt", "y2017_Pt", x))

fpkm_df_16_17 <- inner_join(df_16, df_17, by = 'Gene')
head(fpkm_df_16_17)


# head(t(fpkm_df[,-1]))
mtx_16_17 <- t(fpkm_df_16_17[,-1])
# remove rows with identical values https://stackoverflow.com/a/40317343/310453
mtx_16_17 <- mtx_16_17[ , apply(mtx_16_17, 2, var) != 0]

ptPCA_16_17 <- prcomp(scale(mtx_16_17))

gg_df_2017 <- as.data.frame(ptPCA_2017$x[,1:2]) %>%
  rownames_to_column(var = 'ID') %>%
  left_join(pt_info_2017, by = 'ID') %>%
  mutate(ID = paste0('2017-', ID)) %>%
  select(ID, PC1, PC2, Response)

info_16 <- mutate(pt_info,      ID = paste0("y2016_", ID)) %>%
  mutate(Article = 'Hugo et al, 2016') %>%
  select(ID, Article, Response)
info_17 <- mutate(pt_info_2017, ID = paste0("y2017_", ID)) %>%
  mutate(Article = 'Garcia-Diaz et al, 2017') %>%
  select(ID, Article, Response)
info_16_17 <- rbind(info_16, info_17)

gg_df_16_17 <- as.data.frame(ptPCA_16_17$x[,1:2]) %>%
  rownames_to_column(var = 'ID') %>%
  left_join(info_16_17, by = 'ID')

ggplot(gg_df_16_17) +
  aes(x = PC1, y = PC2, col = Response, shape = Article, size = Article) +
  geom_point() +
  ggtitle('Samples from two studies together') +
  # http://www.sthda.com/english/wiki/ggplot2-point-shapes
  scale_shape_manual(values=c(2, 20)) +
  scale_size_manual(values=c(4, 2)) +
  theme_bw()
ggsave('PCA.16_17.pdf', path = OUT_DIR)

