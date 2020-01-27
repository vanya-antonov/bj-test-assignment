
library(readxl)
library(dplyr)
library(tibble)

library(ggplot2)

library(ComplexHeatmap)


###

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'

###

melanoma_es <- read.delim(paste0(DATA_DIR, 'genesets_2016.GSVA_scores.txt'), row.names = 1) %>%
  as.matrix()

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


# check
all(colnames(melanoma_es) %in% rownames(pt_info))
pt_subset <- pt_info[colnames(melanoma_es), ]
R_pt  <- pt_subset %>% filter(Response == "R")  %>% pull('ID')
NR_pt <- pt_subset %>% filter(Response == "NR") %>% pull('ID')

pdf(paste0(OUT_DIR, 'NR_vs_R_2016.heatmap.pdf'), width = 14, height = 8)
ht_NR <- Heatmap(melanoma_es[, NR_pt],
                 heatmap_legend_param = list(
                   title_position = "topcenter", color_bar = "continuous", legend_direction = "horizontal",
                   title = "GSVA score"),
                 width = unit(6, "cm"),
                 column_title = sprintf("%d NR patients", length(NR_pt)))
ht_R  <- Heatmap(melanoma_es[, R_pt],
                 column_title = sprintf("%d R patients", length(R_pt)),
                 width = unit(6, "cm"),
                 show_heatmap_legend = FALSE)
draw(ht_NR + ht_R,
     heatmap_legend_side = "bottom",
     row_title = sprintf('%d gene sets', nrow(melanoma_es)))
dev.off()




gg_df <- data.frame(
  patient = c(NR_pt, R_pt),
  score = c(
    apply(melanoma_es[, NR_pt], 2, mean),
    apply(melanoma_es[, R_pt], 2, mean)),
  type = c(rep('NR', length(NR_pt)),
           rep('R',  length(R_pt)))
)


THR <- 0.12
gg_df %>%
  arrange(-score) %>%
  mutate(patient = factor(patient, levels = patient)) %>%
  ggplot() +
  aes(x= patient, y = score, fill = type) +
  geom_bar(stat = 'identity') +
  geom_hline(aes(yintercept=THR)) +
  geom_text(aes(nrow(gg_df)-2, THR, label = THR, vjust = -1))
  ylab('Mean GSVA score')
ggsave('NR_vs_R_2016.barplot.pdf', path = OUT_DIR)  
