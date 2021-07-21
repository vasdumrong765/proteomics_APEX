
#########################
## Dumrongprechachan et al 2021 Nature Comms
## dSPNs: NES-H2B-LCK comparison
#########################

## start up library
library(tidyverse)
library(MSstatsTMT)
library(data.table)
#library(ggpubr)
library(ggrepel)
#library(gplots)
library(edgeR)
library(psych)
library(robustbase)
library(expss)
#library(limma)
library(readxl)
#library(RColorBrewer)
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/med_SL_norm.R")
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/plot_label.R")
library(ComplexHeatmap)

#########################
## Load data (shotgun)
#########################

# Read PSMs PD output into R
# S/N export from PD
# Note that if you use read_excel column names will be wrong (correct e.g., `Spectrum.File` not `Spectrum File`)
raw.pd_shotgun <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/04142021_AAVAPEX_unfrac_PSMs.txt", sep="\t", header=TRUE) %>% 
  dplyr::filter(Isolation.Interference.... < 70) %>%
  dplyr::filter(Quan.Info == '') %>%
  dplyr::filter(Percolator.q.Value < 0.01)

# Read PD protein tabs into R
raw.protein.pd_shotgun <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/04142021_AAVAPEX_unfrac_proteins.xlsx")

# Read annotation file
annotation.pd_shotgun <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/PD_set1_annotation.txt", sep="\t", header=TRUE)


#########################
## Load data (high pH reverse phase fractionation)
#########################

raw.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/04142021_AAVAPEX_frac_PSMs.txt", sep="\t", header=TRUE) %>% 
  dplyr::filter(Isolation.Interference.... < 70) %>%
  dplyr::filter(Quan.Info == '') %>%
  dplyr::filter(Percolator.q.Value < 0.01)

raw.protein.pd_frac <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/04142021_AAVAPEX_frac_proteins.xlsx") 

annotation.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set2/PD_set2_annotation.txt", sep="\t", header=TRUE)

#hist(raw.pd_frac$Ion.Inject.Time..ms., breaks=40)

#########################
## MSstatsTMT PSM selection
#########################

# Converting PD output to MSstats format
# Remove Proteins with 1 Feature
# Remember that MSstatsTMT only use unique peptides for protein summarization and quantifications
# Shared peptides are not currently implemented.
# Removed shared peptides

input.pd_shotgun <- PDtoMSstatsTMTFormat(raw.pd_shotgun,
                                      annotation.pd_shotgun,
                                      which.proteinid = 'Master.Protein.Accessions',
                                      useNumProteinsColumn = FALSE,
                                      useUniquePeptide = TRUE,
                                      rmPSM_withMissing_withinRun = FALSE,
                                      rmPSM_withfewMea_withinRun = TRUE, #features with 1-2 measurements
                                      rmProtein_with1Feature = TRUE,
                                      summaryforMultipleRows = sum) %>% filter(!grepl(";", ProteinName))

input.pd_frac <- PDtoMSstatsTMTFormat(raw.pd_frac,
                                 annotation.pd_frac,
                                 which.proteinid = 'Master.Protein.Accessions',
                                 useNumProteinsColumn = FALSE,
                                 useUniquePeptide = TRUE,
                                  rmPSM_withMissing_withinRun = FALSE,
                                 rmPSM_withfewMea_withinRun = TRUE, #features with 1-2 measurements
                                 rmProtein_with1Feature = TRUE,
                                 summaryforMultipleRows = sum) %>% filter(!grepl(";", ProteinName))



##################################################
## MSstatsTMT protein summarization
##################################################

quant.msstats0_shotgun <- proteinSummarization(input.pd_shotgun,
                                       method="msstats",
                                       global_norm=FALSE,
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

quant.msstats0_frac <- proteinSummarization(input.pd_frac,
                                       method="msstats",
                                       global_norm=FALSE,
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

# violin plot --- Fig2
quant.msstats0_frac$Condition <- factor(quant.msstats0_frac$Condition,
                                        levels = c('CTRL', 'H2B', 'NES', 'LCK', 'Norm'))
# set some colors by condition
my_colors <- c('dodgerblue3', 'darkolivegreen4', 'gray60', 'hotpink3','coral3')

protein_violin <- ggplot(data =  quant.msstats0_frac, 
                         mapping = aes(x = BioReplicate, y = Abundance)) + 
  geom_violin(trim = TRUE,
              mapping = aes(fill = Condition)) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = NULL) + 
  stat_summary(fun=median, geom="crossbar", size=0.2, color="black") + 
  #facet_wrap(~ Condition, scale = "free", ncol = 6) +
  #scale_y_continuous(limits = c(-4,12)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_blank()) +
  ggtitle("MSstats Protein Abundance")
protein_violin


##################################################
## Wide data format
##################################################

# wide format protein abundance and set negative abundance to NA
quant.msstats0_shotgun_wide <- quant.msstats0_shotgun %>% dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

quant.msstats0_frac_wide <- quant.msstats0_frac %>% dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

# replace 0 with NA
quant.msstats0_shotgun_wide[quant.msstats0_shotgun_wide < 0] <- NA
quant.msstats0_frac_wide[quant.msstats0_frac_wide < 0] <- NA

# remove NA from the clustering analysis
quant.msstats0_shotgun_temp <- quant.msstats0_shotgun_wide[-1] %>% drop_na()
quant.msstats0_frac_temp <- quant.msstats0_frac_wide[-1] %>% drop_na()

# re-arrange columns
col_order <- c('Norm_1', 'Norm_2', 'Norm_3', 'Norm_4',
               'CTRL_1', 'CTRL_2', 'CTRL_3', 'CTRL_4', 
               'H2B_1', 'H2B_2','H2B_3', 'H2B_4',
               'NES_1', 'NES_2', 'NES_3', 'NES_4',
               'LCK_1', 'LCK_2', 'LCK_3', 'LCK_4')

quant.msstats0_shotgun_temp <- quant.msstats0_shotgun_temp[, col_order]
quant.msstats0_frac_temp <- quant.msstats0_frac_temp[, col_order]


##################################################
## EdgeR multidimensional clusting MDS plot and box plot
##################################################

# create a DGEList object
prot_matrix_shotgun <- DGEList(as.matrix(quant.msstats0_shotgun_temp))
prot_matrix_frac <- DGEList(as.matrix(quant.msstats0_frac_temp))

# set some colors by condition
colors = c(rep('coral3', 4), 
           rep('dodgerblue3', 4), 
           rep('darkolivegreen4',4), 
           rep('gray20',4), 
           rep('hotpink3',4))

# check the clustering
plotMDS(prot_matrix_shotgun, col = colors, main="single sample")
plotMDS(prot_matrix_frac, col = colors, main="fractionated samples")
plotMDS(prot_matrix_frac, pch = 19, col = colors, main="plotMDS - sample similarity")

# box plots
boxplot(quant.msstats0_shotgun_temp, col = colors)
boxplot(quant.msstats0_frac_temp, col = colors)

# plotMDS data dimension 1 boxplot
data_plotMDS <- plotMDS(prot_matrix_frac, pch = 19, col = colors, main="plotMDS - sample similarity")
data_plotMDS_dim1 <- data_plotMDS[['x']] %>% data.frame()
data_plotMDS_dim1$dim2 <- data_plotMDS[['y']]
data_plotMDS_dim1$condition <- c(rep('Norm',4), rep('CTRL', 4), rep('H2B', 4), rep('NES',4), rep('LCK',4))
colnames(data_plotMDS_dim1) <- c('dim1', 'dim2','Condition')

data_plotMDS_dim1_mean <- data_plotMDS_dim1 %>%
  group_by(Condition) %>% summarise(log2mean = mean(dim1))

data_plotMDS_dim1 %>%
  ggplot(aes(Condition, dim1)) +
  geom_boxplot() +
  geom_point(size=2, shape=19, position = position_jitter()) +
  #geom_boxplot(aes(color=factor(PSM))) +
  #facet_wrap(~PSM) +
  ylab('leading log2 FC dim1') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1))


##################################################
## Use fractionated dataset for the following analysis
##################################################

##################################################
## Filter 1 - Cre+ vs Cre- comparison using MSstatsTMT linear mixed model
##################################################

levels(droplevels(quant.msstats0_frac$Condition))

#create a comparison matrix
comparison1 <- matrix(c(-1, 1, 0, 0,
                        -1, 0, 1, 0,
                        -1, 0, 0, 1),
                      nrow=3, ncol=4, byrow=TRUE)
colnames(comparison1) <- c("CTRL","NES", "H2B", "LCK")
row.names(comparison1) <- c("NES-CTRL", "H2B-CTRL", "LCK-CTRL")
comparison1

#perform MSstats group comparison as indicated above
test.contrast1 <- groupComparisonTMT(data = quant.msstats0_frac,
                                     contrast.matrix = comparison1,
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE)


positive_enrich <- test.contrast1 %>% filter(adj.pvalue <= 0.05)
pos_enrich_NES <- positive_enrich %>% filter(Label == 'NES-CTRL') %>% 
  dplyr::select(Protein)
pos_enrich_LCK <- positive_enrich %>% filter(Label == 'LCK-CTRL') %>% 
  dplyr::select(Protein)
pos_enrich_H2B <- positive_enrich %>% filter(Label == 'H2B-CTRL') %>% 
  dplyr::select(Protein)


##################################################
## Filter 2.1 - H2B vs NES comparison
##################################################

H2B_comp2 <- quant.msstats0_frac %>% filter(Condition == 'H2B' | Condition == 'NES') %>%
  filter(Protein %in% pos_enrich_H2B$Protein)

test.contrast2_1 <- groupComparisonTMT(data = H2B_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

##################################################
## Filter 2.2 - LCK vs NES comparison
##################################################

LCK_comp2 <- quant.msstats0_frac %>% filter(Condition == 'LCK' | Condition == 'NES') %>%
  filter(Protein %in% pos_enrich_LCK$Protein)

test.contrast2_2 <- groupComparisonTMT(data = LCK_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

##################################################
## Filter 2.3 - H2B vs LCK comparison
##################################################

H2B_comp3 <- quant.msstats0_frac %>% filter(Condition == 'LCK' | Condition == 'H2B') %>%
  #filter(Protein %in% pos_enrich_LCK$Protein) %>%
  filter(Protein %in% pos_enrich_H2B$Protein)

test.contrast2_3 <- groupComparisonTMT(data = H2B_comp3,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)


rm(H2B_comp2)
rm(LCK_comp2)
rm(H2B_comp3)
rm(pos_enrich_NES)
rm(pos_enrich_LCK)
rm(pos_enrich_H2B)


##################################################
## Use PD output to map protein IDs
##################################################

## find proteins that are quantified by MSstatsTMT
msstatsProt <- intersect(quant.msstats0_frac_wide$Protein, raw.protein.pd_frac$Accession)

## filter rows in pd.prot.raw
## select only useful columns
pd.prot.raw <- dplyr::filter(raw.protein.pd_frac, Accession %in% msstatsProt) %>%
  dplyr::select(
    Accession,
    `Protein FDR Confidence: Combined`,
    Description,
    `Coverage [%]`,
    `# Peptides`,
    `# Unique Peptides`,
    `# PSMs`,
    `MW [kDa]`,
    `Cellular Component`,
    `Gene Symbol`
  )


##################################################
## Merge mapped_ID to test.contrast dataframes
##################################################

mapped_ID <-  supp_table2 %>%
  dplyr::select(Accession:`Gene Symbol`)

## for H2B-NES comparison
test.contrast2_1_mapped <- merge(x=test.contrast2_1, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")
## for LCK-NES comparison
test.contrast2_2_mapped <- merge(x=test.contrast2_2, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")

## for H2B-LCK comparison
test.contrast2_3_mapped <- merge(x=test.contrast2_3, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")


##################################################
## H2B-APEX volcano plot
##################################################

test.contrast2_1_mapped <- plot_label(test.contrast2_1_mapped,
                                      FCvalue = 3,
                                      FDRsig = 0.05,
                                      FDRplot = 0.05)

# H2B-Cre(neg) volcano plot
cre_pos_H2B <- test.contrast1 %>% dplyr::filter(Label == 'H2B-CTRL')
for (i in seq_along(cre_pos_H2B$Protein)){
  cre_pos_H2B$sig[i] <- ifelse(cre_pos_H2B$adj.pvalue[i] < 0.05 &
                                 cre_pos_H2B$log2FC[i] >= 0, 1, 0)}

volcano_df1 <- ggplot(data = cre_pos_H2B,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(sig)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df1

# H2B-NES volcano plot
volcano_df3 <- ggplot(data = test.contrast2_1_mapped,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2","black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  xlim(-3,4) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df3


##################################################
## H2B-APEX rank plot
##################################################

### Rank plot
test.contrast2_ranked_H2B <- rank_label(test.contrast2_1_mapped,
                                   GOterm = 'nucleus',
                                   FDRsig = 0.05,
                                   FDRcutoff = 'pos',
                                   skip = 4)

rank_plot_H2B <- ggplot(data = test.contrast2_ranked_H2B,
                    mapping = aes(x = Rank, y = log2FC)) + 
  geom_point(aes(colour=factor(GOterm)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("gray50", "coral2")) +
  #geom_hline(yintercept = log2(1.5), 
  #                         linetype = "dotted", size = 0.6)
  geom_text_repel(aes(label = Protein_label_select2), 
                  max.overlaps = Inf,
                  force = 3, #direction = "y",
                  #nudge_x = 0.2, nudge_y = 0.2,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.2, #segment line thickness
                  size=4,
                  fontface = 'italic') + #font size
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

rank_plot_H2B


##################################################
## LCK-APEX comparison with Label
##################################################

# LCK-Cre(neg) volcano plot
cre_pos_LCK <- test.contrast1 %>% dplyr::filter(Label == 'LCK-CTRL')
for (i in seq_along(cre_pos_LCK$Protein)){
  # adj.pvalue < 0.05
  cre_pos_LCK$sig[i] <- ifelse(cre_pos_LCK$adj.pvalue[i] < 0.05 &
                                 cre_pos_LCK$log2FC[i] >= 0, 1, 0)}

volcano_df1_LCK <- ggplot(data = cre_pos_LCK,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(sig)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df1_LCK


test.contrast2_3_mapped <- plot_label(test.contrast2_3_mapped,
                                      FCvalue = 3,
                                      FDRsig = 0.05,
                                      FDRplot = 0.05)

# H2B-LCK volcano plot
volcano_df3_LCK <- ggplot(data = test.contrast2_3_mapped,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2","black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  xlim(-3,4) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df3_LCK


##################################################
## LCK-APEX ranked plot
##################################################

## for LCK rank plot
## import annotated version
test.contrast2_ranked_LCK <- rank_label(test.contrast2_3_mapped,
                                    GOterm = 'membrane',
                                    FDRsig = 0.05,
                                    FDRcutoff = 'neg',
                                    skip = 4)

rank_plot_LCK <- ggplot(data = test.contrast2_ranked_LCK,
                    mapping = aes(x = Rank, y = log2FC)) + 
  geom_point(aes(colour=factor(GOterm)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("gray50", "dodgerblue2")) +
  #geom_hline(yintercept = log2(1.5), 
  #                         linetype = "dotted", size = 0.6)
  geom_text_repel(aes(label = Protein_label), 
                  force = 4, #direction = "y",
                  max.overlaps = Inf,
                  #nudge_x = 0.2, nudge_y = 0.2,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.2, #segment line thickness
                  size=4,
                  fontface = 'italic') + #font size
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

rank_plot_LCK


##################################################
## Remarks: without reference channel normalization
##################################################

quant.msstats0_frac_noBATCH <- proteinSummarization(input.pd_frac,
                                                    method="msstats",
                                                    global_norm=FALSE,
                                                    reference_norm=FALSE,
                                                    remove_norm_channel = FALSE,
                                                    remove_empty_channel = TRUE,                  
                                                    MBimpute = FALSE,
                                                    maxQuantileforCensored = NULL)

quant.msstats0_frac_noBATCH_wide <- quant.msstats0_frac_noBATCH %>% 
  dplyr::select(Protein, Abundance, BioReplicate) %>% 
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

# replace 0 with NA
quant.msstats0_frac_noBATCH_wide[quant.msstats0_frac_noBATCH_wide < 0] <- NA

# remove NA from the clustering analysis
quant.msstats0_frac_noBATCH_temp <- quant.msstats0_frac_noBATCH_wide[-1] %>% drop_na()
quant.msstats0_frac_noBATCH_temp <- quant.msstats0_frac_noBATCH_temp[, col_order]

# create a DGEList object
prot_matrix_frac_noBATCH <- DGEList(as.matrix(quant.msstats0_frac_noBATCH_temp))
plotMDS(prot_matrix_frac_noBATCH, pch = 19, col = colors, main="fractionated samples no Batch effect correction")


##################################################
## Remarks: protein overlap between shotgun and fractionation
##################################################

# use setdiff to find non-overlap elements in the first list
prot_setdiff_shotgun <- setdiff(raw.protein.pd_shotgun$Accession, raw.protein.pd_frac$Accession)
prot_setdiff_frac <- setdiff(raw.protein.pd_frac$Accession, raw.protein.pd_shotgun$Accession)

# use intersect to find the intersected elements
prot_intersect <- intersect(raw.protein.pd_frac$Accession, raw.protein.pd_shotgun$Accession)


##################################################
## Remarks: Cre-negative plots Fig S6
##################################################

# LCK-Cre(neg) volcano plot
cre_pos_LCK <- test.contrast1 %>% dplyr::filter(Label == 'LCK-CTRL')
for (i in seq_along(cre_pos_LCK$Protein)){
  # adj.pvalue < 0.05
  cre_pos_LCK$sig[i] <- ifelse(cre_pos_LCK$adj.pvalue[i] < 0.05 &
                                 cre_pos_LCK$log2FC[i] > 0, 1, 0)}

volcano_df1_LCK <- ggplot(data = cre_pos_LCK,
                          aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(sig)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df1_LCK

# NES-Cre(neg) volcano plot
cre_pos_NES <- test.contrast1 %>% dplyr::filter(Label == 'NES-CTRL')
for (i in seq_along(cre_pos_NES$Protein)){
  # adj.pvalue < 0.05
  cre_pos_NES$sig[i] <- ifelse(cre_pos_NES$adj.pvalue[i] < 0.05 &
                                 cre_pos_NES$log2FC[i] > 0, 1, 0)}

volcano_df1_NES <- ggplot(data = cre_pos_NES,
                          aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(sig)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df1_NES


#### LCK-NES comparison volcano plot
# NES-LCK volcano plot
for (i in seq_along(test.contrast2_2_mapped$Protein)){
    #pvalue < 0.05
    test.contrast2_2_mapped$sig[i] <- ifelse(test.contrast2_2_mapped$pvalue[i] < 0.05, 1, 0)
    # get Gene.Symbol for significantly enriched protein.
    test.contrast2_2_mapped$Protein_label[i] <- 
      ifelse(test.contrast2_2_mapped$pvalue[i] < 0.05,
             test.contrast2_2_mapped$`Gene Symbol`[i], "")
    
    # to only display certain foldchange (both + and -)
    test.contrast2_2_mapped$Protein_label_select[i] <- 
      ifelse(
        (test.contrast2_2_mapped$pvalue[i] < 0.015 & test.contrast2_2_mapped$log2FC[i] > log2(1.5)) |
        (test.contrast2_2_mapped$pvalue[i] < 0.1 & test.contrast2_2_mapped$log2FC[i] < log2(1.2)),
             #abs(test.contrast2_2_mapped$log2FC)[i] > log2(1.5), #both +/- FC
             test.contrast2_2_mapped$`Protein_label`[i], "")
}

volcano_df2_NESLCK <- ggplot(data = test.contrast2_2_mapped,
                      aes(x = log2FC, y = (-1)*log10(pvalue))) +
  geom_point(size = 1, show.legend = FALSE) + 
  #scale_colour_manual(values = c("dodgerblue2","black", "coral2")) +
  #geom_hline(yintercept = (-1)*log10(0.05), 
  #           linetype = "dotted", size = 0.6) +
  #xlim(-3,4) +
  geom_text_repel(aes(label = Protein_label_select), 
                  max.overlaps = Inf,
                  force = 5, #direction = "y",
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4, fontface = 'italic') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df2_NESLCK


### annotation overlap
cre_pos_LCK_mapped <- merge(x=cre_pos_LCK, y=mapped_ID, 
                            by.x="Protein", by.y = "Accession")

cre_pos_LCK_mapped_mem <- cre_pos_LCK_mapped %>%
  filter(grepl('membrane',`Cellular Component`))

cre_pos_LCK_mapped_cyto <- cre_pos_LCK_mapped %>%
  filter(grepl('cytosol',`Cellular Component`) | grepl('cytoplasm',`Cellular Component`))

length(intersect(cre_pos_LCK_mapped_mem$Protein, cre_pos_LCK_mapped_cyto$Protein))
length(setdiff(cre_pos_LCK_mapped_mem$Protein, cre_pos_LCK_mapped_cyto$Protein))
length(setdiff(cre_pos_LCK_mapped_cyto$Protein, cre_pos_LCK_mapped_mem$Protein))


### rank plot with membrane annotation
cre_pos_LCK_rank <- test.contrast2_2_mapped
cre_pos_LCK_rank$mem <- ifelse(grepl('membrane',
                                     cre_pos_LCK_rank$`Cellular Component`, 
                                     ignore.case = TRUE),1, 0)
# rank based on log2FC
cre_pos_LCK_rank <- cre_pos_LCK_rank %>% arrange(log2FC)
  #arrange(desc(log2FC))
cre_pos_LCK_rank$Rank <- seq_along(cre_pos_LCK_rank$Protein)

rank_plot_LCKNES <- ggplot(data = cre_pos_LCK_rank,
                        mapping = aes(x = Rank, y = log2FC)) + 
  geom_point(aes(colour=factor(mem)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("gray50", "lightblue1")) +
  #geom_hline(yintercept = log2(1), linetype = "dotted", size = 0.6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
rank_plot_LCKNES