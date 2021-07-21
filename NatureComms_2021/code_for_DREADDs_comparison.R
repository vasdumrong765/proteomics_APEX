
#########################
## Dumrongprechachan et al 2021 Nature Comms
## DREADDs
#########################

## start up library
library(tidyverse)
library(MSstatsTMT)
library(data.table)
library(ggpubr)
library(ggrepel)
library(gplots)
library(edgeR)
library(psych)
#library(reader)
#library(GeneOverlap)
library(robustbase)
library(expss)
library(limma)
library(readxl)
library(RColorBrewer)
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/med_SL_norm.R")
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/plot_label.R")
library(ComplexHeatmap)

#########################
## Load data
#########################

# Read PSMs PD output into R
# S/N export from PD
# Note that if you use read_excel column names will be wrong (correct e.g., `Spectrum.File` not `Spectrum File`)
raw.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set3/04142021_AAVAPEX_DREADDs_PSMs.txt", sep="\t", header=TRUE) %>% 
  dplyr::filter(Isolation.Interference.... < 70) %>%
  dplyr::filter(Quan.Info == '') %>%
  dplyr::filter(Percolator.q.Value < 0.01)

# Read PD protein tabs into R
# for intensity values
raw.protein.pd_frac <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set3/04142021_AAVAPEX_DREADDs_proteins.xlsx")

# Read annotation file
annotation.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set3/PD_DREADD_annotation.txt", sep="\t", header=TRUE)

#########################
## MSstatsTMT PSM selection
#########################

# Converting PD output to MSstats format
# Remove Proteins with 1 Feature
# Remember that MSstatsTMT only use unique peptides for protein summarization and quantifications
# Shared peptides are not currently implemented.
input.pd_frac <- PDtoMSstatsTMTFormat(raw.pd_frac,
                                      annotation.pd_frac,
                                      which.proteinid = 'Master.Protein.Accessions',
                                      useNumProteinsColumn = FALSE,
                                      useUniquePeptide = TRUE,
                                      rmPSM_withMissing_withinRun = FALSE,
                                      rmPSM_withfewMea_withinRun = TRUE, #features with 1-2 measurements
                                      rmProtein_with1Feature = TRUE,
                                      summaryforMultipleRows = sum)

# Removed shared peptides
input.pd_frac <- input.pd_frac %>% filter(!grepl(";", ProteinName))

# preivew
# head(input.pd)


##################################################
## Use MSstats for protein summarization
##################################################

quant.msstats0 <- proteinSummarization(input.pd_frac,
                                       method="msstats",
                                       global_norm=FALSE,
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

# violin plot
quant.msstats0$Condition <- factor(quant.msstats0_frac$Condition,
                                        levels = c('control', 'NES', 'hM3Dq', 'mCherry', 'Norm'))
# set some colors by condition
my_colors <- c('dodgerblue3', 'darkolivegreen4', 'gray60', 'hotpink3','coral3')

protein_violin <- ggplot(data =  quant.msstats0, 
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
## EdgeR multidimensional clusting MDS plot
##################################################

# wide format protein abundance and set negative abundance to NA
quant.msstats0_wide <- quant.msstats0 %>% dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)
quant.msstats0_wide0 <- quant.msstats0_wide

quant.msstats0_wide[quant.msstats0_wide < 0] <- NA

# remove NA from the clustering analysis
quant.msstats0_wide_temp <- quant.msstats0_wide[-1] %>% drop_na()

# re-arrange columns
col_order <- c('Norm_1', 'Norm_2',
               'control_1', 'control_2', 'control_3',
               'NES_1', 'NES_2', 'NES_3',
               'hM3Dq_1', 'hM3Dq_2', 'hM3Dq_3', 'hM3Dq_4', 'hM3Dq_5', 'hM3Dq_6', 'hM3Dq_7',
               'mCherry_1', 'mCherry_2','mCherry_3', 'mCherry_4', 'mCherry_5', 'mCherry_6', 'mCherry_7')
quant.msstats0_wide_temp <- quant.msstats0_wide_temp[, col_order]

# create a DGEList object
prot_matrix <- DGEList(as.matrix(quant.msstats0_wide_temp))

# set some colors by condition
colors = c(rep('coral3', 2), 
           rep('dodgerblue3', 3), 
           rep('darkolivegreen4',3), 
           rep('gray20',7), 
           rep('hotpink3',7))

# check the clustering
plotMDS(prot_matrix, col = colors, main="plotMDS - sample similarity")
boxplot(quant.msstats0_wide_temp, col = colors)


##################################################
## Filter 1 - Cre+ vs Cre- comparison using MSstatsTMT linear mixed model
##################################################

levels(droplevels(quant.msstats0$Condition))

#create a comparison matrix
comparison1 <- matrix(c(-1, 1, 0, 0,
                        -1, 0, 1, 0,
                        -1, 0, 0, 1),
                      nrow=3, ncol=4, byrow=TRUE)
colnames(comparison1) <- c("control","NES", "hM3Dq", "mCherry")
row.names(comparison1) <- c("NES-control", "hM3Dq-control", "mCherry-control")
comparison1


#perform MSstats group comparison as indicated above
test.contrast1 <- groupComparisonTMT(data = quant.msstats0,
                                     contrast.matrix = comparison1,
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE)


##################################################
## Positively enriched protein lists (adj.pvalue <= 0.05 & log2FC > 0 over pull-down control)
##################################################

positive_enrich <- test.contrast1 %>% filter(adj.pvalue <= 0.05 & log2FC > 0)
pos_enrich_NES <- positive_enrich %>% filter(Label == 'NES-control') %>% 
  dplyr::select(Protein)
pos_enrich_hM3Dq <- positive_enrich %>% filter(Label == 'hM3Dq-control') %>% 
  dplyr::select(Protein)
pos_enrich_mCherry <- positive_enrich %>% filter(Label == 'mCherry-control') %>% 
  dplyr::select(Protein)


##################################################
## Filter 2 - H2B vs NES comparison
##################################################

# Filter 2.1 - hM3Dq vs NES comparison
hM3Dq_comp2 <- quant.msstats0 %>% filter(Condition == 'hM3Dq' | Condition == 'NES') %>%
  filter(Protein %in% pos_enrich_hM3Dq$Protein)
levels(droplevels(hM3Dq_comp2$Condition))

test.contrast2_1 <- groupComparisonTMT(data = hM3Dq_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

# Filter 2.2 - mCherry vs NES comparison
mCherry_comp2 <- quant.msstats0 %>% filter(Condition == 'mCherry' | Condition == 'NES') %>%
  filter(Protein %in% pos_enrich_mCherry$Protein)
levels(droplevels(mCherry_comp2$Condition))

test.contrast2_2 <- groupComparisonTMT(data = mCherry_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

##################################################
## H2B enriched protein lists (adj.pvalue < 0.05 & log2FC > 0 NES)
##################################################

# reuse variables!!!! 
pos_enrich_hM3Dq <- test.contrast2_1 %>% filter(adj.pvalue <= 0.05 & log2FC > 0)
pos_enrich_mCherry <- test.contrast2_2 %>% filter(adj.pvalue <= 0.05 & log2FC > 0)

# Union the two dataset
pos_enrich_H2B <- full_join(pos_enrich_hM3Dq, pos_enrich_mCherry, by = 'Protein') %>% dplyr::select(Protein)

# msstats full normalization
quant.msstats0_globnorm <- proteinSummarization(input.pd_frac,
                                       method="msstats",
                                       global_norm=TRUE,
                                       reference_norm=TRUE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)

# wide format protein abundance and set negative abundance to NA
quant.msstats0_globnorm_wide <- quant.msstats0_globnorm %>% 
  filter(Condition == 'hM3Dq' | Condition == 'mCherry') %>%
  filter(Protein %in% pos_enrich_H2B$Protein) %>% 
  dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

# perform additional median substract because MSstatsTMT was unable to align the median at the protein level
quant.msstats0_globnorm_wide[,-1] <- med_norm(quant.msstats0_globnorm_wide[,-1])
quant.msstats0_globnorm_wide[quant.msstats0_globnorm_wide < 0] <- NA

# remove NA from the clustering analysis
quant.msstats0_globnorm_wide_temp <- quant.msstats0_globnorm_wide[,-1] %>% drop_na()

# re-arrange columns
col_order <- c('hM3Dq_1', 'hM3Dq_2', 'hM3Dq_3', 'hM3Dq_4', 'hM3Dq_5', 'hM3Dq_6', 'hM3Dq_7',
               'mCherry_1', 'mCherry_2','mCherry_3', 'mCherry_4', 'mCherry_5', 'mCherry_6', 'mCherry_7')
quant.msstats0_globnorm_wide_temp <- quant.msstats0_globnorm_wide_temp[, col_order]

# create a DGEList object
prot_matrix_DREADD <- DGEList(as.matrix(quant.msstats0_globnorm_wide_temp))

# set some colors by condition
colors = c(#rep('gray20',7), 
           rep('hotpink3',7),
           rep('darkturquoise',7))


# check the clustering
# plotMDS(prot_matrix, pch = 'o', col = colors, main="plotMDS - sample similarity")
plotMDS(prot_matrix_DREADD, pch = 'o', col = colors, main="plotMDS - sample similarity")
boxplot(quant.msstats0_globnorm_wide_temp, col = colors)


##################################################
## hM3Dq vs mCherry comparison
##################################################

comp3 <- quant.msstats0_globnorm %>% filter(Condition == 'hM3Dq' | Condition == 'mCherry') %>%
  filter(Protein %in% pos_enrich_H2B$Protein) %>% dplyr::select(!c(Abundance))

# use med_norm for comparison
comp3_med <- quant.msstats0_globnorm_wide %>%
  pivot_longer(cols = -c(Protein),
               names_to = 'BioReplicate',
               values_to = 'Abundance')

comp3_med <- merge(comp3_med, comp3, by = c('Protein', 'BioReplicate'))

test.contrast3 <- groupComparisonTMT(data = comp3_med,
                                     contrast.matrix = 'pairwise',
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = FALSE,
                                     remove_empty_channel = TRUE)

#violin plot
comp3_med %>% ggplot(mapping = aes(x = BioReplicate, y = Abundance)) + 
  geom_violin(trim = TRUE,
              mapping = aes(fill = Condition)) +
  scale_fill_manual(values = c('lightpink1','lightblue1')) +
  scale_x_discrete(labels = NULL) + 
  stat_summary(fun=median, geom="crossbar", size=0.2, color="black") + 
  #facet_wrap(~ Condition, scale = "free", ncol = 6) +
  #scale_y_continuous(limits = c(-4,12)) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())


##################################################
## comparison MAPPED
##################################################

mapped_ID <-  raw.protein.pd_frac %>%
  dplyr::select(Accession:`Gene Symbol`)

## for H2B-NES comparison
test.contrast1_mapped <- merge(x=test.contrast1, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession") %>% 
  plot_label(FCvalue = 3, FDRsig = 0.05, FDRplot = 0.05, direction = 'pos')

## for H2B-NES comparison
test.contrast2_1_mapped <- merge(x=test.contrast2_1, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")
test.contrast2_2_mapped <- merge(x=test.contrast2_2, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")

test.contrast2_combined <- rbind(test.contrast2_1, test.contrast2_2)
test.contrast2_mapped <- merge(x=test.contrast2_combined, y=mapped_ID, 
                               by.x="Protein", by.y = "Accession")

## for hM3Dq-mCherry comparison
test.contrast3_mapped <- merge(x=test.contrast3, y=mapped_ID, 
                                 by.x="Protein", by.y = "Accession")

##################################################
## Volcano plots
##################################################

# volcano plot over Cre-negative
volcano_df1 <- ggplot(data = test.contrast1_mapped,
                      aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour = factor(sig)), size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6) +
  facet_wrap(~Label, scale = "free", ncol = 6) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

volcano_df1

# volcano plot over APEX NES
test.contrast2_mapped <- plot_label(test.contrast2_mapped,
                                      FCvalue = 3,
                                      FDRsig = 0.05,
                                      FDRplot = 0.05,
                                      direction = 'pos')

volcano_df2_mapped <- ggplot(data = test.contrast2_mapped,
                             aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(sig)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  facet_wrap(~Label, scale = "free") + 
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
volcano_df2_mapped


# volcano plot hM3Dq vs mCherry
test.contrast3_mapped <- plot_label(test.contrast3_mapped,
                                    FCvalue = 1.75,
                                    FDRsig = 0.05,
                                    FDRplot = 0.05,
                                    direction = 'both')

volcano_df3_mapped <- ggplot(data = test.contrast3_mapped,
                             aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2", "black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  geom_text_repel(aes(label = Protein_label_select), 
                  max.overlaps = Inf,
                  force = 5, #direction = "y",
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4, fontface = 'italic') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
  
volcano_df3_mapped


##################################################
## Remarks: non H2B-enriched proteins comparison
##################################################

# wide format protein abundance and set negative abundance to NA
# only keep hM3Dq and mCherry conditions
quant.msstats0_globnorm_wide2 <- quant.msstats0_globnorm %>% 
  filter(Condition == 'hM3Dq' | Condition == 'mCherry') %>%
  dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)


# perform additional median substract because MSstatsTMT was unable to align the median at the protein level
quant.msstats0_globnorm_wide2[,-1] <- med_norm(quant.msstats0_globnorm_wide2[,-1])
quant.msstats0_globnorm_wide2[quant.msstats0_globnorm_wide2 < 0] <- NA


# list of Cre-positive proteins
pos_enrich_Cre <- positive_enrich %>% filter(Label == 'hM3Dq-control' | Label == 'mCherry-control') %>%
  dplyr::select(Protein) %>% unique()

# list of proteins that are NOT H2B enriched and Cre-postive
comp4 <- quant.msstats0_globnorm %>% filter(Condition == 'hM3Dq' | Condition == 'mCherry') %>%
  dplyr::filter(!(Protein %in% pos_enrich_H2B$Protein)) %>% #not H2B-enriched
  dplyr::filter(Protein %in% pos_enrich_Cre$Protein) %>%    #only Cre-positive
dplyr::select(!c(Abundance))

# use med_norm for comparison
comp4_med <- quant.msstats0_globnorm_wide2 %>%
  pivot_longer(cols = -c(Protein),
               names_to = 'BioReplicate',
               values_to = 'Abundance')

comp4_med <- merge(comp4_med, comp4, by = c('Protein', 'BioReplicate'))


test.contrast4 <- groupComparisonTMT(data = comp4_med,
                                     contrast.matrix = 'pairwise',
                                     moderated = TRUE,
                                     adj.method = 'BH',
                                     remove_norm_channel = FALSE,
                                     remove_empty_channel = TRUE)


## for non-H2B hM3Dq-mCherry comparison
test.contrast4_mapped <- merge(x=test.contrast4, y=mapped_ID, 
                               by.x="Protein", by.y = "Accession")

test.contrast4_mapped <- plot_label(test.contrast4_mapped,
                                    FCvalue = 1.75,
                                    FDRsig = 0.05,
                                    FDRplot = 0.05,
                                    direction = 'both')

volcano_df4_mapped <- ggplot(data = test.contrast4_mapped,
                             aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2", "black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  geom_text_repel(aes(label = Protein_label_select), 
                  max.overlaps = Inf,
                  force = 5, #direction = "y",
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4, fontface = 'italic') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

volcano_df4_mapped


##################################################
## Remarks: non H2B-enriched proteins comparison2
##################################################

# H2B protein list was generated by union of mCherry-enriched and hM3Dq-enriched over NES
# how about proteins that are only enriched in either of them but not both?

# hM3Dq only
diff_enrich_H2B_1 <- setdiff(pos_enrich_hM3Dq$Protein, pos_enrich_mCherry$Protein)
# mCherry only
diff_enrich_H2B_2 <- setdiff(pos_enrich_mCherry$Protein, pos_enrich_hM3Dq$Protein)
# intersection
intersect_enrich_H2B <- intersect(pos_enrich_mCherry$Protein, pos_enrich_hM3Dq$Protein)


test.contrast5_remarks <- rbind(test.contrast2_combined, test.contrast3) %>%
  dplyr::filter(Protein %in% pos_enrich_H2B$Protein) %>% 
  dplyr::select(Protein, Label, log2FC, adj.pvalue) %>%
  pivot_wider(names_from = Label,
              names_glue = "{Label}_{.value}",
              values_from = c(log2FC,  adj.pvalue)) %>% 
  merge(y=mapped_ID, by.x="Protein", by.y = "Accession")

data_tmp <- test.contrast5_remarks

for (i in seq_along(data_tmp$Protein)){
  # adj.pvalue < FDR
  data_tmp$sig[i] <- ifelse(data_tmp$`hM3Dq-mCherry_adj.pvalue`[i] < 0.05, 1, 0)
  
  # up or down regulation
  data_tmp$reg[i] <- ifelse(data_tmp$sig[i] == 1, ifelse(data_tmp$`hM3Dq-mCherry_log2FC`[i] > 0, 'up', 'down'), 'none')
  
  # hM3Dq, mCherry, both H2B enriched
  if (data_tmp$Protein[i] %in% diff_enrich_H2B_1) {
    data_tmp$H2B_enrich[i] <- 'hM3Dq_only'
  } else if (data_tmp$Protein[i] %in% diff_enrich_H2B_2) {
    data_tmp$H2B_enrich[i] <- 'mCherry_only'
  } else {data_tmp$H2B_enrich[i] <- 'both'}
}

test.contrast5_remarks <- data_tmp
rm(data_tmp)

# both
test.contrast5_remarks_1 <- test.contrast5_remarks %>%
  dplyr::filter(H2B_enrich == 'both')

# hM3Dq_only
test.contrast5_remarks_2 <- test.contrast5_remarks %>%
  dplyr::filter(H2B_enrich == 'hM3Dq_only')

# mCherry_only
test.contrast5_remarks_3 <- test.contrast5_remarks %>%
  dplyr::filter(H2B_enrich == 'mCherry_only')

scatter_remarks_2 <- ggplot(data=test.contrast5_remarks_2,
                           aes(x=`hM3Dq-mCherry_log2FC`, y=`hM3Dq-NES_log2FC`)) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
scatter_remarks_2  
  

scatter_remarks_3 <- ggplot(data=test.contrast5_remarks_3,
                            aes(x=`hM3Dq-mCherry_log2FC`, y=`mCherry-NES_log2FC`)) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2", "black")) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
scatter_remarks_3


# can we plot them all together!
scatter_remarks_0 <- ggplot(data=test.contrast5_remarks,
                            aes(x=`hM3Dq-NES_log2FC`, y=`mCherry-NES_log2FC`)) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2","black", "coral2")) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
scatter_remarks_0


# Fig 5 - main fig. all H2B_enriched proteins
# scatter plot
scatter_remarks_1 <- ggplot(data=test.contrast5_remarks,
                            aes(x=`hM3Dq-NES_log2FC`, y=`mCherry-NES_log2FC`)) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2","black", "coral2")) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
scatter_remarks_1


##################################################
## Remarks: heatmap
##################################################

#list of differentially regulated proteins
diff_reg_all <- rbind(test.contrast3_mapped, test.contrast4_mapped) %>%
  filter(adj.pvalue < 0.05)

# intensities values
diff_reg_all_wide <- quant.msstats0_globnorm_wide2 %>%
  filter(Protein %in% diff_reg_all$Protein)

# with IRS and med_norm
heatmap_input0 <- diff_reg_all_wide[,-1] %>% drop_na()
Heatmap(as.matrix(heatmap_input0),
        show_row_names = FALSE)


##################################################
## Remarks: diff reg all
##################################################
#list of differentially regulated proteins
diff_reg_all <- rbind(test.contrast3_mapped, test.contrast4_mapped) %>%
  filter(adj.pvalue < 0.05)

# intensities values
diff_reg_all_wide <- quant.msstats0_globnorm_wide2 %>%
  filter(Protein %in% diff_reg_all$Protein)

# all comparison
all_comp_wide <- rbind(test.contrast3_mapped, test.contrast4_mapped)

all_comp_wide <- merge(all_comp_wide, quant.msstats0_globnorm_wide2, by='Protein')

# write.table(diff_reg_all, "/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set3/diff_reg_all_log2FC.txt", sep="\t", row.names = FALSE)


##################################################
## Remarks: JunB interactors
##################################################

JunB_PPI_IDs <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/cytoscape_network_analysis/20210506_JunB_network_manually_combined_nodes_info.xlsx")

JunB_PPI_df <- raw.protein.pd_frac %>%
  filter(Accession %in% JunB_PPI_IDs$Entry) 
colnames(JunB_PPI_df)[29:50] <- annotation.pd_frac %>% filter(Fraction == 1) %>% 
  dplyr::select(BioReplicate) %>% unlist()

JunB_PPI_df <- JunB_PPI_df[,c(2,3,4,5,6,9,11,22,29:50)]

JunB_PPI_df_long <- JunB_PPI_df %>% dplyr::select(Accession, hM3Dq_1:Norm_2) %>%
  pivot_longer(cols = -Accession,
               names_to = 'BioReplicate',
               values_to = 'Intensity') %>%
  mutate(Condition = str_extract(BioReplicate, '^[^_]+')) %>%
  group_by(Condition, Accession) %>%
  summarise(cond_avg = mean(Intensity, na.rm=TRUE)) %>%
  filter(Condition %in% c('hM3Dq', 'mCherry')) %>%
  pivot_wider(names_from = Condition,
              values_from = cond_avg) %>%
  mutate(log2FC = log2(hM3Dq / mCherry)) %>%
  merge(JunB_PPI_IDs, by.x = 'Accession', by.y = 'Entry')


##################################################
## Remarks: JunB interactors + comp3_med
##################################################

#########################
## For low abundance JunB interactors; allowing n=1 unique peptide
## Directly use reporter ion intensties
#########################

# include JunB interactor PSMs that were filtered out by MSstatsTMT
JunB_PPI_IDs <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/cytoscape_network_analysis/20210506_JunB_network_manually_combined_nodes_info.xlsx")

# proteins that were filtered out by MSstatsTMT and find those PSMs
JunB_PPI_extra <- JunB_PPI_IDs[!(JunB_PPI_IDs$Entry %in% input.pd_frac$ProteinName), ]
JunB_PPI_extra <- raw.pd_frac %>% filter(Master.Protein.Accessions %in% JunB_PPI_extra$Entry) %>%
  mutate(ProteinName = Master.Protein.Accessions) %>%
  mutate(PeptideSequence = Annotated.Sequence) %>%
  mutate(PSM = paste(PeptideSequence, Charge, sep='_'))
colnames(JunB_PPI_extra)[31:41] <- c('126', '127N', '127C', '128N', '128C', '129N', '129C', '130N',
                                     '130C', '131N', '131C')
JunB_PPI_extra <- JunB_PPI_extra %>%
  dplyr::select(ProteinName, Spectrum.File, `126`:`131C`)

# manually inspect to remove redundant or empty PSMs
JunB_PPI_extra <- JunB_PPI_extra[-c(2,3,7,10,13),]

# map other columns in the long format
JunB_PPI_extra <- JunB_PPI_extra %>%
  pivot_longer(cols = -c(ProteinName, Spectrum.File),
               names_to = 'Channel',
               values_to = 'Abundance') %>% drop_na() %>%
  left_join(annotation.pd_frac[,-c(2,8,9)], by = c("Spectrum.File" = "Run", "Channel" = "Channel")) %>%
  mutate(Run = paste(Mixture, TechRepMixture, sep = '_')) %>% dplyr::select(-Spectrum.File)
colnames(JunB_PPI_extra)[1] <- c('Protein')

#log2 transform intensity
JunB_PPI_extra$Abundance <- log2(JunB_PPI_extra$Abundance)

#reorder column and only retain hM3Dq and mCherry
JunB_PPI_extra <- JunB_PPI_extra[, colnames(comp3_med)] %>%
  dplyr::filter(Condition %in% c('hM3Dq', 'mCherry'))


#### combining with comp3_med
comp3_med_plus <- rbind(comp3_med, JunB_PPI_extra)

## for hM3Dq-mCherry comparison
test.contrast3_plus <- groupComparisonTMT(data = comp3_med_plus,
                                          contrast.matrix = 'pairwise',
                                          moderated = TRUE,
                                          adj.method = 'BH',
                                          remove_norm_channel = FALSE,
                                          remove_empty_channel = TRUE)

test.contrast3_plus_mapped <- merge(x=test.contrast3_plus, y=mapped_ID, 
                                    by.x="Protein", by.y = "Accession")

# volcano plot hM3Dq vs mCherry
test.contrast3_plus_mapped <- plot_label(test.contrast3_plus_mapped,
                                         FCvalue = 1.75,
                                         FDRsig = 0.05,
                                         FDRplot = 0.05,
                                         direction = 'both')

volcano_df3_plus_mapped <- ggplot(data = test.contrast3_plus_mapped,
                                  aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(aes(colour=factor(reg)), size = 1, show.legend = FALSE) + 
  scale_colour_manual(values = c("dodgerblue2", "black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), 
             linetype = "dotted", size = 0.6) +
  geom_text_repel(aes(label = Protein_label_select), 
                  max.overlaps = Inf,
                  force = 5, #direction = "y",
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4, fontface = 'italic') +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())

volcano_df3_plus_mapped

# box plot of proteins with 1 unique peptides
JunB_PPI_extra_list <- JunB_PPI_IDs[!(JunB_PPI_IDs$Entry %in% input.pd_frac$ProteinName), ]
comp3_med_plus$Condition <- factor(comp3_med_plus$Condition,
                                   levels = c('mCherry', 'hM3Dq'))
comp3_med_plus$Protein <- factor(comp3_med_plus$Protein,
                                 levels = c('P01101', 'P13346', 'P15066', 'P47930',
                                            'P63248', 'Q8BUN5', 'Q91Y86'))

comp3_med_plus %>% dplyr::filter(Protein %in% JunB_PPI_extra_list$Entry) %>%
  ggplot(aes(Condition, Abundance)) +
  geom_boxplot() +
  geom_point(size=1, shape=19, position = position_jitter()) +
  #geom_boxplot(aes(color=factor(Protein))) +
  facet_wrap(~Protein, ncol = 7) +
  theme_bw() + theme(panel.border = element_rect(size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1))