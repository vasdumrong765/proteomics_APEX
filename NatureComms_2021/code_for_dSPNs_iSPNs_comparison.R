#########################
## Dumrongprechachan et al 2021 Nature Comms
## D1 vs A2a comparison
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
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/med_SL_norm.R")
source("/Users/vasdumrong/Box/VasD_projects/Rfunctions/plot_label.R")
library(ComplexHeatmap)


#########################
## Load data
#########################

# Read PSMs PD output into R
# S/N export from PD
# Note that if you use read_excel column names will be wrong (correct e.g., `Spectrum.File` not `Spectrum File`)
raw.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set4/D1_D2P149099-105_limitedsearch_PSMs.txt", sep="\t", header=TRUE) %>% 
  dplyr::filter(Isolation.Interference.... < 70) %>%
  dplyr::filter(Quan.Info == '') %>%
  dplyr::filter(Percolator.q.Value < 0.05)

# Read PD protein tabs into R
# for intensity values
raw.protein.pd_frac <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set4/D1_D2P149099-105_limitedsearch_Proteins.xlsx") #%>%
  #dplyr::filter(`Protein FDR Confidence: Combined` %in% c('Medium', 'High'))

# Read annotation file
annotation.pd_frac <- read.table("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set4/PD_D1D2_annotation_limitedsearch.txt", sep="\t", header=TRUE)


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

##################################################
## Use MSstats for protein summarization using LogSum method
##################################################

# Just 1 set of TMT - no reference normalization
# Will remove contaminants manually, so go ahead and perform global normalization
quant.msstats0 <- proteinSummarization(input.pd_frac,
                                       method="msstats",
                                       global_norm=TRUE,
                                       reference_norm=FALSE,
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE,                  
                                       MBimpute = FALSE,
                                       maxQuantileforCensored = NULL)


##################################################
## EdgeR multidimensional clustering MDS plot before global protein-level normalization
##################################################

# wide format protein abundance and set negative abundance to NA
quant.msstats0_wide <- quant.msstats0 %>% dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

quant.msstats0_wide[quant.msstats0_wide <= 0] <- NA

# remove NA from the clustering analysis
quant.msstats0_wide_temp <- quant.msstats0_wide[-1] %>% drop_na()

# re-arrange columns
col_order <- c('A2aCre_NES_1', 'A2aCre_NES_2', 'A2aCre_NES_3',
               'A2aCre_H2B_1','A2aCre_H2B_2','A2aCre_H2B_3','A2aCre_H2B_4','A2aCre_H2B_5',
               'D1Cre_NES_1','D1Cre_NES_2','D1Cre_NES_3',
               'D1Cre_H2B_1','D1Cre_H2B_2','D1Cre_H2B_3','D1Cre_H2B_4','D1Cre_H2B_5')
quant.msstats0_wide_temp <- quant.msstats0_wide_temp[, col_order]

# create a DGEList object
prot_matrix <- DGEList(as.matrix(quant.msstats0_wide_temp))

# set some colors by condition
colors = c(rep('coral3', 3), 
           rep('dodgerblue3', 5), 
           rep('darkolivegreen4',3), 
           rep('gray20',5))

# check the clustering
plotMDS(prot_matrix, pch = 19, col = colors, main="plotMDS - sample similarity")
plotMDS(prot_matrix, col = colors, main="plotMDS - sample similarity")
boxplot(quant.msstats0_wide_temp, col = colors)


##################################################
## Global protein-level normalization
##################################################

# wide format protein abundance and set negative abundance to NA
quant.msstats0_wide_medSL <- 
  quant.msstats0 %>% dplyr::select(Protein, Abundance, BioReplicate) %>%
  pivot_wider(names_from = BioReplicate,
              values_from = Abundance)

quant.msstats0_wide_medSL_prot <- quant.msstats0_wide_medSL[,1]
quant.msstats0_wide_medSL <- med_norm(log2(quant.msstats0_wide_medSL[,-1]))
#quant.msstats0_wide_medSL <- med_norm(log2(SL_Norm(2^quant.msstats0_wide_medSL[,-1])))
quant.msstats0_wide_medSL[quant.msstats0_wide_medSL <= 0] <- NA

quant.msstats0_wide_medSL <- cbind(quant.msstats0_wide_medSL_prot, quant.msstats0_wide_medSL)

boxplot(quant.msstats0_wide_medSL[,-1], col=colors)

# long format data
quant.msstats0_wide_medSL_long <- pivot_longer(quant.msstats0_wide_medSL,
                                               cols = -Protein,
                                               names_to = 'BioReplicate',
                                               values_to = 'Abundance')
quant.msstats0_wide_medSL_long <- merge(quant.msstats0_wide_medSL_long,
                                        quant.msstats0[,-3],
                                        by=c('Protein','BioReplicate'))


##################################################
## Filter 1 - Cre+ vs Cre- comparison using previously generated list
##################################################

## filter 1 derived from DREADD dataset
data3_cre_neg_comp <- readRDS("/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/reanalysis/analysis_PD_set4/test.contrast1_dataset3.rds")

# strict
contaminant_list <- data3_cre_neg_comp %>%
    dplyr::filter(!(adj.pvalue < 0.05 & log2FC > 0))
contaminant_list <- unique(contaminant_list$Protein)
length(contaminant_list)

keratin_mouse <- read_excel(path="/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/mmu_keratin.xlsx")

## Putative Cre-positive proteome
quant.msstats0_Crepos <- dplyr::filter(quant.msstats0_wide_medSL_long, 
                                       !(Protein %in% contaminant_list)) %>%
  dplyr::filter(!(Protein %in% keratin_mouse$Entry))


##################################################
## Filter 2 - H2B vs NES comparison
##################################################

# Filter 2.1 - D1Cre H2B vs NES comparison
D1Cre_comp2 <- quant.msstats0_Crepos %>% 
  filter(Condition == 'D1Cre_NES' | Condition == 'D1Cre_H2B')
levels(droplevels(D1Cre_comp2$Condition))

test.contrast2_1 <- groupComparisonTMT(data = D1Cre_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

# Filter 2.2 - A2aCre vs NES comparison
A2aCre_comp2 <- quant.msstats0_Crepos %>% filter(Condition == 'A2aCre_NES' | Condition == 'A2aCre_H2B')
levels(droplevels(A2aCre_comp2$Condition))

test.contrast2_2 <- groupComparisonTMT(data = A2aCre_comp2,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

volcano_df2_1 <- ggplot(data = test.contrast2_1,
                        aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df2_2 <- ggplot(data = test.contrast2_2,
                        aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df2_1
volcano_df2_2


##################################################
## Filter 3 - D1Cre vs A2aCre comparison
##################################################

NES_comp3 <- quant.msstats0_Crepos %>% filter(Condition == 'D1Cre_NES' | Condition == 'A2aCre_NES')
levels(droplevels(NES_comp3$Condition))

test.contrast3_1 <- groupComparisonTMT(data = NES_comp3,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)


H2B_comp3 <- quant.msstats0_Crepos %>% filter(Condition == 'D1Cre_H2B' | Condition == 'A2aCre_H2B')
levels(droplevels(H2B_comp3$Condition))

test.contrast3_2 <- groupComparisonTMT(data = H2B_comp3,
                                       contrast.matrix = 'pairwise',
                                       moderated = TRUE,
                                       adj.method = 'BH',
                                       remove_norm_channel = FALSE,
                                       remove_empty_channel = TRUE)

volcano_df3_1 <- ggplot(data = test.contrast3_1,
                        aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df3_2 <- ggplot(data = test.contrast3_2,
                        aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df3_1
volcano_df3_2

##################################################
## Volcano plot filtered by H2B vs NES proteins
##################################################

H2B_enriched_D1Cre <- test.contrast2_1 %>%
  dplyr::filter(log2FC > 0 & adj.pvalue < 0.05) %>% 
  dplyr::select(Protein)
H2B_enriched_A2aCre <- test.contrast2_2 %>%
  dplyr::filter(log2FC > 0 & adj.pvalue < 0.05) %>% 
  dplyr::select(Protein)

H2B_enriched_union <- unique(rbind(H2B_enriched_D1Cre, H2B_enriched_A2aCre))

NES_enriched <- test.contrast2_1 %>% 
  dplyr::filter(!(Protein %in% H2B_enriched_union$Protein))

test.contrast3_1_NES <-  test.contrast3_1 %>%
  dplyr::filter(Protein %in% NES_enriched$Protein)

test.contrast3_2_H2B <- test.contrast3_2 %>%
  dplyr::filter(Protein %in% H2B_enriched_union$Protein)


volcano_df3_1_NES <- ggplot(data = test.contrast3_1_NES,
                            aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "coral2")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df3_2_H2B <- ggplot(data = test.contrast3_2_H2B,
                            aes(x = log2FC, y = (-1)*log10(adj.pvalue))) +
  geom_point(size = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "red")) +
  geom_hline(yintercept = (-1)*log10(0.05), linetype = "dotted", size = 0.6)

volcano_df3_1_NES
volcano_df3_2_H2B

##################################################
## Heatmap after a nominal pvalue cutoff
##################################################

nominal_cutoff_NES <- dplyr::filter(test.contrast3_1_NES, pvalue < 0.30)

heatmap_input_NES <- quant.msstats0_wide_medSL[,c(1, 7:9, 15:17)] %>%
                     dplyr::filter(Protein %in% nominal_cutoff_NES$Protein) %>% drop_na()

Heatmap(as.matrix(heatmap_input_NES[,-1]))
  

nominal_cutoff_H2B <- dplyr::filter(test.contrast3_2_H2B, pvalue < 0.20)

heatmap_input_H2B <- quant.msstats0_wide_medSL[,c(1, 2:6, 10:14)] %>%
  dplyr::filter(Protein %in% nominal_cutoff_H2B$Protein) %>% drop_na()

Heatmap(as.matrix(heatmap_input_H2B[,-1]))


##################################################
## Correlation with transcripts data
##################################################
# List of transcripts that are differentially expressed between D1 and A2a SPNs
transcripts <- read_excel('/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/D1D2_transcripts_Uniprot.xlsx')


## Compare with TRAPseq data
transcripts_Heiman <- read_excel('/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/mapped.xlsx', sheet=1)
colnames(transcripts_Heiman)[1] <- 'yourlist'
#RNAdata_Heiman <- read_excel(path='/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/raw_from_papers.xlsx', sheet=1)
RNAdata_Heiman <- read_excel(path='/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/Heiman_map_nonredundant.xlsx', sheet=1)
RNAdata_Heiman <- merge(RNAdata_Heiman, transcripts_Heiman, by.x='original_list', by.y='yourlist')
RNAdata_Heiman <- RNAdata_Heiman[,c('Entry.x', 'log2FC', 'Input_list')]

# NES proteins
proteinNES_Heiman <- test.contrast3_1_NES %>%
                     dplyr::filter(Protein %in% RNAdata_Heiman$Entry) %>%
                     dplyr::select(Protein, log2FC)
proteinNES_Heiman$log2FC <- proteinNES_Heiman$log2FC * (-1)

# H2B proteins
proteinH2B_Heiman <- test.contrast3_2_H2B %>%
  dplyr::filter(Protein %in% RNAdata_Heiman$Entry) %>%
  dplyr::select(Protein, log2FC)
proteinH2B_Heiman$log2FC <- proteinH2B_Heiman$log2FC * (-1)

protein_Heiman_all <- rbind(proteinH2B_Heiman, proteinNES_Heiman)


corr_Heiman <- merge(protein_Heiman_all, RNAdata_Heiman, by.x='Protein', by.y='Entry.x')

corr_Heiman %>%
  ggplot(aes(x=log2FC.x, y=log2FC.y)) +
  geom_point() +
  xlab('log2FC_Protein(D1Cre - A2aCre)') + 
  ylab('log2FC_TRAPseq(D1 - D2)') +
  stat_cor(method = 'pearson', 
           cor.coef.name = 'R',
           label.x.npc = 'center', label.y.npc = 'bottom') + 
  #geom_smooth(method="lm") +
  geom_text_repel(aes(label = Input_list), 
                  force = 5, #direction = "y",
                  max.overlaps = Inf,
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4) #font size  

corr_Heiman %>%
  ggplot(aes(x=log2FC.x, y=log2FC.y)) +
  geom_point() +
  xlab('log2FC_Protein(D1Cre / A2aCre)') + 
  ylab('log2FC_TRAPseq(D1 / D2), Heiman etal 2008') +
  stat_cor(method = 'spearman', 
           cor.coef.name = 'rho',
           label.x.npc = 'center', label.y.npc = 'bottom') + 
  #geom_smooth(method="lm") +
  geom_text_repel(aes(label = Input_list), 
                  force = 5, #direction = "y",
                  max.overlaps = Inf,
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(0, 'lines'), #draw every segment
                  segment.size = 0.5, #segment line thickness
                  size=4) #font size  


##################################################
## Correlation with transcripts data - Gocke etal 2016
##################################################

RNAdata_GockeS5 <- read_excel(path='/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/raw_from_papers.xlsx', sheet=2)
transcripts_GockeS5 <- read_excel('/Users/vasdumrong/Box/VasD_projects/Project_APEX_AAV/analysis/D1vsD2_transcripts/mapped.xlsx', sheet=2)
colnames(transcripts_GockeS5)[1] <- 'yourlist'


RNAdata_GockeS5 <- merge(RNAdata_GockeS5, transcripts_GockeS5, by.x='Rank', by.y='yourlist')
RNAdata_GockeS5 <- RNAdata_GockeS5[,c('Entry', 'D1-D2', 'Rank')]

# NES proteins
proteinNES_GockeS5 <- test.contrast3_1_NES %>%
  dplyr::filter(Protein %in% RNAdata_GockeS5$Entry) %>%
  dplyr::select(Protein, log2FC)
proteinNES_GockeS5$log2FC <- proteinNES_GockeS5$log2FC * (-1)

# H2B proteins
proteinH2B_GockeS5 <- test.contrast3_2_H2B %>%
  dplyr::filter(Protein %in% RNAdata_GockeS5$Entry) %>%
  dplyr::select(Protein, log2FC)
proteinH2B_GockeS5$log2FC <- proteinH2B_GockeS5$log2FC * (-1)

protein_GockeS5_all <- rbind(proteinH2B_GockeS5, proteinNES_GockeS5)


corr_GockeS5 <- merge(protein_GockeS5_all, RNAdata_GockeS5, by.x='Protein', by.y='Entry')
corr_GockeS5 <- as.data.frame(corr_GockeS5) %>%
  dplyr::mutate(Rank_select = ifelse((log2FC > 0.5 & `D1-D2` > 0.5) | (log2FC < -0.5 & `D1-D2` < -0.5), 
                                     Rank, ''))

corr_GockeS5 %>%
  ggplot(aes(x=log2FC, y=`D1-D2`)) +
  geom_point() +
  #geom_text(aes(label=Rank), hjust=1, yjust=1) + 
  xlab('log2FC_Protein(D1Cre / A2aCre)') + 
  ylab(expression(paste('log2FC_scRNAseq(D1 / D2), Gocke etal 2016'))) +
  stat_cor(method = 'pearson', 
           cor.coef.name = 'R',
           label.x.npc = 'center', label.y.npc = 'bottom') +
  geom_smooth(method="lm", size=0.6, colour='firebrick2') +
  geom_text_repel(aes(label = Rank_select), 
                  force = 5, #direction = "y",
                  max.overlaps = Inf,
                  nudge_x = 0.3, nudge_y = 0.3,
                  min.segment.length = unit(2, 'lines'), #draw every segment
                  segment.size = 0, #segment line thickness
                  size=4) + #font size
  geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted')