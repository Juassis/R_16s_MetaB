###################################################################################################
#Library#
library(phyloseq)
library(microbiome)
library(ggpubr)
library(grid)
library(tidyr)
library(reshape2)
library(reshape)
library(ggrepel)
###################################################################################################
#Phyloseq
Twin_Shime_Longitudinal_Study <- (readRDS("/Users/flb202/Documents/Meu_PC/16s/Shime_20_December/DADA2/Longitudinal_Shime/New_Object_Plot/Twin_Shime_Longitudinal_Study.rds"))
pseq <- Twin_Shime_Longitudinal_Study

#Colors
cols <- c("Fecal" = "#264D59", "AC1" = "#77A515", "TC1" = "#D46C4E", "DC1" = "#43978D", "AC2" = "#77A515", "TC2" = "#D46C4E", "DC2" = "#43978D")
cols4 <- c("Fecal" = "#D46C4E", "AC" = "#77A515", "TC" = "#264D59", "DC" = "#43978D")
cols4_Inv <- c("Fecal" = "#264D59", "AC" = "#77A515", "TC" = "#D46C4E", "DC" = "#43978D")
#

###################################################################################################
####Figure1: A
###AlphaDiversity: Shannon and Richness
##Filter pseq
Shime_Plots_1 <- subset_samples(pseq, Main_Analysis =="YES" )

alpha_diversity_Shime <- estimate_richness(Shime_Plots_1,  measures = c("Shannon", "Observed"))
df <- data.frame(alpha_diversity_Shime, sample_data(Shime_Plots_1))
df$Day <- factor(df$Day, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8 ,9 ,10 ,11 ,12 , 13, 14, 15, 16 ,17 ,18, 19, 20, 21, 22), labels = c("Fecal", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22))
df$Compartment_Shime <- factor(df$Compartment_Shime, levels = c("Fecal", "AC1", "TC1", "DC1", "AC2", "TC2", "DC2"))

#Shannon
p2 <- ggplot(df, aes(x=Day, y=Shannon, color = Compartment_Shime, group = Compartment_Shime, shape = Shime))  +
  theme_pubr(border = TRUE) +
  geom_point(size=6) +
  scale_shape_manual(values=c(15, 20,18)) +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5), 
        axis.text.y = element_text(size = 14, hjust = 1),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none") +
  scale_colour_manual(values = cols) + 
  scale_x_discrete(limits=c("Fecal", 1, 2, " ", " ", 5, 6, 7," "," "," ", " ", 12, 13," ", 15, 16, 17," "," ", 20, 21, 22)) +
  geom_line(aes(linetype=Shime))+
  scale_linetype_manual(values=c("twodash", "solid", "dashed"))+
  ylim(1, 4) +
  labs(hjust = 1, face= "bold", x="Days ", y = "Shannon Diversity") 
p2

#Richness
p <- ggplot(df, aes(x=Day, y=Observed, color = Compartment_Shime, group = Compartment_Shime, shape = Shime))  +
  theme_pubr(border = TRUE) +
  geom_point(size=6) +
  scale_shape_manual(values=c(15, 20,18)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, hjust = 1),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none") +
  scale_colour_manual(values = cols) + 
  scale_x_discrete(limits=c("Fecal", 1, 2, " ", " ", 5, 6, 7," "," "," ", " ", 12, 13," ", 15, 16, 17," "," ", 20, 21, 22)) +
  geom_line(aes(linetype=Shime))+
  scale_linetype_manual(values=c("twodash", "solid", "dashed"))+
  ylim(0, 120) +
  labs(hjust = 0.5, face= "bold", x=" ", y = "Richness") 

p

#Exporting
gA <- ggplotGrob(p)
gB <- ggplotGrob(p2)

#Combine the plots
g = rbind(gA, gB, size = "first")
g = rbind(gA, gB, size = "last")

#Draw it
grid.newpage()
grid.draw(g)

#Save
ggsave("Fig1A_Alpha_Diversity.pdf", units="in", width=10, height=7.7, dpi=300)
dev.off()

#Statistics: pairwise comparision using non-parametric test (Wilcoxon test)#

PY <- levels(df$Compartment) # get the variables
PY.pairs <- combn(seq_along(PY), 2, simplify = FALSE, FUN = function(i)PY[i])
print(PY.pairs)
remove <- c (2, 3, 4, 5, 6, 7,8, 9, 11, 12, 13, 14)
PY.pairs2 <- PY.pairs[-c(2, 3, 4, 5, 6, 7,8, 9, 11, 12, 13, 14)];

#Plot measures
df %>%
  gather(key = metric, value = value, c("Observed", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
  ggplot(aes(x = Compartment_Shime, y = value, color = Compartment_Shime)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Compartment_Shime), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_pubr(border = TRUE) +
  stat_compare_means(comparisons = PY.pairs2, label = "p.signif") + 
  stat_compare_means(label.y = 1, size = 2) +
  scale_colour_manual(values = cols) +
  theme(legend.position="none") 

dev.off()

###################################################################################################
####Figura1:B
###PCOA
##Relative Abundance
Shime_Plots_1 <- transform_sample_counts(Shime_Plots_1, function(otu) otu/sum(otu))

#
expt = prune_taxa(names(sort(taxa_sums(Shime_Plots_1), TRUE)[1:2000]), Shime_Plots_1)
ord<- sqrt(phyloseq::distance(Shime_Plots_1,Type = "samples", "jsd"))
pcoa=ordinate(Shime_Plots_1, "PCoA", distance=ord)

#ggsave("BetaDiversity_JSD_TwinShime.pdf", units="in", width=18, height=7.7, dpi=300)
ord_plot <- (ordplot <- plot_ordination(expt, pcoa, "samples", color="Compartment") + 
               geom_point(size = 6, shape = 16) +
               geom_text_repel(aes(label = Day), nudge_x = 0.04, size = 3.0, segment.alpha = 0.5) +
               theme_pubr(border = TRUE) +
               theme(axis.text=element_text(size=14), 
                     axis.text.x = element_text(size = 12, hjust = 0.5), 
                     axis.title.y = element_text(size = 18),
                     legend.text=element_text(size=14), 
                     legend.title=element_text(size=0),
                     #legend.position="none",
                     axis.title.x = element_text(size = 18), 
                     strip.text.x = element_text(size = 20, face = "bold"))) +
  scale_colour_manual(values = cols4_Inv) #+
labs(x="PCo1 [50.8%]", y = "PCo2 [18.9%]", element_text(face = "bold"))
#labs(x="PCo1 [49.8%]", y = "PCo2 [17.9%]", element_text(face = "bold"))
ord_plot

newSTorder = c("Fecal", "AC", "TC", "DC")
ord_plot$data$Compartment <- as.character(ord_plot$data$Compartment)
ord_plot$data$Compartment <- factor(ord_plot$data$Compartment, levels=newSTorder)
print(ord_plot)
###################################################################################################
####Figura1:C
###OneDay Interval#######################
##Phyloseq
Shime_Plots_1
#JSD Distance

jsd_dist<- sqrt(phyloseq::distance(Shime_Plots_1, "jsd"))
pcoa_jsd <- ordinate(Shime_Plots_1, method = "PCoA", distance = jsd_dist)

#########Only PC1
Shime_Plots_1@sam_data$Compartment
a<-pcoa_jsd$values
b<-as.data.frame(pcoa_jsd$vectors)
samdf<-data.frame(unclass(Shime_Plots_1@sam_data))
row.names(samdf)<-sample_names(Shime_Plots_1)
axis1<-data.frame(b$Axis.1)
colnames(axis1)<-"PC1"
row.names(axis1)<-row.names(b)
info<-merge.data.frame(axis1,samdf, by=0)
info$Day<-as.factor(info$Day)

#
class(jsd_dist)
distance_matrix<-matrix(jsd_dist)
distance_matrix<-as.matrix(jsd_dist)
distance_matrix<-as.data.frame(distance_matrix)

#Making Comparison
samples_1Day <- c("F2.S10D1","S10D1.S12D12","S15D1.S16D1","S16D1.S17D12","S112D12.S113D1","S115D1.S116D1","S116D1.S117D12","S120D1.S121D1","S121D1.S123D12","F2.S10D2","S10D2.S12D22","S15D2.S16D2","S16D2.S17D22","S112D22.S113D2","S115D2.S116D2","S116D2.S117D22","S120D2.S121D2","S121D2.S123D22","F2.S10A1","S10A1.S12A12","S15A1.S16A1","S16A1.S17A12","S112A12.S113A1","S115A1.S116A1","S116A1.S117A12","S120A1.S121A1","S121A1.S123A12","F2.S10A2","S10A2.S12A22","S15A2.S16A2","S16A2.S17A22","S112A22.S113A2","S115A2.S116A2","S116A2.S117A22","S120A2.S121A2","S121A2.S123A22","F2.S10T1","S10T1.S12T12","S15T1.S16T1","S16T1.S17T12","S112T12.S113T1","S115T1.S116T1","S116T1.S117T12","S120T1.S121T1","S121T1.S123T12","F2.S10T2","S10T2.S12T22","S15T2.S16T2","S16T2.S17T22","S112T22.S113T2","S115T2.S116T2","S116T2.S117T22","S120T2.S121T2","S121T2.S123T22")

#
distance_matrix$samples<-row.names(distance_matrix)
distance_matrix_melt<-melt.data.frame(distance_matrix, id.vars = "samples")
distance_matrix_melt$Comparision<-paste(distance_matrix_melt$samples, distance_matrix_melt$variable, sep = ".")
distance_matrix_melt_filt<-subset(distance_matrix_melt, Comparision %in% samples_1Day)
row.names(distance_matrix_melt_filt)<-as.character(distance_matrix_melt_filt$variable)
distance_matrix_melt_filt_info<-merge(distance_matrix_melt_filt, samdf, by=0)

#Factor
distance_matrix_melt_filt_info$Compartment_Shime <- factor(distance_matrix_melt_filt_info$Compartment_Shime, levels = c("Fecal", "AC1", "AC2","TC1", "TC2","DC1", "DC2"))
distance_matrix_melt_filt_info$Compartment <- factor(distance_matrix_melt_filt_info$Compartment, levels = c("Fecal", "AC","TC","DC"))
distance_matrix_melt_filt_info$Shime <- factor(distance_matrix_melt_filt_info$Shime, levels = c("SHIME-1", "SHIME-2"))
distance_matrix_melt_filt_info$Day <- factor(distance_matrix_melt_filt_info$Day, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8 ,9 ,10 ,11 ,12 , 13, 14, 15, 16 ,17 ,18, 19, 20, 21, 22), labels = c(0, "1", "2", 3, 4, 5, "6","7",8 ,9 ,10 , 11, 12,"13", 14, 15,"16", "17", 18,19, 20,"21","22"))

#Plot
plotjsd <- ggplot(distance_matrix_melt_filt_info, aes(Day, value, color=Compartment_Shime, group = Compartment_Shime, shape = Shime))+
  theme_pubr(border = TRUE) +
  geom_point() +
  geom_point(size=6) +
  scale_shape_manual(values=c(20,18)) +
  labs(hjust = 0.5, x="Days", y = "1-Day β-diversity") +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position ="none") +
  scale_colour_manual(values = cols) +
  scale_x_discrete(limits=c("1", "2", " ", " ", " ", "6", "7"," "," "," ", " ", "13", " ","16", "17", " "," ","21", "22")) +
  geom_line(aes(linetype=Shime)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  ylim(0.0, 1) 
  
plotjsd

#Save
ggsave("JSD_DIS_all_1DayInterval_T.pdf", units="in", width=12, height=6, dpi=300)
dev.off()

###################################################################################################
####Fig2:A
###PCA Metabolites
library(MetaboAnalystR)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "Metabolites.tsv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "pdf", 300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "pdf", 300, width=NA)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "pdf", 300, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "pdf", 300, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "pdf", 300, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "pdf", 300, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "pdf", 300, width=NA, 1,2)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_1_", "psd", 300, width=NA, 1,2,0.0,1,0)
shapeVec<-c(16,16,15,17,16)
#cols = Same from 16s
mSet<-UpdateGraphSettings(mSet, cols, shapeVec)
shapeVec<-c(16,16,15,17,16)
mSet<-UpdateGraphSettings(mSet, cols, shapeVec)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_2_", "pdf", 300, width=7, 1,2,0.0,1,0)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_2_", "pdf", 300, width=3.5, 1,2,0.0,1,0)

###################################################################################################
####Supplementary
##Fig2SB
compAC_TC_DC <-c("S10A1.S10A2", "S12A12.S12A22","S15A1.S15A2", "S16A1.S16A2", "S17A12.S17A22", "S112A12.S112A22", "S113A1.S113A2", "S115A1.S115A2", "S115A1.S115A2", "S116A1.S116A2", "S117A12.S117A22", "S120A1.S120A2", "S121A1.S121A2", "S123A12.S123A22","S10T1.S10T2", "S12T12.S12T22","S15T1.S15T2", "S16T1.S16T2", "S17T12.S17T22", "S112T12.S112T22", "S113T1.S113T2", "S115T1.S115T2", "S115T1.S115T2", "S116T1.S116T2", "S117T12.S117T22", "S120T1.S120T2", "S121T1.S121T2", "S123T12.S123T22","S10D1.S10D2", "S12D12.S12D22","S15D1.S15D2", "S16D1.S16D2", "S17D12.S17D22", "S112D12.S112D22", "S113D1.S113D2", "S115D1.S115D2", "S115D1.S115D2", "S116D1.S116D2", "S117D12.S117D22", "S120D1.S120D2", "S121D1.S121D2", "S123D12.S123D22")

ggplot(distance_matrix_melt_filt_info, aes(Reactor, value, color=Compartment_Shime, group = Compartment_Shime))+
  theme_pubr(border = TRUE) +
  geom_boxplot() +
  labs(hjust = 0.5, x="Compartment", y = "Microbiome composition difference between matched samples for the two units,\n shown separately for the three compartments (Jensen-Shannon distance)") +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position ="none") +
  scale_colour_manual(values = cols) +
  ylim(0.0, 1) 
###################################################################################################
###Fig3S
##Triplicates
#Filter all samples with triplicates
sample_variables(pseq)

Shime_Plots_Only_Triplicates <- subset_samples(pseq, Triplicates =="A1" | Triplicates=="A2" | Triplicates=="A3")
Shime_Plots_Only_Triplicates

#Relative Abundance
Shime_Plots_Only_Triplicates_F <- transform_sample_counts(Shime_Plots_Only_Triplicates, function(otu) otu/sum(otu))
jsd_dist<- sqrt(phyloseq::distance(Shime_Plots_Only_Triplicates_F, "jsd"))
pcoa_jsd <- ordinate(Shime_Plots_Only_Triplicates_F, method = "PCoA", distance = jsd_dist)

exptT = prune_taxa(names(sort(taxa_sums(Shime_Plots_Only_Triplicates_F), TRUE)), Shime_Plots_Only_Triplicates_F)
ordT2<- sqrt(phyloseq::distance(Shime_Plots_Only_Triplicates_F,Type = "samples", "jsd"))
pcoa=ordinate(Shime_Plots_Only_Triplicates_F, "PCoA", distance=ordT2)

#
Shime_Plots_Only_Triplicates_F@sam_data$Compartment
a<-pcoa$values
b<-as.data.frame(pcoa$vectors)
samdf<-data.frame(unclass(Shime_Plots_Only_Triplicates_F@sam_data))
row.names(samdf)<-sample_names(Shime_Plots_Only_Triplicates_F)

#
class(jsd_dist)
distance_matrix<-matrix(jsd_dist)
distance_matrix<-as.matrix(jsd_dist)
distance_matrix<-as.data.frame(distance_matrix)
distance_matrix <-(distance_matrix)
distance_matrix

#
compAll_T <-c("F1.F2", "F2.F3", "F3.F1","S12D13.S12D11","S12D11.S12D12","S12D12.S12D13","S12D23.S12D21","S12D21.S12D22","S12D22.S12D23","S17D13.S17D11","S17D11.S17D12","S17D12.S17D13","S17D23.S17D21","S17D21.S17D22","S17D22.S17D23","S112D13.S112D11","S112D11.S112D12","S112D12.S112D13","S112D23.S112D21","S112D21.S112D22","S112D22.S112D23","S117D13.S117D11","S117D11.S117D12","S117D12.S117D13","S117D23.S117D21","S117D21.S117D22","S117D22.S117D23","S123D13.S123D11","S123D11.S123D12","S123D12.S123D13","S123D23.S123D21","S123D21.S123D22","S123D22.S123D23","S12A13.S12A11","S12A11.S12A12","S12A12.S12A13","S12A23.S12A21","S12A21.S12A22","S12A22.S12A23","S17A13.S17A11","S17A11.S17A12","S17A12.S17A13","S17A23.S17A21","S17A21.S17A22","S17A22.S17A23","S112A13.S112A11","S112A11.S112A12","S112A12.S112A13","S112A23.S112A21","S112A21.S112A22","S112A22.S112A23","S117A13.S117A11","S117A11.S117A12","S117A12.S117A13","S117A23.S117A21","S117A21.S117A22","S117A22.S117A23","S123A13.S123A11","S123A11.S123A12","S123A12.S123A13","S123A23.S123A21","S123A21.S123A22","S123A22.S123A23","S12T13.S12T11","S12T11.S12T12","S12T12.S12T13","S12T23.S12T21","S12T21.S12T22","S12T22.S12T23","S17T13.S17T11","S17T11.S17T12","S17T12.S17T13","S17T23.S17T21","S17T21.S17T22","S17T22.S17T23","S112T13.S112T11","S112T11.S112T12","S112T12.S112T13","S112T23.S112T21","S112T21.S112T22","S112T22.S112T23","S117T13.S117T11","S117T11.S117T12","S117T12.S117T13","S117T23.S117T21","S117T21.S117T22","S117T22.S117T23","S123T13.S123T11","S123T11.S123T12","S123T12.S123T13","S123T23.S123T21","S123T21.S123T22","S123T22.S123T23")

#
distance_matrix$samples<-row.names(distance_matrix)
distance_matrix_melt<-melt.data.frame(distance_matrix, id.vars = "samples")
distance_matrix_melt$Comparision<-paste(distance_matrix_melt$samples, distance_matrix_melt$variable, sep = ".")
distance_matrix_melt_filt<-subset(distance_matrix_melt, Comparision %in% compAll_T)
row.names(distance_matrix_melt_filt)<-as.character(distance_matrix_melt_filt$variable)
distance_matrix_melt_filt_info<-merge(distance_matrix_melt_filt, samdf, by=0)

#Factor
distance_matrix_melt_filt_info$Day<-as.factor(distance_matrix_melt_filt_info$Day)
distance_matrix_melt_filt_info$Compartment_Shime <- factor(distance_matrix_melt_filt_info$Compartment_Shime, levels = c("Fecal", "AC1", "AC2","TC1", "TC2","DC1", "DC2"))
distance_matrix_melt_filt_info$Compartment <- factor(distance_matrix_melt_filt_info$Compartment, levels = c("Fecal", "AC","TC","DC"))
distance_matrix_melt_filt_info$Shime <- factor(distance_matrix_melt_filt_info$Shime, levels = c("SHIME-1", "SHIME-2"))

###Fig3S:C
ggsave("Triplicates_Evaluation.pdf", units="in", width=14, height=8, dpi=300)

sup2 <- ggplot(subset(distance_matrix_melt_filt_info, Compartment =="TC" | Compartment == "DC" | Compartment == "AC"), aes(Day, value, color=Compartment_Shime, group = Compartment_Shime)) +
  geom_jitter(shape=16, size = 2, position=position_jitter(0.2), color = "black") +
  theme_pubr(border = TRUE) +
  facet_grid(Shime~Compartment) +
  theme(axis.text=element_text(size=14), 
        axis.text.x = element_text(size = 12, hjust = 0.5), 
        axis.title.y = element_text(size = 18, hjust = 1),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=0),
        legend.position="none",
        axis.title.x = element_text(size = 18), 
        strip.text.x = element_text(size = 20, face = "bold"),
        strip.text.y = element_text(size = 20, face = "bold")) +
  labs(x="Days", y = "Jensen-Shannon Distance (JSD)") +
  ylim(0, 1) 
print(sup2)

newSTorder = c("AC", "TC", "DC")
sup2$data$Compartment <- as.character(sup2$data$Compartment)
sup2$data$Compartment <- factor(sup2$data$Compartment, levels=newSTorder)
print(sup2)
dev.off()

###Fig3S:A,B
#Alpha Diversity Triplicates
alpha_diversity_Shime_T <- estimate_richness(Shime_Plots_Only_Triplicates,  measures = c("Observed", "Shannon"))
map <- sample_data(Shime_Plots_Only_Triplicates)
alpha <- cbind(alpha_diversity_Shime_T, map)
alpha$Day <- factor(alpha$Day)
alpha$Compartment_Shime <- factor(alpha$Compartment_Shime, levels = c("Fecal", "AC1", "TC1", "DC1", "AC2", "TC2", "DC2"))
alpha$Compartment <- factor(alpha$Compartment, levels = c("Fecal", "AC","TC","DC"))
alpha$Shime <- factor(alpha$Shime, levels = c("SHIME-1", "SHIME-2"))

#Plot
ggsave("Triplicates_Evaluation_Shannon.pdf", units="in", width=14, height=8, dpi=300)

ggplot(subset(alpha, Compartment =="TC" | Compartment == "DC" | Compartment == "AC"), aes(Day, Shannon, color=Compartment_Shime, group = Compartment_Shime)) +
  geom_jitter(shape=16, size = 2, position=position_jitter(0.2), color = "black") +
  theme_pubr(border = TRUE) +
  facet_grid(Shime~Compartment) +
  theme(axis.text=element_text(size=14), 
        axis.text.x = element_text(size = 12, hjust = 0.5), 
        axis.title.y = element_text(size = 18, hjust = 1),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=0),
        legend.position="none",
        axis.title.x = element_text(size = 18), 
        strip.text.x = element_text(size = 20, face = "bold"),
        strip.text.y = element_text(size = 20, face = "bold")) +
  labs(hjust = 0.5, face= "bold", x="Days ", y = "Shannon Diversity") +
  ylim(1, 4) 
dev.off()

##Plot
ggsave("Triplicates_Evaluation_Richness.pdf", units="in", width=14, height=8, dpi=300)

ggplot(subset(alpha, Compartment =="TC" | Compartment == "DC" | Compartment == "AC"), aes(Day, Observed, color=Compartment_Shime, group = Compartment_Shime)) +
  geom_jitter(shape=16, size = 2, position=position_jitter(0.2), color = "black") +
  theme_pubr(border = TRUE) +
  facet_grid(Shime~Compartment) +
  theme(axis.text=element_text(size=14), 
        axis.text.x = element_text(size = 12, hjust = 0.5), 
        axis.title.y = element_text(size = 18, hjust = 1),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=0),
        legend.position="none",
        axis.title.x = element_text(size = 18), 
        strip.text.x = element_text(size = 20, face = "bold"),
        strip.text.y = element_text(size = 20, face = "bold")) +
  labs(face= "bold", x="Days ", y = "Richness") +
  ylim(0, 120) 
dev.off()

#################
#Comp All Compartments, Beta diversity

##  Comparison para mudar
comp_all_Shime2_AC2_TC2_DC2_AC1TC1 <-c("S10A2.S10D2", "S12A22.S12D22","S15A2.S15D2", "S16A2.S16D2", "S17A22.S17D22", "S112A22.S112D22", "S113A2.S113D2", "S115A2.S115D2", "S115A2.S115D2", "S116A2.S116D2", "S117A22.S117D22", "S120A2.S120D2", "S121A2.S121D2", "S123A22.S123D22","S10A2.S10T2", "S12A22.S12T22","S15A2.S15T2", "S16A2.S16T2", "S17A22.S17T22", "S112A22.S112T22", "S113A2.S113T2", "S115A2.S115T2", "S115A2.S115T2", "S116A2.S116T2", "S117A22.S117T22", "S120A2.S120T2", "S121A2.S121T2", "S123A22.S123T22","S10T1.S10D1", "S12T12.S12D12","S15T1.S15D1", "S16T1.S16D1", "S17T12.S17D12", "S112T12.S112D12", "S113T1.S113D1", "S115T1.S115D1", "S115T1.S115D1", "S116T1.S116D1", "S117T12.S117D12", "S120T1.S120D1", "S121T1.S121D1", "S123T12.S123D12")


ggplot(distance_matrix_melt_filt_infoB, aes(Day, value, color=Compartment_Shime, group = Compartment_Shime))+
  theme_pubr(border = TRUE) +
  geom_point(size=6) +
  scale_shape_manual(values=c(20,18)) +
  labs(hjust = 0.5, x="Days", y = "Microbiome composition difference between all compartments\n (Jensen-Shannon distance") +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position ="none") +
  scale_colour_manual(values = TC_DC) +
  scale_x_discrete(limits=c(1, 2, " ", " ", 5, 6, 7," "," "," ", " ", 12, 13," ", 15, 16, 17," "," ", 20, 21, 22)) +
  geom_line(aes(linetype=Reactor3)) +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  ylim(0.0, 1) 
####################################################################################################################################
#FigS4
#comp_all_Shime2_AC1_TC1_DC1_AC1TC1

#Plot
ggplot(distance_matrix_melt_filt_infoB, aes(Day, value, color=Compartment_Shime, group = Compartment_Shime))+
  theme_pubr(border = TRUE) +
  geom_point(size=6) +
  scale_shape_manual(values=c(20,18)) +
  labs(hjust = 0.5, x="Days", y = "Microbiome composition difference between all compartments\n (Jensen-Shannon distance") +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position ="none") +
  scale_colour_manual(values = TC_DC) +
  scale_x_discrete(limits=c(1, 2, " ", " ", 5, 6, 7," "," "," ", " ", 12, 13," ", 15, 16, 17," "," ", 20, 21, 22)) +
  geom_line(aes(linetype=Reactor3)) +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  ylim(0.0, 1) 

dev.off()

####################################################################################################################################
###Fig5s:
#comp
comp_AC_TC_DC_Fecal <-c("F2.S10A1", "F2.S10A2", "F2.S12A12", "F2.S12A22","F2.S15A1", "F2.S15A2", "F2.S16A1", "F2.S16A2", "F2.S17A12", "F2.S17A22", "F2.S112A12", "F2.S112A22", "F2.S113A1", "F2.S113A2", "F2.S115A1", "F2.S115A2", "F2.S115A1", "F2.S115A2", "F2.S116A1", "F2.S116A2", "F2.S117A12", "F2.S117A22", "F2.S120A1", "F2.S120A2", "F2.S121A1", "F2.S121A2", "F2.S123A12", "F2.S123A22","F2.S10T1", "F2.S10T2", "F2.S12T12", "F2.S12T22","F2.S15T1", "F2.S15T2", "F2.S16T1", "F2.S16T2", "F2.S17T12", "F2.S17T22", "F2.S112T12", "F2.S112T22", "F2.S113T1", "F2.S113T2", "F2.S115T1", "F2.S115T2", "F2.S115T1", "F2.S115T2", "F2.S116T1", "F2.S116T2", "F2.S117T12", "F2.S117T22", "F2.S120T1", "F2.S120T2", "F2.S121T1", "F2.S121T2", "F2.S123T12", "F2.S123T22","F2.S10D1","F2.S10D2", "F2.S12D12", "F2.S12D22","F2.S15D1", "F2.S15D2", "F2.S16D1", "F2.S16D2", "F2.S17D12", "F2.S17D22", "F2.S112D12", "F2.S112D22", "F2.S113D1", "F2.S113D2", "F2.S115D1", "F2.S115D2", "F2.S115D1", "F2.S115D2", "F2.S116D1", "F2.S116D2", "F2.S117D12", "F2.S117D22", "F2.S120D1", "F2.S120D2", "F2.S121D1", "F2.S121D2", "F2.S123D12", "F2.S123D22")

#Plot
ggplot(distance_matrix_melt_filt_info, aes(Day, value, color=Compartment_Shime, group = Compartment_Shime, shape = Shime))+
  theme_pubr(border = TRUE) +
  geom_point(size=6) +
  scale_shape_manual(values=c(20,18)) +
  labs(hjust = 0.5, x="Days", y = "Microbiome composition difference between all collected samples \n to the fecal inoculum (Jensen-Shannon distance)") +
  theme(axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.x.bottom = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position ="none") +
  scale_colour_manual(values = cols) +
  scale_x_discrete(limits=c(1, 2, " ", " ", 5, 6, 7," "," "," ", " ", 12, 13," ", 15, 16, 17," "," ", 20, 21, 22)) +
  geom_line(aes(linetype=Shime)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  ylim(0.0, 1) 

dev.off()

