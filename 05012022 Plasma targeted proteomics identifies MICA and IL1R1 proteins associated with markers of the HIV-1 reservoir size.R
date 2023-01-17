## Analysis script Plasma targeted proteomics identifies MICA and IL1R1 proteins associated with markers of the HIV-1 reservoir size
# Readme
# Created by: Marc Blaauw (marc.blaauw@radboudumc.nl)
# Created on: 05-01-2023
# Contact at: Vasiliki Matzaraki (vasiliki.matzaraki@radboudumc.nl)
# Principal Investigator: Andre van der Ven 

# Load libraries
library(limma)
library(calibrate)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Hmisc)

# Load datasets
hiv200metadata<-read_delim("hiv200metadata.csv")
HIVreservoirdataset<- read_delim("HIVreservoirdataset.csv")

# Import metadata 200HIV and HIVreservoir dataset
hiv200metadata<-read_delim("hiv200metadata.csv")
HIVreservoirdataset<- read_delim("HIVreservoirdataset.csv")
HIVreservoirdataset<-HIVreservoirdataset[c(2:8)]

## Figure 1 - correlation plot
clinicalplot<-hiv200metadata[c(1,4,5,14,19,20,54)]
correlationfinal<-merge(clinicalplot, HIVreservoirinterest, by= "record_id")
correlationfinal<-na.omit(correlationfinal)
colnames(correlationfinal)[2] <- "Sex" 
colnames(correlationfinal)[3] <- "Age"
colnames(correlationfinal)[4] <- "HIV duration"
colnames(correlationfinal)[5] <- "CD4 nadir"
colnames(correlationfinal)[6] <- "CD4 latest"
colnames(correlationfinal)[7] <- "cART duration"
colnames(correlationfinal)[8] <- "CA HIV RNA"
colnames(correlationfinal)[9] <- "CA HIV DNA" 
correlationfinal2<-correlationfinal[c(2:9)]

cor_mat <- rcorr(as.matrix(correlationfinal2), type = "spearman")
cor_p_mat <-cor_mat$P                 
cor_p_mat[upper.tri(cor_p_mat)] <- NA

cor_r_mat<-cor_mat$r
cor_r_mat[lower.tri(cor_r_mat)] <- NA

cor_p_mat <- as.data.frame(reshape2::melt(cor_p_mat, na.rm = TRUE))
cor_r_mat <-as.data.frame(reshape2::melt(cor_r_mat, na.rm = TRUE))
cor_p_mat$value <- p.adjust(cor_p_mat$value, method = "fdr", n = 28)

myfmt <- function(v) {
  ifelse(v<1e-16, "<1e-16", ifelse(v<0.001, format(v, digits=3, scientific=T), round(v, 3)))
}

tiff(file = "Figure1_correlationplot_proteomics_and_reservoir.tiff", units= "in", height=8 , width=10, res = 600)
ggplot() + 
  geom_tile(aes(x= Var1, y = Var2, fill = value), data = cor_r_mat) + 
  geom_text(aes(x= Var1, y = Var2, label = round(value,3)), data = cor_r_mat, color = "black", size=4) +
  geom_text(aes(x= Var1, y= Var2, label = myfmt(value)), data = cor_p_mat, color = "black", size = 4) + 
  scale_fill_gradient2(low = "dodgerblue", high = "firebrick",midpoint = 0, limit=c(-1,1), name= "Spearman's Rho") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "top",
        legend.direction = "horizontal", 
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1.0, size = 11),
        axis.text.y = element_text(size = 11))
dev.off()

-------------------------------------------------------------------------------------------------------------
### PCA analysis and detection of outliers using the OLINK panel
  
OLINKproteinsnew<-read_xlsx("200HIV_Explore_QC_210_1293.xlsx")
OLINKreservoir<-read_xlsx("OLINKreservoirclinic.xlsx")
colnames(OLINKproteinsnew)[1] <- "record_id"  #change first column name 
OLINKreservoirnew<- OLINKreservoir[c(1:3)]  #only select reservoir measurements
OLINKreservoirnew2<- merge(OLINKproteinsnew, OLINKreservoirnew, by="record_id")

## -------------------- PERFORM PCA ANALYSIS ---------------------------------------------------------------------------------------------------
## --------USING THE 200HIV SAMPLES ----------------------------------------------------------------------------------------------------------
mydata.olink.pca <- prcomp(OLINKreservoirnew2[c(2:1294)], center = TRUE, scale. = TRUE) ### add only the column numbers with the protein measurements

# Plot the PC1 vs PC2
png(file = "PCA_plot_proteomic_analsis_200HIV_with_outliers.png", units= "cm", height=8 , width=12, res = 300)
ggbiplot(mydata.olink.pca,ellipse=FALSE,obs.scale = 0.8, var.scale = 0.5, var.axes = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1))
dev.off()

## check which subjects are outliers
# The standard way to detect outliers is the criterion of being more than 3 standard deviations away from the mean.
U_mydata <- mydata.olink.pca$x
mydata_outliers <- apply(U_mydata, 2, function(x) which( abs(x - mean(x)) > (3 * sd(x)) ))
head(mydata_outliers)

##coloring outliers
Outlierlist<-read_xlsx("Olinkoutlierlist.xlsx")
OLinkreservoiroutlierlist<-merge(OLINKreservoirnew2, Outlierlist, by = "record_id")

png(file = "supplementary_figure1a.png", units= "cm", height=12 , width=12, res = 300)
ggbiplot(mydata.olink.pca,ellipse=FALSE,  labels=rownames(OLINKreservoirnew2$record_id), groups = OLinkreservoiroutlierlist$Outlier, obs.scale = 0.5, var.scale = 0.5, var.axes = FALSE)+
  theme_bw()+
  scale_colour_manual(labels = c("No outlier", "Outlier"),  values =c("black", "red")) +
  geom_point(aes(colour = OLinkreservoiroutlierlist$Outlier), size =0.4) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1))
dev.off()

## PCA analysis after removing outliers
OLINKproteinsnewnooutliers<-read_xlsx("200HIV.xlsx")
colnames(OLINKproteinsnewnooutliers)[1] <- "record_id"
OLINKreservoirnewnooutliers<- merge(OLINKproteinsnewnooutliers, OLINKreservoirnew, by="record_id")
OLINKreservoirnewnooutliers<-OLINKreservoirnewnooutliers[-2]

no_outliers_mydata.pca <- prcomp(OLINKreservoirnewnooutliers[c(2:1294)], center = TRUE, scale. = TRUE)

# Plot the PC1 vs PC2 
png(file = "Supplementary_figure1b.png", units = "cm", height = 12, width = 12, res = 300)
ggbiplot(no_outliers_mydata.pca,ellipse=FALSE,obs.scale = 0.5, var.scale = 0.5, var.axes = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1))
dev.off()
-------------------------------------------------------------------------------------------------------
# Association analysis between plasma proteomics and CA HIV DNA and RNA
# model 1: showed below: RNA or DNA ~ proteome + Age + Sex
# model 2:  RNA or DNA ~ proteome + Age + Sex + CD4 nadir + HIV duration
  
metadataHIV200<-read_xlsx("HIV200_dataset.xlsx")
## filter relevant variables
HIV200<-metadataHIV200[c(1, 4, 5, 34, 41, 42, 79)]
## add reservoirdata
OLINKreservoirdata <- HIV200 %>% left_join(HIVreservoirdataset)

# Association analysis between plasma proteomics and CA HIV RNA (continuous variable)

# Import OLINK data from the 200HIV and then create matrices
OLINKreservoirRNA<-OLINKreservoirdata[c(1, 10, 2, 3)]
OLINKreservoirRNA<-na.omit(OLINKreservoirRNA)
OLINKreservoirRNAwithproteins<-merge(OLINKreservoirRNA, OLINKproteinsnewnooutliers, by="record_id")
OLINKreservoirRNA<-OLINKreservoirRNAreservoirRNAwithproteins[-5]

#creating a matrix with different variables to correct for
AgeR <- OLINKreservoirRNA$AGE
SexR <- OLINKreservoirRNA$GEN
RNA <-OLINKreservoirRNA$LOG_RNA

design_covRNA <- model.matrix(~as.numeric(RNA) + as.numeric(AgeR) + as.numeric(SexR))
colnames(design_covRNA) <-c("Participant", "RNAreservoirsize" , "Age", "Sex")

# Generate the protein expression matrix ----------------------------------------------------------------------------
# Transform NPX data for proteins to represent rows - no group names
OLINKreservoirRNAp <- OLINKreservoirRNA[,c(-1:-4)] #keep only proteins
OLINkreservoirRNAp_tr <- as.matrix(t(OLINkreservoirRNAp))  # Transpose data

# Fitting models for associated genes----------------------------------------------------------------
fit1_covRNA <-lmFit(OLINkreservoirRNAp_tr, design_covRNA)
fit2_covRNA <- eBayes(fit1_covRNA)

DEP_covRNA <- topTable(fit2_covRNA, coef = "RNAreservoirsize", adjust = "BH", number = 1293)

# Save results and maken volcano plot
write.csv(DEP_covRNA, "results diff analyse total OLINK panel for CA HIV RNA Age and Sex corrected.csv")
limma_covR <- read.csv("results diff analyse total OLINK panel for CA HIV RNA Age and Sex corrected.csv", check.names = FALSE)
colnames(limma_covR)[1] <- 'protein'

# Figure2b
tiff(file = "Figure 2b OLINK Volcanoplot RNA.tiff", units= "in", height=10 , width=10, res = 150)
ggplot(limma_covR, aes(x = logFC, y = -log10(P.Value),label = protein)) +
  geom_text_repel(data = subset(limma_cov, P.Value < 0.05 & logFC > 0.02 | logFC <  -0.3), box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  geom_point(color = ifelse(limma_cov$P.Value > 0.05, "black", "red"), size = 2) +
  theme_bw() +
  xlab("logFC") +
  ylab("-log10(P value)") +
  xlab("Effect size")+
  geom_vline(
    xintercept = 0,
    col = "gray50",
    linetype = "dotted",
    size = 1.5
  ) +
  geom_hline(
    yintercept = 1.3,
    col = "gray",
    linetype = "dotted",
    size = 1.5)+
  annotation_custom(text_high,xmin=1,xmax=0.3,ymin=-0.04,ymax=-0.04) + 
  annotation_custom(text_low,xmin=-0.5,xmax=-0.3,ymin=-0.04,ymax=-0.04) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    plot.title = element_text(color = "black", size = 14, face = "bold"),
    plot.caption = element_text(color = "black", face = "italic"))+
  theme( panel.grid.major = element_blank(),text = element_text(size=14),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=14))
dev.off()

---------------------------------------------------------------

---------------------------------------------------------------
  
