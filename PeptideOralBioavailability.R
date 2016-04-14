# Working directory 
setwd("/Users/fabienplisson/Desktop/Github shares/Peptide and Oral Bioavailability 2014")

# Libraries
library(ggplot2)
library(Rmisc)
library(lattice)
library(plyr)
library(gclus)
library(TSP)
library(cluster)
library(seriation)
library(corrgram)
library(grid)
library(reshape2)
library(RColorBrewer)

# Distribution Analysis of Molecular Properties for orally absorbed cyclic peptides 
Data <- read.csv("./MolecularProperties_distribution.csv", header=TRUE)
attach(Data)
colnames(Data) <- c("Compound", "F", "MolWeight", "HBD", "HBA", "NandO", "logPmean", "logPsd", "logSmean", "logSsd", "NPSAmean", "NPSAsd", "PSAmean", "PSAsd", "Rgyrmean", "Rgyrsd", "Gsolvmean", "Gsolvsd", "RotBonds", "Amine","Amidine and Guanidine", "CarboxylicAcid", "Amide")

#Subset dataset with only peptides with measured oral availability in rats (F%)
Data.sub <- subset(Data, is.finite(F)) # is.finite() opposed to is.na(variable)

# Distributions of various PC properties (features) against measured oral availability in rats (F%)
## MOLECULAR WEIGHT
MW <- ggplot(Data.sub, aes(x=Data.sub$MolWeight, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Molecular Weight (Da)", y="Oral Bioavailability in Rats (F%)") + xlim(400, 1600) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=500, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=500, colour="red", linetype="longdash")
# + geom_text(x=500, label="500", y=0, colour="red", text=element_text(size=12))
MWvw <- ggplot(Data.sub, aes(x=Data.sub$MolWeight, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Molecular Weight (Da)", y="Oral Bioavailability in Rats (F%)") + xlim(400, 1600) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=600, xmax=1200, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=c(600, 1200), colour="lightblue", linetype="longdash")

## HYDROGEN BOND DONORS
HBD <- ggplot(Data.sub, aes(x=Data.sub$HBD, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Hydrogen Bond Donors", y="Oral Bioavailability in Rats (F%)") + xlim(0, 20) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=5, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=5, colour="red", linetype="longdash")
HBDvw <- ggplot(Data.sub, aes(x=Data.sub$HBD, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Hydrogen Bond Donors", y="Oral Bioavailability in Rats (F%)") + xlim(0, 20) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=12, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=12, colour="lightblue", linetype="longdash")

## HYDROGEN BOND ACCEPTORS
HBA <- ggplot(Data.sub, aes(x=Data.sub$HBA, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Hydrogen Bond Acceptors", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=10, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=10, colour="red", linetype="longdash")
HBAvw <- ggplot(Data.sub, aes(x=Data.sub$HBA, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Hydrogen Bond Acceptors", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=12, xmax=16, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=c(12,16), colour="lightblue", linetype="longdash")

## SUM OF HBD and HBA
Data.sub <- transform(Data.sub, TotHB=Data.sub$HBD + Data.sub$HBA)
THB <- ggplot(Data.sub, aes(x=Data.sub$TotHB, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Total of Hydrogen Bond Groups", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=20, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=20, colour="red", linetype="longdash")

## SUM OF NITROGEN AND OXYGEN ATOMS
NO <- ggplot(Data.sub, aes(x=Data.sub$NandO, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Sum of N and O", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=10, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=10, colour="red", linetype="longdash")

## OCTANOL-WATER PARTITION COEFFICIENT logP (QlogPo/w)
logP <- ggplot(Data.sub, aes(x=Data.sub$logPmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="logP", y="Oral Bioavailability in Rats (F%)") + xlim(-5, 10) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=0, xmax=5, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=c(0,5), colour="red", linetype="longdash")
logP2 <- ggplot(Data.sub, aes(x=Data.sub$logPmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="logP", y="Oral Bioavailability in Rats (F%)") + xlim(-5, 10) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=0, xmax=5, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=c(0,5), colour="red", linetype="longdash") + geom_errorbarh(aes(xmin=Data.sub$logPmean-Data.sub$logPsd, xmax=Data.sub$logPmean+Data.sub$logPsd), height=1)
logP2vw <- ggplot(Data.sub, aes(x=Data.sub$logPmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="logP", y="Oral Bioavailability in Rats (F%)") + xlim(-5, 10) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-2, xmax=6, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=c(-2,6), colour="lightblue", linetype="longdash") + geom_errorbarh(aes(xmin=Data.sub$logPmean-Data.sub$logPsd, xmax=Data.sub$logPmean+Data.sub$logPsd), height=1)

## SOLUBILITY COEFFICIENT logS (QlogS)
logS <- ggplot(Data.sub, aes(x=Data.sub$logSmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="logS", y="Oral Bioavailability in Rats (F%)") + xlim(-12, 6) + ylim(0,100) + theme(axis.text=element_text(size=12)) 
logS2 <- ggplot(Data.sub, aes(x=Data.sub$logSmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="logS", y="Oral Bioavailability in Rats (F%)") + xlim(-12, 6) + ylim(0,100) + theme(axis.text=element_text(size=12)) + geom_errorbarh(aes(xmin=Data.sub$logSmean-Data.sub$logSsd, xmax=Data.sub$logSmean+Data.sub$logSsd), height=1)

## ROTATABLE BONDS
RB <- ggplot(Data.sub, aes(x=Data.sub$RotBonds, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Rotatable Bonds", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=10, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=10, colour="red", linetype="longdash")
RBvw <- ggplot(Data.sub, aes(x=Data.sub$RotBonds, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Rotatable Bonds", y="Oral Bioavailability in Rats (F%)") + xlim(0, 40) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=15, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=15, colour="lightblue", linetype="longdash")

## POLAR SURFACE AREA
PSA <- ggplot(Data.sub, aes(x=Data.sub$PSAmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="3D Polar Surface Area FISA (A2)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 400) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=140, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=140, colour="red", linetype="longdash")
PSA2 <- ggplot(Data.sub, aes(x=Data.sub$PSAmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="3D Polar Surface Area FISA (A2)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 400) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=-Inf, xmax=140, ymin=-Inf, ymax=Inf, alpha=0.1, fill="red") + geom_vline(xintercept=140, colour="red", linetype="longdash") + geom_errorbarh(aes(xmin=Data.sub$PSAmean-Data.sub$PSAsd, xmax=Data.sub$PSAmean+Data.sub$PSAsd), height=1)
PSA2vw <- ggplot(Data.sub, aes(x=Data.sub$PSAmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="3D Polar Surface Area FISA (A2)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 400) + ylim(0,100) + theme(axis.text=element_text(size=12)) + annotate("rect", xmin=180, xmax=320, ymin=-Inf, ymax=Inf, alpha=0.1, fill="lightblue") + geom_vline(xintercept=c(180, 320), colour="lightblue", linetype="longdash") + geom_errorbarh(aes(xmin=Data.sub$PSAmean-Data.sub$PSAsd, xmax=Data.sub$PSAmean+Data.sub$PSAsd), height=1)

## NON POLAR SURFACE AREA
NPSA <- ggplot(Data.sub, aes(x=Data.sub$NPSAmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="3D Non Polar Surface Area FOSA (A2)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 1000) + ylim(0,100) + theme(axis.text=element_text(size=12))
NPSA2 <- ggplot(Data.sub, aes(x=Data.sub$NPSAmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="3D Non Polar Surface Area FOSA (A2)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 1000) + ylim(0,100) + theme(axis.text=element_text(size=12)) + geom_errorbarh(aes(xmin=Data.sub$NPSAmean-Data.sub$NPSAsd, xmax=Data.sub$NPSAmean+Data.sub$NPSAsd), height=1)

## RADIUS OF GYRATION
RGYR <- ggplot(Data.sub, aes(x=Data.sub$Rgyrmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Radius of gyration (A)", y="Oral Bioavailability in Rats (F%)") + xlim(2, 8) + ylim(0,100) + theme(axis.text=element_text(size=12))
RGYR2 <- ggplot(Data.sub, aes(x=Data.sub$Rgyrmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Radius of gyration (A)", y="Oral Bioavailability in Rats (F%)") + xlim(2, 8) + ylim(0,100) + theme(axis.text=element_text(size=12)) + geom_errorbarh(aes(xmin=Data.sub$Rgyrmean-Data.sub$Rgyrsd, xmax=Data.sub$Rgyrmean+Data.sub$Rgyrsd), height=1)

## ENERGY OF SOLVATION
Gsolv <- ggplot(Data.sub, aes(x=Data.sub$Gsolvmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Solvation Energy (kcal/mol)", y="Oral Bioavailability in Rats (F%)") + xlim(50, 150) + ylim(0,100) + theme(axis.text=element_text(size=12))
Gsolv2 <- ggplot(Data.sub, aes(x=Data.sub$Gsolvmean, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Solvation Energy (kcal/mol)", y="Oral Bioavailability in Rats (F%)") + xlim(50, 150) + ylim(0,100) + theme(axis.text=element_text(size=12)) + geom_errorbarh(aes(xmin=Data.sub$Gsolvmean-Data.sub$Gsolvsd, xmax=Data.sub$Gsolvmean+Data.sub$Gsolvsd), height=1)

## AMIDE GROUPS
AM <- ggplot(Data.sub, aes(x=Data.sub$Amide, y=Data.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Amide groups (Size)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 20) + ylim(0,100) + theme(axis.text=element_text(size=12))

multiplot(MW, HBD, HBA, logP, NO, RB, cols=2)

# Correlation between descriptors
corrgram(Data.sub, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt, main="Correlation between Descriptors")

# Feature engineering: Flexibility with RMSD values
# Scatterplots Mean Backbone RMSD values (conformations in CHCl3 vs H2O)
# and plot against relative difference for specific molecular properties; PSA, NPSA, SASA, Volume, Globularity, R-gyration
myRMSD <- read.csv("./RelativeDifferenceProperty_Flexibility.csv", header=TRUE)
attach(myRMSD)
colnames(myRMSD) <- c("Compound", "Number","MedianRMSD", "MeanRMSD", "StdErrRMSD", "RelDiffSASA", "RelDiffNPSA", "RelDiffPSA", "RelDiffVol","RelDiffGlob", "ReldiffLogP", "RelDiffRgyr", "F")
myRMSD.sub <- subset(myRMSD, is.finite(F))
attach(myRMSD.sub)
colnames(myRMSD.sub) <- c("Compound", "Number", "MedianRMSD", "MeanRMSD", "StdErrRMSD", "RelDiffSASA", "RelDiffNPSA", "RelDiffPSA", "RelDiffVol","RelDiffGlob", "ReldiffLogP", "RelDiffRgyr", "F")

# Can create a specific dataframe for each plot
RMSD.df <- data.frame(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffPSA, z=myRMSD$F)

# Oral bioavailability F%
FRmsd <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$F)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Oral Bioavailability in Rats (F%)") + xlim(0, 7) + ylim(0,100) + theme(axis.text=element_text(size=12)) + geom_errorbarh(aes(xmin=myRMSD.sub$MeanRMSD-myRMSD.sub$StdErrRMSD, xmax=myRMSD.sub$MeanRMSD+myRMSD.sub$StdErrRMSD), height=1)
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/F% RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(FRmsd)
dev.off()

FRmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$F)) + theme_bw() + geom_text(label=myRMSD.sub$Number, size=5) + labs(x="Backbone RMSD value (A)", y="Oral Bioavailability in Rats (F%)")  + xlim(0, 7) + ylim(0,100) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/F% RMSD Sub dataset + Numbers 27072014.pdf", width=8.27, height=5.83)
plot(FRmsd2)
dev.off()

## SASA
Sasa.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffSASA)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - SASA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/SASA RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(Sasa.Rmsd)
dev.off()

Sasa.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffSASA)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - SASA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/SASA RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Sasa.Rmsd2)
dev.off()

Sasa.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffSASA)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - SASA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/SASA RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(Sasa.Rmsd3)
dev.off()


## PSA
Psa.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffPSA)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - PSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/PSA RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(Psa.Rmsd)
dev.off()

Psa.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffPSA)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - PSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/PSA RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Psa.Rmsd2)
dev.off()

Psa.Rmsd3 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffPSA)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - PSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/PSA RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Psa.Rmsd2)
dev.off()

Psa.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffPSA)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - PSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/PSA RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(Psa.Rmsd3)
dev.off()

Psa.Rmsd4 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffPSA)) + theme_bw() + geom_point(aes(colour=myRMSD$F), shape=19, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - PSA (%)") + xlim(0,6) + ylim(0,1) + theme(axis.text=element_text(size=12)) + scale_color_gradient("F%", low="yellow", high="red")

## NPSA
NPsa.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffNPSA)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - NPSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/NPSA RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(NPsa.Rmsd)
dev.off()

NPsa.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffNPSA)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - NPSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/NPSA RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(NPsa.Rmsd2)
dev.off()

NPsa.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffNPSA)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - NPSA (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/NPSA RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(NPsa.Rmsd3)
dev.off()

## Volume
Vol.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffVol)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Volume (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Volume RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(Vol.Rmsd)
dev.off()

Vol.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffVol)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Volume (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Volume RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Vol.Rmsd2)
dev.off()

Vol.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffVol)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Volume (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Volume RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(Vol.Rmsd3)
dev.off()

## Globularity
Glob.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffGlob)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Globularity (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Globularity RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(Glob.Rmsd)
dev.off()

Glob.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffGlob)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Globularity (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Globularity RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Glob.Rmsd2)
dev.off()

Glob.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffGlob)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Globularity (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Globularity RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(Glob.Rmsd3)
dev.off()

## QlogP o/w
LogP.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$ReldiffLogP)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd)
dev.off()

LogP.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$ReldiffLogP)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd2)
dev.off()

LogP.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$ReldiffLogP)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd3)
dev.off()

LogP.Rmsd4 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$ReldiffLogP)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(-1,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP2 RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd4)
dev.off()

LogP.Rmsd5 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$ReldiffLogP)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(-1,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP2 RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd5)
dev.off()

LogP.Rmsd6 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$ReldiffLogP)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - QLogPo/w (%)") + xlim(0, 6) + ylim(-1,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/LogP2 RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(LogP.Rmsd6)
dev.off()

## R-gyration
Rgyr.Rmsd <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffRgyr)) + theme_bw() + geom_point(shape=21, fill="darkgreen", size=3, na.rm=TRUE) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Rgyr (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Rgyr RMSD Full dataset 22072014.pdf", width=8.27, height=5.83)
plot(Rgyr.Rmsd)
dev.off()

Rgyr.Rmsd2 <- ggplot(myRMSD.sub, aes(x=myRMSD.sub$MeanRMSD, y=myRMSD.sub$RelDiffRgyr)) + theme_bw() + geom_point(shape=21, fill="blue", size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Rgyr (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Rgyr RMSD Sub dataset 22072014.pdf", width=8.27, height=5.83)
plot(Rgyr.Rmsd2)
dev.off()

Rgyr.Rmsd3 <- ggplot(myRMSD, aes(x=myRMSD$MeanRMSD, y=myRMSD$RelDiffRgyr)) + theme_bw() + geom_text(label=myRMSD$Number, size=3) + labs(x="Backbone RMSD value (A)", y="Relative Difference - Rgyr (%)") + xlim(0, 6) + ylim(0,1) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Rgyr RMSD Labels 22072014.pdf", width=8.27, height=5.83)
plot(Rgyr.Rmsd3)
dev.off()

# Build Ternary plot from Moments of Inertia to evaluate Peptide geometry space 22072014
# Load dataset
MoI.df <- read.csv("./MomentOfInertia_21072014.csv", header=TRUE)
attach(MoI.df)
colnames(MoI.df) <- c("Compound", "Number","F", "MeanRMSD", "XH2O", "YH2O", "XCHCl3", "YCHCl3") 

#in H2O
p1 <- ggplot(MoI.df, aes(x=MoI.df$XH2O, y=MoI.df$YH2O)) + theme_bw() + geom_point(shape=21, fill="red", size=3, na.rm=TRUE) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Molecular Shape Diversity in H2O 22072014.pdf", width=8.27, height=5.83)
plot(p1)
dev.off()

p1i <- ggplot(MoI.df, aes(x=MoI.df$XH2O, y=MoI.df$YH2O)) + theme_bw() + geom_text(label=MoI.df$Number, size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Molecular Shape Diversity in H2O 22072014 + Numbers.pdf", width=8.27, height=5.83)
plot(p1i)
dev.off()

#in CHCl3
p2 <- ggplot(MoI.df, aes(x=MoI.df$XCHCl3, y=MoI.df$YCHCl3)) + theme_bw() + geom_point(shape=21, fill="blue", size=3, na.rm=TRUE) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Molecular Shape Diversity in CHCl3 22072014.pdf", width=8.27, height=5.83)
plot(p2)
dev.off()

p2i <- ggplot(MoI.df, aes(x=MoI.df$XCHCl3, y=MoI.df$YCHCl3)) + theme_bw() + geom_text(label=MoI.df$Number, size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Molecular Shape Diversity in CHCl3 22072014 + Numbers.pdf", width=8.27, height=5.83)
plot(p2i)
dev.off()

# In CHCl3 with values F% in rats
p3 <- ggplot(MoI.df, aes(x=MoI.df$XCHCl3, y=MoI.df$YCHCl3)) + theme_bw() + geom_point(aes(colour=MoI.df$F), shape=19, size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12)) + scale_color_gradient("F%", low="yellow", high="red")
pdf(file="/Users/f.plisson/Desktop/OrallyAbsorbedPeptides/Figures/Molecular Shape Diversity in CHCl3 22072014.pdf", width=8.27, height=5.83)
plot(p3)
dev.off()

# Tracing individual geometric / conformational change from CHCl3 to H2O environment with arrows
x <- c(MoI.df$XCHCl3)
y <- c(MoI.df$YCHCl3)
xend <- c(MoI.df$XH2O)
yend <- c(MoI.df$YH2O)
p4 <- ggplot() + geom_segment(MoI.df, mapping=aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(angle=25, length=unit(0.10, "inches")), size=0.5, colour="black") + theme_bw() + geom_point(MoI.df, mapping=aes(x=x, y=y), shape=21, fill="blue", size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
p4
p5 <- ggplot() + geom_segment(MoI.df, mapping=aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(angle=25, length=unit(0.15, "inches")), size=0.5, colour="black") + theme_bw() + geom_point(MoI.df, mapping=aes(x=x, y=y), shape=21, fill="blue", size=3) + geom_point(MoI.df, mapping=aes(x=xend, y=yend), shape=21, fill="red", size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12))
p5
p6 <- ggplot() + geom_segment(MoI.df, mapping=aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(angle=25, length=unit(0.10, "inches")), size=0.5, colour="black") + theme_bw() + geom_point(MoI.df, mapping=aes(x=x, y=y, colour=MoI.df$F), shape=19, size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12)) + scale_color_gradient("F%", low="yellow", high="red")
p6

# Just subset with F% values and Molecular Shape Diversity Change
MoI.df2 <- subset(MoI.df, is.finite(F))
x2 <- c(MoI.df2$XCHCl3)
y2 <- c(MoI.df2$YCHCl3)
xend2 <- c(MoI.df2$XH2O)
yend2 <- c(MoI.df2$YH2O)
p7 <- ggplot() + geom_segment(MoI.df2, mapping=aes(x=x2, y=y2, xend=xend2, yend=yend2), arrow=arrow(angle=25, length=unit(0.10, "inches")), size=0.5, colour="black") + theme_bw() + geom_point(MoI.df2, mapping=aes(x=x2, y=y2, colour=MoI.df2$F), shape=19, size=3) + labs(x="I1/I3", y="I2/I3") + xlim(0.00, 1.00) + ylim(0.50,1.00) + theme(axis.text=element_text(size=12)) + scale_color_gradient("F%", low="yellow", high="red")
p7

# Distributions of all 3 polarity components ratios (PSA, NPSA, other) against gradual F%
# in CHCl3 and H2O
# Load the dataset
Ratio.df <- read.csv("./Polarity Ratios and Charges.csv", header=TRUE)
colnames(Ratio.df) <- c("Compound", "Number", "F", "CHCl3npsa", "CHCl3psa", "CHCl3other", "WATERnpsa", "WATERpsa", "WATERother", "Charge" )
Ratio.df2 <- subset(Ratio.df, is.finite(F))
colnames(Ratio.df2) <- c("Compound", "Number", "F", "CHCl3npsa", "CHCl3psa", "CHCl3other", "WATERnpsa", "WATERpsa", "WATERother", "Charge" )

F <- factor(Ratio.df2$F)
npsa <- Ratio.df2$CHCl3npsa
psa <- Ratio.df2$CHCl3psa
other <- Ratio.df2$CHCl3other
npsa2 <- Ratio.df2$WATERnpsa
psa2 <- Ratio.df2$WATERpsa
other2 <- Ratio.df2$WATERother
number <- Ratio.df2$Number

# Representations as Stacked Barplots
data <- data.frame(F, npsa, psa, other)
mx <- melt(data)
ggplot(mx, aes(x = F)) + geom_bar(aes(weight=value, fill = variable), position = 'fill') + scale_fill_manual(values = rev(brewer.pal(6, "Blues"))) + theme_bw() + theme(axis.text=element_text(size=12))

data2 <- data.frame(F, npsa2, psa2, other2)
mx2 <- melt(data2)
ggplot(mx2, aes(x = F)) + geom_bar(aes(weight=value, fill = variable), position = 'fill') + scale_fill_manual(values = rev(brewer.pal(6, "Reds"))) + theme_bw() + theme(axis.text=element_text(size=12))

# To be continued
## Correlation map between features / variables
## Feature selection
## Predictive modeling ?