#!/usr/bin/Rscript
# -*- mode: R -*-
# Run the KS test and plot the CDF plots for all 7 microarray experiments.
# Notice that this script reads text files with the bnumbers of the targets. The names of the files
# are e.g. ryhB_targets_file_all_conds.txt make sure you have these files for all sRNAs.

library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(grDevices)
library(plyr)
library(beeswarm)

tcol = rgb(204,0,255, maxColorValue = 255)
tt <- read.csv("table_of_bnums_gnames.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")

# Read the list of targets, run KS-Test and plot the figure
# Assume that the name of the column of the names is ORF
ks.plots <- function(seff, outfile, tarsfile, sRNA){
tarlist <- read.csv(tarsfile, header=FALSE, stringsAsFactors=FALSE)
#avmedian = quantile(seff[seff$ORF!='', 'AveExpr'], 0.25)#median(seff[seff$ORF!='','AveExpr'])
test_vals = seff[seff$ORF %in% tarlist$V1, 't']
back_vals = seff[seff$ORF!='' & ! seff$ORF %in% tarlist$V1, 't']
seff$istar <- factor((seff$ORF %in% tarlist$V1)+1)
print(ks.test(test_vals, back_vals))
print(wilcox.test(abs(test_vals), abs(back_vals)))
ksres = ks.test(test_vals, back_vals)
#plot(ecdf(test_vals), main="", xlab="Log fold change", ylab="Cumulative frequency", col=tcol, las=1 )
#plot(ecdf(back_vals), add=TRUE, col="black", las=1)
gp <- ggplot(seff[seff$ORF!='',], aes(logFC, AveExpr, color=istar, alpha=istar))+ scale_color_manual(values=c('black', tcol)) +theme_bw() + theme(legend.position="none") 
cairo_ps(filename=paste0(outfile, "_scatter.eps"), width=3.3, height=3.3, family="Helvetica")
par(ps=10, cex=1, mar=c(3, 3, 1, 1)+0.1, mgp=c(1.8, 0.7, 0))
print(gp + geom_point(size=2) + geom_point(size=2, subset=.(istar==2))+ scale_alpha_manual(values=c(0.2, 0.8))+xlab("Log Fold Change") + ylab("Average Expression"))
dev.off()
cairo_ps(filename=paste0(outfile, "_ecdf.eps"), width=3.3, height=3.3, family="Helvetica")
par(ps=10, cex=1, mar=c(3, 3, 1, 1)+0.1, mgp=c(1.8, 0.7, 0))
gp <- ggplot(seff[seff$istar==1,], aes(t, color=istar))+ scale_color_manual(values=c('black', tcol)) +theme_bw() + theme(legend.position="none") 
print(gp + stat_ecdf(geom="point")+stat_ecdf(geom="smooth")+stat_ecdf(data=seff[seff$istar==2,], geom="point")+stat_ecdf(data=seff[seff$istar==2,], geom="smooth") +xlab('Fold change t statistic')+ylab('Distribution'))
#bees <- beeswarm(logFC ~ istar, data=seff[seff$ORF!='',], spacing=2)
#print(ggplot(bees, aes(x,y,color=x.orig)) + geom_point(shape=19) + scale_color_manual(values=c('black', tcol))+xlab('')+ylab('Log Fold Change'))
dev.off()
}


#No. 1: RyhB (Masse et al 2005, GSE3105)

ryhBgse <- getGEO("GSE3105")
eset <- ryhBgse[[1]]
design <- matrix(c(rep(1,12), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,4), rep(0,4), rep(1,4), rep(0,8), rep(1,4)), ncol=4)
colnames(design) <- c("intercept", "ryhB", "Fur", "Fe-SO4")
exprs(eset) <- log2(exprs(eset)+0.01)
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="ryhB", number=nrow(eset))
write.table(reff, "ryhB_Masse_microarray_ryhB_eff.txt", sep="\t")
print("RyhB KS-test (Masse 2005)")
ks.plots(reff, "ryhB_All_Masse", "ryhB_targets_file_all_conds.txt", "RyhB")



# No. 2 GcvB (Sharma et al 2011, GSE26573) Salmonella
# We will test R1 and R2 separately, maybe one has more incluence
# Read the GSE SOFT file and generate a ExpressionSet object
gcvBgse <- getGEO("GSE26573", GSEMatrix=FALSE)
probesets <- Table(GPLList(gcvBgse)[[1]])$GeneName
data.matrix <- do.call('cbind', lapply(GSMList(gcvBgse), function(x) 
                       {tab <- Table(x)
                        mymatch <- match(probesets, tab$GeneName)
                        return(tab$VALUE[mymatch])}))
data.matrix <- apply(data.matrix, 2, function(x) {as.numeric(as.character(x))})
# If a gene appears more than once, take only the first probe

rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(GSMList(gcvBgse))
data.matrix <- data.matrix[!duplicated(row.names(data.matrix)),]
pdata <- data.frame(samples=names(GSMList(gcvBgse)))
rownames(pdata) <- names(GSMList(gcvBgse))
pheno <- as(pdata, "AnnotatedDataFrame")
eset <- new('ExpressionSet', exprs=log2(data.matrix), phenoData=pheno)
# Generate a design matrix
#design <- matrix(c(rep(1,8),rep(1,6), rep(0,2), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,4), rep(0,4)), ncol=4)
#colnames(design) <- c("intercept", "gcvB", "R1", "R2")
design <- matrix(c(rep(1,8), rep(1,2), rep(0,2), rep(1,2), rep(0,2), rep(1,4), rep(0,4)), ncol=3)
colnames(design) <- c("intercept", "R1", "R2")
# Run lmFit of limma
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
# Read the name to bnum table

# For the sRNA gene. REMOVED
#geff <- topTable(fit, coef="gcvB", number=nrow(eset))
#geff$Gene <- row.names(geff)
#geff.bnum <- merge(geff, tt, by="Gene")
#print("GcvB gene itself")
#ks.plots(geff.bnum, "gcvB_gene_log_Sharma.eps", "gcvB_targets_file_log.txt")

#R1 only
geff <- topTable(fit, coef="R1", number=nrow(eset))
geff$Gene <- row.names(geff)
geff.bnum <- merge(geff, tt, by="Gene")
write.table(geff.bnum, "gcvB_Sharma_gcvB_microarray_R1_effect.txt", sep="\t")
print("GcvB R1 region only")
ks.plots(geff.bnum, "gcvB_R1_All_Sharma", "gcvB_targets_file_all_conds.txt", "GcvB-R1")

#R2 only
geff <- topTable(fit, coef="R2", number=nrow(eset))
geff$Gene <- row.names(geff)
geff.bnum <- merge(geff, tt, by="Gene")
write.table(geff.bnum, "gcvB_Sharma_gcvB_microarray_R2_effect.txt", sep="\t")
print("GcvB R2 region only")
ks.plots(geff.bnum, "gcvB_R2_All_Sharma", "gcvB_targets_file_all_conds.txt", "GcvB-R2")


# No. 3 fnrS (Durand 2010 GSE19655)
fnrSgse <- getGEO("GSE19655")
eset <- fnrSgse[[1]]
design <- matrix(c(rep(1,6), 0, 1, 0, 1, 0, 1), ncol=2)
colnames(design) <- c("intercept", "fnrS")
exprs(eset) <- log2(exprs(eset)+0.01)
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="fnrS", number=nrow(eset))
write.table(reff, "fnrS_Durand_microarray_fnrS_effect.txt", sep="\t")
print("FnrS KS-test (Durand 2010)")
ks.plots(reff, "fnrS_All_Durand", "fnrS_targets_file_all_conds.txt", "FnrS")

# No. 4 mgrR (Moon and Gottesman 2009 GSE18935)
mgrRgse <- getGEO("GSE18935")
eset <- mgrRgse[[1]]
design <- matrix(c(rep(1,6), rep(-1,4), rep(1,2), rep(0,2), rep(1,4) ), ncol=3)
colnames(design) <- c("intercept", "mgrR", "pBR")
exprs(eset) <- log2(exprs(eset)+0.01)
exprs(eset) <- exprs(eset)[, 3:8]
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="mgrR", number=nrow(eset))
write.table(reff, "mgrR_Moon_microarray_mgrR_effect.txt", sep="\t")
print("MgrR KS-test (Moon and Gottesman, 2009)")
ks.plots(reff, "mgrR_All_Moon", "mgrR_targets_file_all_conds.txt", "MgrR")

# No. 4 mgrR (Moon and Gottesman 2009 GSE18935) - plasmid with mgrR/empty plasmid
mgrRgse <- getGEO("GSE18935")
eset <- mgrRgse[[1]]
design <- matrix(c(rep(1,4), rep(-1,2), rep(1,2)), ncol=2)
colnames(design) <- c("intercept", "mgrR")
exprs(eset) <- log2(exprs(eset)+0.01)
exprs(eset) <- exprs(eset)[, 5:8]
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="mgrR", number=nrow(eset))
write.table(reff, "mgrR_Moon_microarray_mgrR_plasmid_induction_effect.txt", sep="\t")
print("MgrR KS-test plasmid effect(Moon and Gottesman, 2009)")
ks.plots(reff, "mgrR_All_Moon_plasmid_induction", "mgrR_targets_file_all_conds.txt", "MgrR")

# No. 4 mgrR (Moon and Gottesman 2009 GSE18935)
mgrRgse <- getGEO("GSE18935")
eset <- mgrRgse[[1]]
design <- matrix(c(rep(1,4), rep(0,2), rep(-1,2)), ncol=2)
colnames(design) <- c("intercept", "mgrR")
exprs(eset) <- log2(exprs(eset)+0.01)
exprs(eset) <- exprs(eset)[, 1:4]
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="mgrR", number=nrow(eset))
write.table(reff, "mgrR_Moon_microarray_mgrR_WT_vs_KO_effect.txt", sep="\t")
print("MgrR KS-test WT/KO(Moon and Gottesman, 2009)")
ks.plots(reff, "mgrR_All_Moon_WT_vs_KO", "mgrR_targets_file_all_conds.txt", "MgrR")


#No 5. ArcZ (Papenfort 2009 GSE17771) Salmonella
# There are 3 experiments, test the KO against WT and OE
arcZgse <- getGEO("GSE17771")
eset <- arcZgse[[1]]
design <- matrix(c(rep(1,6), rep(0,3), rep(1,3)), ncol=2)
colnames(design) <- c("intercept", "arcZ")
#exprs(eset) <- log2(exprs(eset)+0.01)
exprs(eset) <- exprs(eset)[, c(4:6, 13:15)]
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="arcZ", number=nrow(eset))
# Add the bnumbers
r2 <- subset(reff, select=-ORF)
r2$Gene <- r2$SYNONYM
r2 <- merge(r2, tt, by="Gene")
write.table(r2, "arcZ_Papenfort_microarray_arcZ_effect.txt", sep="\t")
print("ArcZ KO vs WT KS-test (Papenfort, 2009)")
ks.plots(r2, "arcZ_All_Papenfort", "arcZ_targets_file_all_conds.txt", "ArcZ")

# ArcZ pulse expression (inducible plasmid pBAD)
eset <- arcZgse[[1]]
design <- matrix(c(rep(1,6), rep(0,3), rep(1,3)), ncol=2)
colnames(design) <- c("intercept", "arcZ")
#exprs(eset) <- log2(exprs(eset)+0.01)
exprs(eset) <- exprs(eset)[, c(7:12)]
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="arcZ", number=nrow(eset))
# Add the bnumbers
r2 <- subset(reff, select=-ORF)
r2$Gene <- r2$SYNONYM
r2 <- merge(r2, tt, by="Gene")
write.table(r2, "arcZ_Papenfort_microarray_arcZ_pBAD_effect.txt", sep="\t")
print("ArcZ pBAD vs control KO KS-test (Papenfort, 2009)")
ks.plots(r2, "arcZ_All_pBAD_Papenfort", "arcZ_targets_file_all_conds.txt", "ArcZ")




# No. 6 cyaR (DeLay and Gottesman 2009
# Read the table and build a Eset
cyaR <- read.csv("DeLay_Gottesman_cyaR_microarray.txt", sep="\t")
cyaRtab <- subset(cyaR, select=c("WT1", "CyaR1", "WT2", "CyaR2"))
cyaRtab <- apply(cyaRtab, 2, function(x) {as.numeric(as.character(x))})
rownames(cyaRtab) <- cyaR$ORF
pdata <- data.frame(samples=c("WT1", "CyaR1", "WT2", "CyaR2"))
rownames(pdata) <- c("WT1", "CyaR1", "WT2", "CyaR2")
pheno <- as(pdata, "AnnotatedDataFrame")
eset <-  new('ExpressionSet', exprs=cyaRtab, phenoData=pheno)
exprs(eset) <- log2(exprs(eset) + 0.1)
design <- matrix(c(rep(1,4), 0, 1, 0, 1), ncol=2)
colnames(design) <- c("intercept", "cyaR")
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="cyaR", number=nrow(eset))
# Remove trailing and leading spaces
reff$ORF <- sapply(rownames(reff), function (x) gsub("^\\s+|\\s+$", "", x))
write.table(reff, "cyaR_DeLay_microarray_cyaR_effect.txt", sep="\t")
print("CyaR KS-Test (DeLay 2009)")
ks.plots(reff, "cyaR_All_DeLay", "cyaR_targets_file_all_conds.txt", "CyaR")

# No. 7 Spf (Biesel and Storz 2011, GSE24875)
spfgse <- getGEO("GSE24875")
eset <- spfgse[[1]]
design <- matrix(c(rep(1,6), 0, 1, 0, 1, 0, 1), ncol=2)
colnames(design) <- c("intercept", "spf")
exprs(eset) <- log2(exprs(eset)+0.01)
fit <- lmFit(eset, design=design)
fit <- eBayes(fit)
reff <- topTable(fit, coef="spf", number=nrow(eset))
write.table(reff, "spf_Biesel_microarray_spf_effect.txt", sep="\t")
print("Spf KS-test (Biesel 2011)")
ks.plots(reff, "spf_All_Biesel", "spf_targets_file_all_conds.txt", "Spf")
