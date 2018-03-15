
#Sample script of the preprocessing steps. For detailed description --> Methods of https://www.biorxiv.org/content/early/2018/02/14/265447


#######################################################
library(GEOquery)
getGEOSuppFiles("GSE32591")

untar("/Users/saezlab/Desktop/KD/GSE32591/GSE32591_RAW.tar", exdir = "dat")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels = list.files("dat/", pattern = "[gz]")
sapply(paste("dat", cels, sep = "/"), gunzip)

#Describing experiment
#This is just a text file which describes the chip names, and the source of the biological samples hybridised to them.
#Open a new terminal window and type:

$ ls dat/*.CEL > dat/phenodata.txt

#Open this file in a text editor. 
#Currently this is a single column listing the files just a list of the CEL files. 
#The final file needs 3 tab-delimited columns named ‘Name’ ‘FileName’ and ‘Target’. 
#The FileName and Name columns will be identical in this case.
#This method originates from : http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

setwd("~/Desktop/KD/GSE32591")

library(simpleaffy)



#import CEL files into R
celfiles_ln = read.affy(covdesc = "pheno.txt", path = "dat")




library(RColorBrewer)
cols <- brewer.pal(8, "Set1")


#Boxplot of raw data
pdf("boxplot_intensity_gse32591CEL_nonorm.pdf")
boxplot(celfiles_ln, col=cols)
dev.off()

#Density histograms depicting the distribution of the raw data
pdf("histogram_intensity_gse32591CEL_nonorm.pdf")
hist(celfiles_ln, col=cols)
dev.off()

#RNA degradation plot

ar = AffyRNAdeg(celfiles_ln, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)

library(affyPLM)

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_ln)


#REL and NUSE plots to assess the homogeneity of probe sets

pdf("REL_intensity_gse32591_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()

pdf("NUSE_intensity_gse32591_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()



# Visualising the relationships between the samples using heirarchical clustering

eset = exprs(celfiles_ln)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()


pdf("heatmap_cor_gse32591_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_ln)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()



#Background correction and log2 transformation without normalization using the RMA package.
celnrma_nonorm = affy::rma(celfiles_ln, normalize = FALSE)
save(celnrma_nonorm, file = "celnrma_nonorm.RData")


#Remove probesets without Entrez Gene identifiers
#Probe with the highest IQR was retained as the representative of that given gene in the dataset
celfiles_ln_filtered_nonorm = nsFilter(celnrma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)




#Annotation for Affymetrix Human Genome U133A Array 
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_ln_filtered_nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



ln_gse32591_nonorm = cbind(Symbols, exprs(celfiles_ln_filtered_nonorm$eset))
ln_gse32591_nonorm = ln_gse32591_nonorm[order(ln_gse32591_nonorm[, 1]), ]

rownames(ln_gse32591_nonorm) = ln_gse32591_nonorm[,1]
ln_gse32591_nonorm = ln_gse32591_nonorm[, -1]

colnames_ln32 = colnames(ln_gse32591_nonorm)
rownames_ln32 = rownames(ln_gse32591_nonorm)

ln_gse32591_nonorm = matrix(as.numeric(ln_gse32591_nonorm), ncol = 93, nrow = 12437)

colnames(ln_gse32591_nonorm) = colnames_ln32
rownames(ln_gse32591_nonorm) = rownames_ln32


save(ln_gse32591_nonorm, file = "ln_gse32591_nonorm.RData")
######################################################################################
