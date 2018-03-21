			  ####Full and ordered script of the integrated CKD microarray data#########


#Data sets with GEO numbers:
#Glomerular data: GSE20602; GSE32591; GSE37460; GSE47183; GSE50469
#Tubular data: GSE32591; GSE35487; GSE35488; GSE37455; GSE47184; GSE69438


#1.Downloading and quality control assessment of the data sets.
#This method originates from : http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
#The following script is constituting the steps involved in the preprocessing analysis conducted
# on the public CKD data. Steps: downloading the raw data; quality control assessment; platform-specific merging of the data
#This preproccessing was done for each data set separately. 
#A simpler version of the script can be found in preprocessing_sample.R
#Experiment: GSE93798
#Platform: GPL570


library(GEOquery)
setwd("~/Desktop/KD/GSE93798")

getGEOSuppFiles("GSE93798")

setwd("~/Desktop/KD/GSE93798/GSE93798")

untar("GSE93798/GSE93798_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

library(simpleaffy)
library(affy)
library(hgu133a.db)
library(hgu133acdf)


celfiles_iga = read.affy(covdesc = "phenodataiga.txt", path = "data")

iga_rma_nonorm = rma(celfiles_iga, normalize = FALSE)
save(iga_rma_nonorm, file = "iga_rma_nonorm.RData")

save(iga_rma, file = "iga_rma.RData")

library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("liu17_iga_rawdata_boxplot.pdf")
boxplot(celfiles_iga, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)
  
boxplot(iga_rma_nonorm, col=cols)

# Plot a density vs log intensity histogram for the unnormalised data
pdf("hist_iga_raw.pdf")
hist(celfiles_iga, col=cols)
dev.off()
# Plot a density vs log intensity histogram for the normalised data
hist(iga_rma_nonorm, col = cols)

#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.

pdf("cor_rawdata.pdf")
heatmap.2(cor(exprs(celfiles_iga)), col = colorRampPalette(c("blue","white","red"))(n=200), trace = "none" )
dev.off()


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_iga)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("RLE_iga17.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_iga17.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_iga)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("clusters_iga17_raw.pdf")
plot(clusters)
dev.off()
#Filtering data

#Now we have looked at the data, we can go on to analyse it. 
#The first stage of analysis is to filter out uninformative data such as control probesets and other internal controls as well as
#removing genes with low variance, that would be unlikely to pass statistical tests for differential expression, 
#or are expressed uniformly close to background detection levels. 
#The modifiers to nsFilter below tell nsFilter not to remove probesets without Entrez Gene identifiers, 
#or have duplicated Entrez Gene identifiers.

celfiles.filtered_iga  = nsFilter(iga_rma_nonorm, require.entrez=TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles.filtered_iga, file = "celfiles.filtered_iga.RData")

celfiles.filtered_iga_nonorm  = nsFilter(iga_rma_nonorm, require.entrez=TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles.filtered_iga_nonorm, file = "celfiles.filtered_iga_nonorm.RData")


https://www.biostars.org/p/53870/ --> eventually used this site to annotate gene symbols


#Annotation
library(hgu133plus2.db)
library(annotate)

probes=row.names(celfiles.filtered_iga_nonorm$eset)
Symbols96 = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_ID96 = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
iga_nonorm = cbind(Entrez_ID96, Symbols96, exprs(celfiles.filtered_iga_nonorm$eset))
iga_nonorm = iga_nonorm[order(iga_nonorm[, 1]), ]

rownames(iga_nonorm) = iga_nonorm[,1]
iga_nonorm = iga_nonorm[, -1]
colnames_iganonorm = colnames(iga_nonorm)
rownames_iganonorm = rownames(iga_nonorm)

iga_nonorm = matrix(as.numeric(iga_nonorm), ncol = 42, nrow = 20514)

colnames(iga_nonorm) = colnames_iganonorm
rownames(iga_nonorm) = rownames_iganonorm

save(iga_nonorm, file = "iga_nonorm.RData")

----------------------------------------------------------------------



gene.sym_iga = getSYMBOL(rownames(celfiles.filtered_iga$eset), "hgu133plus2") 


#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
iga = cbind(gene.sym_iga, exprs(celfiles.filtered_iga$eset))
iga = iga[order(iga[, 1]), ]

rownames(iga) = iga[, 1]
iga = iga[, -1]

#RMA-normalised and annotated gene expression data of GSE69438, log2-transformed values (part of RMA)
save(iga, file = "iga.RData")

genesymboliga = rownames(iga)
iga = matrix(as.numeric(iga), ncol = 42, nrow = 20514)

colnames(iga) = colnames(exprs(celfiles.filtered_iga$eset))

rownames(iga) = genesymboliga

save(iga, file = "iga.RData")

#major batch effect in the data --> not included in the analysis


#2. Tissue Transcriptome Driven Identification of Epidermal Growth Factor as a Chronic Kidney Disease Biomarker
# GSE69438
#Platform: GPL570

setwd("/Users/saezlab/Desktop/KD")
library(GEOquery)
getGEOSuppFiles("GSE69438")

untar("/Users/saezlab/Desktop/KD/GSE69438/GSE69438_RAW.tar", exdir = "raw")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels2 = list.files("raw/", pattern = "[gz]")
sapply(paste("raw", cels2, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE69438")

library(simpleaffy)
celfiles2 = read.affy(covdesc = "phenodata.txt", path = "raw")


#no quantile normalization
celfiles2rma_nonorm = rma(celfiles2, normalize = FALSE)
save(celfiles2rma_nonorm, file = "celfiles2rma_nonorm.RData")

save(celfiles2_rma, file = "celfiles2_rma.RData")


#cyclic loess##only used when analysed separately, if integrated with other data set, then cyclic loess normalize when they are already merged.
gse69438_cloess = normalizeBetweenArrays(exprs(celfiles2rma_nonorm), method = "cyclicloess") #normalize expression between patients
save(gse69438_cloess, file = "loess_normalised_gse6938.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_gse69438CEL_nonorm.pdf")
boxplot(celfiles2, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_gse69438CEL_cyclicloess.pdf")
boxplot(gse69438_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse69438CEL_nonorm.pdf")
hist(celfiles2, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse69438CEL_cyclicloess.pdf")
hist(gse69438_cloess, col = cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles2)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse69438CEL_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse69438CEL_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse69438_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

pdf("gse69438_corALL.RData")
heatmap.2(cor(exprs(celfiles2)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles2_filtered = nsFilter(celfiles2_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
celfiles2_filterednonorm = nsFilter(celfiles2rma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles2_filterednonorm, file = "celfiles2_filterednonorm.RData")



#Annotation
library(hgu133plus2.db)
library(annotate)

#Annotation
library(hgu133plus2.db)
library(annotate)

probes=row.names(celfiles2_filterednonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



gse69438_nonorm = cbind(Symbols, exprs(celfiles2_filterednonorm$eset))
gse69438_nonorm = gse69438_nonorm[order(gse69438_nonorm[, 1]), ]

rownames(gse69438_nonorm) = gse69438_nonorm[,1]
gse69438_nonorm = gse69438_nonorm[, -1]
colnames_gse69438_nonorm = colnames(gse69438_nonorm)
rownames_gse69438_nonorm = rownames(gse69438_nonorm)

gse69438_nonorm = matrix(as.numeric(gse69438_nonorm), ncol = 42, nrow = 20514)

colnames(gse69438_nonorm) = colnames_gse69438_nonorm
rownames(gse69438_nonorm) = rownames_gse69438_nonorm

save(gse69438_nonorm, file = "gse69438_nonorm.RData")

#REMOVE CKD SAMPLES (5)

gse69438_nonorm = gse69438_nonorm[, -25]
gse69438_nonorm = gse69438_nonorm[, -21]
gse69438_nonorm = gse69438_nonorm[, -21]
gse69438_nonorm = gse69438_nonorm[, -21]
gse69438_nonorm = gse69438_nonorm[, -21]
gse69438_nonorm = gse69438_nonorm[, -21]


save(gse69438_nonorm, file = "gse69438_nonorm.RData")



read.csv2(GSE69438_INFO.csv, header = TRUE, sep = "/")

----------------------------------------------------------------------------------

gene.sym = getSYMBOL(rownames(celfiles2_filtered$eset), "hgu133plus2") 


#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
rty = cbind(gene.sym, exprs(celfiles2_filtered$eset))
rty = rty[complete.cases(rty), ]
rty = rty[order(rty[, 1]), ]

rownames(rty) = rty[,1]
rty = rty[, -1]
save(rty, file = "rty.RData")


rty1 = rty[order(rty1$gene.sym), ]
rownames(rty1) = rty1$gene.sym
rty1 = rty1[, -1]

#RMA-normalised and annotated gene expression data of GSE69438, log2-transformed values (part of RMA)
save(rty1, file = "rty1.RData")
#######################################################################################################################################

#3. GSE20602, Human Nephrosclerosis Triggers a Hypoxia-Related Glomerulopathy
#Human Kidney biopsies of patients with NSC(14) and Tumor free kidney specimens
#from patients undergoing tumor nephrectomy, #Glomeruli. Platform GPL96

library(GEOquery)
getGEOSuppFiles("GSE20602")

untar("/Users/saezlab/Desktop/KD/GSE20602/GSE20602_RAW.tar", exdir = "NSC")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_nsc = list.files("NSC/", pattern = "[gz]")
sapply(paste("NSC", cels_nsc, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE20602")

library(simpleaffy)
celfiles_nsc = read.affy(covdesc = "phenodata.txt", path = "NSC")

celfiles_nsc_rma = rma(celfiles_nsc)
save(celfiles_nsc_rma, file = "celfiles_nsc_rma.RData")

#no quantile norm.
celfiles_nscrma_nonorm = rma(celfiles_nsc, normalize = FALSE)

nsc_cloess = normalizeBetweenArrays(nsc_nonorm, method = "cyclicloess") #normalize expression between patients
save(nsc_cloess, file = "loess_normalised_nsc.RData")



library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_nscCEL_nonorm.pdf")
boxplot(celfiles_nsc, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_gse69438CEL_cyclicloess.pdf")
boxplot(nsc_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_nscCEL_nonorm.pdf")
hist(celfiles_nsc, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_nscCEL_cyclicloess.pdf")
hist(nsc_cloess, col = cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_nsc)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_nscCEL_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_nscCEL_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = nsc_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()


heatmap.2(cor(exprs(celfiles_nsc)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )

#Remove NSC13, clusters with healthy

celfiles_nsc = celfiles_nsc[,-10]
save(celfiles_nsc, file = "celfiles_nsc.RData")

pdf("celfiles_nsc_allbyall_woNSC13.RData")
heatmap.2(cor(exprs(celfiles_nsc)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

celfiles_nscrma_nonorm = rma(celfiles_nsc, normalize = FALSE)
save(celfiles_nscrma_nonorm, file = "celfiles_nscRMA_nonorm.RData")



#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_nscfilt_nonorm = nsFilter(celfiles_nscrma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_nscfilt_nonorm, file = "celfiles_nscfilt_nonorm.RData")


#Annotation
library(hgu133a.db)
library(annotate)


probes=row.names((celfiles_nscfilt_nonorm$eset))
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))

nsc_nonorm = cbind(Symbols, exprs(celfiles_nscfilt_nonorm$eset))
nsc_nonorm = nsc_nonorm[order(nsc_nonorm[, 1]), ]

rownames(nsc_nonorm) = nsc_nonorm[,1]
nsc_nonorm = nsc_nonorm[, -1]
colnames_nsc_nonorm = colnames(nsc_nonorm)
rownames_nsc_nonorm = rownames(nsc_nonorm)

nsc_nonorm = matrix(as.numeric(nsc_nonorm), ncol = 17, nrow = 12437)

colnames(nsc_nonorm) = colnames_nsc_nonorm
rownames(nsc_nonorm) = rownames_nsc_nonorm

save(nsc_nonorm, file = "nsc_nonorm.RData")


-------------------------------------------------------------
#Annotation
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)
gene.sym_nsc = getSYMBOL(rownames(celfiles_nsc_filtered$eset), "hgu133a") 


#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
nsc = cbind(gene.sym_nsc, exprs(celfiles_nsc_filtered$eset))
nsc = nsc[order(nsc[, 1]), ]

rownames(nsc) = nsc[,1]
nsc = nsc[, -1]
save(nsc, file = "nsc.RData")

nsc1 = matrix(as.numeric(nsc), ncol = 18, nrow = 12437)

colnames(nsc1) = colnames(nsc)
nsc = nsc1
save(nsc, file = "nsc.RData")

###############################################################

#4. GSE69814, Comparison of Glomerular Transcriptome Profiles of 
#Adult-Onset Steroid Sensitive Focal Segmental Glomerulosclerosis
#and Minimal Change Disease
#Tong et al.(2015). Chinese patients, merely. Glomerular Tissue only. #Gene 1.0 ST Array#GPL6244

library(GEOquery)
getGEOSuppFiles("GSE69814")

untar("/Users/saezlab/Desktop/KD/GSE69814/GSE69814_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_tong = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_tong, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE69814")

#for this microarray, we need to use oligo package.


library(oligo)
celFiles <- list.celfiles('myCELs', full.names=TRUE)

celfiles_tong = list.celfiles('/Users/saezlab/Desktop/KD/GSE69814/data' , full.names = TRUE)
raw_tong = read.celfiles(celfiles_tong)
celfiles_tong_rma = oligo::rma(raw_tong, target ='core')


save(celfiles_tong_rma, file = "celfiles_tong_rma.RData")

source("https://bioconductor.org/biocLite.R")
biocLite("pd.hugene.1.0.st.v1")
library(pd.hugene.1.0.st.v1)

library(hugene10sttranscriptcluster.db)

#Extract probe ids, entrez symbols, and entrez ids


celfiles_tong_filtered = nsFilter(celfiles_tong_rma), require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_tong_filtered, file = "celfiles_tong_filtered.RData")


probes=row.names(celfiles_tong_rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
tong=cbind(probes,Symbols,Entrez_IDs, exprs(celfiles_tong_rma))



#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

tong = tong[complete.cases(tong), ]
tong = tong[order(tong[, 2]), ]

save(tong, file = "tong.RData")  #not included  due oligo platform

###################################################################################

#5.GSE50469, The molecular phenotype of endocapillary proliferation 
#in IgA nephropathy and potential modulation 
#by bioactive small molecules, Hodgin et al. (2014)
#ERCB cohort: 16:platform: GPL96
#Toronto: 6, platform: GPL570
#Separate the Two when analyse
#Only Glomerulus


library(GEOquery)
getGEOSuppFiles("GSE50469")

untar("/Users/saezlab/Desktop/KD/GSE50469/GSE50469_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels3 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels3, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE50469")


#5/1.Toronto IgA Data
library(simpleaffy)
celfiles_toronto = read.affy(covdesc = "phenodata_toronto_gpl570.txt", path = "data")

celfiles_toronto_nonorm = affy::rma(celfiles_toronto, normalize = FALSE)

save(celfiles_toronto_rma, file = "celfiles_toronto_rma.RData")

#no quant.norm
celfiles_toronto_nonorm = affy::rma(celfiles_toronto, normalize = FALSE)
save(celfiles_toronto_nonorm, file = "celfiles_toronto_nonorm.RData")

#cyclic loess
tor_cloess = normalizeBetweenArrays(iga_toronto_nonorm, method = "cyclicloess") #normalize expression between patients
save(tor_cloess, file = "loess_normalised_toronto.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_torontoCEL_nonorm.pdf")
boxplot(celfiles_toronto, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_gse50469CEL_cyclicloess.pdf")
boxplot(tor_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse50469CEL_nonorm.pdf")
hist(celfiles_toronto, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse50469CEL_cyclicloess.pdf")
hist(tor_cloess, col = cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_toronto)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse50469CEL_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse50469CEL_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = tor_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

pdf("celfiles_gse50469_toronto.RData")
heatmap.2(cor(exprs(celfiles_toronto)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_toronto_filtered = nsFilter(celfiles_toronto_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_toronto_filtered, file = "celfiles_toronto_filtered.RData")

celfiles_torontofilt_nonorm = nsFilter(celfiles_toronto_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_torontofilt_nonorm, file = "celfiles_torontofilt_nonorm.RData")

#Annotation
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_torontofilt_nonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
iga_toronto_nonorm = cbind(Symbols, exprs(celfiles_torontofilt_nonorm$eset))
iga_toronto_nonorm = iga_toronto_nonorm[order(iga_toronto_nonorm[, 1]), ]

rownames(iga_toronto_nonorm) = iga_toronto_nonorm[,1]
iga_toronto_nonorm = iga_toronto_nonorm[, -1]
colnames_igatoronto = colnames(iga_toronto_nonorm)
rownames_igatoronto = rownames(iga_toronto_nonorm)

iga_toronto_nonorm = matrix(as.numeric(iga_toronto_nonorm), ncol = 6, nrow = 20514)

colnames(iga_toronto_nonorm) = colnames_igatoronto
rownames(iga_toronto_nonorm) = rownames_igatoronto

save(iga_toronto_nonorm, file = "iga_toronto_nonorm.RData")



#5/2. ERCB IgA Data
library(simpleaffy)
celfiles_ercb = read.affy(covdesc = "phenodata_ercb_gpl96.txt", path = "data")

celfiles_ercb_rma = rma(celfiles_ercb)

save(celfiles_toronto_rma, file = "celfiles_toronto_rma.RData")
#no quant.norm
celfiles_ercbrma_nonorm = affy::rma(celfiles_ercb, normalize = FALSE)
save(celfiles_ercbrma_nonorm, file = "celfiles_ercbrma_nonorm.RData")


#cyclic loess
ercb_cloess = normalizeBetweenArrays(iga_ercb_nonorm, method = "cyclicloess") #normalize expression between patients
save(ercb_cloess, file = "loess_normalised_ercb.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_ercbCEL_nonorm.pdf")
boxplot(celfiles_ercb, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_ercbCEL_cyclicloess.pdf")
boxplot(ercb_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_ercbCEL_nonorm.pdf")
hist(celfiles_ercb, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse50469CEL_cyclicloess.pdf")
hist(ercb_cloess, col = cols)
dev.off()

#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_ercb)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse50469ERCBCEL_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse50469ERCBCEL_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = ercb_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

pdf("celfiles_gse50469_ercb.RData")
heatmap.2(cor(exprs(celfiles_ercb)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()



#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_ercb_filtered = nsFilter(celfiles_ercb_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_toronto_filtered, file = "celfiles_toronto_filtered.RData")

celfiles_ercb_filterednonorm = nsFilter(celfiles_ercbrma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_ercb_filterednonorm, file = "celfiles_ercb_filterednonorm.RData")

#Annotation
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)

probes=row.names(celfiles_ercb_filterednonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))

iga_ercb_nonorm = cbind(Symbols, exprs(celfiles_ercb_filterednonorm$eset))
iga_ercb_nonorm = iga_ercb_nonorm[order(iga_ercb_nonorm[, 1]), ]

rownames(iga_ercb_nonorm) = iga_ercb_nonorm[,1]
iga_ercb_nonorm = iga_ercb_nonorm[, -1]

colnames_ercb = colnames(iga_ercb_nonorm)
rownames_ercb = rownames(iga_ercb_nonorm)

iga_ercb_nonorm = matrix(as.numeric(iga_ercb_nonorm), ncol = 16, nrow = 12437)

colnames(iga_ercb_nonorm) = colnames_ercb
rownames(iga_ercb_nonorm) = rownames_ercb

save(iga_ercb_nonorm, file = "iga_ercb_nonorm.RData")


##################################################################################################

#6. GSE32591, Berhier et al. (2012)
#Glomerulus and Tubular
#Lupus Nephritis and Living Donor
#Platform: GPL96, Affymetrix GeneChip Human Genome HG-U133A Custom CDF 

library(GEOquery)
getGEOSuppFiles("GSE32591")

untar("/Users/saezlab/Desktop/KD/GSE32591/GSE32591_RAW.tar", exdir = "dat")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels = list.files("dat/", pattern = "[gz]")
sapply(paste("dat", cels, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

setwd("~/Desktop/KD/GSE32591")

library(simpleaffy)
celfiles_ln = read.affy(covdesc = "pheno.txt", path = "dat")

celfiles_ln_rma = rma(celfiles_ln)

save(celfiles_ln_rma, file = "celfiles_ln_rma.RData")

celnrma_nonorm = affy::rma(celfiles_ln, normalize = FALSE)
save(celnrma_nonorm, file = "celnrma_nonorm.RData")



#cyclic loess
gse32591_cloess = normalizeBetweenArrays(ln_gse32591_nonorm, method = "cyclicloess") #normalize expression between patients
save(gse32591_cloess, file = "loess_normalised_gse32591.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_gse32591CEL_nonorm.pdf")
boxplot(celfiles_ln, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_gse32591CEL_cyclicloess.pdf")
boxplot(gse32591_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse32591CEL_nonorm.pdf")
hist(celfiles_ln, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse32591_cyclicloess.pdf")
hist(gse32591_cloess, col = cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_ln, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_ln)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse32591_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse32591_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse32591_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse32591_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_ln)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#hierarchical clustering using the cyclic loess normalized data
pdf("heatmap_cor_gse32591_cloess.pdf")
heatmap.2(cor(gse32591_cloess), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()



#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_ln_filtered = nsFilter(celfiles_ln_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_ln_filtered, file = "celfiles_ln_filtered.RData")

celfiles_ln_filtered_nonorm = nsFilter(celnrma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_ln_filtered_nonorm, file = "celfiles_ln_filtered_nonorm.RData")


#Annotation: hgu133a 
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_ln_filtered_nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
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


##################################################################################################################

#7. GSE35488, Reich et al. (2010). A molecular signature of proteinuria in glomerulonephritis. PLoS One
#Only Tubulointerstitial compartment of kidney bipsies of patients with IgaN (n=25) and Control (n =6)
#Platform: GPL96



library(GEOquery)
getGEOSuppFiles("GSE35488")

untar("/Users/saezlab/Desktop/KD/GSE35488/GSE35488_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_reich = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_reich, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

setwd("~/Desktop/KD/GSE35488")

library(simpleaffy)
celfiles_reich = read.affy(covdesc = "pheno.txt", path = "data 12.18.30")

celfiles_reich_rma = rma(celfiles_reich)

save(celfiles_reich_rma, file = "celfiles_reich_rma.RData")

#no quant
celfiles_reich_nonorm = affy::rma(celfiles_reich, normalize = FALSE)
save(celfiles_reich_nonorm, file = "celfiles_reich_nonorm.RData")


#cyclic loess
gse35488_cloess = normalizeBetweenArrays(iga_gse35488_nonorm, method = "cyclicloess") #normalize expression between patients
save(gse35488_cloess, file = "loess_normalised_gse35488.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_gse35488_nonorm.pdf")
boxplot(celfiles_reich, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)

#cyclic loess normalized
pdf("boxplot_intensity_gse35488_cyclicloess.pdf")
boxplot(gse35488_cloess, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse35488_nonorm.pdf")
hist(celfiles_reich, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse35488_cyclicloess.pdf")
hist(gse35488_cloess, col = cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_reich, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_reich)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse35488_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse35488_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse35488_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_reich)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#hierarchical clustering using the cyclic loess normalized data
pdf("heatmap_cor_gse35488_cloess.pdf")
heatmap.2(cor(gse35488_cloess), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()




#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_reich_filtered = nsFilter(celfiles_reich_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_reich_filtered, file = "celfiles_reich_filtered.RData")

celfiles_reich_filtered_nonorm = nsFilter(celfiles_reich_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_reich_filtered_nonorm, file = "celfiles_reich_filtered_nonorm.RData")

#Annotation: hgu133a 
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_reich_filtered_nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
iga_gse35488_nonorm = cbind(Symbols, exprs(celfiles_reich_filtered_nonorm$eset))
iga_gse35488_nonorm = iga_gse35488_nonorm[order(iga_gse35488_nonorm[, 1]), ]

rownames(iga_gse35488_nonorm) = iga_gse35488_nonorm[,1]
iga_gse35488_nonorm = iga_gse35488_nonorm[, -1]

col_gse35488 = colnames(iga_gse35488_nonorm)
row_gse35488 = rownames(iga_gse35488_nonorm)

iga_gse35488_nonorm = matrix(as.numeric(iga_gse35488_nonorm), ncol = 31, nrow = 12437)

colnames(iga_gse35488_nonorm) = col_gse35488
rownames(iga_gse35488_nonorm) = row_gse35488

save(iga_gse35488_nonorm, file = "iga_gse35488_nonorm.RData")


#7/2. GSE35487, Reich et al. (2010). A molecular signature of proteinuria in glomerulonephritis. PLoS One
#Only Tubulointerstitial compartment of kidney bipsies of patients with IgaN (n=25) and Control (n =6)
#Platform: GPL96



library(GEOquery)
getGEOSuppFiles("GSE35487")

untar("/Users/saezlab/Desktop/KD/GSE35487/GSE35487_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_reich87 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_reich87, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

setwd("~/Desktop/KD/GSE35487")

library(simpleaffy)
celfiles_reich87 = read.affy(covdesc = "phenodata.txt", path = "data")

save(celfiles_reich87, file = "celfiles_reich.RData")

#no quant
celfiles_reich_nonorm = affy::rma(celfiles_reich, normalize = FALSE)
save(celfiles_reich_nonorm, file = "celfiles_reich_nonorm.RData")


#cyclic loess
gse35488_cloess = normalizeBetweenArrays(iga_gse35488_nonorm, method = "cyclicloess") #normalize expression between patients
save(gse35488_cloess, file = "loess_normalised_gse35488.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values

pdf("boxplot_intensity_gse35487_nonorm.pdf")
boxplot(celfiles_reich87, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)



# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse35487_nonorm.pdf")
hist(celfiles_reich87, col=cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_reich87, log.it = T)
pdf("rnadeg_gse35487.RData")
plotAffyRNAdeg(ar)
dev.off()
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_reich87)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse35487_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse35487_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse35487_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_reich87)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#no quant
celfiles_reich87_nonorm = affy::rma(celfiles_reich87, normalize = FALSE)
save(celfiles_reich87_nonorm, file = "celfiles_reich87_nonorm.RData")

#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_reich87_filtered = nsFilter(celfiles_reich87_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_reich87_filtered, file = "celfiles_reich87_filtered.RData")


#Annotation: hgu133a 
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_reich87_filtered$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
iga_gse35487_nonorm = cbind(Symbols, exprs(celfiles_reich87_filtered$eset))
iga_gse35487_nonorm = iga_gse35487_nonorm[order(iga_gse35487_nonorm[, 1]), ]

rownames(iga_gse35487_nonorm) = iga_gse35487_nonorm[,1]
iga_gse35487_nonorm = iga_gse35487_nonorm[, -1]

col_gse35487 = colnames(iga_gse35487_nonorm)
row_gse35487 = rownames(iga_gse35487_nonorm)

iga_gse35487_nonorm = matrix(as.numeric(iga_gse35487_nonorm), ncol = 31, nrow = 12437)

colnames(iga_gse35487_nonorm) = col_gse35487
rownames(iga_gse35487_nonorm) = row_gse35487

save(iga_gse35487_nonorm, file = "iga_gse35487_nonorm.RData")


###################################################################################################

#8. GSE37455 is part of the super series GSE37463, of which looked bad after limmma DEG.
# Cross-species transcriptional network analysis defines shared inflammatory responses 
#in murine and human lupus nephritis. J Immunol 2012 by Berthier et al.
#Two Platforms: GPL570 and GPL96


library(GEOquery)
getGEOSuppFiles("GSE37455")

untar("/Users/saezlab/Desktop/KD/GSE37455/GSE37455_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_gse37455 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_gse37455, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

# Separate Analysis of samples corresponding to different platforms
#8/1. GPL570


setwd("~/Desktop/KD/GSE37455")
setwd("~/Desktop/KD/GSE37455/data")


library(simpleaffy)
celfiles_570 = read.affy(covdesc = "pheno570.txt", path = "GPL11670_570")

celfiles_570_rma = rma(celfiles_570)

save(celfiles_570_rma, file = "celfiles_570_rma.RData")


#no quant
celfiles570_nonorm = affy::rma(celfiles_570, normalize = FALSE)
save(celfiles570_nonorm, file = "celfiles570_nonorm.RData")


#cyclic loess
gse37455__cloess = normalizeBetweenArrays(gse37455_570_nonorm, method = "cyclicloess") #normalize expression between patients
save(gse37455_570_cloess, file = "loess_normalised_gse37455.RData")



library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse37455_nonorm.pdf")
boxplot(celfiles570_nonorm, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


#cyclic loess normalized
pdf("boxplot_intensity_gse37455_cyclicloess.pdf")
boxplot(gse37455_570_cloess , col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse37455_nonorm.pdf")
hist(celfiles570_nonorm, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse37455_cyclicloess.pdf")
hist(gse37455_570_cloess, col = cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_570, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37455_570.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37455_570.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()



#Remove 1st, 11th and 14th
"GSM920052_H7-Tub-LD1164.CEL"
"GSM920062_H7-Tub-LD1175.CEL" 
"GSM920065_H7-Tub-LD1179.CEL"

celfiles_570 = celfiles_570[, -1]
celfiles_570 = celfiles_570[, -10]
celfiles_570 = celfiles_570[, -12]

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37455_570_without3.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37455_570_without3.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

save(celfiles_570, file = "celfiles_570.RData")
#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37455_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_570)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#no quant
celfiles570_nonorm = affy::rma(celfiles_570, normalize = FALSE)
save(celfiles570_nonorm, file = "celfiles570_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_570_filtered_nonorm = nsFilter(celfiles570_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_570_filtered_nonorm, file = "celfiles_570_filtered_nonorm.RData")


#Annotation: hgu133plus2
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_570_filtered_nonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse37455_570_nonorm = cbind(Symbols, exprs(celfiles_570_filtered_nonorm$eset))
gse37455_570_nonorm = gse37455_570_nonorm[order(gse37455_570_nonorm[, 1]), ]

rownames(gse37455_570_nonorm) = gse37455_570_nonorm[,1]
gse37455_570_nonorm = gse37455_570_nonorm[, -1]

col_gse37455 = colnames(gse37455_570_nonorm)
row_gse37455 = rownames(gse37455_570_nonorm)

gse37455_570_nonorm = matrix(as.numeric(gse37455_570_nonorm), ncol = 15, nrow = 20514)

colnames(gse37455_570_nonorm) = col_gse37455
rownames(gse37455_570_nonorm) = row_gse37455
  
save(gse37455_570_nonorm, file = "gse37455_570_nonorm.RData")   #first patient has been removed




#8/2.  Separate Analysis of samples corresponding to different platforms
#GPL96


setwd("~/Desktop/KD/GSE37455")
setwd("~/Desktop/KD/GSE37455/data")


library(simpleaffy)
celfiles_96 = read.affy(covdesc = "pheno96.txt", path = "GPL14663_96 ")


setwd("~/Desktop/KD/GSE37455/data/GPL14663_96 ")
save(celfiles_96_rma, file = "celfiles_96_rma.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse37455_96_nonorm.pdf")
boxplot(celfiles_96, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


#cyclic loess normalized
pdf("boxplot_intensity_gse37455_96_cyclicloess.pdf")
boxplot(gse37455_96_cloess , col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse37455_96_nonorm.pdf")
hist(celfiles_96, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse37455_96_cyclicloess.pdf")
hist(gse37455_96_cloess, col = cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_96, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37455_96_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37455_96_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()
#first three has to be gone. remove the first three observations

celfiles_96 = celfiles_96[, -1]
celfiles_96 = celfiles_96[, -1]
celfiles_96 = celfiles_96[, -1]


#I removed --> plotAffyRNAdeg(ar["GSM920052_H7-Tub-LD1164.CEL"])
GSM920034_H1-Tub-LD1.CEL
GSM920035_H1-Tub-LD2.CEL
GSM920036_H1-Tub-LD3.CEL

#RNA degradation

ar = AffyRNAdeg(celfiles_96, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37455_96_nonorm_without_first3patients.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37455_96_nonorm_without_first3patients.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()



#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37455_96_nonorm_without_first3patients.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#hierarchical clustering using the cyclic loess normalized data
pdf("heatmap_cor_gse37455_cloess.pdf")
heatmap.2(cor(gse37455_570_cloess), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#Normalization

#1. No normalization
celfiles_96_rma_nonorm = affy::rma(celfiles_96, normalize = FALSE)
save(celfiles_96_rma_nonorm, file = "celfiles_96_rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_96_filtered_nonorm = nsFilter(celfiles_96_rma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_96_filtered_nonorm, file = "celfiles_96_filtered_nonorm.RData")

#Annotation: hgu133a
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)

probes=row.names(celfiles_96_filtered_nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))

#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse37455_96_nonorm = cbind(Symbols, exprs(celfiles_96_filtered_nonorm$eset))
gse37455_96_nonorm = gse37455_96_nonorm[order(gse37455_96_nonorm[, 1]), ]

rownames(gse37455_96_nonorm) = gse37455_96_nonorm[,1]
gse37455_96_nonorm = gse37455_96_nonorm[, -1]

col_gse37455_96 = colnames(gse37455_96_nonorm)
row_gse37455_96 = rownames(gse37455_96_nonorm)


gse37455_96_nonorm = matrix(as.numeric(gse37455_96_nonorm), ncol = 20, nrow = 12437)

colnames(gse37455_96_nonorm) = col_gse37455_96
rownames(gse37455_96_nonorm) = row_gse37455_96

save(gse37455_96_nonorm, file = "gse37455_96_nonorm.RData")


#9. GSE37460 is part of the super series GSE37463, of which looked bad after limmma DEG.
# Cross-species transcriptional network analysis defines shared inflammatory responses 
#in murine and human lupus nephritis. J Immunol 2012 by Berthier et al.
#Two Platforms: GPL570 and GPL96
#Glom Tissue only



library(GEOquery)
getGEOSuppFiles("GSE37460")

untar("/Users/saezlab/Desktop/KD/GSE37460/GSE37460_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_gse37460 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_gse37460, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#9/1  Separate Analysis of samples corresponding to different platforms
#GPL570


setwd("~/Desktop/KD/GSE37460")
setwd("~/Desktop/KD/GSE37460/data")


library(simpleaffy)
celfiles_570 = read.affy(covdesc = "pheno570.txt", path = "GPL11670_570")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse37460_570_nonorm.pdf")
boxplot(celfiles_570, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


#cyclic loess normalized
pdf("boxplot_intensity_gse37455_96_cyclicloess.pdf")
boxplot(gse37455_96_cloess , col=cols)
dev.off()

# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse37460_570_nonorm.pdf")
hist(celfiles_570, col=cols)
dev.off()

# Plot a density vs log intensity histogram for the normalised data
pdf("histogram_intensity_gse37455_96_cyclicloess.pdf")
hist(gse37455_96_cloess, col = cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_570, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37460_570_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37469_570_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()


#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_570_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_570)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#Normalization

#1. No normalization
celfiles_570_rma_nonorm = affy::rma(celfiles_570, normalize = FALSE)
save(celfiles_570_rma_nonorm, file = "celfiles_570_rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_570_filtered = nsFilter(celfiles_570_rma, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_570_filtered, file = "celfiles_570_gse37460_filtered.RData")

#no quant norm data filtering
celfiles_570_filtered_37460nonorm = nsFilter(celfiles_570_rma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_570_filtered_37460nonorm, file = "celfiles_570_gse37460_filtered_37460nonorm.RData")

#Annotation: hgu133plus2
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_570_filtered_37460nonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse37460_570_nonorm = cbind(Symbols, exprs(celfiles_570_filtered_37460nonorm$eset))
gse37460_570_nonorm = gse37460_570_nonorm[order(gse37460_570_nonorm[, 1]), ]

rownames(gse37460_570_nonorm) = gse37460_570_nonorm[,1]
gse37460_570_nonorm = gse37460_570_nonorm[, -1]

col_gse37460_570 = colnames(gse37460_570_nonorm)
row_gse37460_570 = rownames(gse37460_570_nonorm)

gse37460_570_nonorm = matrix(as.numeric(gse37460_570_nonorm), ncol = 18, nrow = 20514)

colnames(gse37460_570_nonorm) = col_gse37460_570
rownames(gse37460_570_nonorm) = row_gse37460_570

save(gse37460_570_nonorm, file = "gse37460_570_nonorm.RData")
-------------------------------------------------
#****
Symbols = sort(Symbols)
rownames(gse37460_570) = Symbols
--------------------------------------------

#9/2.  Separate Analysis of samples corresponding to different platforms
#GPL96


setwd("~/Desktop/KD/GSE37460")
setwd("~/Desktop/KD/GSE37460/data")


library(simpleaffy)
celfiles_96 = read.affy(covdesc = "pheno96.txt", path = "GPL14663_96")



setwd("~/Desktop/KD/GSE37460/data/GPL14663_96")
save(celfiles_96_rma, file = "celfiles_96_gse37460_rma.RData")

celfiles_96_rma37460_nonorm = affy::rma(celfiles_96, normalize = FALSE)
save(celfiles_96_rma37460_nonorm, file = "celfiles_96_rma37460_nonorm.RData")




library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse37460_96_nonorm.pdf")
boxplot(celfiles_96, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse37460_96_nonorm.pdf")
hist(celfiles_96, col=cols)
dev.off()

#RNA degradation

ar = AffyRNAdeg(celfiles_96, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37460_96_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()


#We can also look at the relationships between the samples using heirarchical clustering:

eset = gse35488_cloess
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_cloess.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#Removal of the first 4 patients (Healthy Living Donor)
celfiles_96 = celfiles_96[,-1]
celfiles_96 = celfiles_96[,-1]
celfiles_96 = celfiles_96[,-1]
celfiles_96 = celfiles_96[,-1]

#RNA degradation

ar = AffyRNAdeg(celfiles_96, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse37460_96_nonorm_withoutfirst4.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_withoutfirst4_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()


#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_96)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_nonorm_without4.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


#remove IgA24 and HT13,  they cluster separately than the rest
celfiles_96 = celfiles_96[, -4]
celfiles_96 = celfiles_96[, -29]
save(celfiles_96, file = "celfiles96.RData")


#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_without6_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

pdf("REL_intensity_gse37460_96_nonorm_without6.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_without6.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()


#In total I removed 6 patients

celfiles_96 = celfiles_96[, -40] #IgA9
celfiles_96 = celfiles_96[, -29] #IgA24
celfiles_96 = celfiles_96[, -23] #IgA18
celfiles_96 = celfiles_96[, -21] #IgA16

#total--> removed 10 patients

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_without10_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#AlL-by-All Correlation looks much better.

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

pdf("REL_intensity_gse37460_96_nonorm_without10pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_without10.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

celfiles_96 = celfiles_96[, -18]
celfiles_96 = celfiles_96[, -22]

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

pdf("REL_intensity_gse37460_96_nonorm_without12pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_without12.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_without13_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


celfiles_96 = celfiles_96[, -37]

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

pdf("REL_intensity_gse37460_96_nonorm_without13pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_without13.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_without13_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

celfiles_96 = celfiles_96[, -36]

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96)

pdf("REL_intensity_gse37460_96_nonorm_without14pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse37460_96_without14.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_96_without14_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


save(celfiles_96, file = "celfiles_96.RData")


#Removed samples:
GSM920350_H1-Glom-LD1.CEL
GSM920351_H1-Glom-LD2.CEL	
GSM920352_H1-Glom-LD3.CEL	
GSM920353_H1-Glom-LD4.CEL	

GSM920416_H4-Glom-IgA24.CEL	
GSM920389_H2-Glom-HT13.CEL	

GSM920416_H4-Glom-IgA23.CEL	
GSM920427_H4-Glom-IgA9.CEL	
GSM920409_H4-Glom-IgA18.CEL	
GSM920407_H4-Glom-IgA16.CEL	
GSM920411_H4-Glom-IgA1.CEL
GSM920404_H4-Glom-IgA13.CEL

GSM1046924_H4-Glom-LD7.CEL
GSM1046923_H4-Glom-LD5.CEL #clusters with IgAN


celfiles_96_rma37460_nonorm = affy::rma(celfiles_96, normalize = FALSE)
save(celfiles_96_rma37460_nonorm, file = "celfiles_96_rma37460_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_96_filtered_37460_nonorm = nsFilter(celfiles_96_rma37460_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_96_filtered_37460_nonorm, file = "celfiles_96_filtered_37460_nonorm.RData")
#Annotation: hgu133a
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_96_filtered_37460_nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse37460_96_nonorm = cbind(Symbols, exprs(celfiles_96_filtered_37460_nonorm$eset))
gse37460_96_nonorm = gse37460_96_nonorm[order(gse37460_96_nonorm[, 1]), ]

rownames(gse37460_96_nonorm) = gse37460_96_nonorm[,1]
gse37460_96_nonorm = gse37460_96_nonorm[, -1]

col_37460_96 = colnames(gse37460_96_nonorm)
row_37460_96 = rownames(gse37460_96_nonorm)


gse37460_96_nonorm = matrix(as.numeric(gse37460_96_nonorm), ncol = 37, nrow = 12437)

colnames(gse37460_96_nonorm) = col_37460_96 
rownames(gse37460_96_nonorm) = row_37460_96 


save(gse37460_96_nonorm, file = "gse37460_96_nonorm.RData") #updated

#######################################################################################################

#10. Defining cell-types specificity at the transcriptional level in human disease
#GSE47185 superseries constituting: GSE47183 & GSE47184
#2 different platforms
#1.  GPL570 but the data were analyzed with a custom CDF environment (Hs133P_Hs_ENTREZG.cdf)
#2.  GPL96 (in the study they used https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14663)

#I am going to analyse the two tissue types separately, plus the platform types, as well. 
#SO -Tissue1 --> GPL570 & GPL96
# Tissue2 --> GPL570 $ GPL96

#GSE47183/GLOMERULI, #122 samples in two platforms
library(GEOquery)
getGEOSuppFiles("GSE47183")

untar("/Users/saezlab/Desktop/KD/GSE47183/GSE47183_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_gse47183 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_gse47183, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/




#10/1.  Separate Analysis of samples corresponding to different platforms
#GPL570


setwd("~/Desktop/KD/GSE47183")
setwd("~/Desktop/KD/GSE47183/data")


library(simpleaffy)
celfiles_570_183 = read.affy(covdesc = "phenodata_gpl570.txt", path = "GPL570")



celfiles_570_rma183_nonorm = affy::rma(celfiles_570_183, normalize = FALSE)

save(celfiles_570_rma183_nonorm, file = "celfiles_570_rma183_nonorm.RData")



library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse47183_570_nonorm.pdf")
boxplot(celfiles_570_183, col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse47183_570_nonorm.pdf")
hist(celfiles_570_183, col=cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_570_183, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570_183)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47183_570_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47183_570_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()


#remove the 54th individual   NUSE > 1.10

celfiles_570_183 = celfiles_570_183[, -44]

#Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570_183)

#RNA degradation
ar = AffyRNAdeg(celfiles_570_183, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)

pdf("REL_intensity_gse47183_570_nonorm_wo44.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47183_570_nonorm_wo44.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

Removed sample:
GSM1146312_H7_Glom-RPGN935

save(celfiles_570_183, file = "celfiles_570_183.RData")

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_570_183)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47183_570.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse37460_570_nonorm_wo44.pdf")
heatmap.2(cor(exprs(celfiles_570_183)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#Normalization

#1. No normalization
celfiles_570_rma_nonorm = affy::rma(celfiles_570_183, normalize = FALSE)
save(celfiles_570_rma_nonorm, file = "celfiles_570_rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_570_filtered_183nonorm = nsFilter(celfiles_570_rma_nonorm , require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_570_filtered_183nonorm, file = "celfiles_570_filtered_183nonorm.RData")



#Annotation: hgu133plus2
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_570_filtered_183nonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse47183_570_nonorm = cbind(Symbols, exprs(celfiles_570_filtered_183nonorm$eset))
gse47183_570_nonorm = gse47183_570_nonorm[order(gse47183_570_nonorm[, 1]), ]

rownames(gse47183_570_nonorm) = gse47183_570_nonorm[,1]
gse47183_570_nonorm = gse47183_570_nonorm[, -1]

col_gse47183_570_nonorm = colnames(gse47183_570_nonorm)
row_gse47183_570_nonorm = rownames(gse47183_570_nonorm)

gse47183_570_nonorm = matrix(as.numeric(gse47183_570_nonorm), ncol = 53, nrow = 20514)

colnames(gse47183_570_nonorm) = col_gse47183_570_nonorm
rownames(gse47183_570_nonorm) = row_gse47183_570_nonorm

save(gse47183_570_nonorm, file = "gse47183_570_nonorm.RData")




#10/2.  Separate Analysis of samples corresponding to different platforms
#GPL96


setwd("~/Desktop/KD/GSE47183")
setwd("~/Desktop/KD/GSE47183/data")


library(simpleaffy)
celfiles_96_183 = read.affy(covdesc = "phenodata_gpl96.txt", path = "GPL96")


setwd("~/Desktop/KD/GSE47183/data/GPL96")
save(celfiles_96_rma, file = "celfiles_96_gse47183_rma.RData")




library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse47183_96_nonorm.pdf")
boxplot(celfiles_96_183 , col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse47183_96_nonorm.pdf")
hist(celfiles_96_183, col=cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_96_183, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96_183)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47183_96_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47183_96_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_96_183)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47183_96.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47183_96_nonorm_wo44.pdf")
heatmap.2(cor(exprs(celfiles_96_183)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#Removal

celfiles_96_183 = celfiles_96_183[, -4]
celfiles_96_183 = celfiles_96_183[, -11]
celfiles_96_183 = celfiles_96_183[, -46]


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96_183)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47183_96_nonorm_WO5.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47183_96_nonorm_WO5.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_96_183)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47183_96.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47183_96_nonorm_wo4.pdf")
heatmap.2(cor(exprs(celfiles_96_183)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#remove the remaining TMD samples

celfiles_96_183 = celfiles_96_183[, -11]
celfiles_96_183 = celfiles_96_183[, -11]

# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96_183)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47183_96_nonorm_WO5.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47183_96_nonorm_WO5.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_96_183)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47183_96.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47183_96_nonorm_wo5.pdf")
heatmap.2(cor(exprs(celfiles_96_183)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()


save(celfiles_96_183, file = "celfiles_96_183.RData")


#Normalization

#1. No normalization
celfiles_96_183_rma_nonorm = affy::rma(celfiles_96_183, normalize = FALSE)
save(celfiles_96_183_rma_nonorm, file = "celfiles_96_183_rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_96_filtered_183nonorm = nsFilter(celfiles_96_183_rma_nonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_96_filtered_183nonorm, file = "celfiles_96_filtered_183nonorm.RData")


#Annotation: hgu133a
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_96_filtered_183nonorm$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse47183_96_nonorm = cbind(Symbols, exprs(celfiles_96_filtered_183nonorm$eset))
gse47183_96_nonorm = gse47183_96_nonorm[order(gse47183_96_nonorm[, 1]), ]

rownames(gse47183_96_nonorm) = gse47183_96_nonorm[,1]
gse47183_96_nonorm = gse47183_96_nonorm[, -1]

col_gse47183_96_nonorm = colnames(gse47183_96_nonorm)
row_gse47183_96_nonorm = rownames(gse47183_96_nonorm)

gse47183_96_nonorm = matrix(as.numeric(gse47183_96_nonorm), ncol = 63, nrow = 12437)

colnames(gse47183_96_nonorm) = col_gse47183_96_nonorm
rownames(gse47183_96_nonorm) = row_gse47183_96_nonorm


save(gse47183_96_nonorm, file = "gse47183_96_nonorm.RData")

#Omitted Samples:

GSM1146204_H1-Glom-DN5.CEL
GSM1146212_H1-Glom-TMD.CEL
GSM1146213_H1-Glom-TMD1.CEL
GSM1146214_H1-Glom-TMD5.CEL
GSM1146248_H5-Glom-MGN403.CEL







#10/3,4GSE47184/Tubular, #107 samples in two platforms
library(GEOquery)
getGEOSuppFiles("GSE47184")

untar("/Users/saezlab/Desktop/KD/GSE47184/GSE47184_RAW.tar", exdir = "data")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cels_gse47184 = list.files("data/", pattern = "[gz]")
sapply(paste("data", cels_gse47184, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/




#10/3.  Separate Analysis of samples corresponding to different platforms
#GPL570/Tubular


setwd("~/Desktop/KD/GSE47184")
setwd("~/Desktop/KD/GSE47184/data")


library(simpleaffy)
celfiles_570_184 = read.affy(covdesc = "pheno_gpl570.txt", path = "GPL570")



#no quant
celfiles_570_rma184_nonorm = affy::rma(celfiles_570_184, normalize = FALSE)
save(celfiles_570_rma184_nonorm, file = "celfiles_570_gse47184_rma184_nonorm.RData")


library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse47184_570_nonorm.pdf")
boxplot(celfiles_570_184 , col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse47184_570_nonorm.pdf")
hist(celfiles_570_184, col=cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_570_184, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570_184)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47184_570_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47184_570_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()



#remove samples high IQR and NUSE > 1.05

celfiles_570_184 = celfiles_570_184[, -3]
celfiles_570_184 = celfiles_570_184[, -14]
celfiles_570_184 = celfiles_570_184[, -34]


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_570_184)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47184_570_nonorm_wo3.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47184_570_nonorm_wo3.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47184_570_nonorm_wo5.pdf")
heatmap.2(cor(exprs(celfiles_570_184)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()




#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_570_184)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47184_570.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47184_570_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_570_184)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#Normalization

#1. No normalization
celfiles_570_184_RMAnonorm = affy::rma(celfiles_570_184, normalize = FALSE)
save(celfiles_570_184_RMAnonorm, file = "celfiles_570_184rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.


celfiles_570_filtered_184nonorm = nsFilter(celfiles_570_184_RMAnonorm  , require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_570_filtered_184nonorm, file = "celfiles_570_filtered_184nonorm.RData")



#Annotation: hgu133plus2
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes=row.names(celfiles_570_filtered_184nonorm$eset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse47184_570_nonorm = cbind(Symbols, exprs(celfiles_570_filtered_184nonorm$eset))
gse47184_570_nonorm = gse47184_570_nonorm[order(gse47184_570_nonorm[, 1]), ]

rownames(gse47184_570_nonorm) = gse47184_570_nonorm[,1]
gse47184_570_nonorm = gse47184_570_nonorm[, -1]

col_gse47184_570_nonorm = colnames(gse47184_570_nonorm)
row_gse47184_570_nonorm = rownames(gse47184_570_nonorm)

gse47184_570_nonorm = matrix(as.numeric(gse47184_570_nonorm), ncol = 40, nrow = 20514)

colnames(gse47184_570_nonorm) = col_gse47184_570_nonorm
rownames(gse47184_570_nonorm) = row_gse47184_570_nonorm

save(gse47184_570_nonorm, file = "gse47184_570_nonorm.RData")


REMOVED SAMPLES:
  
  GSM1146389_H7-Tub-DN1110.CEL
  GSM1146401_H7-Tub-MCD1111.CEL
  GSM1146422_H7-Tub-RPGN1130.CEL	


#10/4.  Separate Analysis of samples corresponding to different platforms
#GPL96/ Tubular


setwd("~/Desktop/KD/GSE47184")
setwd("~/Desktop/KD/GSE47184/data")


library(simpleaffy)
celfiles_96_184 = read.affy(covdesc = "pheno_gpl96.txt", path = "GPL96")


setwd("~/Desktop/KD/GSE47184/data/GPL96")
save(celfiles_96_rma, file = "celfiles_96_gse47184_rma.RData")

celfiles_96_rma_184nonorm = affy::rma(celfiles_96_184, normalize = FALSE)
save(celfiles_96_rma_184nonorm, file = "celfiles_96_rma_184nonorm.RData")




library(RColorBrewer)
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised intensity values, log2 transformed and background corrected
pdf("boxplot_intensity_gse47184_96_nonorm.pdf")
boxplot(celfiles_96_184 , col=cols)
dev.off()

# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles

library(affyPLM)


# Plot a density vs log intensity histogram for the unnormalised data
pdf("histogram_intensity_gse47184_96_nonorm.pdf")
hist(celfiles_96_184, col=cols)
dev.off()


#RNA degradation

ar = AffyRNAdeg(celfiles_96_184, log.it = T)
plotAffyRNAdeg(ar)
summaryAffyRNAdeg(ar)



#From these plots we can conclude that there are no major deviations amongst the 12 chips, 
#and normalisation has brought the intensities from all of the chips into distributions with similar characteristics. 
#To take a closer look at the situation on a per-chip level we can use affyPLM. 
#affyPLM allows us to visualise statistical characteristics of the CEL files.


# Perform probe-level metric calculations on the CEL files:
celfiles.qc = fitPLM(celfiles_96_184)

# Create an image of GSM24662.CEL:

image(celfiles.qc, which=1, add.legend=TRUE)

# affyPLM also provides more informative boxplots
# RLE (Relative Log Expression) plots should have
# values close to zero.
pdf("REL_intensity_gse47184_96_nonorm.pdf")
RLE(celfiles.qc, main="RLE")
dev.off()
# We can also use NUSE (Normalised Unscaled Standard Errors).
pdf("NUSE_intensity_gse47184_96_nonorm.pdf")
NUSE(celfiles.qc, main="NUSE")
dev.off()

#We can also look at the relationships between the samples using heirarchical clustering:

eset = exprs(celfiles_570_184)
distance = dist(t(eset), method = "maximum")
clusters = hclust(distance)
pdf("hclust_47184_570.pdf")
plot(clusters)
dev.off()

#hierarchical clustering using the original data
pdf("heatmap_cor_gse47184_96_nonorm.pdf")
heatmap.2(cor(exprs(celfiles_96_184)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#Removal of samples
#removal of cadaveric donor samples
celfiles_96_184 = celfiles_96_184[, -1]
celfiles_96_184 = celfiles_96_184[, -1]
celfiles_96_184 = celfiles_96_184[, -1]
celfiles_96_184 = celfiles_96_184[, -1]

celfiles_96_184 = celfiles_96_184[, -39] #FSGS347

celfiles_96_184 = celfiles_96_184[, -39] #HT11

celfiles_96_184 = celfiles_96_184[, -39] #IgA27

celfiles_96_184 = celfiles_96_184[, -40] #MGN303

celfiles_96_184 = celfiles_96_184[, -39]


REMOVED SAMPLES:
  
GSM1146323_H1-Tub-CD1.CEL	
GSM1146324_H1-Tub-CD2.CEL	
GSM1146325_H1-Tub-CD3.CEL	
GSM1146326_H1-Tub-CD4.CEL	

GSM1146365_H5-Tub-FSGS347.CEL	
GSM1146366_H5-Tub-HT11.CEL	
GSM1146367_H5-Tub-IgA27.CEL	
GSM1146368_H5-Tub-MCD8.CEL	
GSM1146369_H5-Tub-MGN303.CEL	

save(celfiles_96_184, file = "celfiles_96_184.RData")


#hierarchical clustering using the original data
pdf("heatmap_cor_gse47184_96_nonorm_wo9p.pdf")
heatmap.2(cor(exprs(celfiles_96_184)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()



#Normalization

#1. No normalization
celfiles_96_184_RMAnonorm = affy::rma(celfiles_96_184, normalize = FALSE)
save(celfiles_96_184_RMAnonorm, file = "celfiles_96_184rma_nonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_96_filtered = nsFilter(celfiles_96_184_RMAnonorm , require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_96_filtered, file = "celfiles_96_gse47184_filtered.RData")




#Annotation: hgu133a
library(hgu133plus2.db)
library(annotate)
library(hgu133a.db)


probes= row.names(celfiles_96_filtered$eset)
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))



#rty has been filtered for duplicate Entrenz, so it has 50% compared to rty, that has not been filtered for duplicateEntrenz.
gse47184_96_nonorm = cbind(Symbols, exprs(celfiles_96_filtered$eset))
gse47184_96_nonorm = gse47184_96_nonorm[order(gse47184_96_nonorm[, 1]), ]

rownames(gse47184_96_nonorm) = gse47184_96_nonorm[,1]
gse47184_96_nonorm = gse47184_96_nonorm[, -1]

col_gse47184_96_nonorm = colnames(gse47184_96_nonorm)
row_gse47184_96_nonorm = rownames(gse47184_96_nonorm)


gse47184_96_nonorm = matrix(as.numeric(gse47184_96_nonorm), ncol = 55, nrow = 12437)

colnames(gse47184_96_nonorm) = col_gse47184_96_nonorm
rownames(gse47184_96_nonorm) = row_gse47184_96_nonorm

save(gse47184_96_nonorm, file = "gse47184_96_nonorm.RData")

#removing two observervations, IgA and HT, which have only one-one samples

gse47184_96_nonorm = gse47184_96_nonorm[, -44] #twice

save(gse47184_96_nonorm, file = "gse47184_96_nonorm.RData")



#Not included, b/c do not wanna include a third platform in the study
#3. GSE30528, Karolina Woroniecka/ Susztak lab/ Diabetic Nephropathy (22 samples) #Glomeruli. Platform GPL571, but used the one for GPL570

library(GEOquery)
getGEOSuppFiles("GSE30528")

untar("/Users/saezlab/Desktop/KD/GSE30528/GSE30528_RAW.tar", exdir = "DN_Glom")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cel_dnglom = list.files("DN_Glom/", pattern = "[gz]")
sapply(paste("DN_Glom", cel_dnglom, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE30528")

library(simpleaffy)
celfiles_dnglom = read.affy(covdesc = "pheno.txt", path = "DN_Glom")

#no normalization
celfiles_dnglom_rmanonorm = rma(celfiles_dnglom, normalize = FALSE)
save(celfiles_dnglom_rmanonorm, file = "celfiles_dnglom_rmanonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_dnglom_filtered = nsFilter(celfiles_dnglom_rmanonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_dnglom_filtered, file = "celfiles_dnglom_filtered.RData")



#Annotation
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(annotate)


probes=row.names((celfiles_dnglom_filtered$eset))
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))

dnglom_nonorm = cbind(Symbols, exprs(celfiles_dnglom_filtered$eset))
dnglom_nonorm = dnglom_nonorm[order(dnglom_nonorm[, 1]), ]

rownames(dnglom_nonorm) = dnglom_nonorm[,1]
dnglom_nonorm = dnglom_nonorm[, -1]
colnames_dnglom_nonorm = colnames(dnglom_nonorm)
rownames_dnglom_nonorm = rownames(dnglom_nonorm)

dnglom_nonorm = matrix(as.numeric(dnglom_nonorm), ncol = 22, nrow = 12436)

colnames(dnglom_nonorm) = colnames_dnglom_nonorm
rownames(dnglom_nonorm) = rownames_dnglom_nonorm

save(dnglom_nonorm, file = "dnglom_nonorm.RData")



#3. GSE30529, Karolina Woroniecka/ Susztak lab/ Diabetic Nephropathy (22 samples) #Tubulointerstitium. Platform GPL571 !!! third platform. think about the exclusion of this study
#this study was excluded.....
library(GEOquery)
getGEOSuppFiles("GSE30528")

untar("/Users/saezlab/Desktop/KD/GSE30528/GSE30528_RAW.tar", exdir = "DN_Glom")

#use this when untar does not work
system("defaults write org.R-project.R force.LANG en_US.UTF-8")

cel_dnglom = list.files("DN_Glom/", pattern = "[gz]")
sapply(paste("DN_Glom", cel_dnglom, sep = "/"), gunzip)

#Describing experiment --> phenodata --> http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#Loading and Normalising the Data
setwd("~/Desktop/KD/GSE30528")

library(simpleaffy)
celfiles_dnglom = read.affy(covdesc = "pheno.txt", path = "DN_Glom")

#no normalization
celfiles_dnglom_rmanonorm = rma(celfiles_dnglom, normalize = FALSE)
save(celfiles_dnglom_rmanonorm, file = "celfiles_dnglom_rmanonorm.RData")


#Filtering out probesets without Entrenz Gene Identifiers, but no low variance filtering.

celfiles_dnglom_filtered = nsFilter(celfiles_dnglom_rmanonorm, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = FALSE)
save(celfiles_dnglom_filtered, file = "celfiles_dnglom_filtered.RData")



#Annotation
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(annotate)


probes=row.names((celfiles_dnglom_filtered$eset))
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))

dnglom_nonorm = cbind(Symbols, exprs(celfiles_dnglom_filtered$eset))
dnglom_nonorm = dnglom_nonorm[order(dnglom_nonorm[, 1]), ]

rownames(dnglom_nonorm) = dnglom_nonorm[,1]
dnglom_nonorm = dnglom_nonorm[, -1]
colnames_dnglom_nonorm = colnames(dnglom_nonorm)
rownames_dnglom_nonorm = rownames(dnglom_nonorm)

dnglom_nonorm = matrix(as.numeric(dnglom_nonorm), ncol = 22, nrow = 12436)

colnames(dnglom_nonorm) = colnames_dnglom_nonorm
rownames(dnglom_nonorm) = rownames_dnglom_nonorm

save(dnglom_nonorm, file = "dnglom_nonorm.RData")

#################################################################################################################################################################
#Platform-based hierarchical structure creation
#

#1. Glomerulus and GPL96 (Platform1)


Glom_96 = cbind(nsc_nonorm, iga_ercb_nonorm, ber_glom, gse47183_96_nonorm, gse37460_96_nonorm)

#2. Glomerulus and GPL570 (Platform2)

Glom_570 = cbind(iga_toronto_nonorm, gse47183_570_nonorm, gse37460_570_nonorm)



#3. Tubulointerstitium and GPL96 (Platform1)

Tub_96 = cbind(ber_tub, gse47184_96_nonorm, gse37455_96_nonorm, iga_gse35488_nonorm, iga_gse35487_nonorm)
save(Tub_96, file = "Tub_96.RData")

#4. Tubulointerstitium and GPL570 (Platform2)


Tub_570 = cbind(gse47184_570_nonorm, gse37455_570_nonorm, gse69438_nonorm)
save(Tub_570, file = "Tub_570.RData")

save(Glom_96, Glom_570, Tub_96, Tub_570, file = "Platform_Hierarchy_new.RData")
#########################################################################################################################################################
#Batch effect mitigation procedure/Platform 96

rm(list = ls())
library(scater)
library(scran)
library(biomaRt)
library(LSD)
library(limma)
library(gplots)
library(YuGene)
library(RColorBrewer)


##INPUT FILE FOR THE FOLLOWING SCRIPT:
https://drive.google.com/drive/folders/1mBkawIQ0zxbeepJLXyZvvMGJpiDdZ0M2 (Platform_Hierarchy_new.RData)
#use Platform_Hiearachy_new.RData


##############importing the data
load("/Users/admin/Documents/data/ckd/microarray/new/Platform_Hierarchy_new.RData")

glom_96 = as.data.frame(Glom_96)
colnames(glom_96) = sapply(colnames(glom_96), function(x) gsub(" ", "", x, fixed = TRUE))
glom_96 = glom_96[,order(colnames(glom_96))]
glom_96 = glom_96[order(rownames(glom_96)),]
glom_96 = as.matrix(glom_96)
info_glom_96 = read.table("/Users/admin/Documents/data/ckd/microarray/new/Glom_96_Info.csv", sep = ";", header = T)
info_glom_96$Accession = sapply(info_glom_96$Accession, function(x) gsub(" ", "", x, fixed = TRUE))
info_glom_96 = info_glom_96[order(info_glom_96$Accession),]
info_glom_96$Group = as.character(info_glom_96$Group)
if(!(all(colnames(glom_96) == info_glom_96$Accession))) { message("glom_96: names don't match") }



###remove NSC and two patients from FSGS
ww = which(info_glom_96$Group == "NSC")
glom_96 = glom_96[,-ww]
info_glom_96 = info_glom_96[-ww,]
glom_96_final = glom_96


#######################tumor nephrectomy###################
ww = which(info_glom_96$Group == "Tumor Nephrectomy")
norm_expr = normalizeBetweenArrays(glom_96[,ww], method = "cyclicloess")
design = cbind(rep(1,18), c(rep(0,14), rep(1, 4)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]

top = top[order(top$adj.P.Val),]
go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 10)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(glom_96[-www_temp,ww]), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[2]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_tumourNephrectomy.pdf")
(plot(go,go2))
dev.off()
go = go[which.min(go2)]
www = which(rownames(glom_96) %in% rownames(top[1:1000,]))

##pca plot
infoFrame1 = data.frame(Cell = colnames(glom_96), info_glom_96)
rownames(infoFrame1) = colnames(glom_96)
infoFrame_glom_96 = new("AnnotatedDataFrame", data = infoFrame1)
glom_96_y = newSCESet(countData = glom_96, phenoData = infoFrame_glom_96)
glom_96_y = calculateQCMetrics(glom_96_y)
pdf("/Users/admin/Documents/data/ckd/microarray/new/pca_glom_96_tumourNephrectomy.pdf")
print(plotPCA(glom_96_y[-www,ww], ntop = 10000000000, ncomponents = 3, colour_by = "Study", exprs_values = "counts"))
dev.off()

glom_96_final[www,ww] = NA
###########################################################






#######################iga###################
ww = which(info_glom_96$Group == "IgaN")
norm_expr = normalizeBetweenArrays(glom_96[,ww], method = "cyclicloess")
design = cbind(rep(1,18), c(rep(0,14), rep(1, 4)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]

top = top[order(top$adj.P.Val),]
go = seq(1,length(which(top$adj.P.Val < 0.00000000000000005)), by = 1)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(glom_96[-www_temp,ww]), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[2]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_IgAN.pdf")
(plot(go,go2))
dev.off()
www = which(rownames(glom_96) %in% rownames(top[1:1000,]))


##pca plot
infoFrame1 = data.frame(Cell = colnames(glom_96), info_glom_96)
rownames(infoFrame1) = colnames(glom_96)
infoFrame_glom_96 = new("AnnotatedDataFrame", data = infoFrame1)
glom_96_y = newSCESet(countData = glom_96, phenoData = infoFrame_glom_96)
glom_96_y = calculateQCMetrics(glom_96_y)


###Sup.Fig.1D######################
pdf("/Users/admin/Documents/data/ckd/microarray/new/pca_glom_96_IgAN_after.pdf")
print(plotPCA(glom_96_y[-www,ww],  ntop = 10000000000, ncomponents = 2, colour_by = "Study", exprs_values = "counts"))
dev.off()

###Sup.Fig.1A#################
pdf("/Users/admin/Documents/data/ckd/microarray/new/pca_glom_96_IgAN_before.pdf")
print(plotPCA(glom_96_y[,ww],  ntop = 10000000000, ncomponents = 2, colour_by = "Study", exprs_values = "counts"))
dev.off()

###Sup.Fig1C############
pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_IgAN_pvalue.pdf")
plot(all,go2 * 100, type = "l", lwd = 5, xlim = c(15,45), ylim = c(10.5,16))
dev.off()

pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_IgAN_volcano.pdf")
limma::volcanoplot(fit, coef = 2, highlight = 50, names = rownames(top))
dev.off()

####Sup.Fig.1B############
pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_IgAN_ma.pdf")
plot(top$AveExpr, top$logFC, col= "#00000020")
dev.off()


glom_96_final[www,ww] = NA
###########################################################



#######################HLD###################
ww = which(info_glom_96$Group == "Healthy Living Donor")
norm_expr = normalizeBetweenArrays(glom_96[,ww], method = "cyclicloess")
design = cbind(rep(1,17), c(rep(0,3), rep(1,14)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]

top = top[order(top$adj.P.Val),]
go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 10)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(glom_96[-www_temp,ww]), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[2]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/pcaTEST_glom_96_healthyLivingDonor.pdf")
(plot(go,go2))
dev.off()
go = go[which.min(go2)]
www = which(rownames(glom_96) %in% rownames(top[1:3000,]))

##pca plot
infoFrame1 = data.frame(Cell = colnames(glom_96), info_glom_96)
rownames(infoFrame1) = colnames(glom_96)
infoFrame_glom_96 = new("AnnotatedDataFrame", data = infoFrame1)
glom_96_y = newSCESet(countData = glom_96, phenoData = infoFrame_glom_96)
glom_96_y = calculateQCMetrics(glom_96_y)
pdf("/Users/admin/Documents/data/ckd/microarray/new/pca_glom_96_healthyLivingDonor.pdf")
print(plotPCA(glom_96_y[-www,ww], ntop = 10000000000, ncomponents = 3, colour_by = "Study", exprs_values = "counts"))
dev.off()


glom_96_final[www,ww] = NA
#############################################################


###save##
processed_glom_96 = list(data = glom_96, info = info_glom_96, data_NA = glom_96_final)
save(processed_glom_96, file = "/Users/admin/Documents/data/ckd/microarray/new/process_glom_96.RData")
#########################################
#Batch effect mitigation procedure/Platform 570


rm(list = ls())
library(scater)
library(scran)
library(biomaRt)
library(LSD)
library(limma)
library(gplots)
library(YuGene)
library(RColorBrewer)



##############importing the data#########
load("/Users/admin/Documents/data/ckd/microarray/new/Platform_Hierarchy_new.RData")

glom_570 = as.data.frame(Glom_570)
colnames(glom_570) = sapply(colnames(glom_570), function(x) gsub(" ", "", x, fixed = TRUE))
glom_570 = glom_570[,order(colnames(glom_570))]
glom_570 = glom_570[order(rownames(glom_570)),]
glom_570 = as.matrix(glom_570)
info_glom_570 = read.table("/Users/admin/Documents/data/ckd/microarray/new/glom_570_Info.csv", sep = ";", header = T)
info_glom_570$Accession = sapply(info_glom_570$Accession, function(x) gsub(" ", "", x, fixed = TRUE))
info_glom_570 = info_glom_570[order(info_glom_570$Accession),]
info_glom_570$Group = as.character(info_glom_570$Group)
if(!(all(colnames(glom_570) == info_glom_570$Accession))) { message("glom_570: names don't match") }
#########################################




###remove TN and one patient from IgA####
no1 = which(info_glom_570$Group == "Tumor Nephrectomy")
no3 = which(info_glom_570$Group == "IgAN") #remove all because small number compared to other study
no = unique(c(no1,no3))
glom_570 = glom_570[,-no]
info_glom_570 = info_glom_570[-no,]

ww = which(info_glom_570$Group == "Healthy Living Donor")
info_glom_570$Group[ww] = "Healthy Living Donor 570"
###########################################





temp = glom_570
temp = YuGene(temp)
glom_570_y = matrix(0, ncol = ncol(temp), nrow = nrow(temp))
for (i in 1:ncol(temp)) {
	glom_570_y[,i] = temp[,i]
}
rownames(glom_570_y) = rownames(glom_570)
infoFrame1 = data.frame(Cell = colnames(glom_570), info_glom_570)
rownames(infoFrame1) = colnames(glom_570)
infoFrame_glom_570 = new("AnnotatedDataFrame", data = infoFrame1)
glom_570_y = newSCESet(countData = glom_570_y, phenoData = infoFrame_glom_570)
glom_570_y = calculateQCMetrics(glom_570_y)


processed_glom_570 = list(data = glom_570, info = info_glom_570)
save(processed_glom_570, file = "/Users/admin/Documents/data/ckd/microarray/new/process_glom_570.RData")
####################################################################################################################################################

#Merging process_glom_96 and process_glom_570 with their correspinding meta data: info_glom_96 and info_glom_570

#resultant object: 
#process_glom_all.RData")
#processed_glom_all$info -> meta data
#processed_glom_all$final -> expression data

###################################################################################################################

#PCA; Differential expression; MetaDe for p-value aggregation.




rm(list = ls())
library(scater)
library(scran)
library(biomaRt)
library(LSD)
library(limma)
library(gplots)


setwd("~/Desktop/KD/Platform_Hierarchy")
load("~/Desktop/KD/Platform_Hierarchy/Process_All_Glom/process_glom_all.RData")

##INPUT FILE FOR THE FOLLOWING SCRIPT:
https://drive.google.com/drive/folders/1mBkawIQ0zxbeepJLXyZvvMGJpiDdZ0M2 (process_glom_all.RData)
#process_glom_all.RData
#processed_glom_all$info -> meta data
#processed_glom_all$final -> expression data


Glom_All_Info = processed_glom_all$info
save(Glom_All_Info, file = "Glom_All_Info.RData")

#omitting NA, 2587 genes remained.
norm = na.omit(processed_glom_all$final)
save(norm, file = "norm_noNA.RData")

a = norm
a = a[order(rownames(a)),]
a = a[,order(colnames(a))]
info = processed_glom_all$info
info = info[order(info$Accession),]
rownames(info) = info$Accession

i = data.frame(Cell = colnames(norm), info)
rownames(i) = colnames(a)
infoFrame = new("AnnotatedDataFrame", data = info)
sc_glom = newSCESet(countData = a, phenoData = infoFrame)
sc_glom = calculateQCMetrics(sc_glom)

save(sc_glom, file = "sc_glom.RData")


#PCA by Study after removal of the platform batch-effect inducing genes
pdf("pca_byStudy_Glom_noNA.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Study", exprs_values = "counts", ntop = 2587 )
dev.off()

#PCA by Group after removal of the platform batch-effect inducing genes
pdf("pca_byGroup_Glom_noNA.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Group", exprs_values = "counts", ntop = 2587 )
dev.off()

#PCA by platform after removal of the platform batch-effect inducing genes
pdf("pca_byPlatform_Glom_nonNA.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Platform", exprs_values = "counts", ntop = 2587)
dev.off()

#####################################################################################


#For this analysis, we used the version from which the healthy living donor samples were removed and the tumor nephrectomy samples were left in the data set.

#Remove all the Healthy Living Donor Samples, because a lot of genes have been removed from this cohort.

write.csv2(processed_glom_all$final, file = "processed_glom_all_final.csv")

###healthy living donor samples were manually removed from the csv file containing the expression data.

#matching to info data

rownames(processed_glom_all_final) = processed_glom_all_final$X
processed_glom_all_final = processed_glom_all_final[, -1]
save(processed_glom_all_final, file = "processed_glom_all_final_noLD.RData")

write.csv2(info, file = "info_glom.csv")
info_glom <- read.csv("~/Desktop/KD/Platform_Hierarchy/info_glom.csv", sep=";")
rownames(info_glom) = info_glom[,1]
info_glom = info_glom[,-1]

save(info_glom, file = "info_glom.Rdata")


trans_norm = t(norm)
write.csv2(trans_norm, file = "trans_norm.csv")


save(norm, file = "norm.RData") --> #right format
  
  
ty = list(processed_glom_all$final, norm)
ty1 = Reduce(intersect, lapply(ty, colnames))
norm_a = processed_glom_all$final[,ty1 ]


#6289 remaining genes; HLD removed
norm = na.omit(norm_a)
save(norm, file = "norm_noNA.RData")


a = norm
a = a[order(rownames(a)),]
a = a[,order(colnames(a))]
info_glom = info_glom[order(info_glom$Accession),]
rownames(info_glom) = info_glom$Accession


i = data.frame(Cell = colnames(norm), info_glom)
rownames(i) = colnames(a)
infoFrame = new("AnnotatedDataFrame", data = info_glom)
sc_glom = newSCESet(countData = a, phenoData = infoFrame)
sc_glom = calculateQCMetrics(sc_glom)

save(sc_glom, file = "sc_glom.RData")
save(info_glom, infoFrame, file = "info_glom.RData")


#PCA by Study after removal of the platform batch-effect inducing genes
pdf("pca_byStudy_Glom_noHLD.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Study", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by Group after removal of the platform batch-effect inducing genes
pdf("pca_byGroup_Glom_noHLD.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Group", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by platform after removal of the platform batch-effect inducing genes
pdf("pca_byPlatform_Glom_noHLD.pdf")
plotPCA(sc_glom, ncomponents = 3, colour_by = "Platform", exprs_values = "counts", ntop = 6289)
dev.off()




#group-specific plots

#1.1.PCA of gene expression of Tumor Nephrectomy separated by Study
pdf("pca_TN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Tumor Nephrectomy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#1.2.PCA of gene expression of Tumor Nephrectomy separated by Study
pdf("pca_TN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Tumor Nephrectomy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#2.1.PCA of gene expression of DN separated by Study
pdf("pca_DN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Diabetic Nephropathy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#2.2.PCA of gene expression of DN separated by Study
pdf("pca_TN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Diabetic Nephropathy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#3.1.PCA of gene expression of FSGS separated by Study
pdf("pca_FSGS_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "FSGS"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#3.2.PCA of gene expression of FSGS separated by Study
pdf("pca_FSGS_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "FSGS"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#4.1.PCA of gene expression of FSGS_MCD separated by Study
pdf("pca_FSGS_MCD_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "FSGS_MCD"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#4.2.PCA of gene expression of FSGS_MCD separated by Study
pdf("pca_FSGS_MCD_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "FSGS_MCD"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#5.1.PCA of gene expression of HT separated by Study
pdf("pca_HT_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Hypertensive Nephropathy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#5.2.PCA of gene expression of HT separated by Study
pdf("pca_HT_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Hypertensive Nephropathy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#6.1.PCA of gene expression of IgAN separated by Study
pdf("pca_IgAN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "IgAN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#6.2.PCA of gene expression of IgAN separated by Study
pdf("pca_IgAN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "IgAN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#7.1.PCA of gene expression of LN separated by Study
pdf("pca_LN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Lupus Nephritis"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#7.2.PCA of gene expression of LN separated by Study
pdf("pca_TN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "Lupus Nephritis"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#8.1.PCA of gene expression of MCD separated by Study
pdf("pca_MCD_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "MCD"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#8.2.PCA of gene expression of MCD separated by Study
pdf("pca_MCD_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "MCD"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#9.1.PCA of gene expression of MGN separated by Study
pdf("pca_MGN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "MGN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#9.2.PCA of gene expression of MGN separated by Study
pdf("pca_MGN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "MGN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#10.1.PCA of gene expression of RPGN separated by Study
pdf("pca_RPGN_glom_byStudy.pdf")
plotPCA(sc_glom[, infoFrame$Group == "RPGN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#10.2.PCA of gene expression of RPGN separated by Study
pdf("pca_RPGN_glom_byPlatform.pdf")
plotPCA(sc_glom[, infoFrame$Group == "RPGN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#histogram before norm

pdf("hist_glom_all_nonorm.pdf")
hist(counts(sc_glom))
dev.off()

pdf("bozplot_glom_all_nonorm.pdf")
boxplot(counts(sc_glom))
dev.off()

pdf("cor_heat_glom_all_nonorm.pdf")
heatmap.2(cor(counts(sc_glom)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()

#normalization with Cyclic Loess
norm_glom_cl = normalizeBetweenArrays(counts(sc_glom), method = "cyclicloess")
save(norm_glom_cl, file = "norm_glom_cl.RData")

#scater again on cyclic loess normalized data

a = norm_glom_cl
a = a[order(rownames(a)),]
a = a[,order(colnames(a))]
info_glom = info_glom[order(info_glom$Accession),]
rownames(info_glom) = info_glom$Accession


i = data.frame(Cell = colnames(norm_glom_cl), info_glom)
rownames(i) = colnames(a)
infoFrame = new("AnnotatedDataFrame", data = info_glom)
sc_glom_cl = newSCESet(countData = a, phenoData = infoFrame)
sc_glom_cl = calculateQCMetrics(sc_glom_cl)

save(sc_glom_cl, infoFrame, file = "sc_glom_cl.RData")



#PCA by Study after removal of the platform batch-effect inducing genes
pdf("pca_byStudy_Glom_noHLD_cloess.pdf")
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Study", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by Group after removal of the platform batch-effect inducing genes
pdf("pca_byGroup_Glom_noHLD_cloess.pdf")
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Group", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by platform after removal of the platform batch-effect inducing genes
pdf("pca_byPlatform_Glom_noHLD_cloess.pdf")
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Platform", exprs_values = "counts", ntop = 6289)
dev.off()

#histogram before norm> dev.off()

hist(counts(sc_glom_cl))

pdf("hist_glom_cloess.pdf")
hist(counts(sc_glom_cl))
dev.off()

pdf("boxplot_glom_cloess.pdf")
boxplot(counts(sc_glom_cl))
dev.off()

pdf("cor_heat_glom_cloess.pdf")
heatmap.2(cor(counts(sc_glom_cl)), col = colorRampPalette(c("blue","white","hotpink"))(n=200),trace = "none" )
dev.off()


#########Separating the Conditions into distinct objects###################
1.
DN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Diabetic Nephropathy"])

dn_tn = cbind(DN_Glom, TN_Glom)

2. 
FSGS_Glom = counts(sc_glom_cl[, infoFrame$Group == "FSGS"])

fsgs_tn = cbind(FSGS_Glom, TN_Glom)


3. 
FSGS_MCD_Glom = counts(sc_glom_cl[, infoFrame$Group == "FSGS_MCD"])

fsgs_mcd_tn = cbind(FSGS_MCD_Glom, TN_Glom)

fsgs_fsgsmcd = cbind(FSGS_Glom, FSGS_MCD_Glom)



4. 
HT_Glom = counts(sc_glom_cl[, infoFrame$Group == "Hypertensive Nephropathy"])

ht_tn =  cbind(HT_Glom, TN_Glom)


5. 
IgAN_Glom = counts(sc_glom_cl[, infoFrame$Group == "IgAN"])

igan_tn =  cbind(IgAN_Glom, TN_Glom)


6. 
LN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Lupus Nephritis"])

ln_tn =  cbind(LN_Glom, TN_Glom)


7. 
MCD_Glom = counts(sc_glom_cl[, infoFrame$Group == "MCD"])

mcd_tn =  cbind(MCD_Glom, TN_Glom)

mcd_fsgsmcd = cbind(MCD_Glom, FSGS_MCD_Glom)

mcd_fsgs = cbind(FSGS_Glom, MCD_Glom)


8. 
MGN_Glom = counts(sc_glom_cl[, infoFrame$Group == "MGN"])

mgn_tn =  cbind(MGN_Glom, TN_Glom)


9. 
RPGN_Glom = counts(sc_glom_cl[, infoFrame$Group == "RPGN"])

rpgn_tn =  cbind(RPGN_Glom, TN_Glom)


10. 
TN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Tumor Nephrectomy"])




save(DN_Glom, FSGS_Glom, FSGS_MCD_Glom, HT_Glom, IgAN_Glom, LN_Glom, MCD_Glom, MGN_Glom, RPGN_Glom, TN_Glom, file = "ckd_glom_sepMat.RData")

save(dn_tn, fsgs_tn, fsgs_mcd_tn, ht_tn, igan_tn, ln_tn, mcd_tn, mgn_tn, rpgn_tn, fsgs_fsgsmcd, mcd_fsgs, mcd_fsgsmcd,  file = "combined_disease_control.RData")



#Hierarchical Clustering of Normalized Gene Expression of CKD Entities

load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/ckd_glom_sepMat.RData")

source("~/Desktop/KD/Luz/lib_enrichment_scores.r")
library(gplots)
load("~/Desktop/KD/Luz/geneset2.RData")



gg = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,DN_Glom, HT_Glom, RPGN_Glom, TN_Glom)

creating a scaffold matrix that will later be filled with the average of each gene/disease
glom = matrix(9, ncol = 6289, nrow = 10) 



# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {glom[i,] = rowMeans(gg[[i]])
}


rownames(glom) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN")
colnames(glom) = rownames(FSGS_Glom)

glom = t(glom)



#cyclic loess normalized, average gex.
save(glom, file = "glom.RData")

#Scaling and Recentering gene expression using R scale method (Luz)
glom_scaled = gene_expression_statistic(glom, method = "scale", rnaseq = FALSE)

heatmap.2(cor(glom_scaled, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none",
          scale = "none", notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), cexCol = 1.5)




#Limma analysis for differential expression between Disease and Tumor Nephrectomy Samples

#1. Diabetic Nephropathy vs Tumor Nephrectomy

#intersection


#get the fold changes; the data was already log transformed (RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_dn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Diabetic Nephropathy")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_dn, file = "DN_Glom_TN_FC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_cloess_DNglom.pdf")
hist(fc_dn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoDN = infoFrame[infoFrame$Group == "Diabetic Nephropathy"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

idntn = combine(infoDN, infoTN)
save(idntn, file = "idntn_info_DN_TN.RData")

diabetic_vs_tn = data.frame(idntn$Accession, idntn$Study, idntn$Group, idntn$Platform, idntn$Tissue)
colnames(diabetic_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(diabetic_vs_tn, file = "diabetic_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(diabetic_vs_tn, file = "diabetic_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(dn_tn, file = "dn_tn.csv")


#basic limma analysis for differential expression
design_dn = cbind(rep(1,length(idntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Diabetic Nephropathy"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_dn = lmFit(dn_tn, design_dn)
fit2_dn = eBayes(fit_dn)
top_dn = topTable(fit2_dn, adjust="BH", number = 50000000, coef = 2)

save(top_dn, fit2_dn,  file = "DN_TN_Glom_limma.RData")

saveRDS(top_dn, "top_dn.rds")

write.csv2(top_dn, file = "top_dn.csv")


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = dn_tn
ww = which(rownames(aa) %in% rownames(top_dn) [(which(top_dn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_DN_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

DNTN_glom_decide = decideTests(fit2_dn, method = "separate", adjust.method = "BH", p.value = 0.05)
save(DNTN_glom_decide, file = "DNTN_glom_decide.RData")

pdf("maplot_limma_DN_Glom.pdf")
limma::plotMA(fit2_dn, status = DNTN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Diabetic Nephropathy and Tumor Nephrectomy")
dev.off()



#2. FSGS vs Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_fsgs, file = "FSGS_Glom_TN_FC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_FSGS_TN_glom.pdf")
hist(fc_dn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS = infoFrame[infoFrame$Group == "FSGS"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgstn = combine(infoFSGS, infoTN)
save(ifsgstn, file = "ifsgstn_info_FSGS_TN.RData")

fsgs_vs_tn = data.frame(ifsgstn$Accession, ifsgstn$Study, ifsgstn$Group, ifsgstn$Platform, ifsgstn$Tissue)
colnames(fsgs_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(fsgs_vs_tn, file = "fsgs_vs_tn_info.RData")
#info data with FSGS & TN for Hyojin
write.csv2(fsgs_vs_tn, file = "fsgs_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(fsgs_tn, file = "fsgs_tn.csv")


#basic limma analysis for differential expression
design_fsgs = cbind(rep(1,length(ifsgstn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "FSGS"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_fsgs = lmFit(fsgs_tn, design_fsgs)
fit2_fsgs = eBayes(fit_fsgs)
top_fsgs = topTable(fit2_fsgs, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs, fit2_fsgs,  file = "FSGS_TN_Glom_limma.RData")
saveRDS(top_fsgs, "top_fsgs.rds")


write.csv2(top_fsgs, file = "top_fsgs.csv")


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = fsgs_tn
ww = which(rownames(aa) %in% rownames(top_fsgs) [(which(top_fsgs$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_FSGS_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

FSGSTN_glom_decide = decideTests(fit2_fsgs, method = "separate", adjust.method = "BH", p.value = 0.05)
save(FSGSTN_glom_decide, file = "FSGSTN_glom_decide.RData")

pdf("maplot_limma_FSGS_TN_Glom.pdf")
limma::plotMA(fit2_fsgs, status = FSGSTN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Focal Segmental Glomerulosclerosis and Tumor Nephrectomy")
dev.off()


#3. FSGS_MCD vs Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs_mcd = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_fsgs_mcd, file = "FSGS_MCD_Glom_TN_FC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_FSGS_MCD_TN_glom.pdf")
hist(fc_fsgs_mcd, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS_MCD = infoFrame[infoFrame$Group == "FSGS_MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgsmcdtn = combine(infoFSGS_MCD, infoTN)
save(ifsgsmcdtn, file = "ifsgsmcdtn_info_FSGS_TN.RData")


fsgs_mcd_vs_tn = data.frame(ifsgsmcdtn$Accession, ifsgsmcdtn$Study, ifsgsmcdtn$Group, ifsgsmcdtn$Platform, ifsgsmcdtn$Tissue)
colnames(fsgs_mcd_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(fsgs_mcd_vs_tn, file = "fsgs_mcd_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(fsgs_mcd_vs_tn, file = "fsgs_mcd_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(fsgs_mcd_tn, file = "fsgs_mcd_tn.csv")


#basic limma analysis for differential expression
design_fsgs_mcd = cbind(rep(1,length(ifsgsmcdtn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_fsgs_mcd = lmFit(fsgs_mcd_tn, design_fsgs_mcd)
fit2_fsgs_mcd = eBayes(fit_fsgs_mcd)
top_fsgs_mcd = topTable(fit2_fsgs_mcd, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs_mcd, fit2_fsgs_mcd,  file = "FSGS_MCD_TN_Glom_limma.RData")
saveRDS(top_fsgs_mcd, "top_fsgs_mcd.rds")


write.csv2(top_fsgs_mcd, file = "top_fsgs_mcd.csv")


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = fsgs_mcd_tn
ww = which(rownames(aa) %in% rownames(top_fsgs_mcd) [(which(top_fsgs_mcd$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_FSGS_MCD_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

FSGS_MCD_TN_glom_decide = decideTests(fit2_fsgs_mcd, method = "separate", adjust.method = "BH", p.value = 0.05)
save(FSGS_MCD_TN_glom_decide, file = "FSGS_MCD_TN_glom_decide.RData")

pdf("maplot_limma_FSGS_MCD_TN_Glom.pdf")
limma::plotMA(fit2_fsgs_mcd, status = FSGS_MCD_TN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Focal Segmental Glomerulosclerosis + MCD and Tumor Nephrectomy")
dev.off()



#3. FSGS vs FSGS_MCD
#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs_mcd01 = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) ) 

save(fc_fsgs_mcd01, file = "FSGS_MCD01_GlomFC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_FSGS_MCD01_glom.pdf")
hist(fc_fsgs_mcd01, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS_MCD = infoFrame[infoFrame$Group == "FSGS_MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgsmcd = combine(infoFSGS_MCD, infoFSGS)
save(ifsgsmcd, file = "info_FSGSMCD.RData")


#basic limma analysis for differential expression
design_fsgsmcd01 = cbind(rep(1,length(ifsgsmcd$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "FSGS"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"]))))
fit_fsgs_mcd01 = lmFit(fsgs_fsgsmcd, design_fsgsmcd01)
fit2_fsgs_fsgsmcd = eBayes(fit_fsgs_mcd01)
top_fsgs_fsgsmcd = topTable(fit2_fsgs_fsgsmcd, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs_fsgsmcd, fit2_fsgs_fsgsmcd,  file = "FSGS_FSGSMCD_Glom_limma.RData")

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = fsgs_fsgsmcd
ww = which(rownames(aa) %in% rownames(top_fsgs_fsgsmcd) [(which(top_fsgs_fsgsmcd$adj.P.Val < 0.05))])

#no diff. exp. genes between FSGS and FSGS_MCD


#3. MCD vs  FSGS_MCD

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_fsgsmcd = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) ) 

save(fc_mcd_fsgsmcd, file = "MCD_FSGSMCD_FC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_MCD_FSGSMCD_glom.pdf")
hist(fc_mcd_fsgsmcd, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

imcd_fsgsmcd = combine(infoMCD, infoFSGS_MCD)
save(imcd_fsgsmcd, file = "imcd_fsgsmcd.RData")


#basic limma analysis for differential expression
design_mcd_fsgsmcd = cbind(rep(1,length(imcd_fsgsmcd$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"]))))
fit_mcd_fsgsmcd = lmFit(mcd_fsgsmcd, design_mcd_fsgsmcd)
fit2_mcd_fsgsmcd = eBayes(fit_mcd_fsgsmcd)
top_mcd_fsgsmcd = topTable(fit2_mcd_fsgsmcd, adjust="BH", number = 50000000, coef = 2)

save(top_mcd_fsgsmcd, fit2_mcd_fsgsmcd,  file = "MCD_FSSGSMCD_Glom_limma.RData")

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_fsgsmcd
ww = which(rownames(aa) %in% rownames(top_mcd_fsgsmcd) [(which(top_mcd_fsgsmcd$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_MCD_FSGSMCD_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_FSGSMCD_glom_decide = decideTests(fit2_mcd_fsgsmcd, method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_FSGSMCD_glom_decide, file = "MCD_FSGSMCD_glom_decide.RData")

pdf("maplot_limma_MCD_FSGSMCD_Glom.pdf")
limma::plotMA(fit2_mcd_fsgsmcd, status = MCD_FSGSMCD_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Minimal Change & FSGS with Minimal Change")
dev.off()



#3. MCD vs  FSGS

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_fsgs02 = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) ) 

save(fc_mcd_fsgs02, file = "fc_mcd_fsgs02.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_MCD_FSGS02_glom.pdf")
hist(fc_mcd_fsgs02, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

imcdfsgs = combine(infoMCD, infoFSGS)
save(imcdfsgs, file = "imcdfsgs.RData")


#basic limma analysis for differential expression
design_imcdfsgs = cbind(rep(1,length(imcdfsgs$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS"]))))
fit_imcdfsgs = lmFit(mcd_fsgs, design_imcdfsgs)
fit2_imcdfsgs = eBayes(fit_imcdfsgs)
top_imcdfsgs = topTable(fit2_imcdfsgs, adjust="BH", number = 50000000, coef = 2)

save(top_imcdfsgs , fit2_imcdfsgs,  file = "MCD_FSGS_Glom_limma.RData")

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_fsgs
ww = which(rownames(aa) %in% rownames(top_imcdfsgs) [(which(top_imcdfsgs$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_MCD_FSGS_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_FSGS_glom_decide = decideTests(fit2_imcdfsgs , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_FSGS_glom_decide, file = "MCD_FSGS_glom_decide.RData")

pdf("maplot_limma_MCD_FSGS_Glom.pdf")
limma::plotMA(fit2_imcdfsgs, status = MCD_FSGS_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Minimal Change & Focal Segmental Glomerulosclerosis")
dev.off()



#3. Hypertensive Nephropathy vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_ht_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Hypertensive Nephropathy")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_ht_tn, file = "fc_ht_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_HT_TN_glom.pdf")
hist(fc_ht_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoHT = infoFrame[infoFrame$Group == "Hypertensive Nephropathy"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ihttn = combine(infoHT, infoTN)
save(ihttn, file = "ihttn.RData")


ht_vs_tn = data.frame(ihttn$Accession, ihttn$Study, ihttn$Group, ihttn$Platform, ihttn$Tissue)
colnames(ht_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(ht_vs_tn, file = "ht_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(ht_vs_tn, file = "ht_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(ht_tn, file = "ht_tn.csv")




#basic limma analysis for differential expression
design_ihttn = cbind(rep(1,length(ihttn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Hypertensive Nephropathy"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_ihttn = lmFit(ht_tn, design_ihttn)
fit2_ihttn = eBayes(fit_ihttn)
top_ihttn = topTable(fit2_ihttn, adjust="BH", number = 50000000, coef = 2)

save(top_ihttn , fit2_ihttn,  file = "ihttn_Glom_limma.RData")
saveRDS(top_ihttn, "top_ht.rds")


write.csv2(top_ihttn, file = "top_ihttn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = ht_tn
ww = which(rownames(aa) %in% rownames(top_ihttn) [(which(top_ihttn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_HT_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("cyan3","black","lightsalmon1"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

HT_TN_glom_decide = decideTests(fit2_ihttn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(HT_TN_glom_decide, file = "HT_TN_glom_decide.RData")

pdf("maplot_limma_HT_TN_Glom.pdf")
limma::plotMA(fit2_ihttn, status = HT_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Hypertensive Nephropathy & Tumor Nephrectomy")
dev.off()


#3. IgAN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_igan_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "IgAN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_igan_tn, file = "fc_igan_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_IgAN_TN_glom.pdf")
hist(fc_igan_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoIgAN = infoFrame[infoFrame$Group == "IgAN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

igantn = combine(infoIgAN, infoTN)
save(igantn, file = "igantn.RData")


igan_vs_tn = data.frame(igantn$Accession, igantn$Study, igantn$Group, igantn$Platform, igantn$Tissue)
colnames(igan_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(igan_vs_tn, file = "igan_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(igan_vs_tn, file = "igan_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(igan_tn, file = "igan_tn.csv")



#basic limma analysis for differential expression
design_igantn = cbind(rep(1,length(igantn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "IgAN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_igantn = lmFit(igan_tn, design_igantn)
fit2_igantn = eBayes(fit_igantn)
top_igantn = topTable(fit2_igantn, adjust="BH", number = 50000000, coef = 2)

save(top_igantn , fit2_igantn,  file = "igantn_Glom_limma.RData")
saveRDS(top_igantn, "top_igan.rds")


write.csv2(top_igantn, file = "top_igantn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = igan_tn
ww = which(rownames(aa) %in% rownames(top_igantn) [(which(top_igantn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_IgAN_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

IgAN_TN_glom_decide = decideTests(fit2_igantn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(IgAN_TN_glom_decide, file = "IgAN_TN_glom_decide.RData")

pdf("maplot_limma_IgAN_TN_Glom.pdf")
limma::plotMA(fit2_igantn, status = IgAN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between IgA Nephropathy & Tumor Nephrectomy")
dev.off()


#3. LN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_ln_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Lupus Nephritis")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_ln_tn, file = "fc_ln_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_LN_TN_glom.pdf")
hist(fc_ln_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoLN = infoFrame[infoFrame$Group == "Lupus Nephritis"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

lntn = combine(infoLN, infoTN)
save(lntn, file = "lntn.RData")

ln_vs_tn = data.frame(lntn$Accession, lntn$Study, lntn$Group, lntn$Platform, lntn$Tissue)
colnames(ln_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(ln_vs_tn, file = "ln_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(ln_vs_tn, file = "ln_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(ln_tn, file = "ln_tn.csv")



#basic limma analysis for differential expression
design_lntn = cbind(rep(1,length(lntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Lupus Nephritis"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_lntn = lmFit(ln_tn, design_lntn)
fit2_lntn = eBayes(fit_lntn)
top_lntn = topTable(fit2_lntn, adjust="BH", number = 50000000, coef = 2)

save(top_lntn , fit2_lntn,  file = "lntn_Glom_limma.RData")

saveRDS(top_lntn, "top_ln.rds")


write.csv2(top_lntn, file = "top_lntn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = ln_tn
ww = which(rownames(aa) %in% rownames(top_lntn) [(which(top_lntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_LN_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

LN_TN_glom_decide = decideTests(fit2_lntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(LN_TN_glom_decide, file = "LN_TN_glom_decide.RData")

pdf("maplot_limma_LN_TN_Glom.pdf")
limma::plotMA(fit2_lntn, status = LN_TN_glom_decide[,2], main = "DEGs between Lupus Nephritis & Tumor Nephrectomy")
dev.off()


#3. MCD vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_mcd_tn, file = "fc_mcd_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_MCD_TN_glom.pdf")
hist(fc_mcd_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

mcdtn = combine(infoMCD, infoTN)
save(mcdtn, file = "mcdtn.RData")

mcd_vs_tn = data.frame(mcdtn$Accession, mcdtn$Study, mcdtn$Group, mcdtn$Platform, mcdtn$Tissue)
colnames(mcd_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(mcd_vs_tn, file = "mcd_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(mcd_vs_tn, file = "mcd_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(mcd_tn, file = "mcd_tn.csv")



#basic limma analysis for differential expression
design_mcdtn = cbind(rep(1,length(mcdtn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_mcdtn = lmFit(mcd_tn, design_mcdtn)
fit2_mcdtn = eBayes(fit_mcdtn)
top_mcdtn = topTable(fit2_mcdtn, adjust="BH", number = 50000000, coef = 2)

save(top_mcdtn , fit2_mcdtn,  file = "MCD_TN_Glom_limma.RData")

 saveRDS(top_mcdtn, "top_mcd.rds")


write.csv2(top_mcdtn, file = "top_mcdtn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_tn
ww = which(rownames(aa) %in% rownames(top_mcdtn) [(which(top_mcdtn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_MCD_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_TN_glom_decide = decideTests(fit2_mcdtn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_TN_glom_decide, file = "MCD_TN_glom_decide.RData")

pdf("maplot_limma_MCD_TN_Glom.pdf")
limma::plotMA(fit2_mcdtn, status = MCD_TN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Minimal Change & Tumor Nephrectomy")
dev.off()


#3. MGN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mgn_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MGN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_mgn_tn, file = "fc_mgn_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_MGN_TN_glom.pdf")
hist(fc_mgn_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMGN = infoFrame[infoFrame$Group == "MGN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

mgntn = combine(infoMGN, infoTN)
save(mgntn, file = "mgntn.RData")


mgn_vs_tn = data.frame(mgntn$Accession, mgntn$Study, mgntn$Group, mgntn$Platform, mgntn$Tissue)
colnames(mgn_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(mgn_vs_tn, file = "mgn_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(mgn_vs_tn, file = "mgn_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(mgn_tn, file = "mgn_tn.csv")




#basic limma analysis for differential expression
design_mgntn = cbind(rep(1,length(mgntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "MGN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_mgntn = lmFit(mgn_tn, design_mgntn)
fit2_mgntn = eBayes(fit_mgntn)
top_mgntn = topTable(fit2_mgntn, adjust="BH", number = 50000000, coef = 2)

save(top_mgntn , fit2_mgntn,  file = "mgntn_Glom_limma.RData")

saveRDS(top_mgntn, "top_mgn.rds")


write.csv2(top_mgntn, file = "top_mgntn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mgn_tn
ww = which(rownames(aa) %in% rownames(top_mgntn) [(which(top_mgntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_MGN_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MGN_TN_glom_decide = decideTests(fit2_mgntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MGN_TN_glom_decide, file = "MGN_TN_glom_decide.RData")

pdf("maplot_limma_MGN_TN_Glom.pdf")
limma::plotMA(fit2_mgntn, status = MGN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Membranous Nephropathy & Tumor Nephrectomy")
dev.off()


3. RPGN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_rpgn_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "RPGN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_rpgn_tn, file = "fc_rpgn_tn.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf("foldChange_byGroup_RPGN_TN_glom.pdf")
hist(fc_rpgn_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoRPGN = infoFrame[infoFrame$Group == "RPGN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

rpgntn = combine(infoRPGN, infoTN)
save(rpgntn, file = "rpgntn.RData")

rpgn_vs_tn = data.frame(rpgntn$Accession, rpgntn$Study, rpgntn$Group, rpgntn$Platform, rpgntn$Tissue)
colnames(rpgn_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(rpgn_vs_tn, file = "rpgn_vs_tn_info.RData")
#info data with DN & TN for Hyojin
write.csv2(rpgn_vs_tn, file = "rpgn_vs_tn.csv")
#expression of DN and TN for Hyojin
write.csv2(rpgn_tn, file = "rpgn_tn.csv")



#basic limma analysis for differential expression
design_rpgntn = cbind(rep(1,length(rpgntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "RPGN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_rpgntn = lmFit(rpgn_tn, design_rpgntn)
fit2_rpgntn = eBayes(fit_rpgntn)
top_rpgntn = topTable(fit2_rpgntn, adjust="BH", number = 50000000, coef = 2)

save(top_rpgntn , fit2_rpgntn,  file = "rpgntn_Glom_limma.RData")

saveRDS(top_rpgntn, "top_rpgn.rds")


write.csv2(top_rpgntn, file = "top_rpgntn.csv")


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = rpgn_tn
ww = which(rownames(aa) %in% rownames(top_rpgntn) [(which(top_rpgntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf("heatmap_RPGN_TN_Glom_limma.pdf")
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

RPGN_TN_glom_decide = decideTests(fit2_rpgntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(RPGN_TN_glom_decide, file = "RPGN_TN_glom_decide.RData")

pdf("maplot_limma_RPGN_TN_Glom.pdf")
limma::plotMA(fit2_rpgntn, status = RPGN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Rapidly Progressive Glomerulosclerosis & Tumor Nephrectomy")
dev.off()

-----------------------------------------------------------------------------------------------------------------------------
  
  
  #Gene x Adjusted-Pvalue
  #1. DN
  
  dn_limma_adjp = top_dn$adj.P.Val    
DN_p = cbind(rownames(top_dn), top_dn$adj.P.Val)
rownames(DN_p) = rownames(top_dn)
DN_p = as.data.frame(DN_p)
DN_p = DN_p[, -1]
DN_p = matrix(as.numeric(DN_p), ncol = 1, nrow = 6289)
DN_p[,1] = dn_limma_adjp
rownames(DN_p) = rownames(top_dn)

save(DN_p, file = "DN_p.RData")


#Gene x Adjusted-Pvalue
#2. FSGS

fsgs_limma_adjp = top_fsgs$adj.P.Val    
FSGS_p = cbind(rownames(top_fsgs), top_fsgs$adj.P.Val)
rownames(FSGS_p) = rownames(top_fsgs)
FSGS_p = as.data.frame(FSGS_p)
FSGS_p = FSGS_p[, -1]
FSGS_p = matrix(as.numeric(FSGS_p), ncol = 1, nrow = 6289)
FSGS_p[,1] = fsgs_limma_adjp
rownames(FSGS_p) = rownames(top_fsgs)

save(FSGS_p, file = "FSGS_p.RData")


#Gene x Adjusted-Pvalue
#3. FSGS_MCD

fsgsmcd_limma_adjp = top_fsgs_mcd$adj.P.Val    
FSGSMCD_p = cbind(rownames(top_fsgs_mcd), top_fsgs_mcd$adj.P.Val)
rownames(FSGSMCD_p) = rownames(top_fsgs_mcd)
FSGSMCD_p = as.data.frame(FSGSMCD_p)
FSGSMCD_p = FSGSMCD_p[, -1]
FSGSMCD_p = matrix(as.numeric(FSGSMCD_p), ncol = 1, nrow = 6289)
FSGSMCD_p[,1] = fsgsmcd_limma_adjp
rownames(FSGSMCD_p) = rownames(top_fsgs_mcd)

save(FSGSMCD_p, file = "FSGSMCD_p.RData")


#Gene x Adjusted-Pvalue
#4. MCD

mcd_limma_adjp = top_mcdtn$adj.P.Val    
MCD_p = cbind(rownames(top_mcdtn), top_mcdtn$adj.P.Val)
rownames(MCD_p) = rownames(top_mcdtn)
MCD_p = as.data.frame(MCD_p)
MCD_p = MCD_p[, -1]
MCD_p = matrix(as.numeric(MCD_p), ncol = 1, nrow = 6289)
MCD_p[,1] = mcd_limma_adjp
rownames(MCD_p) = rownames(top_mcdtn)

save(MCD_p, file = "MCD_p.RData") 


#Gene x Adjusted-Pvalue
#5. HT

ht_limma_adjp = top_ihttn$adj.P.Val    
HT_p = cbind(rownames(top_ihttn), top_ihttn$adj.P.Val)
rownames(HT_p) = rownames(top_ihttn)
HT_p = as.data.frame(HT_p)
HT_p = HT_p[, -1]
HT_p = matrix(as.numeric(HT_p), ncol = 1, nrow = 6289)
HT_p[,1] = ht_limma_adjp
rownames(HT_p) = rownames(top_ihttn)

save(HT_p, file = "HT_p.RData")


#Gene x Adjusted-Pvalue
#6. IgAN

igan_limma_adjp = top_igantn$adj.P.Val    
IgAN_p = cbind(rownames(top_igantn), top_igantn$adj.P.Val)
rownames(IgAN_p) = rownames(top_igantn)
IgAN_p = as.data.frame(IgAN_p)
IgAN_p = IgAN_p[, -1]
IgAN_p = matrix(as.numeric(IgAN_p), ncol = 1, nrow = 6289)
IgAN_p[,1] = igan_limma_adjp
rownames(IgAN_p) = rownames(top_igantn)

save(IgAN_p, file = "IgAN_p.RData")


#Gene x Adjusted-Pvalue
#7. LN

ln_limma_adjp = top_lntn$adj.P.Val    
LN_p = cbind(rownames(top_lntn), top_lntn$adj.P.Val)
rownames(LN_p) = rownames(top_lntn)
LN_p = as.data.frame(LN_p)
LN_p = LN_p[, -1]
LN_p = matrix(as.numeric(LN_p), ncol = 1, nrow = 6289)
LN_p[,1] = ln_limma_adjp
rownames(LN_p) = rownames(top_lntn)

save(LN_p, file = "LN_p.RData")


#Gene x Adjusted-Pvalue
#8. MGN

mgn_limma_adjp = top_mgntn$adj.P.Val    
MGN_p = cbind(rownames(top_mgntn), top_mgntn$adj.P.Val)
rownames(MGN_p) = rownames(top_mgntn)
MGN_p = as.data.frame(MGN_p)
MGN_p = MGN_p[, -1]
MGN_p = matrix(as.numeric(MGN_p), ncol = 1, nrow = 6289)
MGN_p[,1] = mgn_limma_adjp
rownames(MGN_p) = rownames(top_mgntn)

save(MGN_p, file = "MGN_p.RData")


#Gene x Adjusted-Pvalue
#9. RPGN

rpgn_limma_adjp = top_rpgntn$adj.P.Val    
RPGN_p = cbind(rownames(top_rpgntn), top_rpgntn$adj.P.Val)
rownames(RPGN_p) = rownames(top_rpgntn)
RPGN_p = as.data.frame(RPGN_p)
RPGN_p = RPGN_p[, -1]
RPGN_p = matrix(as.numeric(RPGN_p), ncol = 1, nrow = 6289)
RPGN_p[,1] = rpgn_limma_adjp
rownames(RPGN_p) = rownames(top_rpgntn)

save(RPGN_p, file = "RPGN_p.RData")



#putting the adjusted-pvalues into a matrix
x = matrix(9, nrow = 6289, ncol = 9)
rownames(x) = GENE
#Gene.RData --> genes in alphabetical order.

#1.DN

DN_p = DN_p[order(rownames(DN_p)), ]
DN_p = as.data.frame(DN_p)
head(DN_p)

x[,1] = DN_p[,1]

#2.FSGS

FSGS_p = FSGS_p[order(rownames(FSGS_p)), ]
FSGS_p = as.data.frame(FSGS_p)
head(FSGS_p)

x[,2] = FSGS_p[,1]


#3.FSGS_MCD

FSGSMCD_p = FSGSMCD_p[order(rownames(FSGSMCD_p)), ]
FSGSMCD_p = as.data.frame(FSGSMCD_p)
head(FSGSMCD_p)

x[,3] = FSGSMCD_p[,1]


#4.MCD

MCD_p = MCD_p[order(rownames(MCD_p)), ]
MCD_p = as.data.frame(MCD_p)
head(MCD_p)

x[,4] = MCD_p[,1]


#5.HT

HT_p = HT_p[order(rownames(HT_p)), ]
HT_p = as.data.frame(HT_p)
head(HT_p)

x[,5] = HT_p[,1]



#6.IgAN

IgAN_p = IgAN_p[order(rownames(IgAN_p)), ]
IgAN_p = as.data.frame(IgAN_p)
head(IgAN_p)

x[,6] = IgAN_p[,1]


#7.LN

LN_p = LN_p[order(rownames(LN_p)), ]
LN_p = as.data.frame(LN_p)
head(LN_p)

x[,7] = LN_p[,1]


#8.MGN

MGN_p = MGN_p[order(rownames(MGN_p)), ]
MGN_p = as.data.frame(MGN_p)
head(MGN_p)

x[,8] = MGN_p[,1]


#9.RPGN

RPGN_p = RPGN_p[order(rownames(RPGN_p)), ]
RPGN_p = as.data.frame(RPGN_p)
head(RPGN_p)

x[,9] = RPGN_p[,1]

rownames(x) = rownames(DN_p)

colnames(x) = c("DN", "FSGS", "FSGS-MCD", "MCD", "HT", "IgAN", "LN", "MGN", "RPGN")

save(x, file = "x.RData")


##############################################################################################################################
#MetaDE: R package containig various meta-analytic approaches to aggregate p-values and effect sizes.
#Aim: Extraction of DE genes that have small p-values in ALL studies.

library(MetaDE)

p = list()
p$p = x
p$bp = NULL

res = MetaDE::MetaDE.pvalue(p, meta.method = c("maxP"))
save(res, file = "res.RData") #1st MA result using maxP method


count.DEnumber(res, p.cut = c(0.01, 0.05), q.cut = c(0.01, 0.05))

w = which(rownames(res$ind.p) %in% rownames(x) [which(res$meta.analysis$FDR < 0.01)])

ckd_maxP_001 = x[w, ]    


k = which(rownames(res$ind.p) %in% rownames(x) [which(res$meta.analysis$FDR < 0.05)])

ckd_maxP_005 = x[k, ] 

ckd_glom_maxP_001 = rownames(ckd_maxP_001)

ckd_glom_maxP_005 = rownames(ckd_maxP_005)

save(ckd_glom_maxP_001, ckd_glom_maxP_005, file = "ckd_glom_maxP_gene_names.RData")


ckd_glom_maxP = list()
ckd_glom_maxP$fdr01 = ckd_glom_maxP_001
ckd_glom_maxP$HLDfdr01 = ckd_glomHLD_maxP_001  


save(ckd_glom_maxP, file = "CKD_Glom_maxP_List.RData")




#Supplementary Heatmaps for maxP (fdr < 0.01) genes' expression
#Relative to Tumor Nephrectomy

load("/Users/saezlab/Downloads/CKD_Glom_maxP_List.RData")
load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/RelativeTN/ckd_glom.RData")

load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/RelativeTN/ckd_glom.RData")

maxP_tn = as.data.frame(ckd_glom_maxP$fdr01, row.names = ckd_glom_maxP$fdr01)

gexMaxP_TN = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,
                  DN_Glom, HT_Glom, RPGN_Glom, TN_Glom, maxP_tn)

op = Reduce(intersect, lapply(gexMaxP_TN, rownames))

FSGS_Glom_maxpTN = FSGS_Glom[op, ]
FSGS_MCD_Glom_maxpTN = FSGS_MCD_Glom[op, ]
MCD_Glom_maxpTN = MCD_Glom[op, ]
IgAN_Glom_maxpTN = IgAN_Glom[op, ]
LN_Glom_maxpTN =  LN_Glom[op, ]
MGN_Glom_maxpTN = MGN_Glom[op, ]
DN_Glom_maxpTN = DN_Glom[op, ]
HT_Glom_maxpTN = HT_Glom[op, ]
RPGN_Glom_maxpTN = RPGN_Glom[op, ]
TN_Glom_maxpTN = TN_Glom[op, ]
  
Glom_MP = list(FSGS_Glom_maxpTN, FSGS_MCD_Glom_maxpTN, MCD_Glom_maxpTN, IgAN_Glom_maxpTN, LN_Glom_maxpTN,
                MGN_Glom_maxpTN, DN_Glom_maxpTN, HT_Glom_maxpTN, RPGN_Glom_maxpTN, TN_Glom_maxpTN)




#creating a scaffold matrix that will later be filled with the average of each gene/disease
glom_maxPTN = matrix(9, ncol = 1790, nrow = 10) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {glom_maxPTN[i,] = rowMeans(Glom_MP[[i]])
}


rownames(glom_maxPTN) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN")
colnames(glom_maxPTN) = rownames(FSGS_Glom_maxpTN)

glom_maxPTN = t(glom_maxPTN)

save(glom_maxPTN, file = "glom_maxPTN.RData")



#Spearman Correlation of maxP Genes in Subtypes
heatmap.2(cor(glom_maxPTN_norm, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none",
           scale = "none", notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), cexCol = 1.5)


##Relative to Healthy Living Donor

save(DN_Glom_hld, FSGS_Glom_hld, FSGS_MCD_Glom_hld, HT_Glom_hld, IgAN_Glom_hld, LN_Glom_hld, MCD_Glom_hld, MGN_Glom_hld, RPGN_Glom_hld, HLD_Glom_hld, file = "ckd_glom_HLD_sepMat.RData")

load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/With_HLD/Glom_Rel_HLD/ckd_glom_HLD_sepMat.RData")


maxP_hld = as.data.frame(ckd_glom_maxP$HLDfdr01, row.names = ckd_glom_maxP$HLDfdr01)

gexMaxP_HLD = list(FSGS_Glom_hld, FSGS_MCD_Glom_hld, MCD_Glom_hld, IgAN_Glom_hld, LN_Glom_hld, MGN_Glom_hld,
                  DN_Glom_hld, HT_Glom_hld, RPGN_Glom_hld, HLD_Glom_hld, maxP_hld)

hop = Reduce(intersect, lapply(gexMaxP_HLD, rownames))

FSGS_Glom_hld_maxpHLD = FSGS_Glom_hld[hop, ]
FSGS_MCD_Glom_hld_maxpHLD = FSGS_MCD_Glom_hld[hop, ]
MCD_Glom_hld_maxpHLD = MCD_Glom_hld[hop, ]
IgAN_Glom_hld_maxpHLD = IgAN_Glom_hld[hop, ]
LN_Glom_hld_maxpHLD =  LN_Glom_hld[hop, ]
MGN_Glom_hld_maxpHLD = MGN_Glom_hld[hop, ]
DN_Glom_hld_maxpHLD = DN_Glom_hld[hop, ]
HT_Glom_hld_maxpHLD = HT_Glom_hld[hop, ]
RPGN_Glom_hld_maxpHLD = RPGN_Glom_hld[hop, ]
HLD_Glom_maxpHLD = HLD_Glom_hld[hop, ]

Glom_MP_hld = list(FSGS_Glom_hld_maxpHLD, FSGS_MCD_Glom_hld_maxpHLD, MCD_Glom_hld_maxpHLD, IgAN_Glom_hld_maxpHLD, LN_Glom_hld_maxpHLD,
                   MGN_Glom_hld_maxpHLD, DN_Glom_hld_maxpHLD, HT_Glom_hld_maxpHLD, RPGN_Glom_hld_maxpHLD, HLD_Glom_maxpHLD)


save(Glom_MP_hld, gexMaxP_HLD, file = "Glom_MP_HLD.RData")


#creating a scaffold matrix that will later be filled with the average of each gene/disease
glom_maxPHLD = matrix(9, ncol = 517, nrow = 10) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {glom_maxPHLD[i,] = rowMeans(Glom_MP_hld[[i]])
}


rownames(glom_maxPHLD) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "HLD")
colnames(glom_maxPHLD) = rownames(FSGS_Glom_hld_maxpHLD)

glom_maxPHLD = t(glom_maxPHLD)



source("~/Desktop/KD/Luz/lib_enrichment_scores.r")
library(gplots)
load("~/Desktop/KD/Luz/geneset2.RData")

glom_maxPHLD_norm = gene_expression_statistic(glom_maxPHLD, method = "scale", rnaseq = FALSE)

save(glom_maxPHLD, glom_maxPHLD_norm, file = "glom_maxPHLD.RData")


heatmap.2(glom_maxPHLD_norm, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "hotpink"))(n=100), 
          Rowv = NA, dendrogram = "col",
          notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), labRow = rownames(glom_maxPHLD_norm), labCol = colnames(glom_maxPHLD_norm), cexRow = 1, cexCol = 1.5)


#Spearman Correlation of maxP Genes in Subtypes
heatmap.2(cor(glom_maxPHLD_norm, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none",
          scale = "none", notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), cexCol = 1.5)



###Union of relative TN and  relative HLD maxP genes

set.union <- function(a, b) {
  u <- a
  for (i in 1:length(b)) {
    if (!(b[i] %in% u)) {
      u <- append(u, b[i])
    }
  }
  return(u)
}

#2137 genes
maxP_U = set.union(a= ckd_glom_maxP$fdr01, b = ckd_glom_maxP$HLDfdr01)
save(maxP_U, file = "maxP_U.RData")

maxP_U = as.data.frame(maxP_U, row.names = maxP_U)

gexMaxP_U = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,
                  DN_Glom, HT_Glom, RPGN_Glom, TN_Glom, HLD_Glom_hld, maxP_U)

up = Reduce(intersect, lapply(gexMaxP_U, rownames))

FSGS_Glom_maxpU = FSGS_Glom[up, ]
FSGS_MCD_Glom_maxpU = FSGS_MCD_Glom[up, ]
MCD_Glom_maxpU = MCD_Glom[up, ]
IgAN_Glom_maxpU = IgAN_Glom[up, ]
LN_Glom_maxpU =  LN_Glom[up, ]
MGN_Glom_maxpU = MGN_Glom[up, ]
DN_Glom_maxpU = DN_Glom[up, ]
HT_Glom_maxpU = HT_Glom[up, ]
RPGN_Glom_maxpU = RPGN_Glom[up, ]
TN_Glom_maxpU = TN_Glom[up, ]
HLD_Glom_maxpU = HLD_Glom_hld[up, ]

Glom_U = list(FSGS_Glom_maxpU, FSGS_MCD_Glom_maxpU, MCD_Glom_maxpU, IgAN_Glom_maxpU, LN_Glom_maxpU,
               MGN_Glom_maxpU, DN_Glom_maxpU, HT_Glom_maxpU, RPGN_Glom_maxpU, TN_Glom_maxpU, HLD_Glom_maxpU)




#creating a scaffold matrix that will later be filled with the average of each gene/disease
glom_maxPU = matrix(9, ncol = 1017, nrow = 11) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:11) {glom_maxPU[i,] = rowMeans(Glom_U[[i]])
}


rownames(glom_maxPU) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN", "HLD")
colnames(glom_maxPU) = rownames(FSGS_Glom_maxpU)

glom_maxPU = t(glom_maxPU)

save(glom_maxPU, file = "glom_maxPU.RData")


source("~/Desktop/KD/Luz/lib_enrichment_scores.r")
library(gplots)
load("~/Desktop/KD/Luz/geneset2.RData")

glom_maxPU_norm = gene_expression_statistic(glom_maxPU, method = "scale", rnaseq = FALSE)


save(Glom_U, gexMaxP_U, glom_maxPU_norm, file = "Glom_maxP_Union.RData")


heatmap.2(glom_maxPU_norm, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "hotpink"))(n=100), 
          Rowv = NA, dendrogram = "col",
          notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), labRow = rownames(glom_maxPU_norm), labCol = colnames(glom_maxPU_norm), cexRow = 1, cexCol = 1.5)


#Spearman Correlation of maxP Genes in Subtypes
heatmap.2(cor(glom_maxPU_norm, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none",
          scale = "none", notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), cexCol = 1.5)


################################################################################################################################################################
#Effect Size Computation
#This part and its outcome were not included in the CKD manuscript.


#Converting means scores with pooled sd to an effect size (sd mean difference)


#1. Diabetic Nephropathy
#mean value of each gene in Diabetic Nephropathy
dn_mean = data.frame(rowMeans(DN_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_dn = apply(DN_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_dn = sqrt((sd_dn^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
dn_mes2 = data.frame(id = 1:6289, m1 = dn_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_dn, n.t = 13, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_dn_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = dn_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(dn_mean, tn_mean, sd_dn, sd_tn, sd_pooled_dn, dn_mes2, eff_dn_mes2, file = "effsize_DN_vs_TN.RData")



#2. FSGS
#mean value of each gene in FSGS
fsgs_mean = data.frame(rowMeans(FSGS_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_fsgs = apply(FSGS_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_fsgs = sqrt((sd_fsgs^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
fsgs_mes2 = data.frame(id = 1:6289, m1 = fsgs_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_fsgs, n.t = 23, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_fsgs_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = fsgs_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(fsgs_mean, tn_mean, sd_fsgs, sd_tn, sd_pooled_fsgs, fsgs_mes2, eff_fsgs_mes2, file = "effsize_FSGS_vs_TN.RData")




#3. FSGS_MCD
#mean value of each gene in FSGS
fsgsmcd_mean = data.frame(rowMeans(FSGS_MCD_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_fsgsmcd = apply(FSGS_MCD_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_fsgsmcd = sqrt((sd_fsgsmcd^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
fsgsmcd_mes2 = data.frame(id = 1:6289, m1 = fsgsmcd_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_fsgsmcd, n.t = 6, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_fsgsmcd_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = fsgsmcd_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(fsgsmcd_mean, tn_mean, sd_fsgsmcd, sd_tn, sd_pooled_fsgsmcd, fsgsmcd_mes2, eff_fsgsmcd_mes2, file = "effsize_FSGSMCD_vs_TN.RData")


#4. MCD

#mean value of each gene in FSGS
mcd_mean = data.frame(rowMeans(MCD_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_mcd = apply(MCD_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_mcd = sqrt((sd_mcd^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
mcd_mes2 = data.frame(id = 1:6289, m1 = mcd_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_mcd, n.t = 15, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_mcd_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = mcd_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(mcd_mean, tn_mean, sd_mcd, sd_tn, sd_pooled_mcd, mcd_mes2, eff_mcd_mes2, file = "effsize_MCD_vs_TN.RData")


#5. HT
#mean value of each gene in FSGS
ht_mean = data.frame(rowMeans(HT_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_ht = apply(HT_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_ht = sqrt((sd_ht^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
ht_mes2 = data.frame(id = 1:6289, m1 = ht_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_ht, n.t = 14, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_ht_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = ht_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(ht_mean, tn_mean, sd_ht, sd_tn, sd_pooled_ht, ht_mes2, eff_ht_mes2, file = "effsize_HT_vs_TN.RData")



#6. IgAN
#mean value of each gene in FSGS
igan_mean = data.frame(rowMeans(IgAN_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_igan = apply(IgAN_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_igan = sqrt((sd_igan^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
igan_mes2 = data.frame(id = 1:6289, m1 = igan_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_igan, n.t = 36, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_igan_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = igan_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(igan_mean, tn_mean, sd_igan, sd_tn, sd_pooled_igan, igan_mes2, eff_igan_mes2, file = "effsize_IgAN_vs_TN.RData")




#7. LN
#mean value of each gene in FSGS
ln_mean = data.frame(rowMeans(LN_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_ln = apply(LN_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_ln = sqrt((sd_ln^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
ln_mes2 = data.frame(id = 1:6289, m1 = ln_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_ln, n.t = 32, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_ln_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = ln_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(ln_mean, tn_mean, sd_ln, sd_tn, sd_pooled_ln, ln_mes2, eff_ln_mes2, file = "effsize_LN_vs_TN.RData")



#8. MGN
#mean value of each gene in FSGS
mgn_mean = data.frame(rowMeans(MGN_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_mgn = apply(MGN_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_mgn = sqrt((sd_mgn^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
mgn_mes2 = data.frame(id = 1:6289, m1 = mgn_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_mgn, n.t = 20, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_mgn_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = mgn_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(mgn_mean, tn_mean, sd_mgn, sd_tn, sd_pooled_mgn, mgn_mes2, eff_mgn_mes2, file = "effsize_MGN_vs_TN.RData")




#9. RPGN
#mean value of each gene in FSGS
rpgn_mean = data.frame(rowMeans(RPGN_Glom))

#mean value of each gene in Tumor Nephrectomy
tn_mean = data.frame(rowMeans(TN_Glom))


#standard deviation of each gene in Diabetic Nephropathy
sd_rpgn = apply(RPGN_Glom, 1, function(x) sd(x))

##standard deviation of each gene in Tumor Nephrectomy
sd_tn = apply(TN_Glom, 1, function(x) sd(x))

#pooled standard deviation of DN and TN
sd_pooled_rpgn = sqrt((sd_rpgn^2 + sd_tn^2)/2)

#data frame with the above computed things + the sample size of each group
rpgn_mes2 = data.frame(id = 1:6289, m1 = rpgn_mean[,1],m2 = tn_mean[,1], sd.pooled = sd_pooled_rpgn, n.t = 22, n.c =18)


#Effect size computation with other things, including the variance of the effect size
eff_rpgn_mes2 = compute.es::mes2(m.1 = m1, m.2 = m2, n.1 = n.t, n.2 = n.c,  s.pooled = sd.pooled, id = id, data = rpgn_mes2)

---------------------------------------------------------------------------------
  #2. t-staistics-based computation of effect sizes (t-statistic from limma)
  
  dat = data.frame(id = 1:7739, t = top_dn$t, n.t = 7, n.c =18)
uff = tes(t=t, n.1=n.t, n.2=n.c, level=95, dig=2, id=id, data=dat)
-----------------------------------------------------------------------------------
  
  save(rpgn_mean, tn_mean, sd_rpgn, sd_tn, sd_pooled_rpgn, rpgn_mes2, eff_rpgn_mes2, file = "effsize_RPGN_vs_TN.RData")



ckd_glom_g = matrix(9, ncol=9, nrow = 6289)
colnames(ckd_glom_g) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "Lupus Nephritis", "MGN", "DN", "HT", "RPGN")
rownames(ckd_glom_g) = rownames(x)

ckd_glom_g[,1] = eff_fsgs_mes2$g
ckd_glom_g[,2] = eff_fsgsmcd_mes2$g
ckd_glom_g[,3] = eff_mcd_mes2$g
ckd_glom_g[,4] = eff_igan_mes2$g
ckd_glom_g[,5] = eff_ln_mes2$g
ckd_glom_g[,6] = eff_mgn_mes2$g
ckd_glom_g[,7] = eff_dn_mes2$g
ckd_glom_g[,8] = eff_ht_mes2$g
ckd_glom_g[,9] = eff_rpgn_mes2$g

save(ckd_glom_g, file = "ckd_glom_hedgesg.RData")


#Ranked list of Gene Expression Effect Size Values for Each Disease--> input for GOrilla


#1. DN
dn_hedges = ckd_glom_g[,1]
dn_hedges = unlist(dn_hedges)
dn_hedges = matrix(as.numeric(dn_hedges), ncol = 1, nrow = 6289)
rownames(dn_hedges) = rownames(ckd_glom_g)
colnames(dn_hedges) = c("DN")
dn_hedges_ranked = dn_hedges[order(dn_hedges, decreasing = TRUE), ]
save(dn_hedges, dn_hedges_ranked, file = "dn_hedges_ranked.RData")
write.csv2(dn_hedges_ranked, file = "dn_hedges_ranked.csv")


#2. FSGS
fsgs_hedges = ckd_glom_g[,2]
fsgs_hedges = unlist(fsgs_hedges)
fsgs_hedges = matrix(as.numeric(fsgs_hedges), ncol = 1, nrow = 6289)
rownames(fsgs_hedges) = rownames(ckd_glom_g)
colnames(fsgs_hedges) = c("FSGS")
fsgs_hedges_ranked = fsgs_hedges[order(fsgs_hedges, decreasing = TRUE), ]
save(fsgs_hedges, fsgs_hedges_ranked, file = "fsgs_hedges_ranked.RData")
write.csv2(fsgs_hedges_ranked, file = "fsgs_hedges_ranked.csv")


#3. FSGS & MCD
fsgsmcd_hedges = ckd_glom_g[,3]
fsgsmcd_hedges = unlist(fsgsmcd_hedges)
fsgsmcd_hedges = matrix(as.numeric(fsgsmcd_hedges), ncol = 1, nrow = 6289)
rownames(fsgsmcd_hedges) = rownames(ckd_glom_g)
colnames(fsgsmcd_hedges) = c("FSGS & MCD")
fsgsmcd_hedges_ranked = fsgsmcd_hedges[order(fsgsmcd_hedges, decreasing = TRUE), ]
save(fsgsmcd_hedges, fsgsmcd_hedges_ranked, file = "fsgsmcd_hedges_ranked.RData")
write.csv2(fsgsmcd_hedges_ranked, file = "fsgsmcd_hedges_ranked.csv")


#4. MCD
mcd_hedges = ckd_glom_g[,4]
mcd_hedges = unlist(mcd_hedges)
mcd_hedges = matrix(as.numeric(mcd_hedges), ncol = 1, nrow = 6289)
rownames(mcd_hedges) = rownames(ckd_glom_g)
colnames(mcd_hedges) = c("MCD")
mcd_hedges_ranked = mcd_hedges[order(mcd_hedges, decreasing = TRUE), ]
save(mcd_hedges, mcd_hedges_ranked, file = "mcd_hedges_ranked.RData")
write.csv2(mcd_hedges_ranked, file = "mcd_hedges_ranked.csv")

#5. HT
ht_hedges = ckd_glom_g[,5]
ht_hedges = unlist(ht_hedges)
ht_hedges = matrix(as.numeric(ht_hedges), ncol = 1, nrow = 6289)
rownames(ht_hedges) = rownames(ckd_glom_g)
colnames(ht_hedges) = c("HT")
ht_hedges_ranked = ht_hedges[order(ht_hedges, decreasing = TRUE), ]
save(ht_hedges, ht_hedges_ranked, file = "ht_hedges_ranked.RData")
write.csv2(ht_hedges_ranked, file = "ht_hedges_ranked.csv")


#6. IgAN
igan_hedges = ckd_glom_g[,6]
igan_hedges = unlist(igan_hedges)
igan_hedges = matrix(as.numeric(igan_hedges), ncol = 1, nrow = 6289)
rownames(igan_hedges) = rownames(ckd_glom_g)
colnames(igan_hedges) = c("IgAN")
igan_hedges_ranked = igan_hedges[order(igan_hedges, decreasing = TRUE), ]
save(igan_hedges, igan_hedges_ranked, file = "igan_hedges_ranked.RData")
write.csv2(igan_hedges_ranked, file = "igan_hedges_ranked.csv")

#7. LN
ln_hedges = ckd_glom_g[,7]
ln_hedges = unlist(ln_hedges)
ln_hedges = matrix(as.numeric(ln_hedges), ncol = 1, nrow = 6289)
rownames(ln_hedges) = rownames(ckd_glom_g)
colnames(ln_hedges) = c("LN")
ln_hedges_ranked = ln_hedges[order(ln_hedges, decreasing = TRUE), ]
save(ln_hedges, ln_hedges_ranked, file = "ln_hedges_ranked.RData")
write.csv2(ln_hedges_ranked, file = "ln_hedges_ranked.csv")

#8. MGN
mgn_hedges = ckd_glom_g[,8]
mgn_hedges = unlist(mgn_hedges)
mgn_hedges = matrix(as.numeric(mgn_hedges), ncol = 1, nrow = 6289)
rownames(mgn_hedges) = rownames(ckd_glom_g)
colnames(mgn_hedges) = c("MGN")
mgn_hedges_ranked = mgn_hedges[order(mgn_hedges, decreasing = TRUE), ]
save(mgn_hedges, mgn_hedges_ranked, file = "mgn_hedges_ranked.RData")
write.csv2(mgn_hedges_ranked, file = "mgn_hedges_ranked.csv")


#9. RPGN
rpgn_hedges = ckd_glom_g[,9]
rpgn_hedges = unlist(rpgn_hedges)
rpgn_hedges = matrix(as.numeric(rpgn_hedges), ncol = 1, nrow = 6289)
rownames(rpgn_hedges) = rownames(ckd_glom_g)
colnames(rpgn_hedges) = c("RPGN")
rpgn_hedges_ranked = rpgn_hedges[order(rpgn_hedges, decreasing = TRUE), ]
save(rpgn_hedges, rpgn_hedges_ranked, file = "rpgn_hedges_ranked.RData")
write.csv2(rpgn_hedges_ranked, file = "rpgn_hedges_ranked.csv")

#############################################################################################################################

#bind the logFC values together


#Gene x logFC values
#1. DN

dn_limma_logFC = top_dn$logFC   
DN_l = cbind(rownames(top_dn), top_dn$logFC)
rownames(DN_l) = rownames(top_dn)
DN_l = as.data.frame(DN_l)
DN_l = DN_l[, -1]
DN_l = matrix(as.numeric(DN_l), ncol = 1, nrow = 6289)
DN_l[,1] = dn_limma_logFC
rownames(DN_l) = rownames(top_dn)

save(DN_l, file = "DN_l.RData")


#Gene x logFC
#2. FSGS

fsgs_limma_logFC = top_fsgs$logFC   
FSGS_l = cbind(rownames(top_fsgs), top_fsgs$logFC)
rownames(FSGS_l) = rownames(top_fsgs)
FSGS_l = as.data.frame(FSGS_l)
FSGS_l = FSGS_l[, -1]
FSGS_l = matrix(as.numeric(FSGS_l), ncol = 1, nrow = 6289)
FSGS_l[,1] = fsgs_limma_logFC
rownames(FSGS_l) = rownames(top_fsgs)

save(FSGS_l, file = "FSGS_l.RData")


#Gene x logFC
#3. FSGS_MCD

fsgsmcd_limma_logFC = top_fsgs_mcd$logFC    
FSGSMCD_l = cbind(rownames(top_fsgs_mcd), top_fsgs_mcd$logFC)
rownames(FSGSMCD_l) = rownames(top_fsgs_mcd)
FSGSMCD_l = as.data.frame(FSGSMCD_l)
FSGSMCD_l = FSGSMCD_l[, -1]
FSGSMCD_l = matrix(as.numeric(FSGSMCD_l), ncol = 1, nrow = 6289)
FSGSMCD_l[,1] = fsgsmcd_limma_logFC
rownames(FSGSMCD_l) = rownames(top_fsgs_mcd)

save(FSGSMCD_l, file = "FSGSMCD_l.RData")


#Gene x logFC
#4. MCD

mcd_limma_logFC = top_mcdtn$logFC    
MCD_l = cbind(rownames(top_mcdtn), top_mcdtn$logFC)
rownames(MCD_l) = rownames(top_mcdtn)
MCD_l = as.data.frame(MCD_l)
MCD_l = MCD_l[, -1]
MCD_l = matrix(as.numeric(MCD_l), ncol = 1, nrow = 6289)
MCD_l[,1] = mcd_limma_logFC
rownames(MCD_l) = rownames(top_mcdtn)

save(MCD_l, file = "MCD_l.RData") 


#Gene x logFC
#5. HT

ht_limma_logFC = top_ihttn$logFC    
HT_l = cbind(rownames(top_ihttn), top_ihttn$logFC)
rownames(HT_l) = rownames(top_ihttn)
HT_l = as.data.frame(HT_l)
HT_l = HT_l[, -1]
HT_l = matrix(as.numeric(HT_l), ncol = 1, nrow = 6289)
HT_l[,1] = ht_limma_logFC
rownames(HT_l) = rownames(top_ihttn)

save(HT_l, file = "HT_l.RData")


#Gene x logFC
#6. IgAN

igan_limma_logFC = top_igantn$logFC 
IgAN_l = cbind(rownames(top_igantn), top_igantn$logFC)
rownames(IgAN_l) = rownames(top_igantn)
IgAN_l = as.data.frame(IgAN_l)
IgAN_l = IgAN_l[, -1]
IgAN_l = matrix(as.numeric(IgAN_l), ncol = 1, nrow = 6289)
IgAN_l[,1] = igan_limma_logFC
rownames(IgAN_l) = rownames(top_igantn)

save(IgAN_l, file = "IgAN_l.RData")


#Gene x logFC
#7. LN

ln_limma_logFC = top_lntn$logFC
LN_l = cbind(rownames(top_lntn), top_lntn$logFC)
rownames(LN_l) = rownames(top_lntn)
LN_l = as.data.frame(LN_l)
LN_l = LN_l[, -1]
LN_l = matrix(as.numeric(LN_l), ncol = 1, nrow = 6289)
LN_l[,1] = ln_limma_logFC
rownames(LN_l) = rownames(top_lntn)

save(LN_l, file = "LN_l.RData")


#Gene x logFC
#8. MGN

mgn_limma_logFC = top_mgntn$logFC  
MGN_l = cbind(rownames(top_mgntn), top_mgntn$logFC)
rownames(MGN_l) = rownames(top_mgntn)
MGN_l = as.data.frame(MGN_l)
MGN_l = MGN_l[, -1]
MGN_l = matrix(as.numeric(MGN_l), ncol = 1, nrow = 6289)
MGN_l[,1] = mgn_limma_logFC
rownames(MGN_l) = rownames(top_mgntn)

save(MGN_l, file = "MGN_l.RData")


#Gene x logFC
#9. RPGN

rpgn_limma_logFC = top_rpgntn$logFC    
RPGN_l = cbind(rownames(top_rpgntn), top_rpgntn$logFC)
rownames(RPGN_l) = rownames(top_rpgntn)
RPGN_l = as.data.frame(RPGN_l)
RPGN_l = RPGN_l[, -1]
RPGN_l = matrix(as.numeric(RPGN_l), ncol = 1, nrow = 6289)
RPGN_l[,1] = rpgn_limma_logFC
rownames(RPGN_l) = rownames(top_rpgntn)

save(RPGN_l, file = "RPGN_l.RData")



#putting the logFC values into a matrix
ckd_glom_logFC = matrix(9, nrow = 6289, ncol = 9)
rownames(ckd_glom_logFC) = rownames(ckd_glom_g)
#Gene.RData --> genes in alphabetical order.

#1.DN

DN_l = DN_l[order(rownames(DN_l)), ]
DN_l = as.data.frame(DN_l)
head(DN_l)

ckd_glom_logFC[,1] = DN_l[,1]

#2.FSGS

FSGS_l = FSGS_l[order(rownames(FSGS_l)), ]
FSGS_l = as.data.frame(FSGS_l)
head(FSGS_l)

ckd_glom_logFC[,2] = FSGS_l[,1]


#3.FSGS_MCD

FSGSMCD_l = FSGSMCD_l[order(rownames(FSGSMCD_l)), ]
FSGSMCD_l = as.data.frame(FSGSMCD_l)
head(FSGSMCD_l)

ckd_glom_logFC[,3] = FSGSMCD_l[,1]


#4.MCD

MCD_l = MCD_l[order(rownames(MCD_l)), ]
MCD_l = as.data.frame(MCD_l)
head(MCD_l)

ckd_glom_logFC[,4] = MCD_l[,1]


#5.HT

HT_l = HT_l[order(rownames(HT_l)), ]
HT_l = as.data.frame(HT_l)
head(HT_l)

ckd_glom_logFC[,5] = HT_l[,1]



#6.IgAN

IgAN_l = IgAN_l[order(rownames(IgAN_l)), ]
IgAN_l = as.data.frame(IgAN_l)
head(IgAN_l)

ckd_glom_logFC[,6] = IgAN_l[,1]


#7.LN

LN_l = LN_l[order(rownames(LN_l)), ]
LN_l = as.data.frame(LN_l)
head(LN_l)

ckd_glom_logFC[,7] = LN_l[,1]


#8.MGN

MGN_l = MGN_l[order(rownames(MGN_l)), ]
MGN_l = as.data.frame(MGN_l)
head(MGN_l)

ckd_glom_logFC[,8] = MGN_l[,1]


#9.RPGN

RPGN_l = RPGN_l[order(rownames(RPGN_l)), ]
RPGN_l = as.data.frame(RPGN_l)
head(RPGN_l)

ckd_glom_logFC[,9] = RPGN_l[,1]


colnames(ckd_glom_logFC) = c("DN", "FSGS", "FSGS-MCD", "MCD", "HT", "IgAN", "LN", "MGN", "RPGN")

save(ckd_glom_logFC, file = "ckd_glom_logFC.RData")

######################################################################################################




rm(list = ls())
library(gplots)
library(YuGene)
library(RColorBrewer)
library(scater)
library(scran)
library(destiny)


ckd_colors = c("#a6cee3","#1f78b4","black","#33a02c","#fb9a99","#e31a1c","#999999","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
ckd_names = c("FSGS/MCD","MCD","IgA Nephropathy","FSGS","Tumor Nephrectomy","Healthy Living Donor","Lupus Nephritis","Diabetic Nephropathy","TMD","Hypertensive Nephropathy","MGN","RPGN")
ckd_names = c(ckd_names[6],ckd_names[5],ckd_names[4],ckd_names[1],ckd_names[2],ckd_names[3],ckd_names[7],ckd_names[11],ckd_names[8],ckd_names[10],ckd_names[9],ckd_names[12])
ckd_colors = c(ckd_colors[6],ckd_colors[5],ckd_colors[4],ckd_colors[1],ckd_colors[2],ckd_colors[3],ckd_colors[7],ckd_colors[11],ckd_colors[8],ckd_colors[10],ckd_colors[9],ckd_colors[12])

#order: [1]HLD, [2]Tumor Nephrectomy, [3]FSGS, [4]FSGS/MCD, [5]MCD, [6]IgAN, [7]Lupus Nephritis, [8]MGN, [9]DN, [10]Hypertensive Neph., [11]TMD, [12]RPGN

load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/process_glom_all.RData")


gex = processed_glom_all$final
info = processed_glom_all$info
ww1 = which(info$Group == "Healthy Living Donor")
ww2 = which(info$Group == "Healthy Living Donor 570")
info$Group[c(ww1,ww2)] = "HealthyLivingDonor"
ww3 = which(info$Group == "Tumor Nephrectomy")
ww = unique(sort(c(ww1,ww2,ww3)))


###CKD_numbers###Figure1B###########################
ckd_numbers = table(info$Group, info$Platform)
ckd_numbers = rbind(ckd_numbers[4,], ckd_numbers[11,], ckd_numbers[2,], ckd_numbers[3,], ckd_numbers[8,], ckd_numbers[6,], ckd_numbers[7,], ckd_numbers[9,], ckd_numbers[1,], ckd_numbers[5,], ckd_numbers[10,])
rownames(ckd_numbers) = ckd_names[-11]
pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/number_of_patients_glom_platform.pdf")
heatmap.2(ckd_numbers, scale = "none", dendrogram = "none", Rowv = NA, Colv = NA, col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Reds")))(n=100), trace = "none", density.info = "none", margins = c(15,15), cellnote = ckd_numbers, notecol = "black", breaks = seq(0,40,length.out = 101))
dev.off()

ckd_numbers = table(info$Group, info$Study)
ckd_numbers = rbind(ckd_numbers[4,], ckd_numbers[11,], ckd_numbers[2,], ckd_numbers[3,], ckd_numbers[8,], ckd_numbers[6,], ckd_numbers[7,], ckd_numbers[9,], ckd_numbers[1,], ckd_numbers[5,], ckd_numbers[10,])
rownames(ckd_numbers) = ckd_names[-11]
pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/number_of_patients_glom_study.pdf")
heatmap.2(ckd_numbers, scale = "none", dendrogram = "none", Rowv = NA, Colv = NA, col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Reds")))(n=100), trace = "none", density.info = "none", margins = c(15,15), cellnote = ckd_numbers, notecol = "black", breaks = seq(0,40,length.out = 101))
dev.off()
##################





########Correlation Heatmap############Figure1C###########################
ckd_colors2 = c(ckd_colors[9],ckd_colors[3],ckd_colors[4],ckd_colors[1],ckd_colors[10],ckd_colors[6],ckd_colors[7],ckd_colors[5],ckd_colors[8],ckd_colors[12],ckd_colors[2])
pal_study = c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0")
pal_platform = c("#377eb8","#e41a1c")



aa = cor(gex, method = "spearman", use = "pairwise.complete.obs")



aa2 = aa
rownames(aa2) = info$Group
colnames(aa2) = info$Group
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(rownames(aa2) == (unique(rownames(aa2)))[i])
  colors[wtemp] = ckd_colors2[i]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/corr_gex_glom_cluster_group.pdf")
heatmap.2(aa2, scale = "none", dendrogram = "row", col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Blues")))(n=100), trace = "none", RowSideColors = colors, ColSideColors = colors, density.info = "none", labRow = "", labCol = "", margins = c(25,25))
legend("topright",legend = ckd_names[-11], col = ckd_colors[-11], lty= 1, lwd = 10, bty = "n")
dev.off()





aa2 = aa
rownames(aa2) = info$Study
colnames(aa2) = info$Study
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors2 = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(colnames(aa2) == (unique(rownames(aa2)))[i])
  colors2[wtemp] = pal_study[i]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/corr_gex_glom_cluster_study.pdf")
heatmap.2(aa2, scale = "none", dendrogram = "row", col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Blues")))(n=100), trace = "none", RowSideColors = colors2, ColSideColors = colors2, density.info = "none", labRow = "", labCol = "", margins = c(25,25))
legend("topright",legend = unique(rownames(aa2)), col = pal_study, lty= 1, lwd = 10, bty = "n")
dev.off()





aa2 = aa
rownames(aa2) = info$Platform
colnames(aa2) = info$Platform
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors3 = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(colnames(aa2) == (unique(rownames(aa2)))[i])
  colors3[wtemp] = pal_platform[i]
}
pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/corr_gex_glom_cluster_platform.pdf")
heatmap.2(aa2, scale = "none", dendrogram = "row", col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Blues")))(n=100), trace = "none", RowSideColors = colors3, ColSideColors = colors3, density.info = "none", labRow = "", labCol = "", margins = c(25,25))
legend("topright",legend = unique(rownames(aa2)), col = pal_platform, lty= 1, lwd = 10, bty = "n")
dev.off()
###################################################



################ % variance explained#############SupplementaryFigure1E###########################
www_all = which((apply(gex, 1, function(x) any(is.na(x)))) == TRUE)
temp = gex[-www_all,]
temp = YuGene(temp)
gex_y = matrix(0, ncol = ncol(temp), nrow = nrow(temp))
for (i in 1:ncol(temp)) {
  gex_y[,i] = temp[,i]
}
rownames(gex_y) = rownames(gex[-www_all,])

infoFrame1 = data.frame(Cell = colnames(gex), info)
rownames(infoFrame1) = colnames(gex)
infoFrame1 = new("AnnotatedDataFrame", data = infoFrame1)

gex_y = newSCESet(countData = gex_y, phenoData = infoFrame1)
gex_y = calculateQCMetrics(gex_y)

pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure1/variance_gex_glom.pdf")
gex_y = plotExplanatoryVariables(gex_y, exprs_values = "counts", variables = c("Group", "Study", "Platform"), return_object=TRUE)
dev.off()
############################################################################################


###################diffusion map Figure2B & SupplementaryFigure2##############################
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

load("/Users/saezlab/Downloads/CKD_Glom_maxP_List (1).RData")
ckd_glom_maxP = ckd_glom_maxP$fdr01
gex1 = gex[which(rownames(gex) %in% ckd_glom_maxP),]


www_all = which((apply(gex1, 1, function(x) any(is.na(x)))) == TRUE)

temp = gex1[-www_all,]
temp = YuGene(temp)
gex1_y = matrix(0, ncol = ncol(temp), nrow = nrow(temp))
for (i in 1:ncol(temp)) {
  gex1_y[,i] = temp[,i]
}
rownames(gex1_y) = rownames(gex1[-www_all,])

colnames(gex1_y) = info$Group
gex1_y = gex1_y[,order(colnames(gex1_y))]
info2 = info
info2 = info2[,order(colnames(info2))]
colnames(gex1_y) = info2$Accession

infoFrame1 = data.frame(Cell = colnames(gex1_y), info2)
rownames(infoFrame1) = colnames(gex1)
infoFrame1 = new("AnnotatedDataFrame", data = infoFrame1)

gex1_y = newSCESet(countData = gex1_y, phenoData = infoFrame1)
gex1_y = calculateQCMetrics(gex1_y)
dm = DiffusionMap(gex1_y, sigma = "local")
ev = eigenvectors(dm)


pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure2/tn_only/diffusion_gex_glom.pdf")
scatterplot3d(ev[,c(1,2,3)], box = F, grid = F, color = colorsT, pch = 16, type = "p")
addgrids3d(ev[,c(1,2,3)], grid = c("xy","xz","yz"), col.grid = "#99999950")
dev.off()


pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure2/tn_only/diffusion_gex_glom12.pdf")
plot(ev[,c(1,2)], col = colorsT, pch = 16, type = "p")
dev.off()

pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure2/tn_only/diffusion_gex_glom13.pdf")
plot(ev[,c(1,3)], col = colorsT, pch = 16, type = "p")
dev.off()


pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure2/tn_only/diffusion_gex_glom23.pdf")
plot(ev[,c(2,3)], col = colorsT, pch = 16, type = "p")
dev.off()
###################################################################

#######################Healthy vs Tumor; not included in the manuscript##########################
rm(list = ls())
library(scater)
library(scran)
library(biomaRt)
library(LSD)
library(limma)
library(gplots)
library(YuGene)
library(RColorBrewer)



load("/Users/admin/Documents/data/ckd/microarray/new/process_glom_all.RData")
info = processed_glom_all$info
gex = processed_glom_all$final



who = (apply(gex, 1, function(x) any(is.na(x))))
who2 = which(who != TRUE)
who = which(who == TRUE)
gex = gex[who2,]


hld = which(info$Group == "Healthy Living Donor")
tn = which(info$Group == "Tumor Nephrectomy")
hlde = gex[,hld]
tne = gex[,tn]

norm_expr = normalizeBetweenArrays(cbind(hlde,tne), method = "cyclicloess")
hlde = rowMeans(norm_expr[,1:17])
tne = rowMeans(norm_expr[,18:35])
design = cbind(rep(1,35), c(rep(1,17), rep(0,18)))

fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]
top = top[order(top$adj.P.Val),]

ww = (which(top$adj.P.Val < 0.05))
ww_up = (which(top$adj.P.Val < 0.05 & top$logFC > 0))
ww_down = (which(top$adj.P.Val < 0.05 & top$logFC < 0))



pdf("/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure3/maplot.pdf")
plot(top[-ww,]$AveExpr, top[-ww,]$logFC, col = "#00000050", ylim = c(-3,3), xlim = c(3,14), pch = 20)
points(top[ww_up,]$AveExpr, top[ww_up,]$logFC, col = "#FF69B450", pch = 20)
points(top[ww_down,]$AveExpr, top[ww_down,]$logFC, col = "#0000FF50", pch = 20)
legend("topright", c("Tumor Nephrectomy","No Change","Healthy Living Donor"), fill = c("hotpink","black","blue"), bty = "n")
dev.off()



write(rownames(top[ww_up,]), file = "/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure3/maplot_up_tnf.txt", ncolumns = 1)
write(rownames(top[ww_down,]), file = "/Users/admin/Documents/data/ckd/microarray/new/ferenc/final_plots/Figure3/maplot_up_healthy.txt", ncolumns = 1)

#top[ww_up[1:10],]
#top[ww_down[1:10],]

#########################################################################




load("/Users/saezlab/Documents/Schubert Model Matrix/model.RData")  #model.RData: PROGENy model matrix found in github --> saezlab/progeny/R/model.r

load("ckd_glom_sepMat.RData") #created in Process_All_Glom.r   Part --> #########Separating the Conditions into distinct objects###################



#Scaling relative to tumor nephrectomy

meancontrolALL = apply(TN_Glom, 1, mean)
meancontrolALL = matrix(as.numeric(meancontrolALL), ncol = 1, nrow = 6289)
sdcontrolALL = apply(TN_Glom ,1, sd)
sdcontrolALL = matrix(as.numeric(sdcontrolALL), ncol = 1, nrow = 6289)

save(meancontrolALL, sdcontrolALL, file = "Mean_SD_TumorN_Glom_ALL.RData")



1. #FSGS Relative to Control
for (i in 1:length(FSGS_Glom)) {
  FSGS_Glom_relativeTN = (FSGS_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


2. #FSGS-MCd Relative to Control
for (i in 1:length(FSGS_MCD_Glom)) {
  FSGS_MCD_Glom_relativeTN = (FSGS_MCD_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}



3. #MCD Relative to Control
for (i in 1:length(MCD_Glom)) {
  MCD_Glom_relativeTN = (MCD_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


4. #IgAN Relative to Control
for (i in 1:length(IgAN_Glom)) {
  IgAN_Glom_relativeTN = (IgAN_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}



5. #LN Relative to Control
for (i in 1:length(LN_Glom)) {
  LN_Glom_relativeTN = (LN_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


6. #MGN Relative to Control
for (i in 1:length(MGN_Glom)) {
  MGN_Glom_relativeTN = (MGN_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


7. #DN Relative to Control
for (i in 1:length(DN_Glom)) {
  DN_Glom_relativeTN = (DN_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


8. #HT Relative to Control
for (i in 1:length(HT_Glom)) {
  HT_Glom_relativeTN = (HT_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}


9. #RPGN Relative to Control
for (i in 1:length(RPGN_Glom)) {
  RPGN_Glom_relativeTN = (RPGN_Glom[1:6289, ] - meancontrolALL[1:6289, ]) / sdcontrolALL[1:6289, ]
}



save(DN_Glom_relativeTN, FSGS_MCD_Glom_relativeTN, FSGS_Glom_relativeTN, IgAN_Glom_relativeTN, LN_Glom_relativeTN, MGN_Glom_relativeTN,
     MCD_Glom_relativeTN, RPGN_Glom_relativeTN, HT_Glom_relativeTN, file = "CKD_Glom_Relative_to_TN_ALL_Genes.RData")


#List, RowMeans-->

ckd_glom_relTN = list(FSGS_Glom_relativeTN, FSGS_MCD_Glom_relativeTN, MCD_Glom_relativeTN, IgAN_Glom_relativeTN, LN_Glom_relativeTN, MGN_Glom_relativeTN,
                      DN_Glom_relativeTN, HT_Glom_relativeTN, RPGN_Glom_relativeTN)

save(ckd_glom_relTN, file = "ckd_glom_relTN.RData")


#creating a scaffold matrix that will later be filled with the average of each gene/disease
ckdglom_RTN = matrix(9, ncol = 6289, nrow = 9) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:9) {ckdglom_RTN[i,] = rowMeans(ckd_glom_relTN[[i]])
}


rownames(ckdglom_RTN) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN")
colnames(ckdglom_RTN) = rownames(FSGS_Glom_relativeTN)

ckdglom_RTN = t(ckdglom_RTN)

save(ckdglom_RTN, file = "ckdglom_RTN.RData")


save(ckdglom_RTN, file = "ckdglom_RTN.rds")


oad("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/RelativeTN/ckd_glom_prog11_RelTN.RData")

Part 8. # Computing the original/observed PROGENy scores.

load("/Users/saezlab/Documents/Schubert Model Matrix/model.RData")

#PROGENy model matrix used can be found:
https://drive.google.com/drive/folders/1mBkawIQ0zxbeepJLXyZvvMGJpiDdZ0M2 (model.RData)  or github.com/saezlab/progeny/R/model.r

8/1 # Getting the overlapping genes between the gene expression matrix (ckd_glom_g: gene x sample matrix containig the average expression of each gene/condition)) and model from saezlab/progeny/R/model.r (GitHub).

list_glom_prog11 = list(ckd_glom_g, model)

id11 = Reduce(intersect, lapply(list_glom_prog11, rownames))

glom_p11 = ckd_glom_g[id11, ]

model11_glom = model[id11, ]


glom_p11  = t(glom_p11 )


8/2. #PROGENy multiplication #model11_glom: PROGENy model matrix found in github --> saezlab/progeny/R/model.r


glom_prog11 = glom_p11  %*% model11_glom
######################################################################################################################################

Part 1.# Matching the gene expression matricies with the PROGENy model matrix

#All expression values represent change from tumor nephrectomy. 

load("/Users/saezlab/Desktop/KD/Platform_Hierarchy/Process_All_Glom/RelativeTN/CKD_Glom_Relative_to_TN_ALL_Genes.RData")
gex_prog11 = list(FSGS_Glom_relativeTN, FSGS_MCD_Glom_relativeTN, MCD_Glom_relativeTN, IgAN_Glom_relativeTN,
                  LN_Glom_relativeTN, MGN_Glom_relativeTN, DN_Glom_relativeTN, HT_Glom_relativeTN, RPGN_Glom_relativeTN, model)

list_gex_p11 = Reduce(intersect, lapply(gex_prog11, rownames))

FSGS_Glom_relativeTN_m = FSGS_Glom_relativeTN[list_gex_p11, ]
FSGS_MCD_Glom_relativeTN_m = FSGS_MCD_Glom_relativeTN[list_gex_p11, ]
MCD_Glom_relativeTN_m = MCD_Glom_relativeTN[list_gex_p11, ]
IgAN_Glom_relativeTN_m = IgAN_Glom_relativeTN[list_gex_p11, ]
LN_Glom_relativeTN_m = LN_Glom_relativeTN[list_gex_p11, ]
MGN_Glom_relativeTN_m = MGN_Glom_relativeTN[list_gex_p11, ]
DN_Glom_relativeTN_m = DN_Glom_relativeTN[list_gex_p11, ]
HT_Glom_relativeTN_m = HT_Glom_relativeTN[list_gex_p11, ]
RPGN_Glom_relativeTN_m = RPGN_Glom_relativeTN[list_gex_p11, ]
model_m = model[list_gex_p11, ]

save(FSGS_Glom_relativeTN_m, FSGS_MCD_Glom_relativeTN_m, MCD_Glom_relativeTN_m, IgAN_Glom_relativeTN_m,
     LN_Glom_relativeTN_m, MGN_Glom_relativeTN_m, DN_Glom_relativeTN_m, 
     HT_Glom_relativeTN_m, RPGN_Glom_relativeTN_m, model_m, file = "Glom_RelTN_MatchedP11.RData")

##############################################################################################################

Part 2. #merging the result of Part 1 into a single matrix:


Glom_Matched11 = cbind(FSGS_Glom_relativeTN_m, FSGS_MCD_Glom_relativeTN_m, MCD_Glom_relativeTN_m, IgAN_Glom_relativeTN_m,
                       LN_Glom_relativeTN_m, MGN_Glom_relativeTN_m, DN_Glom_relativeTN_m, 
                       HT_Glom_relativeTN_m, RPGN_Glom_relativeTN_m)
save(Glom_Matched11, file = "Glom_Matched11.RData")

#################################################################################################################

Part 3. #vector with as many repetitions as many samples we have for each condition. 
vec= c(rep("FSGS_Glom_relativeTN_m", 23), rep("FSGS_MCD_Glom_relativeTN_m", 6), 
       rep("MCD_Glom_relativeTN_m", 15), rep("IgAN_Glom_relativeTN_m", 36), 
       rep("LN_Glom_relativeTN_m", 32), rep("MGN_Glom_relativeTN_m", 20), rep("DN_Glom_relativeTN_m", 13), rep("HT_Glom_relativeTN_m", 14),
       rep("RPGN_Glom_relativeTN_m", 22))

save(vec, file = "vec.RData")

##############################################################################################################

Part 4. #permuting the coloumns 10000 times in the gene expression matrix; generation 10000 matrices contaning these randomized expression values:

  
shuffle_glom = as.matrix(lapply(1:10000, function(x) Glom_Matched11[,sample(ncol(Glom_Matched11))]))

###############################################################################################################

Part 5. #taking the average of each gene's expression in a way that takes the number of samples/condition into account.

newfun = function(num, shuffle_glom, vec) {
  
  unique_vec = unique(vec)
  av_shuffle = matrix(9, ncol = length(unique_vec), nrow = nrow(shuffle_glom[[num]]))
  
  for (i in 1:length(unique_vec)) {
    av_shuffle[,i] = rowMeans(shuffle_glom[[num]][,which(vec == unique_vec[i])])
  }
  return(av_shuffle)
}

something = lapply(1:10000, newfun, shuffle_glom = shuffle_glom, vec = vec)

################################################################################################################ 

Part 6.

6/1.
# creating the transpose of the generated 10000 matrices and storing them in variable trans_something.


transpose = function(i, something) {
  
  t_something = t(something[[i]])
  return(t_something) }


trans_something = lapply(1:10000, transpose, something)   



6/2.
# matrix multiplication to get single PROGENy score/pathway/disease.


#model_m: PROGENy model matrix found in github --> saezlab/progeny/R/model.r

winnie = function(i, trans_something) {
  
  progeny_shuffle = trans_something[[i]] %*% model_m
  return(progeny_shuffle) }


prog_shuf = lapply(1:10000, winnie, trans_something)


#################################################################################################################

Part 7. 

#compute the empirical cumulative distribution function (ECDF) in a pathway specific manner.
#1. EGFR
e = function(i, prog_shuf) {
  
  egfr_prog_shuf = prog_shuf[[i]][ , 1]    #substitute with the desired number of coloumn (pathway)
  
  return(egfr_prog_shuf) }

egfr_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_egfr_glom = ecdf(unlist(egfr_dist))

save(ecdf_egfr_glom, file = "ecdf_egfr_glom.RData")


#2. Hypoxia
e = function(i, prog_shuf) {
  
  hypoxia_prog_shuf = prog_shuf[[i]][ , 2]    #substitute with the desired number of coloumn (pathway)
  
  return(hypoxia_prog_shuf) }

hypoxia_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_hypoxia_glom = ecdf(unlist(hypoxia_dist))

save(ecdf_hypoxia_glom, file = "ecdf_hypoxia_glom.RData")


#3. JAK-STAT
e = function(i, prog_shuf) {
  
  jakstat_prog_shuf = prog_shuf[[i]][ , 3]    #substitute with the desired number of coloumn (pathway)
  
  return(jakstat_prog_shuf) }

jakstat_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_jakstat_glom = ecdf(unlist(jakstat_dist))

save(ecdf_jakstat_glom, file = "ecdf_jakstat_glom.RData")


#4. MAPK
e = function(i, prog_shuf) {
  
  mapk_prog_shuf = prog_shuf[[i]][ , 4]    #substitute with the desired number of coloumn (pathway)
  
  return(mapk_prog_shuf) }

mapk_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_mapk_glom = ecdf(unlist(mapk_dist))

save(ecdf_mapk_glom, file = "ecdf_mapk_glom.RData")


#5. NFkB
e = function(i, prog_shuf) {
  
  nfkb_prog_shuf = prog_shuf[[i]][ , 5]    #substitute with the desired number of coloumn (pathway)
  
  return(nfkb_prog_shuf) }

nfkb_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_nfkb_glom = ecdf(unlist(nfkb_dist))

save(ecdf_nfkb_glom, file = "ecdf_nfkb_glom.RData")


#6. PI3K
e = function(i, prog_shuf) {
  
  pi3k_prog_shuf = prog_shuf[[i]][ , 6]    #substitute with the desired number of coloumn (pathway)
  
  return(pi3k_prog_shuf) }

pi3k_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_pi3k_glom = ecdf(unlist(pi3k_dist))

save(ecdf_pi3k_glom, file = "ecdf_pi3k_glom.RData")


#7. TGFb
e = function(i, prog_shuf) {
  
  tgfb_prog_shuf = prog_shuf[[i]][ , 7]    #substitute with the desired number of coloumn (pathway)
  
  return(tgfb_prog_shuf) }

tgfb_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_tgfb_glom = ecdf(unlist(tgfb_dist))

save(ecdf_tgfb_glom, file = "ecdf_tgfb_glom.RData")

#8. TNFa
e = function(i, prog_shuf) {
  
  tnfa_prog_shuf = prog_shuf[[i]][ , 8]    #substitute with the desired number of coloumn (pathway)
  
  return(tnfa_prog_shuf) }

tnfa_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_tnfa_glom = ecdf(unlist(tnfa_dist))

save(ecdf_tnfa_glom, file = "ecdf_tnfa_glom.RData")


#9. Trail
e = function(i, prog_shuf) {
  
  trail_prog_shuf = prog_shuf[[i]][ , 9]    #substitute with the desired number of coloumn (pathway)
  
  return(trail_prog_shuf) }

trail_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_trail_glom = ecdf(unlist(trail_dist))

save(ecdf_trail_glom, file = "ecdf_trail_glom.RData")


#10. VEGF
e = function(i, prog_shuf) {
  
  vegf_prog_shuf = prog_shuf[[i]][ , 10]    #substitute with the desired number of coloumn (pathway)
  
  return(vegf_prog_shuf) }

vegf_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_vegf_glom = ecdf(unlist(vegf_dist))

save(ecdf_vegf_glom, file = "ecdf_vegf_glom.RData")


#11. p53
e = function(i, prog_shuf) {
  
  p53_prog_shuf = prog_shuf[[i]][ , 11]    #substitute with the desired number of coloumn (pathway)
  
  return(p53_prog_shuf) }

p53_dist = lapply(1:10000, e, prog_shuf = prog_shuf)



ecdf_p53_glom = ecdf(unlist(p53_dist))

save(ecdf_p53_glom, file = "ecdf_p53_glom.RData")



Part 9. #P-value calculation


#Formulas:

#pval = 1 - pathway_spec_ecdf(the ##real## progeny score you want to test)  #this formula for positive PROGENy scores

#pval = pathway_spec_ecdf(the ##real## progeny score you want to test)      #this formula for negative PROGENy scores

  
  oldsport = function(x,y){
    ifelse(x < 0, y(x), 1 - y(x))
  }


1.#EGFR-spec. pvalues

egfr_glomP = oldsport(x = glom_prog11[1:9, 1], y = ecdf_egfr_glom)

2.#Hypoxia-spec. pvalues

hypoxia_glomP = oldsport(x = glom_prog11[1:9, 2], y = ecdf_hypoxia_glom)

3.#JS-spec. pvalues

jakstat_glomP = oldsport(x = glom_prog11[1:9, 3], y = ecdf_jakstat_glom)

4.# MAPK-spec. pvalues

mapk_glomP = oldsport(x = glom_prog11[1:9, 4], y = ecdf_mapk_glom)

5.#NFKB-spec. pvalues

nfkb_glomP = oldsport(x = glom_prog11[1:9, 5], y = ecdf_nfkb_glom)


6.#PI3K-spec. pvalues

pi3k_glomP = oldsport(x = glom_prog11[1:9, 6], y = ecdf_pi3k_glom)

7.#TGFB-spec. pvalues

tgfb_glomP = oldsport(x = glom_prog11[1:9, 7], y = ecdf_tgfb_glom)


8.#TNFa-spec. pvalues

tnfa_glomP = oldsport(x = glom_prog11[1:9, 8], y = ecdf_tnfa_glom)


9.#Trail-spec. pvalues

trail_glomP = oldsport(x = glom_prog11[1:9, 9], y = ecdf_trail_glom)


10.#VEGF-spec. pvalues

vegf_glomP = oldsport(x = glom_prog11[1:9, 10], y = ecdf_vegf_glom)


11.#P53-spec. pvalues

p53_glomP = oldsport(x = glom_prog11[1:9, 11], y = ecdf_p53_glom)



# creating a matrix filled with the p-values for all the pathways in all conditions. 
scaffold_glom_pval = matrix(9, ncol = ncol(glom_prog11), nrow = nrow(glom_prog11[1:9, ]))


#put it generally -->
scaffold_glom_pval[,1] = egfr_glomP
scaffold_glom_pval[,2] = hypoxia_glomP
scaffold_glom_pval[,3] = jakstat_glomP
scaffold_glom_pval[,4] = mapk_glomP
scaffold_glom_pval[,5] = nfkb_glomP
scaffold_glom_pval[,6] = pi3k_glomP
scaffold_glom_pval[,7] = tgfb_glomP
scaffold_glom_pval[,8] = tnfa_glomP
scaffold_glom_pval[,9] = trail_glomP
scaffold_glom_pval[,10] = vegf_glomP
scaffold_glom_pval[,11] = p53_glomP

rownames(scaffold_glom_pval) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN")
colnames(scaffold_glom_pval) = colnames(glom_prog11)

GlomP = scaffold_glom_pval
save(GlomP, file = "GlomP.RData")


Part 10. #FDR-adjustment for p-values --> 

#Computing FDR-adjusted p-values. This is just an example. All this has to be done for all the samples and put it to a matrix.


GlomFDR = matrix(9, ncol = ncol(glom_prog11), nrow = nrow(glom_prog11[1:9, ]))


1.#FSGS
fsgs_glomP = as.vector(GlomP[1, ])
fsgs_glomFDR = p.adjust(fsgs_glomP, method = "BH")    #do it for all subtypes. We use the Benjamini-Hochberg - adjustment

2.#FSGS-MCD
fsgs_mcd_glomP = as.vector(GlomP[2, ])
fsgs_mcd_glomFDR = p.adjust(fsgs_mcd_glomP, method = "BH")    

3.#MCD
mcd_glomP = as.vector(GlomP[3, ])
mcd_glomFDR = p.adjust(mcd_glomP, method = "BH")    

4.#IgAN
igan_glomP = as.vector(GlomP[4, ])
igan_glomFDR = p.adjust(igan_glomP, method = "BH")    

5.#LN
ln_glomP = as.vector(GlomP[5, ])
ln_glomFDR = p.adjust(ln_glomP, method = "BH")    


6.#MGN
mgn_glomP = as.vector(GlomP[6, ])
mgn_glomFDR = p.adjust(mgn_glomP, method = "BH")    


7.#DN
dn_glomP = as.vector(GlomP[7, ])
dn_glomFDR = p.adjust(dn_glomP, method = "BH")    

8.#HT
ht_glomP = as.vector(GlomP[8, ])
ht_glomFDR = p.adjust(ht_glomP, method = "BH")    

9.#RPGN
rpgn_glomP = as.vector(GlomP[9, ])
rpgn_glomFDR = p.adjust(rpgn_glomP, method = "BH")    


GlomFDR[1, ] = fsgs_glomFDR
GlomFDR[2, ] = fsgs_mcd_glomFDR
GlomFDR[3, ] = mcd_glomFDR
GlomFDR[4, ] = igan_glomFDR
GlomFDR[5, ] = ln_glomFDR
GlomFDR[6, ] = mgn_glomFDR
GlomFDR[7, ] = dn_glomFDR
GlomFDR[8, ] = ht_glomFDR
GlomFDR[9, ] = rpgn_glomFDR

rownames(GlomFDR) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN")
colnames(GlomFDR) = c("EGFR", "Hypoxia", "JAK-STAT", "MAPK", "NFKb", "PI3K", "TGFb",
                      "TNFa", "Trail", "VEGF", "p53")

colnames(glom_prog11) = c("EGFR", "Hypoxia", "JAK-STAT", "MAPK", "NFKb", "PI3K", "TGFb",
                      "TNFa", "Trail", "VEGF", "p53")


save(GlomFDR, glom_prog11, file = "GlomFDR.RData")


Part 11. #FRD q-value representation by asterisks --> 



fdr_star_glom = GlomFDR
fdr_star_glom[which(fdr_star_glom >= 0.05)] = NA


fdr_star_glom = mystars_fdr_glom <- ifelse(fdr_star_glom < .001, " **** ", 
                                           ifelse(fdr_star_glom < .005, " *** ",
                                                  ifelse(fdr_star_glom < .01, " ** " ,
                                                         ifelse(fdr_star_glom < .05, " * ",
                                                                 " ")))) 
                                                                 
################Figure4A###########################################################################################                                                                 
library(gplots)
heatmap.2(glom_prog11,
          cellnote = fdr_star_glom,  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=200), 
          scale = "col",
          dendrogram="none",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA",
          Rowv = "NA")
#################################################################################################################################



load("ckd_glom_sepMat.RData") #created in Process_All_Glom.r   Part --> #########Separating the Conditions into distinct objects###################

ckd_glom = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,
                      DN_Glom, HT_Glom, RPGN_Glom, TN_Glom)

save(ckd_glom, file = "ckd_glom.RData")


#creating a scaffold matrix that will later be filled with the average of each gene/disease
ckdglom_gex = matrix(9, ncol = 6289, nrow = 10) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {ckdglom_gex[i,] = rowMeans(ckd_glom[[i]])
}


rownames(ckdglom_gex) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN")
colnames(ckdglom_gex) = rownames(FSGS_Glom)

ckdglom_gex = t(ckdglom_gex)

save(ckdglom_gex, file = "ckdglom_gex.RData")

#####TF Activity Estimation by DoROthEA###Figure3####################################################################################

#Feeding cyclic loess norm values into DoROthEA ###the one I USED!!!!CKD PAPER!!#####

#needed for this script:
#lib_enrichment_scores.r -> https://github.com/saezlab/DoRothEA/src/lib_enrichment_scores.r  or https://drive.google.com/drive/folders/1mBkawIQ0zxbeepJLXyZvvMGJpiDdZ0M2 (lib_enrichment_scores.R)

#geneset2.RData ->  https://drive.google.com/drive/folders/1mBkawIQ0zxbeepJLXyZvvMGJpiDdZ0M2 (geneset2.RData)

#Data needed can be found:  

source("~/Desktop/KD/Luz/lib_enrichment_scores.r")
library(gplots)
load("~/Desktop/KD/Luz/geneset2.RData")

aa_norm = gene_expression_statistic(ckdglom_gex, method = "scale", rnaseq = FALSE)


Glom_gex = SLEA(aa_norm, geneset2, method = "VIPER", filter_E = T)
Glom_gex_Factors = names(Glom_gex$regulons)
which_in_matrixgex = which((rownames(aa_norm)) %in% Glom_gex_Factors)
Glom_gex_transcription_factors = aa_norm[which_in_matrixgex,]

save(Glom_gex, file = "Glom_gex.RData")
save(Glom_gex_transcription_factors, file = "Glom_gex_transcription_factors.RData" )

fwhichglom = which(rownames(Glom_gex$ES) %in% rownames(Glom_gex_transcription_factors) )


Glom_gex_transcription_factors = Glom_gex_transcription_factors[order(rownames(Glom_gex_transcription_factors)), ]
Glom_gex$ES = Glom_gex$ES[order(rownames(Glom_gex$ES)), ]
Glom_gex$FDR = Glom_gex$FDR[order(rownames(Glom_gex$FDR)), ]

##Correlation of TF Activity Scores and Factor Expression################

corr = rep(9, ncol(Glom_gex_transcription_factors))
for (i in 1:ncol(Glom_gex_transcription_factors)) {
  temp = Glom_gex$ES
  corr[i] = cor(Glom_gex_transcription_factors[,i], temp[fwhichglom, i], method = "spearman")
}

corr2 = rep(9, nrow(Glom_gex_transcription_factors))
for(i in 1:length(corr2)) {
  temp = Glom_gex$ES
  temp = temp[fwhichglom,]
  corr2[i] = cor(Glom_gex_transcription_factors[i,], temp[i, ], method = "spearman")
}

save(corr, corr2, file = "correlate_TF_Expression.RData")

temp = Glom_gex$FDR
finalwhich2 = apply(temp[fwhichglom, ],1, function(x) any(x < 0.05))
finalwhich2 = which(finalwhich2 == TRUE)


tf_activity = Glom_gex$ES
tf_activity = tf_activity[fwhichglom,]
tf_activity = tf_activity[finalwhich2,]


correlate = round(corr2[finalwhich2], 3)
save(correlate, file = "factor_expression_corr.RData")


fdr = Glom_gex$FDR
fdr = fdr[fwhichglom,]
fdr = fdr[finalwhich2,]


tf_activity_print = tf_activity
tf_activity_print[which(fdr >= 0.05)] = NA

final_names = rep("w", length(rownames(tf_activity)))
for (i in 1:length(final_names)) {
  final_names[i] = paste0((rownames(tf_activity))[i], " / ", correlate[i])
}

final_col_names = rep("w", ncol(tf_activity))
for (i in 1:length(final_col_names)) {
  final_col_names[i] = paste0((colnames(tf_activity))[i])
}

#rownames(tf_activity) = final_names


#FRD q-value representation by astrisks --> 

Glomerulus

fdr_star_glom = fdr

fdr_star_glom[which(fdr_star_glom >= 0.05)] = NA


fdr_star_glom = mystars_fdr_glom <- ifelse(fdr_star_glom < .001, " **** ", 
                                           ifelse(fdr_star_glom < .005, " *** ",
                                                  ifelse(fdr_star_glom < .01, " ** " ,
                                                         ifelse(fdr_star_glom < .05, " * ",
                                                                " "))))                                                    



heatmap.2(tf_activity,
          cellnote = fdr_star_glom,  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,10),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=200), 
          scale = "none",
          dendrogram="none",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA",
          Rowv = "NA",
          labRow = final_names, labCol = final_col_names, cexRow = 1, cexCol = 1.5)






pdf("tf_activity_glom_gex_scale.pdf", height = 12, width = 10)
heatmap.2(tf_activity, cellnote = fdr_star_glom, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "hotpink"))(n=100), Colv = NA, Rowv = NA, dendrogram = "none", cellnote = round(tf_activity_print,2), notecol = "black", density.info = "none", notecex = 1, margins=c(12,10), labRow = final_names, labCol = final_col_names, cexRow = 1, cexCol = 1.5)
dev.off()


#Checking whether the factors present in tf_activity_effs_print are present in the gene expression data:

#overlap btw tf_glom and ckd_glom_g (gene expression values represented by  scaled Hedges'g effect size values)

list_tf = list(tf_activity, aa_norm)
tf_id = Reduce(intersect, lapply(list_tf, rownames))
ckd_glom_gex_tf = aa_norm[tf_id, ]
save(ckd_glom_gex_tf, file = "ckd_glom_gex_tf.RData")


pdf("tf_gex_tf.pdf")
heatmap.2(ckd_glom_gex_tf,
          cellnote = fdr_star_glom,  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=100), 
          dendrogram="none",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA", Rowv = "NA")
           # turn off column
dev.off()

############Not included######################################
#Bar Plot of the Correlation Coefficients

corr_TF_glom = matrix(as.numeric(correlate), ncol = 9, nrow = 1)
rownames(corr_TF_glom) = c("Spearman Correlation Coefficient")
colnames(corr_TF_glom)= c("FOS", "FOXA2", "FOXM1", "GATA4", "IRF1", "REST", "STAT1", "USF2", "VDR")

barplot(corr_TF_glom[1, ], col = c("chartreuse", "chocolate1", "cornflowerblue",
                              "deeppink1", "darkgreen", "khaki1", 
                              "darkorchid1", "darkgoldenrod1", "thistle1"), ylab = "Spearman Correlation Coefficient", ylim = c(-1,1)


#Hierarchical Clustering of KD Types based on gene expression Spearman's Corr coeff. of significant TFS
heatmap.2(cor(ckd_glom_gex_tf, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none")

#######################################################################
###################################################################################################################










