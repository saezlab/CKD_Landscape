rm(list = ls())
library(scater)
library(scran)
library(biomaRt)
library(LSD)
library(limma)
library(gplots)


setwd("~/Desktop/KD/Platform_Hierarchy")
load("~/Desktop/KD/Platform_Hierarchy/Process_All_Glom/process_glom_all.RData")
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


#overlap with DorotheA output

lll = list(ckd_glom_logFC, tf_activity_effs_print)
id_lll = Reduce(intersect,lapply(lll, rownames))
ckd_glom_logFC_D = ckd_glom_logFC[id_lll, ]

pdf("tf_gex_eff.pdf")
heatmap.2(ckd_glom_logFC_D,
          cellnote = fdr_star_glom,  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=200), 
          scale = "row",
          dendrogram="row",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA") # turn off column
dev.off()
######################################################################################################
