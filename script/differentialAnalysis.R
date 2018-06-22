
# install necessary packages if not present yet 
list.of.packages <- c("scater", "scran","biomaRt", "LSD","limma", "gplots","MetaDE")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


path = "/Users/francescoceccarelli/Desktop/CKD Landscape code/"
load(paste0(path, "Data/process_glom_all.RData"))
dir.create(paste0(path, "RDS_folder"))

#processed_glom_all$info -> meta data
#processed_glom_all$final -> expression data

info = processed_glom_all$info
save(info, file = paste0(path, "Data/Glom_All_Info.RData"))

exp = processed_glom_all$final


#For this analysis, we used the version from which the healthy living donor samples were removed and the tumor nephrectomy samples were left in the data set.
#Remove all the Healthy Living Donor Samples, because a lot of genes have been removed from this cohort.
info2 <- info[!(grepl(pattern = "Healthy Living", x = info$Group)), ]
dim(info2)
sampels_toUse <- info2$Accession
exp2 <- exp[, sampels_toUse]
dim(exp2)
exp2 <- na.omit(exp2)
exp <- exp2
info <- info2
dim(exp)

rownames(exp) <- toupper(rownames(exp))

save(exp, info, file = paste0(path, "Data/exp_info.RData"))


a = exp
a = a[order(rownames(a)),]
a = a[,order(colnames(a))]
info_glom = info[order(info$Accession),]
rownames(info_glom) = info_glom$Accession


i = data.frame(Cell = colnames(exp), info_glom)
rownames(i) = colnames(a)
infoFrame = new("AnnotatedDataFrame", data = info_glom)
sc_glom = newSCESet(countData = a, phenoData = infoFrame)
sc_glom = calculateQCMetrics(sc_glom)

save(sc_glom, file = paste0(path, "Data/sc_glom.RData"))
save(info_glom, infoFrame, file = paste0(path, "Data/info_glom.RData"))


#PCA by Study after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byStudy_Glom_noHLD.pdf"))
plotPCA(sc_glom, ncomponents = 3, colour_by = "Study", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by Group after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byGroup_Glom_noHLD.pdf"))
plotPCA(sc_glom, ncomponents = 3, colour_by = "Group", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by platform after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byPlatform_Glom_noHLD.pdf"))
plotPCA(sc_glom, ncomponents = 3, colour_by = "Platform", exprs_values = "counts", ntop = 6289)
dev.off()


#group-specific plots
#1.1.PCA of gene expression of Tumor Nephrectomy separated by Study
pdf(paste0(path, "Plot/pca_TN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Tumor Nephrectomy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#1.2.PCA of gene expression of Tumor Nephrectomy separated by Study
pdf(paste0(path, "Plot/pca_TN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Tumor Nephrectomy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#2.1.PCA of gene expression of DN separated by Study
pdf(paste0(path, "Plot/pca_DN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Diabetic Nephropathy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#2.2.PCA of gene expression of DN separated by Study
pdf(paste0(path, "Plot/pca_TN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Diabetic Nephropathy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#3.1.PCA of gene expression of FSGS separated by Study
pdf(paste0(path, "Plot/pca_FSGS_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "FSGS"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#3.2.PCA of gene expression of FSGS separated by Study
pdf(paste0(path, "Plot/pca_FSGS_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "FSGS"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#4.1.PCA of gene expression of FSGS_MCD separated by Study
pdf(paste0(path, "Plot/pca_FSGS_MCD_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "FSGS_MCD"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#4.2.PCA of gene expression of FSGS_MCD separated by Study
pdf(paste0(path, "Plot/pca_FSGS_MCD_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "FSGS_MCD"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#5.1.PCA of gene expression of HT separated by Study
pdf(paste0(path, "Plot/pca_HT_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Hypertensive Nephropathy"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#5.2.PCA of gene expression of HT separated by Study
pdf(paste0(path, "Plot/pca_HT_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Hypertensive Nephropathy"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#6.1.PCA of gene expression of IgAN separated by Study
pdf(paste0(path, "Plot/pca_IgAN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "IgAN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#6.2.PCA of gene expression of IgAN separated by Study
pdf(paste0(path, "Plot/pca_IgAN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "IgAN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#7.1.PCA of gene expression of LN separated by Study
pdf(paste0(path, "Plot/pca_LN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Lupus Nephritis"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#7.2.PCA of gene expression of LN separated by Study
pdf(paste0(path, "Plot/pca_TN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "Lupus Nephritis"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#8.1.PCA of gene expression of MCD separated by Study
pdf(paste0(path, "Plot/pca_MCD_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "MCD"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#8.2.PCA of gene expression of MCD separated by Study
pdf(paste0(path, "Plot/pca_MCD_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "MCD"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#9.1.PCA of gene expression of MGN separated by Study
pdf(paste0(path, "Plot/pca_MGN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "MGN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#9.2.PCA of gene expression of MGN separated by Study
pdf(paste0(path, "Plot/pca_MGN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "MGN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()

#10.1.PCA of gene expression of RPGN separated by Study
pdf(paste0(path, "Plot/pca_RPGN_glom_byStudy.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "RPGN"], ncomponents = 3, colour_by = "Study", exprs_value = "counts")
dev.off()

#10.2.PCA of gene expression of RPGN separated by Study
pdf(paste0(path, "Plot/pca_RPGN_glom_byPlatform.pdf"))
plotPCA(sc_glom[, infoFrame$Group == "RPGN"], ncomponents = 3, colour_by = "Platform", exprs_value = "counts")
dev.off()


#histogram before norm

pdf(paste0(path, "Plot/hist_glom_all_nonorm.pdf"))
hist(counts(sc_glom))
dev.off()

pdf(paste0(path, "Plot/boxplot_glom_all_nonorm.pdf"))
boxplot(counts(sc_glom))
dev.off()

pdf(paste0(path, "Plot/cor_heat_glom_all_nonorm.pdf"))
heatmap.2(cor(counts(sc_glom)), col = colorRampPalette(c("blue","white","red"))(n=200),trace = "none" )
dev.off()
#################################

#normalization with Cyclic Loess
norm_glom_cl = normalizeBetweenArrays(counts(sc_glom), method = "cyclicloess")
save(norm_glom_cl, file = paste0(path, "Data/norm_glom_cl.RData"))

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

save(sc_glom_cl, infoFrame, file = paste0(path, "Data/sc_glom_cl.RData"))



#PCA by Study after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byStudy_Glom_noHLD_cloess.pdf"))
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Study", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by Group after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byGroup_Glom_noHLD_cloess.pdf"))
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Group", exprs_values = "counts", ntop = 6289 )
dev.off()

#PCA by platform after removal of the platform batch-effect inducing genes
pdf(paste0(path, "Plot/pca_byPlatform_Glom_noHLD_cloess.pdf"))
plotPCA(sc_glom_cl, ncomponents = 3, colour_by = "Platform", exprs_values = "counts", ntop = 6289)
dev.off()

#histogram before norm> dev.off()

#hist(counts(sc_glom_cl))

pdf(paste0(path, "Plot/hist_glom_cloess.pdf"))
hist(counts(sc_glom_cl))
dev.off()

pdf(paste0(path, "Plot/boxplot_glom_cloess.pdf"))
boxplot(counts(sc_glom_cl))
dev.off()

pdf(paste0(path, "Plot/cor_heat_glom_cloess.pdf"))
heatmap.2(cor(counts(sc_glom_cl)), col = colorRampPalette(c("blue","white","hotpink"))(n=200),trace = "none" )
dev.off()


####Separate the conditions into distinct objects


TN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Tumor Nephrectomy"])


DN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Diabetic Nephropathy"])

dn_tn = cbind(DN_Glom, TN_Glom)

FSGS_Glom = counts(sc_glom_cl[, infoFrame$Group == "FSGS"])

fsgs_tn = cbind(FSGS_Glom, TN_Glom)


FSGS_MCD_Glom = counts(sc_glom_cl[, infoFrame$Group == "FSGS_MCD"])

fsgs_mcd_tn = cbind(FSGS_MCD_Glom, TN_Glom)

fsgs_fsgsmcd = cbind(FSGS_Glom, FSGS_MCD_Glom)


HT_Glom = counts(sc_glom_cl[, infoFrame$Group == "Hypertensive Nephropathy"])

ht_tn =  cbind(HT_Glom, TN_Glom)


IgAN_Glom = counts(sc_glom_cl[, infoFrame$Group == "IgAN"])

igan_tn =  cbind(IgAN_Glom, TN_Glom)


LN_Glom = counts(sc_glom_cl[, infoFrame$Group == "Lupus Nephritis"])

ln_tn =  cbind(LN_Glom, TN_Glom)


MCD_Glom = counts(sc_glom_cl[, infoFrame$Group == "MCD"])

mcd_tn =  cbind(MCD_Glom, TN_Glom)

mcd_fsgsmcd = cbind(MCD_Glom, FSGS_MCD_Glom)

mcd_fsgs = cbind(FSGS_Glom, MCD_Glom)


MGN_Glom = counts(sc_glom_cl[, infoFrame$Group == "MGN"])

mgn_tn =  cbind(MGN_Glom, TN_Glom)


RPGN_Glom = counts(sc_glom_cl[, infoFrame$Group == "RPGN"])

rpgn_tn =  cbind(RPGN_Glom, TN_Glom)




save(DN_Glom, FSGS_Glom, FSGS_MCD_Glom, HT_Glom, IgAN_Glom, LN_Glom, MCD_Glom, MGN_Glom, RPGN_Glom, TN_Glom, file = paste0(path, "Data/ckd_glom_sepMat.RData"))

save(dn_tn, fsgs_tn, fsgs_mcd_tn, ht_tn, igan_tn, ln_tn, mcd_tn, mgn_tn, rpgn_tn, fsgs_fsgsmcd, mcd_fsgs, mcd_fsgsmcd,  file = paste0(path, "Data/combined_disease_control.RData"))



#Hierarchical Clustering of Normalized Gene Expression of CKD Entities

load(paste0(path, "Data/ckd_glom_sepMat.RData"))

source(paste0(path, "Data/lib_enrichment_scores.r"))

load(paste0(path, "Data/geneset2.RData"))


gg = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,DN_Glom, HT_Glom, RPGN_Glom, TN_Glom)

#creating a scaffold matrix that will later be filled with the average of each gene/disease
glom = matrix(9, ncol = 6289, nrow = 10) 


# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {glom[i,] = rowMeans(gg[[i]])
}


rownames(glom) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN")
colnames(glom) = rownames(FSGS_Glom)

glom = t(glom)


#cyclic loess normalized, average gex.
save(glom, file = paste0(path, "Data/glom.RData"))

#Scaling and Recentering gene expression using R scale method (Luz)
glom_scaled = gene_expression_statistic(glom, method = "scale", rnaseq = FALSE)

heatmap.2(cor(glom_scaled, method = "spearman"), col= colorRampPalette(c("blue","white","hotpink"))(n=200), trace = "none",
          scale = "none", notecol = "black", density.info = "none", notecex = 1,
          margins=c(12,10), cexCol = 1.5)

#####################################################################


#Limma analysis for differential expression between Disease and Tumor Nephrectomy Samples

#1. Diabetic Nephropathy vs Tumor Nephrectomy

#intersection


#get the fold changes; the data was already log transformed (RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_dn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Diabetic Nephropathy")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_dn, file = "Data/DN_Glom_TN_FC.RData")
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_cloess_DNglom.pdf"))
hist(fc_dn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoDN = infoFrame[infoFrame$Group == "Diabetic Nephropathy"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

idntn = combine(infoDN, infoTN)
save(idntn, file = paste0(path, "Data/idntn_info_DN_TN.RData"))

diabetic_vs_tn = data.frame(idntn$Accession, idntn$Study, idntn$Group, idntn$Platform, idntn$Tissue)
colnames(diabetic_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(diabetic_vs_tn, file = paste0(path, "Data/diabetic_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(diabetic_vs_tn, file = paste0(path, "Data/diabetic_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(dn_tn, file = paste0(path, "Data/dn_tn.csv"))


#basic limma analysis for differential expression
design_dn = cbind(rep(1,length(idntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Diabetic Nephropathy"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_dn = lmFit(dn_tn, design_dn)
fit2_dn = eBayes(fit_dn)
top_dn = topTable(fit2_dn, adjust="BH", number = 50000000, coef = 2)

save(top_dn, fit2_dn,  file = paste0(path, "Data/DN_TN_Glom_limma.RData"))

saveRDS(top_dn, paste0(path, "RDS_folder/top_dn.rds"))

write.csv2(top_dn, file = paste0(path, "Data/top_dn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = dn_tn
ww = which(rownames(aa) %in% rownames(top_dn) [(which(top_dn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_DN_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

DNTN_glom_decide = decideTests(fit2_dn, method = "separate", adjust.method = "BH", p.value = 0.05)
save(DNTN_glom_decide, file = paste0(path, "Data/DNTN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_DN_Glom.pdf"))
limma::plotMA(fit2_dn, status = DNTN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Diabetic Nephropathy and Tumor Nephrectomy")
dev.off()



#2. FSGS vs Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_fsgs, file = paste0(path, "Data/FSGS_Glom_TN_FC.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_FSGS_TN_glom.pdf"))
hist(fc_dn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS = infoFrame[infoFrame$Group == "FSGS"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgstn = combine(infoFSGS, infoTN)
save(ifsgstn, file = paste0(path, "Data/ifsgstn_info_FSGS_TN.RData"))

fsgs_vs_tn = data.frame(ifsgstn$Accession, ifsgstn$Study, ifsgstn$Group, ifsgstn$Platform, ifsgstn$Tissue)
colnames(fsgs_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(fsgs_vs_tn, file = paste0(path, "Data/fsgs_vs_tn_info.RData"))
#info data with FSGS & TN for Hyojin
write.csv2(fsgs_vs_tn, file = paste0(path, "Data/fsgs_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(fsgs_tn, file = paste0(path, "Data/fsgs_tn.csv"))


#basic limma analysis for differential expression
design_fsgs = cbind(rep(1,length(ifsgstn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "FSGS"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_fsgs = lmFit(fsgs_tn, design_fsgs)
fit2_fsgs = eBayes(fit_fsgs)
top_fsgs = topTable(fit2_fsgs, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs, fit2_fsgs,  file = paste0(path, "Data/FSGS_TN_Glom_limma.RData"))
saveRDS(top_fsgs, paste0(path, "RDS_folder/top_fsgs.rds"))


write.csv2(top_fsgs, file = paste0(path, "Data/top_fsgs.csv"))


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = fsgs_tn
ww = which(rownames(aa) %in% rownames(top_fsgs) [(which(top_fsgs$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_FSGS_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

FSGSTN_glom_decide = decideTests(fit2_fsgs, method = "separate", adjust.method = "BH", p.value = 0.05)
save(FSGSTN_glom_decide, file = paste0(path, "Data/FSGSTN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_FSGS_TN_Glom.pdf"))
limma::plotMA(fit2_fsgs, status = FSGSTN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Focal Segmental Glomerulosclerosis and Tumor Nephrectomy")
dev.off()


#3. FSGS_MCD vs Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs_mcd = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_fsgs_mcd, file = paste0(path, "Data/FSGS_MCD_Glom_TN_FC.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_FSGS_MCD_TN_glom.pdf"))
hist(fc_fsgs_mcd, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS_MCD = infoFrame[infoFrame$Group == "FSGS_MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgsmcdtn = combine(infoFSGS_MCD, infoTN)
save(ifsgsmcdtn, file = paste0(path, "Data/ifsgsmcdtn_info_FSGS_TN.RData"))


fsgs_mcd_vs_tn = data.frame(ifsgsmcdtn$Accession, ifsgsmcdtn$Study, ifsgsmcdtn$Group, ifsgsmcdtn$Platform, ifsgsmcdtn$Tissue)
colnames(fsgs_mcd_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(fsgs_mcd_vs_tn, file = paste0(path, "Data/fsgs_mcd_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(fsgs_mcd_vs_tn, file = paste0(path, "Data/fsgs_mcd_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(fsgs_mcd_tn, file = paste0(path, "Data/fsgs_mcd_tn.csv"))


#basic limma analysis for differential expression
design_fsgs_mcd = cbind(rep(1,length(ifsgsmcdtn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_fsgs_mcd = lmFit(fsgs_mcd_tn, design_fsgs_mcd)
fit2_fsgs_mcd = eBayes(fit_fsgs_mcd)
top_fsgs_mcd = topTable(fit2_fsgs_mcd, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs_mcd, fit2_fsgs_mcd,  file = paste0(path, "Data/FSGS_MCD_TN_Glom_limma.RData"))
saveRDS(top_fsgs_mcd, paste0(path, "RDS_folder/top_fsgs_mcd.rds"))


write.csv2(top_fsgs_mcd, file = paste0(path, "Data/top_fsgs_mcd.csv"))


#select genes that are differentially expressed with an adjusted p value of 0.05
aa = fsgs_mcd_tn
ww = which(rownames(aa) %in% rownames(top_fsgs_mcd) [(which(top_fsgs_mcd$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_FSGS_MCD_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

FSGS_MCD_TN_glom_decide = decideTests(fit2_fsgs_mcd, method = "separate", adjust.method = "BH", p.value = 0.05)
save(FSGS_MCD_TN_glom_decide, file = paste0(path, "Data/FSGS_MCD_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_FSGS_MCD_TN_Glom.pdf"))
limma::plotMA(fit2_fsgs_mcd, status = FSGS_MCD_TN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Focal Segmental Glomerulosclerosis + MCD and Tumor Nephrectomy")
dev.off()



#3. FSGS vs FSGS_MCD
#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_fsgs_mcd01 = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) ) 

save(fc_fsgs_mcd01, file = paste0(path, "Data/FSGS_MCD01_GlomFC.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_FSGS_MCD01_glom.pdf"))
hist(fc_fsgs_mcd01, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoFSGS_MCD = infoFrame[infoFrame$Group == "FSGS_MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ifsgsmcd = combine(infoFSGS_MCD, infoFSGS)
save(ifsgsmcd, file = paste0(path, "Data/info_FSGSMCD.RData"))


#basic limma analysis for differential expression
design_fsgsmcd01 = cbind(rep(1,length(ifsgsmcd$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "FSGS"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"]))))
fit_fsgs_mcd01 = lmFit(fsgs_fsgsmcd, design_fsgsmcd01)
fit2_fsgs_fsgsmcd = eBayes(fit_fsgs_mcd01)
top_fsgs_fsgsmcd = topTable(fit2_fsgs_fsgsmcd, adjust="BH", number = 50000000, coef = 2)

save(top_fsgs_fsgsmcd, fit2_fsgs_fsgsmcd,  file = paste0(path, "Data/FSGS_FSGSMCD_Glom_limma.RData"))

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = fsgs_fsgsmcd
ww = which(rownames(aa) %in% rownames(top_fsgs_fsgsmcd) [(which(top_fsgs_fsgsmcd$adj.P.Val < 0.05))])

#no diff. exp. genes between FSGS and FSGS_MCD


#3. MCD vs  FSGS_MCD

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_fsgsmcd = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS_MCD")])) ) 

save(fc_mcd_fsgsmcd, file = paste0(path, "Data/MCD_FSGSMCD_FC.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_MCD_FSGSMCD_glom.pdf"))
hist(fc_mcd_fsgsmcd, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

imcd_fsgsmcd = combine(infoMCD, infoFSGS_MCD)
save(imcd_fsgsmcd, file = paste0(path, "Data/imcd_fsgsmcd.RData"))


#basic limma analysis for differential expression
design_mcd_fsgsmcd = cbind(rep(1,length(imcd_fsgsmcd$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS_MCD"]))))
fit_mcd_fsgsmcd = lmFit(mcd_fsgsmcd, design_mcd_fsgsmcd)
fit2_mcd_fsgsmcd = eBayes(fit_mcd_fsgsmcd)
top_mcd_fsgsmcd = topTable(fit2_mcd_fsgsmcd, adjust="BH", number = 50000000, coef = 2)

save(top_mcd_fsgsmcd, fit2_mcd_fsgsmcd,  file = paste0(path, "Data/MCD_FSSGSMCD_Glom_limma.RData"))

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_fsgsmcd
ww = which(rownames(aa) %in% rownames(top_mcd_fsgsmcd) [(which(top_mcd_fsgsmcd$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_MCD_FSGSMCD_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_FSGSMCD_glom_decide = decideTests(fit2_mcd_fsgsmcd, method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_FSGSMCD_glom_decide, file = paste0(path, "Data/MCD_FSGSMCD_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_MCD_FSGSMCD_Glom.pdf"))
limma::plotMA(fit2_mcd_fsgsmcd, status = MCD_FSGSMCD_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Minimal Change & FSGS with Minimal Change")
dev.off()



#3. MCD vs  FSGS

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_fsgs02 = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "FSGS")])) ) 

save(fc_mcd_fsgs02, file = paste0(path, "Data/fc_mcd_fsgs02.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_MCD_FSGS02_glom.pdf"))
hist(fc_mcd_fsgs02, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

imcdfsgs = combine(infoMCD, infoFSGS)
save(imcdfsgs, file = paste0(path, "Data/imcdfsgs.RData"))


#basic limma analysis for differential expression
design_imcdfsgs = cbind(rep(1,length(imcdfsgs$Group)), c(rep(0,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(1, length(infoFrame$Group[infoFrame$Group == "FSGS"]))))
fit_imcdfsgs = lmFit(mcd_fsgs, design_imcdfsgs)
fit2_imcdfsgs = eBayes(fit_imcdfsgs)
top_imcdfsgs = topTable(fit2_imcdfsgs, adjust="BH", number = 50000000, coef = 2)

save(top_imcdfsgs , fit2_imcdfsgs,  file = paste0(path, "Data/MCD_FSGS_Glom_limma.RData"))

#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_fsgs
ww = which(rownames(aa) %in% rownames(top_imcdfsgs) [(which(top_imcdfsgs$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_MCD_FSGS_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","black","yellow"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_FSGS_glom_decide = decideTests(fit2_imcdfsgs , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_FSGS_glom_decide, file = paste0(path, "Data/MCD_FSGS_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_MCD_FSGS_Glom.pdf"))
limma::plotMA(fit2_imcdfsgs, status = MCD_FSGS_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Minimal Change & Focal Segmental Glomerulosclerosis")
dev.off()



#3. Hypertensive Nephropathy vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_ht_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Hypertensive Nephropathy")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_ht_tn, file = paste0(path, "Data/fc_ht_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_HT_TN_glom.pdf"))
hist(fc_ht_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoHT = infoFrame[infoFrame$Group == "Hypertensive Nephropathy"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

ihttn = combine(infoHT, infoTN)
save(ihttn, file = paste0(path, "Data/ihttn.RData"))


ht_vs_tn = data.frame(ihttn$Accession, ihttn$Study, ihttn$Group, ihttn$Platform, ihttn$Tissue)
colnames(ht_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(ht_vs_tn, file = paste0(path, "Data/ht_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(ht_vs_tn, file = paste0(path, "Data/ht_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(ht_tn, file = paste0(path, "Data/ht_tn.csv"))




#basic limma analysis for differential expression
design_ihttn = cbind(rep(1,length(ihttn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Hypertensive Nephropathy"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_ihttn = lmFit(ht_tn, design_ihttn)
fit2_ihttn = eBayes(fit_ihttn)
top_ihttn = topTable(fit2_ihttn, adjust="BH", number = 50000000, coef = 2)

save(top_ihttn , fit2_ihttn,  file = paste0(path, "Data/ihttn_Glom_limma.RData"))
saveRDS(top_ihttn, paste0(path, "RDS_folder/top_ht.rds"))


write.csv2(top_ihttn, file = paste0(path, "Data/top_ihttn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = ht_tn
ww = which(rownames(aa) %in% rownames(top_ihttn) [(which(top_ihttn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_HT_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("cyan3","black","lightsalmon1"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

HT_TN_glom_decide = decideTests(fit2_ihttn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(HT_TN_glom_decide, file = paste0(path, "Data/HT_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_HT_TN_Glom.pdf"))
limma::plotMA(fit2_ihttn, status = HT_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Hypertensive Nephropathy & Tumor Nephrectomy")
dev.off()


#3. IgAN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_igan_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "IgAN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_igan_tn, file = paste0(path, "Data/fc_igan_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_IgAN_TN_glom.pdf"))
hist(fc_igan_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoIgAN = infoFrame[infoFrame$Group == "IgAN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

igantn = combine(infoIgAN, infoTN)
save(igantn, file = paste0(path, "Data/igantn.RData"))


igan_vs_tn = data.frame(igantn$Accession, igantn$Study, igantn$Group, igantn$Platform, igantn$Tissue)
colnames(igan_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(igan_vs_tn, file = paste0(path, "Data/igan_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(igan_vs_tn, file = paste0(path, "Data/igan_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(igan_tn, file = paste0(path, "Data/igan_tn.csv"))



#basic limma analysis for differential expression
design_igantn = cbind(rep(1,length(igantn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "IgAN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_igantn = lmFit(igan_tn, design_igantn)
fit2_igantn = eBayes(fit_igantn)
top_igantn = topTable(fit2_igantn, adjust="BH", number = 50000000, coef = 2)

save(top_igantn , fit2_igantn,  file = paste0(path, "Data/igantn_Glom_limma.RData"))
saveRDS(top_igantn, paste0(path, "RDS_folder/top_igan.rds"))


write.csv2(top_igantn, file = paste0(path, "Data/top_igantn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = igan_tn
ww = which(rownames(aa) %in% rownames(top_igantn) [(which(top_igantn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_IgAN_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

IgAN_TN_glom_decide = decideTests(fit2_igantn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(IgAN_TN_glom_decide, file = paste0(path, "Data/IgAN_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_IgAN_TN_Glom.pdf"))
limma::plotMA(fit2_igantn, status = IgAN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between IgA Nephropathy & Tumor Nephrectomy")
dev.off()


#3. LN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_ln_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Lupus Nephritis")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_ln_tn, file = paste0(path, "Data/fc_ln_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Data/foldChange_byGroup_LN_TN_glom.pdf"))
hist(fc_ln_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoLN = infoFrame[infoFrame$Group == "Lupus Nephritis"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

lntn = combine(infoLN, infoTN)
save(lntn, file = paste0(path, "Data/lntn.RData"))

ln_vs_tn = data.frame(lntn$Accession, lntn$Study, lntn$Group, lntn$Platform, lntn$Tissue)
colnames(ln_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(ln_vs_tn, file = paste0(path, "Data/ln_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(ln_vs_tn, file = paste0(path, "Data/ln_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(ln_tn, file = paste0(path, "Data/ln_tn.csv"))



#basic limma analysis for differential expression
design_lntn = cbind(rep(1,length(lntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "Lupus Nephritis"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_lntn = lmFit(ln_tn, design_lntn)
fit2_lntn = eBayes(fit_lntn)
top_lntn = topTable(fit2_lntn, adjust="BH", number = 50000000, coef = 2)

save(top_lntn , fit2_lntn,  file = paste0(path, "Data/lntn_Glom_limma.RData"))

saveRDS(top_lntn, paste0(path, "RDS_folder/top_ln.rds"))


write.csv2(top_lntn, file = paste0(path, "Data/top_lntn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = ln_tn
ww = which(rownames(aa) %in% rownames(top_lntn) [(which(top_lntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_LN_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

LN_TN_glom_decide = decideTests(fit2_lntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(LN_TN_glom_decide, file = paste0(path, "Data/LN_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_LN_TN_Glom.pdf"))
    
limma::plotMA(fit2_lntn, status = LN_TN_glom_decide[,2], main = "DEGs between Lupus Nephritis & Tumor Nephrectomy")
dev.off()


#3. MCD vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mcd_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MCD")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_mcd_tn, file = paste0(path, "Data/fc_mcd_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_MCD_TN_glom.pdf"))
hist(fc_mcd_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMCD = infoFrame[infoFrame$Group == "MCD"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

mcdtn = combine(infoMCD, infoTN)
save(mcdtn, file = paste0(path, "Data/mcdtn.RData"))

mcd_vs_tn = data.frame(mcdtn$Accession, mcdtn$Study, mcdtn$Group, mcdtn$Platform, mcdtn$Tissue)
colnames(mcd_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(mcd_vs_tn, file = paste0(path, "Data/mcd_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(mcd_vs_tn, file = paste0(path, "Data/mcd_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(mcd_tn, file = paste0(path, "Data/mcd_tn.csv"))



#basic limma analysis for differential expression
design_mcdtn = cbind(rep(1,length(mcdtn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "MCD"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_mcdtn = lmFit(mcd_tn, design_mcdtn)
fit2_mcdtn = eBayes(fit_mcdtn)
top_mcdtn = topTable(fit2_mcdtn, adjust="BH", number = 50000000, coef = 2)

save(top_mcdtn , fit2_mcdtn,  file = paste0(path, "Data/MCD_TN_Glom_limma.RData"))

 saveRDS(top_mcdtn, paste0(path, "RDS_folder/top_mcd.rds"))


write.csv2(top_mcdtn, file = paste0(path, "Data/top_mcdtn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mcd_tn
ww = which(rownames(aa) %in% rownames(top_mcdtn) [(which(top_mcdtn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_MCD_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MCD_TN_glom_decide = decideTests(fit2_mcdtn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MCD_TN_glom_decide, file = paste0(path, "Data/MCD_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_MCD_TN_Glom.pdf"))
limma::plotMA(fit2_mcdtn, status = MCD_TN_glom_decide[,2], col = c("blue", "hotpink"), main = "DEGs between Minimal Change & Tumor Nephrectomy")
dev.off()


#3. MGN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_mgn_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "MGN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_mgn_tn, file = paste0(path, "Data/fc_mgn_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Plot/foldChange_byGroup_MGN_TN_glom.pdf"))
hist(fc_mgn_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoMGN = infoFrame[infoFrame$Group == "MGN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

mgntn = combine(infoMGN, infoTN)
save(mgntn, file = paste0(path, "Data/mgntn.RData"))


mgn_vs_tn = data.frame(mgntn$Accession, mgntn$Study, mgntn$Group, mgntn$Platform, mgntn$Tissue)
colnames(mgn_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(mgn_vs_tn, file = paste0(path, "Data/mgn_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(mgn_vs_tn, file = paste0(path, "Data/mgn_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(mgn_tn, file = paste0(path, "Data/mgn_tn.csv"))




#basic limma analysis for differential expression
design_mgntn = cbind(rep(1,length(mgntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "MGN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_mgntn = lmFit(mgn_tn, design_mgntn)
fit2_mgntn = eBayes(fit_mgntn)
top_mgntn = topTable(fit2_mgntn, adjust="BH", number = 50000000, coef = 2)

save(top_mgntn , fit2_mgntn,  file = paste0(path, "Data/mgntn_Glom_limma.RData"))

saveRDS(top_mgntn, paste0(path, "RDS_folder/top_mgn.rds"))


write.csv2(top_mgntn, file = paste0(path, "Data/top_mgntn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = mgn_tn
ww = which(rownames(aa) %in% rownames(top_mgntn) [(which(top_mgntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_MGN_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

MGN_TN_glom_decide = decideTests(fit2_mgntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(MGN_TN_glom_decide, file = paste0(path, "Data/MGN_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_MGN_TN_Glom.pdf"))
limma::plotMA(fit2_mgntn, status = MGN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Membranous Nephropathy & Tumor Nephrectomy")
dev.off()


#3. RPGN vs  Tumor Nephrectomy

#intersection


#get the fold chagen..the data was already logged (by RMA package) so we just do a subtraction...you can make an MA plot or make a histogram, see below
fc_rpgn_tn = ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "RPGN")])) )  - ( rowMeans(counts(sc_glom_cl[,which(infoFrame$Group == "Tumor Nephrectomy")])) ) 

save(fc_rpgn_tn, file = paste0(path, "Data/fc_rpgn_tn.RData"))
#historgram of fold change values..should be centered on 0, otherwise, something is strange either technical (likely) or biological (unlikely)
pdf(paste0(path, "Data/foldChange_byGroup_RPGN_TN_glom.pdf"))
hist(fc_rpgn_tn, breaks = "scott")
lines(rep(0,10000), 0:9999, col = "red", lty = 2, lwd = 2)
dev.off()

#Combine Annotated Data Frame

infoRPGN = infoFrame[infoFrame$Group == "RPGN"]
infoTN = infoFrame[infoFrame$Group == "Tumor Nephrectomy"]

rpgntn = combine(infoRPGN, infoTN)
save(rpgntn, file = paste0(path, "Data/rpgntn.RData"))

rpgn_vs_tn = data.frame(rpgntn$Accession, rpgntn$Study, rpgntn$Group, rpgntn$Platform, rpgntn$Tissue)
colnames(rpgn_vs_tn) = c("Accession", "Study", "Group", "Platform", "Tissue")
save(rpgn_vs_tn, file = paste0(path, "Data/rpgn_vs_tn_info.RData"))
#info data with DN & TN for Hyojin
write.csv2(rpgn_vs_tn, file = paste0(path, "Data/rpgn_vs_tn.csv"))
#expression of DN and TN for Hyojin
write.csv2(rpgn_tn, file = paste0(path, "Data/rpgn_tn.csv"))



#basic limma analysis for differential expression
design_rpgntn = cbind(rep(1,length(rpgntn$Group)), c(rep(1,length(infoFrame$Group[infoFrame$Group == "RPGN"])), rep(0, length(infoFrame$Group[infoFrame$Group == "Tumor Nephrectomy"]))))
fit_rpgntn = lmFit(rpgn_tn, design_rpgntn)
fit2_rpgntn = eBayes(fit_rpgntn)
top_rpgntn = topTable(fit2_rpgntn, adjust="BH", number = 50000000, coef = 2)

save(top_rpgntn , fit2_rpgntn,  file = paste0(path, "Data/rpgntn_Glom_limma.RData"))

saveRDS(top_rpgntn, paste0(path, "RDS_folder/top_rpgn.rds"))


write.csv2(top_rpgntn, file = paste0(path, "Data/top_rpgntn.csv"))


#select genes that are differentially expressed with an adjusted p value of 0;05
aa = rpgn_tn
ww = which(rownames(aa) %in% rownames(top_rpgntn) [(which(top_rpgntn$adj.P.Val < 0.05))])

##heatmap of the expression of those genes across different samples, scaled by row. 
pdf(paste0(path, "Plot/heatmap_RPGN_TN_Glom_limma.pdf"))
heatmap.2(aa[ww,], trace = "none", col = colorRampPalette(c("blue","white","hotpink"))(n=1000), dendrogram="col", scale = "row", labRow = "")
dev.off()

RPGN_TN_glom_decide = decideTests(fit2_rpgntn , method = "separate", adjust.method = "BH", p.value = 0.05)
save(RPGN_TN_glom_decide, file = paste0(path, "Data/RPGN_TN_glom_decide.RData"))

pdf(paste0(path, "Plot/maplot_limma_RPGN_TN_Glom.pdf"))
limma::plotMA(fit2_rpgntn, status = RPGN_TN_glom_decide[,2], col = c("hotpink", "blue"), main = "DEGs between Rapidly Progressive Glomerulosclerosis & Tumor Nephrectomy")
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------
  
  
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

save(DN_p, file = paste0(path, "Data/DN_p.RData"))


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

save(FSGS_p, file = paste0(path, "Data/FSGS_p.RData"))


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

save(FSGSMCD_p, file = paste0(path, "Data/FSGSMCD_p.RData"))


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

save(MCD_p, file = paste0(path, "Data/MCD_p.RData")) 


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

save(HT_p, file = paste0(path, "Data/HT_p.RData"))


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

save(IgAN_p, file = paste0(path, "Data/IgAN_p.RData"))


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

save(LN_p, file = paste0(path, "Data/LN_p.RData"))


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

save(MGN_p, file = paste0(path, "Data/MGN_p.RData"))


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

save(RPGN_p, file = paste0(path, "Data/RPGN_p.RData"))



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

save(x, file = paste0(path, "Data/x.RData"))


##############################################################################################################################
#MetaDE: R package containig various meta-analytic approaches to aggregate p-values and effect sizes.
#Aim: Extraction of DE genes that have small p-values in ALL studies.



p = list()
p$p = x
p$bp = NULL

res = MetaDE::MetaDE.pvalue(p, meta.method = c("maxP"))
save(res, file = paste0(path, "Data/res.RData")) #1st MA result using maxP method


count.DEnumber(res, p.cut = c(0.01, 0.05), q.cut = c(0.01, 0.05))

w = which(rownames(res$ind.p) %in% rownames(x) [which(res$meta.analysis$FDR < 0.01)])

ckd_maxP_001 = x[w, ]    


k = which(rownames(res$ind.p) %in% rownames(x) [which(res$meta.analysis$FDR < 0.05)])

ckd_maxP_005 = x[k, ] 

ckd_glom_maxP_001 = rownames(ckd_maxP_001)

ckd_glom_maxP_005 = rownames(ckd_maxP_005)

save(ckd_glom_maxP_001, ckd_glom_maxP_005, file = paste0(path, "Data/ckd_glom_maxP_gene_names.RData"))


ckd_glom_maxP = list()
ckd_glom_maxP$fdr01 = ckd_glom_maxP_001

#Missing for now 
ckd_glom_maxP$fdr05 = ckd_glom_maxP_005  


save(ckd_glom_maxP, file = paste0(path, "Data/CKD_Glom_maxP_List.RData"))

################################################################################################################################
###################################################################THE END##########################################
####################The following things below are not needed##########
#############################################################