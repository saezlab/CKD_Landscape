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
