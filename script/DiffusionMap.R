
list.of.packages <- c("gplots", "YuGene","RColorBrewer", "scater","scran", "destiny", "scatterplot3d")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


ckd_colors = c("#a6cee3","#1f78b4","black","#33a02c","#fb9a99","#e31a1c","#999999","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
ckd_names = c("FSGS-MCD","MCD","IgA Nephropathy","FSGS","Tumor Nephrectomy","Healthy Living Donor","Lupus Nephritis","Diabetic Nephropathy","TMD","Hypertensive Nephropathy","MGN","RPGN")
ckd_names = c(ckd_names[6],ckd_names[5],ckd_names[4],ckd_names[1],ckd_names[2],ckd_names[3],ckd_names[7],ckd_names[11],ckd_names[8],ckd_names[10],ckd_names[9],ckd_names[12])
ckd_colors = c(ckd_colors[6],ckd_colors[5],ckd_colors[4],ckd_colors[1],ckd_colors[2],ckd_colors[3],ckd_colors[7],ckd_colors[11],ckd_colors[8],ckd_colors[10],ckd_colors[9],ckd_colors[12])

#order: [1]HLD, [2]Tumor Nephrectomy, [3]FSGS, [4]FSGS/MCD, [5]MCD, [6]IgAN, [7]Lupus Nephritis, [8]MGN, [9]DN, [10]Hypertensive Neph., [11]TMD, [12]RPGN

path = "/Users/francescoceccarelli/Desktop/CKD Landscape code/"

load(paste0(path,"Data/process_glom_all.RData"))


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
pdf(paste0(path, "Plot/number_of_patients_glom_platform.pdf"))
heatmap.2(ckd_numbers, scale = "none", dendrogram = "none", Rowv = NA, Colv = NA, col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Reds")))(n=100), trace = "none", density.info = "none", margins = c(15,15), cellnote = ckd_numbers, notecol = "black", breaks = seq(0,40,length.out = 101))
dev.off()

ckd_numbers = table(info$Group, info$Study)
ckd_numbers = rbind(ckd_numbers[4,], ckd_numbers[11,], ckd_numbers[2,], ckd_numbers[3,], ckd_numbers[8,], ckd_numbers[6,], ckd_numbers[7,], ckd_numbers[9,], ckd_numbers[1,], ckd_numbers[5,], ckd_numbers[10,])
rownames(ckd_numbers) = ckd_names[-11]
pdf(paste0(path,"Plot/number_of_patients_glom_study.pdf"))
heatmap.2(ckd_numbers, scale = "none", dendrogram = "none", Rowv = NA, Colv = NA, col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Reds")))(n=100), trace = "none", density.info = "none", margins = c(15,15), cellnote = ckd_numbers, notecol = "black", breaks = seq(0,40,length.out = 101))
dev.off()
##################


########Correlation Heatmap############Figure1C###########################
ckd_colors2 = c(ckd_colors[9],ckd_colors[3],ckd_colors[4],ckd_colors[1],ckd_colors[10],ckd_colors[6],ckd_colors[7],ckd_colors[5],ckd_colors[8],ckd_colors[12],ckd_colors[2])
pal_study = c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0")
pal_platform = c("#377eb8","#e41a1c")

aa = cor(gex, method = "spearman", use = "pairwise.complete.obs")

#by Group
aa2 = aa
rownames(aa2) = info$Group
colnames(aa2) = info$Group
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(rownames(aa2) == (unique(rownames(aa2)))[i])
  colors[wtemp] = ckd_colors2[i]
}
pdf(paste0(path,"Plot/corr_gex_glom_cluster_group.pdf"))
heatmap.2(aa2, scale = "none", dendrogram = "row", col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Blues")))(n=100), trace = "none", RowSideColors = colors, ColSideColors = colors, density.info = "none", labRow = "", labCol = "", margins = c(25,25))
legend("topright",legend = ckd_names[-11], col = ckd_colors[-11], lty= 1, lwd = 10, bty = "n")
dev.off()



#by Study
aa2 = aa
rownames(aa2) = info$Study
colnames(aa2) = info$Study
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors2 = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(colnames(aa2) == (unique(rownames(aa2)))[i])
  colors2[wtemp] = pal_study[i]
}
pdf(paste0(path,"Plot/corr_gex_glom_cluster_study.pdf"))
heatmap.2(aa2, scale = "none", dendrogram = "row", col = colorRampPalette(c("white",brewer.pal(n = 9, name = "Blues")))(n=100), trace = "none", RowSideColors = colors2, ColSideColors = colors2, density.info = "none", labRow = "", labCol = "", margins = c(25,25))
legend("topright",legend = unique(rownames(aa2)), col = pal_study, lty= 1, lwd = 10, bty = "n")
dev.off()




#by Platform
aa2 = aa
rownames(aa2) = info$Platform
colnames(aa2) = info$Platform
aa2 = aa2[order(rownames(aa2)), order(rownames(aa2))]
colors3 = rep("", length = length(info$Group))
for (i in 1:length(unique(rownames(aa2)))) {
  wtemp = which(colnames(aa2) == (unique(rownames(aa2)))[i])
  colors3[wtemp] = pal_platform[i]
}
pdf(paste0(path,"Plot/corr_gex_glom_cluster_platform.pdf"))
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

pdf(paste0(path,"Plot/variance_gex_glom.pdf"))
gex_y = plotExplanatoryVariables(gex_y, exprs_values = "counts", variables = c("Group", "Study", "Platform"), return_object = TRUE)
dev.off()
############################################################################################


###################diffusion map##############################
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

# Union of genes that are significant in both TN and HLD
load(paste0(path, "Data/CKD_Glom_maxP_List.RData"))
ckd_glom_maxP = union(ckd_glom_maxP$fdr01, ckd_glom_maxP$fdr05)
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

dm = DiffusionMap(gex1_y, sigma = "local")
ev = eigenvectors(dm)


# addalpha()
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
colorsT = addalpha(colors, 0.6)

pdf(paste0(path, "Plot/diffusion_gex_glom_union.pdf"))
scatterplot3d(ev[,c(1,3,2)], box = F, grid = F, color = colorsT, pch = 16, type = "p")
addgrids3d(ev[,c(1,3,2)], grid = c("xy","xz","yz"), col.grid = "#99999950")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom12_union.pdf"))
plot(ev[,c(1,2)], col = colorsT, pch = 16, type = "p")
dev.off()

pdf(paste0(path, "Plot/diffusion_gex_glom13_union.pdf"))
plot(ev[,c(1,3)], col = colorsT, pch = 16, type = "p")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom23_union.pdf"))
plot(ev[,c(2,3)], col = colorsT, pch = 16, type = "p")
dev.off()
########################################################################






###################diffusion map##############################
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')


# using the genes that are significant in HLD only
load(paste0(path, "Data/CKD_Glom_maxP_List.RData"))
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



pdf(paste0(path, "Plot/diffusion_gex_glom_hldonly.pdf"))
scatterplot3d(ev[,c(1,3,2)], box = F, grid = F, color = colorsT, pch = 16, type = "p")
addgrids3d(ev[,c(1,3,2)], grid = c("xy","xz","yz"), col.grid = "#99999950")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom12_hldonly.pdf"))
plot(ev[,c(1,2)], col = colorsT, pch = 16, type = "p")
dev.off()

pdf(paste0(path, "Plot/diffusion_gex_glom13_hldonly.pdf"))
plot(ev[,c(1,3)], col = colorsT, pch = 16, type = "p")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom23_hldonly.pdf"))
plot(ev[,c(2,3)], col = colorsT, pch = 16, type = "p")
dev.off()
###################################################################

###################diffusion map Figure2B & SupplementaryFigure2##############################
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

# using the genes that are significant in Tumor Nephrectomy only
load(paste0(path, "Data/CKD_Glom_maxP_List.RData"))
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



pdf(paste0(path, "Plot/diffusion_gex_glom_tnonly.pdf"))
scatterplot3d(ev[,c(1,2,3)], box = F, grid = F, color = colorsT, pch = 16, type = "p")
addgrids3d(ev[,c(1,2,3)], grid = c("xy","xz","yz"), col.grid = "#99999950")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom12_tnonly.pdf"))
plot(ev[,c(1,2)], col = colorsT, pch = 16, type = "p")
dev.off()

pdf(paste0(path, "Plot/diffusion_gex_glom13_tnonly.pdf"))
plot(ev[,c(1,3)], col = colorsT, pch = 16, type = "p")
dev.off()


pdf(paste0(path, "Plot/diffusion_gex_glom23_tnonly.pdf"))
plot(ev[,c(2,3)], col = colorsT, pch = 16, type = "p")
dev.off()
###################################################################


