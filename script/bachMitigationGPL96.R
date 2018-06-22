

list.of.packages <- c("scater", "scran","biomaRt", "LSD","limma", "gplots","YuGene","RColorBrewer","scatterplot3d")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

##############importing the data
path = "/Users/francescoceccarelli/Desktop/CKD Landscape code/"
dir.create(paste0(path, "Plot"))

load(paste0(path,"Data/Platform_Hierarchy_new.rdata"))

glom_96 = as.data.frame(Glom_96)
colnames(glom_96) = sapply(colnames(glom_96), function(x) gsub(" ", "", x, fixed = TRUE))
glom_96 = glom_96[,order(colnames(glom_96))]
glom_96 = glom_96[order(rownames(glom_96)),]
glom_96 = as.matrix(glom_96)
info_glom_96 = read.table(paste0(path,"Data/Glom_96_Info.csv"), sep = ";", header = T)
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

pdf(paste0(path,"Plot/pcaTEST_glom_96_tumourNephrectomy.pdf"))
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
pdf(paste0(path,"Plot/pca_glom_96_tumourNephrectomy.pdf"))
print(plotPCA(glom_96_y[-www,ww], ntop = 10000000000, ncomponents = 3, colour_by = "Study", exprs_values = "counts"))
dev.off()

glom_96_final[www,ww] = NA
###########################################################


#######################IgAN###################
ww = which(info_glom_96$Group == "IgAN")
norm_expr = normalizeBetweenArrays(glom_96[,ww], method = "cyclicloess")
design = cbind(rep(1,36), c(rep(0,16), rep(1, 20)))
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
pdf(paste0(path,"Plot/pcaTEST_glom_96_IgAN.pdf"))
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






pdf(paste0(path,"Plot/pca_glom_96_IgAN_after.pdf"))
plotPCA(glom_96_y[-www,ww],  ntop = 10000000000, ncomponents = 2, colour_by = "Study", exprs_values = "counts")
dev.off()

###Sup.Fig.1A#################
pdf(paste0(path,"Plot/pca_glom_96_IgAN_before.pdf"))
plotPCA(glom_96_y[,ww],  ntop = 10000000000, ncomponents = 2, colour_by = "Study", exprs_values = "counts")
dev.off()

###Sup.Fig1C############
#pdf(paste0(path,"Reproduction/pcaTEST_glom_96_IgAN_pvalue.pdf"))
#plot(, go2 * 100, type = "l", lwd = 5, xlim = c(15,45), ylim = c(10.5,16))
#dev.off()

pdf(paste0(path,"Plot/pcaTEST_glom_96_IgAN_volcano.pdf"))
limma::volcanoplot(fit2, coef = 2, highlight = 50, names = rownames(top))
dev.off()

####Sup.Fig.1B############
pdf(paste0(path,"Plot/pcaTEST_glom_96_IgAN_ma.pdf"))
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

pdf(paste0(path,"Plot/pcaTEST_glom_96_healthyLivingDonor.pdf"))
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

pdf(paste0(path,"Plot/pca_glom_96_healthyLivingDonor.pdf"))
print(plotPCA(glom_96_y[-www,ww], ntop = 10000000000, ncomponents = 3, colour_by = "Study", exprs_values = "counts"))
dev.off()


glom_96_final[www,ww] = NA
#############################################################


###save##
processed_glom_96 = list(data = glom_96, info = info_glom_96, data_NA = glom_96_final)

save(processed_glom_96, file = paste0(path,"Data/process_glom_96.RData"))
#########################################
