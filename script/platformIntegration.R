
# install necessary packages if not present yet 
list.of.packages <- c("scater", "scran","biomaRt", "LSD","limma", "gplots","YuGene","RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


##############importing the data#########
path = "/Users/francescoceccarelli/Desktop/CKD Landscape code/"
 
load(paste0(path,"Data/process_glom_570.RData"))
load(paste0(path, "Data/process_glom_96.RData"))

dat_96 = processed_glom_96$data
info_96 = processed_glom_96$info
dat_96_NA = processed_glom_96$data_NA

dat_570 = processed_glom_570$data
info_570 = processed_glom_570$info
#ww = which(info_570$Group == "Healthy Living Donor")
#info_570$Group[ww] = "Healthy Living Donor 570"


common_genes = intersect(rownames(dat_570), rownames(dat_96))
dat_96 = dat_96[which(rownames(dat_96) %in% common_genes),]
dat_96_NA = dat_96_NA[which(rownames(dat_96_NA) %in% common_genes),]
dat_570 = dat_570[which(rownames(dat_570) %in% common_genes),]
dat_96 = dat_96[order(rownames(dat_96)),]
dat_96_NA = dat_96_NA[order(rownames(dat_96_NA)),]
dat_570 = dat_570[order(rownames(dat_570)),]

dat_96 = dat_96[,order(colnames(dat_96))]
dat_96_NA = dat_96_NA[,order(colnames(dat_96_NA))]
dat_570 = dat_570[,order(colnames(dat_570))]
info_96 = info_96[order(info_96$Accession),]
info_570 = info_570[order(info_570$Accession),]


dat = cbind(dat_96, dat_570)
dat_NA = cbind(dat_96_NA, dat_570)
info = rbind(info_96, info_570)

glom_all = dat
info_glom_all = info
info_glom_all$Group = as.character(info_glom_all$Group)
if(!(all(colnames(glom_all) == info_glom_all$Accession))) { message("glom_all: names don't match") }

infoFrame1 = data.frame(Cell = colnames(glom_all), info_glom_all)
rownames(infoFrame1) = colnames(glom_all)
infoFrame_glom_all = new("AnnotatedDataFrame", data = infoFrame1)
glom_all_y = newSCESet(countData = glom_all, phenoData = infoFrame_glom_all)
glom_all_y = calculateQCMetrics(glom_all_y)


glom_all_NA = dat_NA
glom_all_final = dat_NA
info_glom_all_NA = info
info_glom_all_NA$Group = as.character(info_glom_all_NA$Group)
if(!(all(colnames(glom_all_NA) == info_glom_all_NA$Accession))) { message("glom_all_NA: names don't match") }

who = (apply(glom_all_NA, 1, function(x) any(is.na(x))))
who = which(who == TRUE)

glom_all_NA = glom_all_NA[-who,] 

infoFrame1 = data.frame(Cell = colnames(glom_all_NA), info_glom_all_NA)
rownames(infoFrame1) = colnames(glom_all_NA)
infoFrame_glom_all_NA = new("AnnotatedDataFrame", data = infoFrame1)
glom_all_NA_y = newSCESet(countData = glom_all_NA, phenoData = infoFrame_glom_all_NA)
glom_all_NA_y = calculateQCMetrics(glom_all_NA_y)

#########################################
#expr = glom_all[,ww]
#exprS = ComBat(scale(expr), batch = as.character((info_glom_all$Platform)[ww]))
#exprS = scale(exprS, center = FALSE , scale=1/(attr(scale(expr), "scaled:scale")))
#exprS = scale(exprS, center = -1 * (attr(scale(expr), "scaled:center")), scale=FALSE)
#glom_all[,ww] = exprS
#exprS[exprS < 0] = 0
###############################DN##############################

ww = which(info_glom_all$Group == "Diabetic Nephropathy")
norm_expr = normalizeBetweenArrays(glom_all[,ww], method = "cyclicloess")
design = cbind(rep(1,13), c(rep(0,6), rep(1,7)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]
top = top[order(top$adj.P.Val),]



go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 50)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(counts(glom_all_y[-www_temp,ww])), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[1]
}

www = which(rownames(glom_all) %in% rownames(top[1:4000,]))

##pca plot
pdf(paste0(path,"Plot/pca_glom_all_dn.pdf"))
print(plotPCA(glom_all_y[-www,ww], ntop = 1000000, ncomponents = 3, colour_by = "Platform", exprs_values = "counts"))
dev.off()


ww_dn = ww
www_dn = www
glom_all_final[www,ww] = NA
###########################################################






##############################FSGS##############################
ww = which(info_glom_all$Group == "FSGS")
norm_expr = normalizeBetweenArrays(glom_all[,ww], method = "cyclicloess")
design = cbind(rep(1,23), c(rep(0,13), rep(1,10)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]
top = top[order(top$adj.P.Val),]


go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 50)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(counts(glom_all_y[-www_temp,ww])), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[1]
}

www = which(rownames(glom_all) %in% rownames(top[1:4000,]))

pdf(paste0(path, "Plot/pca_glom_all_fsgs.pdf"))
print(plotPCA(glom_all_y[-www,ww], ncomponents = 3, colour_by = "Platform", exprs_values = "counts"))
dev.off()

ww_fsgs = ww
www_fsgs = www
glom_all_final[www,ww] = NA
###########################################################



##############################mcd##############################
ww = which(info_glom_all$Group == "MCD")
norm_expr = normalizeBetweenArrays(glom_all[,ww], method = "cyclicloess")
design = cbind(rep(1,15), c(rep(0,10), rep(1,5)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]
top = top[order(top$adj.P.Val),]


go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 50)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(counts(glom_all_y[-www_temp,ww])), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[1]
}

www = which(rownames(glom_all) %in% rownames(top[1:4000,]))

pdf(paste0(path,"Plot/pca_glom_all_mcd.pdf"))
print(plotPCA(glom_all_y[-www,ww], ncomponents = 3, colour_by = "Platform", exprs_values = "counts"))
dev.off()

ww_mcd = ww
www_mcd = www
glom_all_final[www,ww] = NA
###########################################################



##############################hld##############################
ww = which(info_glom_all$Group == "Healthy Living Donor")
ww = c(ww, which(info_glom_all$Group == "Healthy Living Donor 570"))
norm_expr = normalizeBetweenArrays(glom_all[,ww], method = "cyclicloess")
design = cbind(rep(1,35), c(rep(0,17), rep(1,18)))
fit = lmFit(norm_expr, design)
fit2 = eBayes(fit)
top = topTable(fit2, adjust="BH", number = 50000000, coef = 2)
top = top[order(rownames(top)),]
top = top[order(top$adj.P.Val),]


go = seq(1,length(which(top$adj.P.Val < 0.01)), by = 50)
go2 = rep(0, length(go))
for(i in 1:length(go)) {
	name_temp = rownames(top[1:go[i],])
	www_temp = which(rownames(norm_expr) %in% name_temp)
	pp_temp = prcomp(t(counts(glom_all_y[-www_temp,ww])), scale = T, center = T)
	eigs = (pp_temp$sdev)^2
	eigs = (eigs)/(sum(eigs))
	go2[i] = eigs[1]
}

www = which(rownames(glom_all) %in% rownames(top[1:8000,]))

pdf(paste0(path, "Plot/pca_glom_all_hld.pdf"))
print(plotPCA(glom_all_y[-www,ww], ncomponents = 3, colour_by = "Platform", exprs_values = "counts"))
dev.off()

ww_mcd = ww
www_mcd = www
glom_all_final[www,ww] = NA
###########################################################


####Final Form of the Data: processed_glom_all#######
processed_glom_all = list(data = glom_all, info = info_glom_all, data_NA = glom_all_NA, final = glom_all_final)
save(processed_glom_all, file = paste0(path, "Data/process_glom_all.RData"))
####
