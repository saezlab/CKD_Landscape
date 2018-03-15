#Batch Effect Mitigation of Data generated of GPL570


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
