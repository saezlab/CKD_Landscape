
path = "/Users/francescoceccarelli/Desktop/CKD Landscape code/"
load(paste0(path, "Data/ckd_glom_sepMat.RData"))  #########Separating the Conditions into distinct objects###################




# download DoRothEA v1 from https://github.com/saezlab/DoRothEA/releases/tag/version1 
# load in function lib_enrichment_scores.r located in the src folder 
# to estimate Transcription Factor activity
source(paste0(path, "Data/lib_enrichment_scores.r"))
# load in CTFRs_v122016.rdata (TF regulons)
# located in data/regulons from DoRothEA v1 folder
load(paste0(path, "Data/geneset2.RData"))

library(gplots)

ckd_glom = list(FSGS_Glom, FSGS_MCD_Glom, MCD_Glom, IgAN_Glom, LN_Glom, MGN_Glom,
                      DN_Glom, HT_Glom, RPGN_Glom, TN_Glom)

save(ckd_glom, file = paste0(path, "Data/ckd_glom.RData"))


#creating a scaffold matrix that will later be filled with the average of each gene/disease
ckdglom_gex = matrix(9, ncol = 6289, nrow = 10) 

# for loop function to make to average of all genes/disease;using the list of matrices created earlier(kd_quantumnormprog14)

for (i in 1:10) {ckdglom_gex[i,] = rowMeans(ckd_glom[[i]])
}


rownames(ckdglom_gex) = c("FSGS", "FSGS-MCD", "MCD", "IgAN", "LN", "MGN", "DN", "HT", "RPGN", "TN")
colnames(ckdglom_gex) = rownames(FSGS_Glom)

ckdglom_gex = t(ckdglom_gex)

save(ckdglom_gex, file = paste0(path, "Data/ckdglom_gex.RData"))

#####TF Activity Estimation by DoROthEA###Figure3####################################################################################

#Feeding cyclic loess norm values into DoROthEA ###the one I USED!!!!CKD PAPER!!#####

#Transcription Factor Activity on 


aa_norm = gene_expression_statistic(ckdglom_gex, method = "scale", rnaseq = FALSE)


Glom_gex = SLEA(aa_norm, geneset2, method = "VIPER", filter_E = T)
Glom_gex_Factors = names(Glom_gex$regulons)
which_in_matrixgex = which((rownames(aa_norm)) %in% Glom_gex_Factors)
Glom_gex_transcription_factors = aa_norm[which_in_matrixgex,]

save(Glom_gex, file = paste0(path, "Data/Glom_gex.RData"))

save(Glom_gex_transcription_factors, file =  paste0(path, "Data/Glom_gex_transcription_factors.RData"))

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

save(corr, corr2, file =  paste0(path, "Data/correlate_TF_Expression.RData"))

temp = Glom_gex$FDR
finalwhich2 = apply(temp[fwhichglom, ],1, function(x) any(x < 0.05))
finalwhich2 = which(finalwhich2 == TRUE)


tf_activity = Glom_gex$ES
tf_activity = tf_activity[fwhichglom,]
tf_activity = tf_activity[finalwhich2,]


correlate = round(corr2[finalwhich2], 3)
save(correlate, file =  paste0(path, "Data/factor_expression_corr.RData"))


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

#Glomerulus

fdr_star_glom = fdr

fdr_star_glom[which(fdr_star_glom >= 0.05)] = NA


fdr_star_glom = mystars_fdr_glom <- ifelse(fdr_star_glom < .001, " **** ", 
                                           ifelse(fdr_star_glom < .005, " *** ",
                                                  ifelse(fdr_star_glom < .01, " ** " ,
                                                         ifelse(fdr_star_glom < .05, " * ",
                                                                " "))))                                                    

pdf(paste0(path, "Plot/tf_activity_glom_gex_scale.pdf"))

heatmap.2(tf_activity,
          cellnote = fdr_star_glom,  # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,10),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=100), 
          scale = "none",
          dendrogram="none",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA",
          Rowv = "NA",
          labRow = final_names, labCol = final_col_names, cexRow = 1, cexCol = 1.5)
dev.off()





#Checking whether the factors present in tf_activity_effs_print are present in the gene expression data:

#overlap btw tf_glom and ckd_glom_g (gene expression values represented by  scaled Hedges'g effect size values)

list_tf = list(tf_activity, aa_norm)
tf_id = Reduce(intersect, lapply(list_tf, rownames))
ckd_glom_gex_tf = aa_norm[tf_id, ]
save(ckd_glom_gex_tf, file =  paste0(path, "Data/ckd_glom_gex_tf.RData"))


pdf(paste0(path, "Plot/tf_gex_tf.pdf"))
heatmap.2(ckd_glom_gex_tf, # same data set for cell labels
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),    # widens margins around plot
          col = colorRampPalette(c("blue","white","hotpink"))(n=100), 
          dendrogram="none",
          symm = F,
          symkey = F,
          symbreaks = T,
          Colv = "NA", Rowv = "NA",labRow = final_names, labCol = final_col_names, cexRow = 1, cexCol = 1.5)
           # turn off column
dev.off()
################################################################################
