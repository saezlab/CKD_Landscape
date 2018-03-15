#I. Aim of this script:
#We aim to get insight into a pathway's influence in 
#a certain disease entity. We can get a simulated empirical distribution by generating 10000 resampled PROGENy matricies.
#The key step is shuffling the original gene expression values (standardised over the sd of tumor nephrectomy), hence randomly allocating different expression values to different disease subtypes. This way the original distribution stays the same, but the values are randomly assigned. We then multiply this 1000 series of gene expression matrices with the PROGENy model matrix and attain one single PROGENY score/pathway in a disease (we have that 1000x).
#Mapping the original PROGENy score/Pathway/Disease onto the generated ECDF,
#we would be able to infer the probability of the score being located within or out of range under the distribution curve.
#This way we can obtain a p-value indicating a  significance of the pathway in that particular disease.

#Steps:
 #We merge all the accessions of all the disease subtypes into one single matrix. We do this separately for the different tissues. 
 #We shuffle the coloumns in a way that different expression values are assigned to different disease labels. It is done by 10000 times. So, we have 1000 matrices containing random expression values (based on the real expression values). 
 #We average the accessions taking into account the number of replicates/disease subtype. #function Newfun.
 #Run PROGENy model matrix on all the newly created random matricies.
 #Compute ECDF pathway-specifically
 #Compute Pathway-Specific P-values
 #pval = 1 - pathway_spec_ecdf(the ##real## progeny score you want to test) #this formula for positive PROGENy scores
 #pval = pathway_spec_ecdf(the ##real## progeny score you want to test)    #this formula for negative PROGENy scores
 #FDR-adjustment (BH) of p-values
 #Visualisation

rm(list = ls())


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
###################################################################
