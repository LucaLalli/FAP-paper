# Prepare Cytokines ####
data_cytokines <- data[, c("Id", "Time", name_cytokines)]

for (i in name_cytokines) {
  data_cytokines[,i] <- as.numeric(data_cytokines[,i])
}

# Removing all IDs with missing values for Friedman test
tmp <- rowSums(is.na(data[, name_cytokines]))
Id_removes <- data$Id[tmp>0]

data_cytokines_tmp <- subset(data_cytokines, !(data_cytokines$Id %in% Id_removes))    

# Scaling values between -1 and 1
min = -1
max = 1
data_cyt_norm_tmp = data_cytokines[, -c(1,2)]
data_cyt_norm_tmp = as.data.frame(apply(data_cyt_norm_tmp, 2, rescale, c(min,max)))

data_cyt_norm_tmp = cbind(data_cytokines[, c(1,2)],data_cyt_norm_tmp)

# Normalizing baseline by computing delta values
data_cyt_norm = NULL
for(i in unique(data_cytokines$Id)){
  data_cyt_norm = rbind(data_cyt_norm, data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                                           data_cyt_norm_tmp$Time == 0,-c(1,2)] - 
                          data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                              data_cyt_norm_tmp$Time == 0, -c(1,2)],
                        data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                            data_cyt_norm_tmp$Time == 3,-c(1,2)] - 
                          data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                              data_cyt_norm_tmp$Time == 0, -c(1,2)],
                        data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                            data_cyt_norm_tmp$Time == 6,-c(1,2)] - 
                          data_cyt_norm_tmp[data_cyt_norm_tmp$Id == i &
                                              data_cyt_norm_tmp$Time == 0, -c(1,2)])
}


data_cyt_norm = cbind(data_cytokines[, c(1,2)], data_cyt_norm)
rm(data_cyt_norm_tmp)

# Friedman test to select only curves that differ significantly
df_friedman <- as.data.frame(matrix(, nrow = length(name_cytokines), ncol = 2))
colnames(df_friedman) <- c("Cytokine", "Pvalue")
df_friedman$Cytokine <- name_cytokines

for (i in name_cytokines) {
  frm <- formula(paste0("`", i, "`", "~", "Time |Id"))
  df_friedman$Pvalue[df_friedman$Cytokine == i] <- round(friedman_test(data_cytokines_tmp, frm)$p, 4)
}

name_cytokines <- df_friedman$Cytokine[df_friedman$Pvalue < 0.05]

data_cytokines_conn <- as.data.frame(matrix(NA, nrow = length(name_ptzs)*length(name_cytokines)*3, ncol = 3))
colnames(data_cytokines_conn) <- c("ID", "Observations", "Time")

# Creation of new dataset in long format for CONNECTOR
k <- 1
vet_newid <- c()
for (i in name_cytokines) {
  tmp <- paste0(data_cytokines$Id, "_", i)
  vet_newid <- c(vet_newid, tmp)
  k <- k+1
}

data_cytokines_conn$ID <- vet_newid
rm(k, vet_newid, tmp)

data_cytokines_conn$Time <- rep(name_timepoints_num, nrow(data_cytokines_conn)/length(name_timepoints_num))

k <- 1
vet_newvars <- c()
for (i in name_cytokines) {
  tmp <- data_cyt_norm[, i]
  vet_newvars <- c(vet_newvars, tmp)
  k <- k+1
}
data_cytokines_conn$Observations <- vet_newvars
rm(k, vet_newvars, tmp)

data_cytokines_conn$Observations <- as.numeric(data_cytokines_conn$Observations)

# Removing only the variables containing a missing value for any of the IDs
vet_removeID <- data_cytokines_conn$ID[is.na(data_cytokines_conn$Observations)]
data_cytokines_conn <- subset(data_cytokines_conn, !(data_cytokines_conn$ID %in% vet_removeID))

# Annotation file

cytokines_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_cytokines_conn)/length(name_timepoints_num), ncol = 3))
colnames(cytokines_annotation) <- c("ID", "IdSample", "Cytokine")
cytokines_annotation$ID <- unique(data_cytokines_conn$ID)
cytokines_annotation$IdSample <- gsub("_.*", "", cytokines_annotation$ID)
cytokines_annotation$Cytokine <- gsub(".*_", "", cytokines_annotation$ID)


# Prepare Emocromo ####
data_emocromo <- data[, c("Id", "Time", name_emocromos)]

for (i in name_emocromos) {
  data_emocromo[,i] <- as.numeric(data_emocromo[,i])
}

# Removing all IDs with missing values for Friedman test
tmp <- rowSums(is.na(data[, name_emocromos]))
Id_removes <- data$Id[tmp>0]

data_emocromo_tmp <- subset(data_emocromo, !(data_emocromo$Id %in% Id_removes))

# Scaling values between -1 and 1
min = -1
max = 1
data_emo_norm_tmp = data_emocromo[,-c(1,2)]
data_emo_norm_tmp = as.data.frame(apply(data_emo_norm_tmp, 2, rescale, c(min,max)))

data_emo_norm_tmp = cbind(data_emocromo[,c(1,2)],data_emo_norm_tmp)

# Normalizing baseline by computing delta values
data_emo_norm = NULL
for(i in unique(data_emocromo$Id)){
  data_emo_norm = rbind(data_emo_norm, data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                                           data_emo_norm_tmp$Time == 0,-c(1,2)] - 
                          data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                              data_emo_norm_tmp$Time == 0, -c(1,2)],
                        data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                            data_emo_norm_tmp$Time == 3,-c(1,2)] - 
                          data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                              data_emo_norm_tmp$Time == 0, -c(1,2)],
                        data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                            data_emo_norm_tmp$Time == 6,-c(1,2)] - 
                          data_emo_norm_tmp[data_emo_norm_tmp$Id == i &
                                              data_emo_norm_tmp$Time == 0, -c(1,2)])
}


data_emo_norm = cbind(data_emocromo[,c(1,2)],data_emo_norm)
rm(data_emo_norm_tmp)


# Friedman test to select only curves that differ significantly
df_friedman <- as.data.frame(matrix(, nrow = length(name_emocromos), ncol = 2))
colnames(df_friedman) <- c("Emocromo", "Pvalue")
df_friedman$Emocromo <- name_emocromos

for (i in name_emocromos) {
  frm <- formula(paste0("`", i, "`", "~", "Time |Id"))
  df_friedman$Pvalue[df_friedman$Emocromo == i] <- round(friedman_test(data_emocromo_tmp, frm)$p, 4)
}

name_emocromos <- df_friedman$Emocromo[df_friedman$Pvalue < 0.05]

# Filter for Neutrophils and NLR
name_emocromos = name_emocromos[which(grepl('Neutrophils', name_emocromos) | grepl('NLR',name_emocromos))]

# Creation dataset for CONNECTOR
data_emocromo_conn <- as.data.frame(matrix(NA, nrow = length(name_ptzs)*length(name_emocromos)*3, ncol = 3))
colnames(data_emocromo_conn) <- c("ID", "Observations", "Time")

k <- 1
vet_newid <- c()
for (i in name_emocromos) {
  tmp <- paste0(data_emocromo$Id, "_", i)
  vet_newid <- c(vet_newid, tmp)
  k <- k+1
}

data_emocromo_conn$ID <- vet_newid
rm(k, vet_newid, tmp)

data_emocromo_conn$Time <- rep(name_timepoints_num, nrow(data_emocromo_conn)/length(name_timepoints_num))

k <- 1
vet_newvars <- c()
for (i in name_emocromos) {
  tmp <- data_emo_norm[, i]
  vet_newvars <- c(vet_newvars, tmp)
  k <- k+1
}
data_emocromo_conn$Observations <- vet_newvars
rm(k, vet_newvars, tmp)

data_emocromo_conn$Observations <- as.numeric(data_emocromo_conn$Observations)

vet_removeID <- data_emocromo_conn$ID[is.na(data_emocromo_conn$Observations)]

data_emocromo_conn <- subset(data_emocromo_conn, !(data_emocromo_conn$ID %in% vet_removeID))

# Annotation file
emocromo_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_emocromo_conn)/length(name_timepoints_num), ncol = 3))
colnames(emocromo_annotation) <- c("ID", "IdSample", "Granulocytes")
emocromo_annotation$ID <- unique(data_emocromo_conn$ID)
emocromo_annotation$IdSample <- gsub("_.*", "", emocromo_annotation$ID)
emocromo_annotation$Granulocytes <- gsub(".*_", "", emocromo_annotation$ID)


# Prepare Myeloids ####

data_myeloids <- data[, c("Id", "Time", name_myeloids)]

for (i in name_myeloids) {
  data_myeloids[,i] <- as.numeric(data_myeloids[,i])
}

# Removing all IDs with missing values for Friedman test
tmp <- rowSums(is.na(data[, name_myeloids]))
Id_removes <- data$Id[tmp>0]

data_myeloids_tmp <- subset(data_myeloids, !(data_myeloids$Id %in% Id_removes))    

# Scaling values between -1 and 1
min = -1
max = 1
data_myel_norm_tmp = data_myeloids[,-c(1,2)]
data_myel_norm_tmp = as.data.frame(apply(data_myel_norm_tmp, 2, rescale, c(min,max)))

data_myel_norm_tmp = cbind(data_myeloids[,c(1,2)],data_myel_norm_tmp)

# Normalizing baseline by computing delta values
data_myel_norm = NULL
for(i in unique(data_myeloids$Id)){
  data_myel_norm = rbind(data_myel_norm, data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                                              data_myel_norm_tmp$Time == 0,-c(1,2)] - 
                           data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                                data_myel_norm_tmp$Time == 0, -c(1,2)],
                         data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                              data_myel_norm_tmp$Time == 3,-c(1,2)] - 
                           data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                                data_myel_norm_tmp$Time == 0, -c(1,2)],
                         data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                              data_myel_norm_tmp$Time == 6,-c(1,2)] - 
                           data_myel_norm_tmp[data_myel_norm_tmp$Id == i &
                                                data_myel_norm_tmp$Time == 0, -c(1,2)])
}


data_myel_norm = cbind(data_myeloids[,c(1,2)],data_myel_norm)
rm(data_myel_norm_tmp)

# Friedman test to select only curves that differ significantly 

df_friedman <- as.data.frame(matrix(, nrow = length(name_myeloids), ncol = 2))
colnames(df_friedman) <- c("Myeloids", "Pvalue")
df_friedman$Myeloids <- name_myeloids

for (i in name_myeloids) {
  frm <- formula(paste0("`", i, "`", "~", "Time |Id"))
  df_friedman$Pvalue[df_friedman$Myeloids == i] <- round(friedman_test(data_myeloids_tmp, frm)$p, 4)
}

name_myeloids <- df_friedman$Myeloids[df_friedman$Pvalue < 0.05]

data_myeloids_conn <- as.data.frame(matrix(NA, nrow = length(name_ptzs)*length(name_myeloids)*3, ncol = 3))
colnames(data_myeloids_conn) <- c("ID", "Observations", "Time")

# Creation of new dataset in long format for CONNECTOR
k <- 1
vet_newid <- c()
for (i in name_myeloids) {
  tmp <- paste0(data_myeloids$Id, "_", i)
  vet_newid <- c(vet_newid, tmp)
  k <- k+1
}

data_myeloids_conn$ID <- vet_newid
rm(k, vet_newid, tmp)

data_myeloids_conn$Time <- rep(name_timepoints_num, nrow(data_myeloids_conn)/length(name_timepoints_num))

k <- 1
vet_newvars <- c()
for (i in name_myeloids) {
  tmp <- data_myel_norm[, i]
  vet_newvars <- c(vet_newvars, tmp)
  k <- k+1
}
data_myeloids_conn$Observations <- vet_newvars
rm(k, vet_newvars, tmp)

data_myeloids_conn$Observations <- as.numeric(data_myeloids_conn$Observations)

# Removing only the variables containing a missing value for any of the IDs
vet_removeID <- data_myeloids_conn$ID[is.na(data_myeloids_conn$Observations)]

data_myeloids_conn <- subset(data_myeloids_conn, !(data_myeloids_conn$ID %in% vet_removeID))

# Annotation file
myeloids_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_myeloids_conn)/length(name_timepoints_num), ncol = 3))
colnames(myeloids_annotation) <- c("ID", "IdSample", "Myeloids")
myeloids_annotation$ID <- unique(data_myeloids_conn$ID)
myeloids_annotation$IdSample <- gsub("_.*", "", myeloids_annotation$ID)
myeloids_annotation$Myeloids <- gsub(".*_", "", myeloids_annotation$ID)


# Prepare Lymphoids ####
data_lymphoids <- data[, c("Id", "Time", name_lymphoids)]

for (i in name_lymphoids) {
  data_lymphoids[,i] <- as.numeric(data_lymphoids[,i])
}

# Removing all IDs with missing values for Friedman test
tmp <- rowSums(is.na(data[, name_lymphoids]))
Id_removes <- data$Id[tmp>0]

data_lymphoids_tmp <- subset(data_lymphoids, !(data_lymphoids$Id %in% Id_removes))    

# Scaling values between -1 and 1
min = -1
max = 1
data_lymph_norm_tmp = data_lymphoids[,-c(1,2)]
data_lymph_norm_tmp = as.data.frame(apply(data_lymph_norm_tmp, 2, rescale, c(min,max)))

data_lymph_norm_tmp = cbind(data_lymphoids[,c(1,2)],data_lymph_norm_tmp)

# Normalizing baseline by computing delta values
data_lymph_norm = NULL
for(i in unique(data_lymphoids$Id)){
  data_lymph_norm = rbind(data_lymph_norm, data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                                 data_lymph_norm_tmp$Time == 0,-c(1,2)] - 
                            data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                  data_lymph_norm_tmp$Time == 0, -c(1,2)],
                          data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                data_lymph_norm_tmp$Time == 3,-c(1,2)] - 
                            data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                  data_lymph_norm_tmp$Time == 0, -c(1,2)],
                          data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                data_lymph_norm_tmp$Time == 6,-c(1,2)] - 
                            data_lymph_norm_tmp[data_lymph_norm_tmp$Id == i &
                                                  data_lymph_norm_tmp$Time == 0, -c(1,2)])
}

data_lymph_norm = cbind(data_lymphoids[,c(1,2)],data_lymph_norm)
rm(data_lymph_norm_tmp)


# Friedman test to select only curves that differ significantly 
df_friedman <- as.data.frame(matrix(, nrow = length(name_lymphoids), ncol = 2))
colnames(df_friedman) <- c("lymphoids", "Pvalue")
df_friedman$lymphoids <- name_lymphoids

for (i in name_lymphoids) {
  frm <- formula(paste0("`", i, "`", "~", "Time |Id"))
  df_friedman$Pvalue[df_friedman$lymphoids == i] <- round(friedman_test(data_lymphoids_tmp, frm)$p, 4)
}

name_lymphoids <- df_friedman$lymphoids[df_friedman$Pvalue < 0.05]

data_lymphoids_conn <- as.data.frame(matrix(NA, nrow = length(name_ptzs)*length(name_lymphoids)*3, ncol = 3))
colnames(data_lymphoids_conn) <- c("ID", "Observations", "Time")

# Creation of new dataset in long format for CONNECTOR
k <- 1
vet_newid <- c()
for (i in name_lymphoids) {
  tmp <- paste0(data_lymphoids$Id, "_", i)
  vet_newid <- c(vet_newid, tmp)
  k <- k+1
}

data_lymphoids_conn$ID <- vet_newid
rm(k, vet_newid, tmp)

data_lymphoids_conn$Time <- rep(name_timepoints_num, nrow(data_lymphoids_conn)/length(name_timepoints_num))

k <- 1
vet_newvars <- c()
for (i in name_lymphoids) {
  tmp <- data_lymph_norm[, i]
  vet_newvars <- c(vet_newvars, tmp)
  k <- k+1
}
data_lymphoids_conn$Observations <- vet_newvars
rm(k, vet_newvars, tmp)

data_lymphoids_conn$Observations <- as.numeric(data_lymphoids_conn$Observations)

# Removing only the variables containing a missing value for any of the IDs
vet_removeID <- data_lymphoids_conn$ID[is.na(data_lymphoids_conn$Observations)]

data_lymphoids_conn <- subset(data_lymphoids_conn, !(data_lymphoids_conn$ID %in% vet_removeID))

# Annotation file
lymphoids_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_lymphoids_conn)/length(name_timepoints_num), ncol = 3))
colnames(lymphoids_annotation) <- c("ID", "IdSample", "lymphoids")
lymphoids_annotation$ID <- unique(data_lymphoids_conn$ID)
lymphoids_annotation$IdSample <- gsub("_.*", "", lymphoids_annotation$ID)
lymphoids_annotation$lymphoids <- gsub(".*_", "", lymphoids_annotation$ID)


# Merge all populations ####
data_tot_norm = cbind(data_cyt_norm, data_emo_norm[,c(3:ncol(data_emo_norm))], data_myel_norm[,c(3:ncol(data_myel_norm))], data_lymph_norm[,c(3:ncol(data_lymph_norm))])

name_tot = colnames(data_tot_norm[,3:ncol(data_tot_norm)])
name_tot = name_tot[!name_tot %in% c('Leukocytes','Monocytes','Lymphocytes','PLR')]

df_friedman_tot <- as.data.frame(matrix(, nrow = length(name_tot), ncol = 2))
colnames(df_friedman_tot) <- c("Variables", "Pvalue")
df_friedman_tot$Variables <- name_tot

for (i in name_tot) {
  vet_removeID <- data_tot_norm$Id[is.na(data_tot_norm[,i])]
  
  data_tmp <- subset(data_tot_norm, !(data_tot_norm$Id %in% vet_removeID))
  frm <- formula(paste0("`", i, "`", "~", "Time |Id"))
  df_friedman_tot$Pvalue[df_friedman_tot$Variables == i] <- round(friedman_test(data_tmp, frm)$p, 4)
}

name_tot <- df_friedman_tot$Variables[df_friedman_tot$Pvalue < 0.05]

data_tot_conn <- as.data.frame(matrix(NA, nrow = length(name_ptzs)*length(name_tot)*3, ncol = 3))
colnames(data_tot_conn) <- c("ID", "Observations", "Time")


## Creation of merged dataset in long format for CONNECTOR

k <- 1
vet_newid <- c()
for (i in name_tot) {
  tmp <- paste0(data_tot_norm$Id, "_", i)
  vet_newid <- c(vet_newid, tmp)
  k <- k+1
}

data_tot_conn$ID <- vet_newid
rm(k, vet_newid, tmp)

data_tot_conn$Time <- rep(name_timepoints_num, nrow(data_tot_conn)/length(name_timepoints_num))

k <- 1
vet_newvars <- c()
for (i in name_tot) {
  tmp <- data_tot_norm[, i]
  vet_newvars <- c(vet_newvars, tmp)
  k <- k+1
}
data_tot_conn$Observations <- vet_newvars
rm(k, vet_newvars, tmp)


data_tot_conn$Observations <- as.numeric(data_tot_conn$Observations)

## Removing only the variables containing a missing value for any of the IDs
vet_removeID <- data_tot_conn$ID[is.na(data_tot_conn$Observations)]
data_tot_conn <- subset(data_tot_conn, !(data_tot_conn$ID %in% vet_removeID))

# Annotation file
tot_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_tot_conn)/length(name_timepoints_num), ncol = 3))
colnames(tot_annotation) <- c("ID", "IdSample", "Variables")
tot_annotation$ID <- unique(data_tot_conn$ID)
tot_annotation$IdSample <- gsub("_.*", "", tot_annotation$ID)
tot_annotation$Variables <- gsub(".*_", "", tot_annotation$ID)

# CONNECTOR pipeline ####
tot_List <- DataFrameImport(TimeSeriesDataFrame = data_tot_conn, AnnotationFrame = tot_annotation)
CurvesPlot <- PlotTimeSeries(data = tot_List, feature = "Variables")
CurvesPlot

## Choice of the number of splines p
CrossLogLike <- BasisDimension.Choice(data = tot_List, p = 2:3 )
CrossLogLike$CrossLogLikePlot

p <- 3
ClusteringList <- ClusterAnalysis(data = tot_List, G = 2:5, p = p, runs = 100, Cores = 8)
indexes <- IndexesPlot.Extrapolation(ClusteringList) 
indexes$Plot

ConsMatrix <- ConsMatrix.Extrapolation(stability.list = ClusteringList)
ConsMatrix$G4$ConsensusPlot

CONNECTORList.FCM.opt <- MostProbableClustering.Extrapolation(stability.list = ClusteringList, G = 4)
FCMplots <- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.opt, feature = "Variables", labels = c("Time","Value"), title = "FCM model")

# Create matrix modulations Patients x variables
df_ID_cluster <-  FCMplots$plotsCluster$ALL$layers[[1]]$data[, c("ID", "Cluster", 'Info')]
df_ID_IDSample <- CurvesPlot$data$LabCurv[,c('ID','IdSample')]

df_ID_cluster = unique(df_ID_cluster)
rownames(df_ID_cluster) = rownames(df_ID_IDSample)
df_ID_cluster$IdSample = df_ID_IDSample$IdSample

df_merged_tot = as.data.frame(matrix(NA, nrow = length(name_ptzs), ncol = length(unique(df_ID_cluster$Info))))
colnames(df_merged_tot) = unique(df_ID_cluster$Info)
rownames(df_merged_tot) = name_ptzs
markers = unique(df_ID_cluster$Info)
samples = rownames(df_merged_tot)

for(j in samples){
  for(i in markers){
    if(length(df_ID_cluster$Cluster[df_ID_cluster$Info == i & df_ID_cluster$IdSample == j]) != 0){
      df_merged_tot[j,i] = df_ID_cluster$Cluster[df_ID_cluster$Info == i & df_ID_cluster$IdSample == j]
    }
  }
}

# Exclusion two patients with complete missing time points
df_merged_tot = df_merged_tot[-which(rownames(df_merged_tot) == 'FAP3' | rownames(df_merged_tot) == 'FAP35') ,]


# Random forest imputations ####
df_imp = as.data.frame(lapply(df_imp, as.factor))
set.seed(4478432)
# prova miss var
list_imp <- list()
for (i in 1:200) {
  imp <- missForest(df_imp, ntree = 2000)$ximp
  list_imp[[i]] = imp
  print(i)
}

imput_df <- data.frame(matrix(NA, nrow = nrow(list_imp[[1]]), ncol = ncol(list_imp[[1]])))
colnames(imput_df) <- colnames(list_imp[[1]])

for (i in 1:nrow(imput_df)) {
  for(j in 1:ncol(imput_df)){
    cell_values <- sapply(list_imp, function(df) df[i, j])
    value <- names(which.max(table(cell_values)))
    imput_df[i,j] <- value
  }
}

rownames(imput_df) = colnames(df_merged_tot)
imput_df = as.data.frame(t(imput_df))

dist_mat = as.dist(mat_dist_f_GL(imput_df))
hierarchical_cluster <- hclust(dist_mat, method = 'ward.D2')
plot(hierarchical_cluster)
clust_hc_3 = cutree(hierarchical_cluster, k=3)
clust_hc_2 = cutree(hierarchical_cluster, k=2)
clust_hc_4 = cutree(hierarchical_cluster, k=4)
table(clust_hc_4)

# Metacluster identifications ####
clust = as.data.frame(cbind(clust_hc_2))
colnames(clust) = 'MetaClust'

clust$MetaClust = factor(clust$MetaClust, levels = 1:2, labels = c('X','Y'), order = T)
clust$IdSample <- rownames(clust)












