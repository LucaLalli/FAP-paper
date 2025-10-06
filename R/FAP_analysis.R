pacman::p_load(Hmisc, rms, ggplot2, dplyr, scales, tidyr, tidyverse, flextable, ComplexHeatmap, umap, officer, readxl, gridExtra,
               cowplot, fda, flexclust,ggplot2, MASS, Matrix, plyr,ggplotify, RColorBrewer, readxl, reshape2, splines, statmod, zoo,
               sfsmisc,shinyWidgets, viridis, dashboardthemes, shinybusy, shinydashboard, shinyjs, tidyr, shinyFiles,devtools, connector, BiocParallel,
               rstatix, ggsci, ggh4x, klaR, factoextra, GSDA, FactoMineR, Amelia,
               mice, missForest, naniar,cluster, purrr, doParallel, foreach,ADAPTS, lme4, sjPlot, ggpubr, rstatix, TSclust, TSdist, ggsankey, remotes, ggdendro, dendextend)

# library(devtools)
# install_github("qBioTurin/connector", ref="master",dependencies=TRUE)

knitr::opts_knit$set(root.dir = "C:/Users/Lenovo/OneDrive - Università degli Studi di Milano/Desktop/Collaborazione INT/FAP")
knitr::opts_chunk$set(include = F)
data_path = file.path("C:/Users/Lenovo/OneDrive - Università degli Studi di Milano/Desktop/Collaborazione INT/FAP")
path_final = file.path("C:/Users/Lenovo/OneDrive - Università degli Studi di Milano/Desktop/Collaborazione INT/FAP/Results_final")
#source('C:/Users/Lenovo/OneDrive - Università degli Studi di Milano/Desktop/Collaborazione INT/Funzioni/distance_matrix_categ.R')
options(scipen = 9999999)


Data_all_raw <- as.data.frame(read_excel("//fileserversvc/biomimmunol/Projects/2023 - Connector FAP - Daveri/Connector etc/Tabella Fenotipo e tutti i dati_V2.xlsx", sheet = "ALL DATA"))


name_ptzs <- colnames(Data_all_raw)[grepl("FAP", colnames(Data_all_raw))]
name_ptzs <- gsub(" ", "", name_ptzs)
name_timepoints <- c("PRE", "3M", "6M")
name_timepoints_num <- c(0, 3, 6)


vet <- c()
for (i in name_ptzs) {
   for (k in name_timepoints) {
      vet <- c(vet, paste0(i,"_", k))
   }
}
colnames(Data_all_raw) <- c("Type", "Markers", vet)
Data_all_raw <- Data_all_raw[-1,]
Data_all_raw <- Data_all_raw[-112,]
Data_all_raw$Type[Data_all_raw$Markers == "Fecal calprotectin (µg/g)"] <- "Fecal calprotectin"

Data_all_raw$Markers[grepl("(103/µl)", Data_all_raw$Markers)] <- gsub(" .*", "", Data_all_raw$Markers[grepl("(103/µl)", Data_all_raw$Markers)])

Data_all_raw$Type <- na.locf(Data_all_raw$Type)

name_types <- unique(Data_all_raw$Type)
name_myeloids <- Data_all_raw$Markers[Data_all_raw$Type == "Myeloid"]
name_lymphoids <- Data_all_raw$Markers[Data_all_raw$Type == "Lymphoid"]
name_cytokines <- Data_all_raw$Markers[Data_all_raw$Type == "Cytokines"]
name_acids <- Data_all_raw$Markers[Data_all_raw$Type == "Esterified Fatty Acid"]
name_emocromos <- Data_all_raw$Markers[Data_all_raw$Type == "Emocromo"]
name_calprotect <- "Fecal calprotectin"
name_all_pops <- c(name_myeloids, name_lymphoids, name_cytokines, name_acids, name_emocromos, name_calprotect)
length(name_all_pops) == nrow(Data_all_raw)

for (i in 3:92) {
   Data_all_raw[,i ] <- as.numeric(Data_all_raw[,i ])
}


data <- as.data.frame(t(Data_all_raw[,-1]))
colnames(data) <- data[1,]
data <- data[-1,]
data$Id <- rownames(data)
data <- data[,c(ncol(data),1:(ncol(data)-1))]
data$Time <- NA
rownames(data) <- NULL
data$Time[grepl("_PRE", data$Id)] <- name_timepoints[1]
data$Time[grepl("_3M", data$Id)] <- name_timepoints[2]
data$Time[grepl("_6M", data$Id)] <- name_timepoints[3]
data$Id <- sub("_.*", "", data$Id)
data <- data[,c(1, ncol(data),2:(ncol(data)-1))]
colnames(data)[2] <- "Time_factor"
data$Time <- NA
data$Time[data$Time_factor == "PRE"] <- 0
data$Time[data$Time_factor == "3M"] <- 3
data$Time[data$Time_factor == "6M"] <- 6
data <- data[,c(1, 121, 3:120)]



##########
# Time series file
##########

data_cytokines <- data[, c("Id", "Time", name_cytokines)]

for (i in name_cytokines) {
   data_cytokines[,i] <- as.numeric(data_cytokines[,i])
}

## Removing all IDs with missing values for Friedman test
tmp <- rowSums(is.na(data[, name_cytokines]))
Id_removes <- data$Id[tmp>0]

data_cytokines_tmp <- subset(data_cytokines, !(data_cytokines$Id %in% Id_removes))    

## Scaling values between -1 and 1
min = -1
max = 1
data_cyt_norm_tmp = data_cytokines[,-c(1,2)]
data_cyt_norm_tmp = as.data.frame(apply(data_cyt_norm_tmp, 2, rescale, c(min,max)))

data_cyt_norm_tmp = cbind(data_cytokines[,c(1,2)],data_cyt_norm_tmp)

## Normalizing baseline by computing delta values

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


data_cyt_norm = cbind(data_cytokines[,c(1,2)],data_cyt_norm)
rm(data_cyt_norm_tmp)




# name_ptzs <- unique(data_cytokines$Id)

## Friedman test to select only curves that differ significantly 

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

## Creation of new dataset in long format for CONNECTOR

## Adding variable name to ID for CONNECTOR

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

## Checking
# data_cytokines_conn$Observations[data_cytokines_conn$ID == "FAP24_IL-18" & data_cytokines_conn$Time == 3] ==
#   data_cyt_norm$`IL-18`[data_cyt_norm$Id == "FAP24" & data_cyt_norm$Time == 3]


data_cytokines_conn$Observations <- as.numeric(data_cytokines_conn$Observations)

## Removing only the variables containing a missing value for any of the IDs

vet_removeID <- data_cytokines_conn$ID[is.na(data_cytokines_conn$Observations)]

data_cytokines_conn <- subset(data_cytokines_conn, !(data_cytokines_conn$ID %in% vet_removeID))

##########
# Annotation file
##########

cytokines_annotation <- as.data.frame(matrix(NA, nrow = nrow(data_cytokines_conn)/length(name_timepoints_num), ncol = 3))
colnames(cytokines_annotation) <- c("ID", "IdSample", "Cytokine")
cytokines_annotation$ID <- unique(data_cytokines_conn$ID)
cytokines_annotation$IdSample <- gsub("_.*", "", cytokines_annotation$ID)
cytokines_annotation$Cytokine <- gsub(".*_", "", cytokines_annotation$ID)


Cytokine_List <- connector::DataFrameImport(TimeSeriesDataFrame = data_cytokines_conn, AnnotationFrame = cytokines_annotation)

CurvesPlot <- connector::PlotTimeSeries(data = Cytokine_List, feature = "Cytokine")
CurvesPlot