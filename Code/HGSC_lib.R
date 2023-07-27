libs4HGSC<-c("beanplot","cowplot","Seurat","EBImage","survival","rms","mixtools","MASS","ggplot2",
             "nnet","ppcor","ROCR","tsne","gplots","ggpubr","EnhancedVolcano","plyr","reshape2","plotrix","stats", "Matrix",
             "Rtsne","lmerTest","devtools","gplots","heatmap3","e1071","openxlsx","RColorBrewer","heatmap3","UpSetR")

v<-lapply(libs4HGSC,function(x) library(x,character.only = T))

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
