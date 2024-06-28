#### Results Section 7 ###
# Figure 7. Perturb-Seq meta-analyses

# Figure 7. High content CRISPR screen in ovarian cancer cells in monoculture 
# and co-culture with targeting NK and T cells.
# 7a. Overview of experimental design (made outside of R).
# 7b. CRISPR-based differential fitness.
# 7c. Perturb-Seq UMAP.
# 7d. Perturb-Seq-based mTIL regulators
# 7e. mTIL DEGs under different perturbations.
# 7f. Differential expression of KO signatures in monoculture and co-culture
# 7g. KOs altering the response to NK cells.
# 7h. Perturb-Seq UMAPs with cell colored based on gene KO signatures.

#' Figure 7 Wrapper Function
#'
#' This function calls code to reproduce main text Figures 7b-h 
#'
#' @param r list object containing the Perturb-Seq data with NK co-culture 
#' @param rslts differentially expressed genes results for Perturb-seq 
#' @param sig gene signature 
#' @param fitnessR fitness scores associated with CRISPR guide sequencing data
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure7_perturbOC<-function(r,rslts,sig,fitnessR){
  # load data if missing
  if(missing(r)){
    r<-readRDS("Data/PerturbSeq_TYKnuNK.rds")
    rslts<-readRDS("Results/PerturbSeq_TYKnuNK_DEGs.rds")
    sig<-readRDS("Results/HGSC_mTIL.rds")
    fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  }

  # wrangle data for appropriate visualization 
  r$scores1c<-cap.mat(r$scores1,cap = 0.05,MARGIN = 2)
  r$condL<-r$cond
  r$condL[r$cond=="Mono"]<-"Monoculture"
  r$condL[r$cond=="1to1"]<-"Co-culture, 1:1"
  r$condL[r$cond=="25to1"]<-"Co-culture, 2.5:1"

  # Plot 7b. CRISPR-based differential fitness.
  HGSC_Fig7B(fitnessR)
  # Plot 7c. Perturb-Seq UMAP.
  HGSC_Fig7C(r = r)
  # Plot 7d. Perturb-Seq-based mTIL regulators
  HGSC_Fig7D(rslts = rslts)
  # Plot 7e. mTIL DEGs under different perturbations.
  HGSC_Fig7E(rslts = rslts,sig = sig)
  # Plot 7f. KOs altering the response to NK cells.
  HGSC_Fig7F(rslts = rslts)
  # Plot 7g. Differential expression of KO signatures in monoculture and co-culture
  HGSC_Fig7G(r = r,rslts = rslts)
  # Plot 7h. Perturb-Seq UMAPs with cell colored based on gene KO signatures.
  HGSC_Fig7H(r = r)
  return()
}

#' Figure 7b. CRISPR-based differential fitness.
#'
#' This is 2D volcano plot where the x and y axis confer fitness z-scores for
#' NK and CD8 T cells, respectively. Statistically significant genes and genes
#' of interest are highlighted. 
#' 
#' @param fitnessR results of differential fitness
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7B<-function(fitnessR){
  # load data if missing 
  if(missing(fitnessR)){
    fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  }
  
  # plot to disk 
  pdf(get.file("Figures/Fig7B.pdf"))
  p<-create_volcano_plot(x = fitnessR$Tcell.z, y = fitnessR$NK.z,
                         dot_names = rownames(fitnessR), 
                         x_label = "T cell selection", 
                         y_label = "NK selection",
                         quadrant_labels = c("Resistance", 
                                             "T cell response\nNK resistance",
                                             "Response", 
                                             "NK response\nT cell resistance"),
                         zcat = 2)
  print(p)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7c. Perturb-Seq UMAP.
#'
#' UMAP embedding of cells from perturb-seq experiment in ovarian cancer cells
#' in mono-culture and co-culture with NK cells. 
#' 
#' @param r list object containing the Perturb-Seq data with NK co-culture 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7C<-function(r){
  # plot to disk 
  pdf(get.file("Figures/Fig7C.pdf"))
  print(umap.ggplot(r$umap,labels = add.n.of.samples(r$condL)))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7d. Perturb-Seq-based MTIL regulators
#'
#' Barplot of impact of perturbation on MTIL expression in comparison to non-
#' targeting controls.  
#' 
#' @param rslts list object containing statistical results from th
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7D<-function(rslts){
  # wrangle results for plotting 
  P<-rslts$sum$sum
  P<-P[P$N<2,]
  z<-P$Z;names(z)<-rownames(P)
  z[abs(z)>8]<-8*sign(z[abs(z)>8])
  
  # plot to disk 
  pdf(get.file("Figures/Fig7D.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "mTIL differential expression",
          xlab = "Gene KO",main = "Perturbation impact on mTIL program",
          col = ifelse(z>2,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Activator","Repressor","NS"), pch = 15,
         col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}


#' Figure 7e. MTIL DEGs under different perturbations.
#'
#' Heatmap of significant gene perturbational effect on MTIL expression. Results
#' of differentially expressed genes between condition
#' 
#' @param rslts results of differentially expressed genes. 
#' @param sig MTIL signature genes to add to visualization 
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7E<-function(rslts,sig){
  # select out the significant data 
  sig1<-get.top.elements(rslts$sum$sum[,1:2],7,min.ci = 0.01)
  p.hot<-c(paste0(sig1$Repressor,".up"),paste0(sig1$Activators,".down"))
  genes<-get.abundant(unlist(intersect.list1(rslts$sig[p.hot],sig$mTIL.up)),5)
  ids<-c(unique(get.strsplit(p.hot[grepl("up",p.hot)],".",1)),"NTC",unique(get.strsplit(p.hot[grepl("down",p.hot)],".",1)))
  prt<-unlist(sig1)
  
  # make matrix for heatmap
  X<-cbind.data.frame(Mono = rslts$Mono$deg.ttest[match(genes,rownames(rslts$Mono$deg.ttest)),prt],
                      Co1 = rslts$Co1$deg.ttest[genes,prt],
                      Co2 = rslts$Co2$deg.ttest[genes,prt])
  rownames(X)<-genes
  X[is.na(X)]<-0
  X[abs(X)>4]<-4*sign(X[abs(X)>4])
  b<-is.element(get.strsplit(colnames(X),".",2),sig1$Repressor)

  # hierarchical clustering for visualization
  hc1 <- hclust(dist(t(X[,b]),method = 'euclidean'), method="complete")
  hc2 <- hclust(dist(t(X[,!b]),method = 'euclidean'), method="complete")
  prt1<-hc1$labels[hc1$order]
  prt1<-c(prt1[grepl("Co",prt1)],prt1[!grepl("Co",prt1)])
  prt2<-hc2$labels[hc2$order]
  prt2<-c(prt2[grepl("Co",prt2)],prt2[!grepl("Co",prt2)])
  prtA<-c(prt1,prt2)
  X1<-X[,prtA]
  b1<-is.element(get.strsplit(prtA,".",2),sig1$Repressor)
  X.lab<-cbind.data.frame(mTIL = ifelse(b1,"Repressor","Activator"),
                          Condition = get.strsplit(prtA,".",1))
  
  # plot to disk 
  pdf(get.file("Figures/Fig7E.pdf"))
  blue_white_red <- colorRampPalette(c("dodgerblue4", "white", "red3"))(100)
  row.names(X.lab) <- colnames(X1)
  pheatmap(t(X1), cluster_cols = T, cluster_rows = F, color = blue_white_red, treeheight_col = 0,
           annotation_row = X.lab, annotation_colors = 
           list(mTIL = c(Repressor = "darkred", Activator = "dodgerblue4"),
             Condition = c(Co2 = "gold", Co1 = "purple4", Mono = "turquoise4") 
           ), filename = get.file("Figures/Fig7E.pdf"), width = 6, height = 8, 
           labels_row = unlist(lapply(colnames(X1), 
                                      function(x) {strsplit(x, split = "\\.")[[1]][2]})))
  return()

}

#' Figure 7f. DEGs under different culture conditions
#'
#' Barplots of differential expression of KO signatures in monoculture and 
#' co-culture KOs altering the response to NK cells.
#' 
#' @param rslts results of differentially expressed genes.  
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7F<-function(rslts){
  # fetch differentially expressed genes
  z<-rslts$prt.deNK1[,"zscores"]
  names(z)<-rownames(rslts$prt.deNK1)
  z<-sort(z)
  z[abs(z)>30]<-30*sign(z[abs(z)>30])
  
  # plot to disk
  pdf(get.file("Figures/Fig7F.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "Differential expression,\nmono- vs. co-culture",
          xlab = "Perturbation signature",main = "",cex.axis = 0.1,
          col = ifelse(z>4,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Higher","Lower","NS"), pch = 15,col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7g. KOs altering the response to NK cells.
#'
#' Beanplots of of gene KOs that mimick response to NK cells 
#' 
#' @param rslts results of differentially expressed genes.  
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7G<-function(r,rslts){
  sig1<-get.top.elements(rslts$prt.deNK1[,c("BH.more","BH.less")],7,min.ci = 0.01)
  par(mfrow=c(2,3),oma = c(2, 1, 0, 5),xpd = T)
  X<-NULL
  for(x in unlist(sig1)){
    X1<-cbind.data.frame(scores = r$scores1[r$b.ctrl,x],cond = r$cond1[r$b.ctrl],
                         prt = x)
    X<-rbind(X,X1)
  }
  X1<-X[is.element(X$prt,sig1$BH.more),]
  X2<-X[!is.element(X$prt,sig1$BH.more),]

  pdf(get.file("Figures/Fig7G.pdf"))
  violin.split(scores = X1$scores, 
               treatment = ifelse(X1$cond=="Mono","Monoculture","Co-culture"),
               conditions = X1$prt,
               ylab = "Pertrbation signature OE",
               xlab = "Perturbation",
               main = "Genetic pertrubation mimicking NK effects")
  violin.split(scores = X2$scores, 
               treatment = ifelse(X2$cond=="Mono","Monoculture","Co-culture"),
               conditions = X2$prt,
               ylab = "Pertrbation signature OE",
               xlab = "Perturbation",
               main = "Genetic pertrubation counteracting NK effects")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7h. Perturb-Seq UMAPs with cell colored based on gene KO signatures.
#'
#' UMAPs of select gene KO perturbational signatures. 
#' 
#' @param r list object of Perturb-seq data. 
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7H<-function(r){
  # boxplots
  f1<-function(x){
    call.boxplot(r$scores1[,x],paste0(r$targets,r$cond1),cex = 0.5,
                 labels = paste0(r$cond1,
                                 ifelse(r$targets==x,paste(", ",x,"KO"),"")),
                 xlab = "Perturbations",
                 ylab = paste0(x," signature OE"),
                 main = x)}
  
  # 
  X<-cbind.data.frame(Conditions = r$condL,
                      r$scores1c[,c("ACTR8.up","IRF1","MED12.up","STAT1")])
  colnames(X)<-gsub(".up","",colnames(X))
  
  # plot to disk
  pdf(get.file("Figures/Fig7H.pdf"))
  print(umap.ggplot(r$umap,labels = X[,-1],size = 0.1))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}
