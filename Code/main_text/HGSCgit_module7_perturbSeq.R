#### Results Section 7 ###
# Figure 7. Perturb-Seq meta-analyses

# Figure 7. High content CRISPR screen in ovarian cancer cells in monoculture 
# and co-culture with targeting NK and T cells.
# 7a. Overview of experimental design (made outside of R).
# 7b. CRISPR-based differential fitness.
# 7c. Perturb-Seq UMAP.
# 7d. Perturb-Seq-based Mtil regulators
# 7e. Mtil DEGs under different perturbations.
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
  if(missing(r)){
    r<-readRDS(get.file("Data/PerturbSeq_TYKnu.rds"))
    rslts1<-readRDS(get.file("Results/PerturbSeq_TYKnu_DEP.rds"))
    fitnessR<-readRDS(get.file("Results/SourceData_Fig7b.rds"))
    sig<-readRDS(get.file("Results/ST_Mtil.sig.rds"))
  }

  # wrangle data for appropriate visualization 
  r$scores1c<-cap.mat(r$scores1,cap = 0.05,MARGIN = 2)
  r$condL<-r$cond
  r$condL[r$cond=="Mono"]<-"Monoculture"
  r$condL[r$cond=="1to1"]<-"Co-culture, 1:1"
  r$condL[r$cond=="25to1"]<-"Co-culture, 2.5:1"

  # Plot 7b. CRISPR-based differential fitness.
  HGSC_Fig7b(fitnessR)
  # Plot 7c. Perturb-Seq UMAP.
  HGSC_Fig7c(r = r)
  # Plot 7d. Perturb-Seq-based Mtil regulators
  HGSC_Fig7d(rslts = rslts)
  # Plot 7e. Mtil DEGs under different perturbations.
  HGSC_Fig7e()
  # Plot 7f. KOs altering the response to NK cells.
  HGSC_Fig7f(rslts = rslts)
  # Plot 7g. Differential expression of KO signatures in monoculture and co-culture
  HGSC_Fig7g(r = r,rslts = rslts)
  # Plot 7h. Perturb-Seq UMAPs with cell colored based on gene KO signatures.
  HGSC_Fig7h(r = r)
  return()
}

#' Figure 7b. CRISPR-based differential fitness.
#' This is 2D volcano plot where the x and y axis confer fitness z-scores for
#' NK and CD8 T cells, respectively. Statistically significant genes and genes
#' of interest are highlighted. 
#' @param fitnessR results of differential fitness
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7b<-function(fitnessR){
  # load data if missing 
  if(missing(fitnessR)){
    fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  }
  
  # plot to disk 
  pdf(get.file("Figures/Fig7b.pdf"))
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
#' UMAP embedding of cells from perturb-seq experiment in ovarian cancer cells
#' in mono-culture and co-culture with NK cells. 
#' @param r list object containing the Perturb-Seq data with NK co-culture 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7c<-function(r){
  pdf(get.file("Figures/Fig7c.pdf"))
  print(umap.ggplot(r$umap,labels = add.n.of.samples(r$condL)))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7d. Perturb-Seq-based Mtil regulators
#' Barplot of impact of perturbation on Mtil expression in comparison to non-
#' targeting controls.  
#' @param rslts list object containing statistical results from th
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig7d<-function(rslts){
  # wrangle results for plotting 
  P<-rslts$prt_Mtil_reg
  P<-P[P$N<2,]
  z<-P$Z;names(z)<-rownames(P)
  z[abs(z)>8]<-8*sign(z[abs(z)>8])
  
  # plot to disk 
  pdf(get.file("Figures/Fig7d.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "Mtil differential expression",
          xlab = "Gene KO",main = "Perturbation impact on Mtil program",
          col = ifelse(z>2,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Activator","Repressor","NS"), pch = 15,
         col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}


#' Figure 7e. Mtil DEGs under different perturbations.
#' Heatmap of significant gene perturbational effect on Mtil expression. Results
#' of differentially expressed genes between condition
#' @param rslts results of differentially expressed genes. 
#' @param sig Mtil signature genes to add to visualization 
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7e<-function(){
  rslts<-readRDS(get.file("Results/SourceData_Fig7e.rds"))
  pdf(get.file("Figures/Fig7e.pdf"))
  blue_white_red <- colorRampPalette(c("dodgerblue4", "white", "red3"))(100)
  pheatmap(t(rslts$DEGs), cluster_cols = T, cluster_rows = F, color = blue_white_red, treeheight_col = 0,
           annotation_row = rslts$labels, annotation_colors = 
           list(Mtil = c(Repressor = "darkred", Activator = "dodgerblue4"),
             Condition = c(Co2 = "gold", Co1 = "purple4", Mono = "turquoise4") 
           ), filename = get.file("Figures/Fig7E.pdf"), width = 6, height = 8, 
           labels_row = unlist(lapply(colnames(rslts$DEGs), 
                                      function(x) {strsplit(x, split = "\\.")[[1]][2]})))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()

}

#' Figure 7f. DEGs under different culture conditions
#' Barplots of differential expression of KO signatures in monoculture and 
#' co-culture KOs altering the response to NK cells.
#' @param rslts results of differentially expressed genes.  
#' @return this function returns nothing, but prints a plot to device.
HGSC_Fig7f<-function(rslts){
  # fetch differentially expressed genes
  z<-rslts$prt_NKresponse_reg[,"zscores"]
  names(z)<-rownames(rslts$prt_NKresponse_reg)
  z<-sort(z)
  z[abs(z)>30]<-30*sign(z[abs(z)>30])
  
  # plot to disk
  pdf(get.file("Figures/Fig7f.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "Differential expression,\nmono- vs. co-culture",
          xlab = "Perturbation signature",main = "",cex.axis = 0.1,
          col = ifelse(z>4,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Higher","Lower","NS"), pch = 15,col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 7g. KOs altering the response to NK cells.
#' Beanplots of of gene KOs that mimick response to NK cells 
#' @param rslts results of differentially expressed genes.  
#' @return Fig7g.pdf in Figures directory.
HGSC_Fig7g<-function(r,rslts){
  sig1<-get.top.elements(rslts$prt_NKresponse_reg[,c("BH.more","BH.less")],7,min.ci = 0.01)
  par(mfrow=c(2,3),oma = c(2, 1, 0, 5),xpd = T)
  X<-NULL
  for(x in unlist(sig1)){
    X1<-cbind.data.frame(scores = r$scores1[r$b.ctrl,x],cond = r$cond[r$b.ctrl],
                         prt = x)
    X<-rbind(X,X1)
  }
  X1<-X[is.element(X$prt,sig1$BH.more),]
  X2<-X[!is.element(X$prt,sig1$BH.more),]

  pdf(get.file("Figures/Fig7g.pdf"))
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
#' UMAPs of select gene KO perturbational signatures. 
#' @param r HGSC Perturb-seq data. 
#' @return Fig7h.pdf in Figures directory.
HGSC_Fig7h<-function(r){
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
