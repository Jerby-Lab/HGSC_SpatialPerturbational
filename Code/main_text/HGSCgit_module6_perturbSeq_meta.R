#### Results Section 6 ###
# Figure 6. Perturb-Seq meta-analyses

# Figure 6a. Perturb-seq hits in K562
# Figure 6b. Perturb-seq hits in RPE1
# Figure 6c. UMAPs of Perturb-seq hits in K5621
# Figure 6d. UMAPs of Perturb-seq hits in RPE1

#' Figure 6 Wrapper Function
#'
#' This function calls code to reproduce main text Figures 6a-d.  
#'
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure6_perturbMeta<-function(rslts,rslts1,rslts2,rslts3){
  # load missing results 
  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds"))
    rslts1<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRa.rds"))
    rslts2<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRi.rds"))
    rslts3<-readRDS(get.file("Results/HGSC_mTIL_RPE1_CRISPRi.rds"))
    HGSC_Figure6.regenerate_perturbSeqMetaA(rslts,rslts1,rslts2,rslts3)
  }
  
  #1. Reproduce Figure 6a: Perturb-seq hits in K562 cell line
  pdf(get.file("Figures/Fig6A.pdf"))
  perturbSeq_sig.reg.barplot(rslts1)
  perturbSeq_sig.reg.barplot(rslts2)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  #2. Reproduce Figure 6b: Perturb-seq hits in RPE1 cell line
  pdf(get.file("Figures/Fig6B.pdf"))
  perturbSeq_sig.reg.barplot(rslts3,0.2)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  #3. UMAPs of select Perturb-seq hits in K562 cell line
  pdf(get.file("Figures/Fig6C.pdf"))
  perturbSeq_OE.umaps(rslts = rslts1$OE,prt = c("IRF1","CEBPA","CEBPE"),name = rslts1$name)
  perturbSeq_OE.umaps(rslts = rslts2$OE,prt = c("PTPN1","KDM1A","CHMP6"),name = rslts2$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  #4. UMAPs of select Perturb-seq hits in RPE1 cell line 
  pdf(get.file("Figures/Fig6D.pdf"))
  perturbSeq_OE.umaps(rslts = rslts3$OE,prt = c("PSMB5","THAP1"),name = rslts3$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

#' Plot Barplots of Perturb-seq meta-analysis results. 
#'
#' Here, the function uses the results from comparing MTIL overall expression 
#' in single cells in perturbed vs. non-targeting controls to make barplots.
#' 
#' @param rslts list object of qunatitative results comparing MTIL OE with 
#' perturbational status in a perturb-seq experiment 
#' @cex.names (default = 0.5) parameter for controlling size of axis labels. 
#' @return this function returns nothing, but prints a plot to device. 
perturbSeq_sig.reg.barplot<-function(rslts,cex.names = 0.5){
  # sort result by significance
  X<-rslts$OE$sum[order(rslts$OE$sum[,1],decreasing = T),]
  
  # make plot
  par(mfrow=c(1,1),oma = c(10, 5, 0, 0),xpd = F)
  b<-abs(X[,1])>5|X[,1]<(-3)
  barplot(height = t(abs(subset(X[b,],select = "OE.ttest.zscore"))),beside = T,las = 2,
          cex.names = cex.names,ylab = "-log10(p-value)",
          main = gsub("_"," ",rslts$name),col = ifelse(X[b,1]>0,"darkgrey","lightblue"))
  legend("top",legend = c("Activate","Repress"), pch = 15,
         col = c("darkgrey","lightblue"))
}

#' Plot UMAP embeddings of cells in Perturb-seq Meta-analysis 
#'
#' UMAP embeddings of perturb-seq experiments for select gene targets and non-
#' targeting control cells. 
#' 
#' @param rslts list object of qunatitative results comparing MTIL OE with 
#' perturbational status in a perturb-seq experiment 
#' @param prt selected gene names in the form of a list 
#' @param fileName path to desired file name 
#' @param name string corresponding to the name of the experiment
#' @return this function returns nothing, but prints a plot to file if fileName
#' provided, or to device otherwise.  
perturbSeq_OE.umaps<-function(rslts,prt,fileName,name){
  # setup results
  X1<-rslts$sum[prt,]
  prt<-prt[order(-abs(X1[,1]))]
  if(missing(prt)){prt<-unlist(rslts$OE.sig)}
  names(prt)<-NULL
  
  # make subplots
  l2<-NULL
  for(x in prt){
    r2<-rslts[[x]]
    r2$targets[r2$targets!=x]<-"zNTC"
    l1<-list(umap.ggplot(r2$umap,labels = add.n.of.samples(r2$targets,sep = " "),
                         size = 0.5,main = paste0(x," Z = ",round(X1[x,"OE.ttest.zscore"],2)),
                         labels.name = "sgRNA")+theme(legend.position="bottom"),
             umap.ggplot(r2$umap,labels = r2$scores1[,"hot100"],size = 0.5,main = name,
                         labels.name = "Program OE")+theme(legend.position="bottom"))
    l2<-c(l2,l1)
  }
  
  # plot
  if(!missing(fileName)){pdf(fileName)}
  call.multiplot(l2,cols = 3,nplots = 6)
  if(!missing(fileName)){dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)}
  return()
}

