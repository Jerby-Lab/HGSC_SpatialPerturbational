#### Results Section 6 ###
# Figure 6. Perturb-Seq meta-analyses

# Figure 6a. Perturb-seq hits in K562
# Figure 6b. Perturb-seq hits in RPE1
# Figure 6c. UMAPs of Perturb-seq hits in K5621
# Figure 6d. UMAPs of Perturb-seq hits in RPE1

#' Figure 6 Wrapper Function
#' This function calls code to reproduce main text Figures 6a-d.  
#' @param rslts1 results from mapping MTIL to K562 CRISPRa data 
#' @param rslts2 results from mapping MTIL to K562 CRISPRi data
#' @param rslts3 results from mapping MTIL to RPE1 CRISPRi data
#' @return Figure 3 files in Figures folder. 
HGSC_Figure6_perturbMeta<-function(rslts1,rslts2,rslts3){
  # load missing results 
  if(missing(rslts1)){
    rslts1<-readRDS(get.file("Results/PerturbSeq_Mtil_K562_CRISPRa.rds"))
    rslts2<-readRDS(get.file("Results/PerturbSeq_Mtil_K562_CRISPRi.rds"))
    rslts3<-readRDS(get.file("Results/PerturbSeq_Mtil_RPE1_CRISPRi.rds"))
  }
  
  #1. Reproduce Figure 6a: Perturb-seq hits in K562 cell line
  pdf(get.file("Figures/Fig6a.pdf"),width = 15,height = 8)
  perturbSeq_sig.reg.barplot(rslts1)
  perturbSeq_sig.reg.barplot(rslts2)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  #2. Reproduce Figure 6b: Perturb-seq hits in RPE1 cell line
  pdf(get.file("Figures/Fig6b.pdf"),width = 15,height = 8)
  perturbSeq_sig.reg.barplot(rslts3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  #3. UMAPs of select Perturb-seq hits in K562 cell line
  pdf(get.file("Figures/Fig6c.pdf"),width = 12,height = 8)
  perturbSeq_OE.umaps(rslts = rslts1,prt = c("IRF1","CEBPA","CEBPE"),name = rslts1$name)
  perturbSeq_OE.umaps(rslts = rslts2,prt = c("PTPN1","KDM1A","CHMP6"),name = rslts2$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  #4. UMAPs of select Perturb-seq hits in RPE1 cell line 
  pdf(get.file("Figures/Fig6d.pdf"),width = 12,height = 8)
  perturbSeq_OE.umaps(rslts = rslts3,prt = c("PSMB5","THAP1"),name = rslts3$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

#' Plot Barplots of Perturb-seq meta-analysis results. 
#' Here, the function uses the results from comparing MTIL overall expression 
#' in single cells in perturbed vs. non-targeting controls to make barplots.
#' @param rslts results from mapping MTIL to Perturb-seq data
#' @cex.names (default = 0.7) parameter for controlling size of axis labels. 
#' @return prints a plot to device. 
perturbSeq_sig.reg.barplot<-function(rslts,cex.names = 0.7){
  # sort result by significance
  X<-rslts$sum[order(rslts$sum[,1],decreasing = T),]
  
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
#' UMAP embeddings of perturb-seq experiments for select gene targets and non-
#' targeting control cells. 
#' @param rslts results from mapping MTIL to Perturb-seq data
#' @param prt selected gene names
#' @param name the name of the experiment
#' @return plots to device.  
perturbSeq_OE.umaps<-function(rslts,prt,name){
  # setup results
  X1<-rslts$sum[prt,]
  prt<-prt[order(-abs(X1[,1]))]
  
  # make subplots
  l2<-NULL
  for(x in prt){
    r2<-rslts[[x]]
    r2$targets[r2$targets!=x]<-"NTC"
    l1<-list(umap.ggplot(r2$umap,labels = add.n.of.samples(r2$targets,sep = " "),reorder.flag = F,
                         size = 0.5,main = paste0(x," Z = ",round(X1[x,"OE.ttest.zscore"],2)),
                         labels.name = "sgRNA")+theme(legend.position="bottom"),
             umap.ggplot(r2$umap,labels = r2$scores1[,"hot100"],size = 0.5,main = name,
                         labels.name = "Program OE")+theme(legend.position="bottom"))
    if(sort(c(x,"NTC"))[1]=="NTC"){
      l1[[1]]<-l1[[1]]+scale_color_manual(values = c("#00BFC4","#F8766D"))
    }
    l2<-c(l2,l1)
  }
  
  # plot
  call.multiplot(l2,cols = 3,nplots = 6)
  return()
}

