#### Results Section 5 ###
# Figure 5. Perturb-Seq meta-analyses
# Supplementary Table 6a. mTIL program regulators.

# Figure 5a. Perturb-seq hits in K5621
# Figure 5b. Perturb-seq hits in RPE1
# Figure 5c. UMAPs of Perturb-seq hits in K5621
# Figure 5d. UMAPs of Perturb-seq hits in RPE1

HGSC_Figure5_perturbMeta<-function(rslts,rslts1,rslts2,rslts3){

  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds"))
    rslts1<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRa.rds"))
    rslts2<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRi.rds"))
    rslts3<-readRDS(get.file("Results/HGSC_mTIL_RPE1_CRISPRi.rds"))
    HGSC_Figure6.regenerate_perturbSeqMetaA(rslts,rslts1,rslts2,rslts3)
  }

  pdf(get.file("Figures/Fig5a.pdf"))
  perturbSeq_sig.reg.barplot(rslts1)
  perturbSeq_sig.reg.barplot(rslts2)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  pdf(get.file("Figures/Fig5b.pdf"))
  perturbSeq_sig.reg.barplot(rslts3,0.2)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  pdf(get.file("Figures/Fig5c.pdf"))
  perturbSeq_OE.umaps(rslts = rslts1$OE,prt = c("IRF1","CEBPA","CEBPE"),name = rslts1$name)
  perturbSeq_OE.umaps(rslts = rslts2$OE,prt = c("PTPN1","KDM1A","CHMP6"),name = rslts2$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  pdf(get.file("Figures/Fig5d.pdf"))
  perturbSeq_OE.umaps(rslts = rslts3$OE,prt = c("PSMB5","THAP1"),name = rslts3$name)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

}

perturbSeq_sig.reg.barplot<-function(rslts,cex.names = 0.5){
  X<-rslts$OE$sum[order(rslts$OE$sum[,1],decreasing = T),]
  par(mfrow=c(1,1),oma = c(10, 5, 0, 0),xpd = F)
  b<-abs(X[,1])>5|X[,1]<(-3)
  barplot(height = t(abs(subset(X[b,],select = "OE.ttest.zscore"))),beside = T,las = 2,
          cex.names = cex.names,ylab = "-log10(p-value)",
          main = gsub("_"," ",rslts$name),col = ifelse(X[b,1]>0,"darkgrey","lightblue"))
  legend("top",legend = c("Activate","Repress"), pch = 15,
         col = c("darkgrey","lightblue"))

}

perturbSeq_OE.umaps<-function(rslts,prt,fileName,name){
  X1<-rslts$sum[prt,]
  prt<-prt[order(-abs(X1[,1]))]
  if(missing(prt)){prt<-unlist(rslts$OE.sig)}
  names(prt)<-NULL

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

  if(!missing(fileName)){pdf(fileName)}
  call.multiplot(l2,cols = 3,nplots = 6)
  if(!missing(fileName)){dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)}
  return()
}
