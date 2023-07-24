#### Results Section 5 ###
# Figure 5. Copy Number Alterations (CNAs) mapping to mTIL and TIL levels in SMI spatial data and TCGA.

HGSC_Figure5_CNAs<-function(r1,rslts1,rslts2){
  
  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_mTIL_CNA.rds"))
    rslts1<-readRDS(get.file("Results/HGSC_CNAs.vs.mTIL_SMI.rds"))
    rslts2<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  }
  
  # HGSC_Fig5A(r1 = r1)
  HGSC_Fig5B(r1 = r1,rslts = rslts1)
  HGSC_Fig5C(rslts = rslts1)
  HGSC_Fig5E(rslts = rslts2)
  
  return()
}

HGSC_Fig5B<-function(r1,rslts){
  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_mTIL_CNA.rds"))
    rslts<-readRDS(get.file("Results/HGSC_mTIL_vsCNAs.rds"))
  }
  
  f1<-function(x){
    idx<-order(r1$cnv[,x])
    call.boxplot(r1$scores[idx,"mTIL.up"],r1$cnv[idx,x],order.flag = F,
               xlab = paste(x,"Copy Number"),ylab = "mTIL program (OE)",
               main = paste(x,"(",call.format.pval(rslts$HLM[x,"mTIL.up.P"]),")"))
  }
  
  l1<-lapply(c("TCF7L2","IFNGR2","AXL","IFNAR1","ACTA2","RUNX1"),f1)
  
  pdf(get.file("Figures/Fig5B.pdf"))
  call.multiplot(l1,nplots = 6,cols = 3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

HGSC_Fig5C<-function(rslts){
  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_mTIL_vsCNAs.rds"))
  }
  
  z<-rslts$Z[,"mTIL.up.Z"]
  names(z)<-rownames(rslts$Z)
  z<-z[p.adjust(10^-abs(z),method = "BH")<0.1]
  z<-sort(z)
  
  pdf(get.file("Figures/Fig5C.pdf"))
  barplot(abs(z),las=2,cex.names = 0.5,xlab = "Genes CNA",
          col = ifelse(z<0,"lightblue","darkred"),
          ylab = "Z-score",main = "CNAs correlated with the mTIL program")
  legend(legend = c("Negative","Positive"), pch = 15, 
         "topright",col = c("lightblue","darkred"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig5E<-function(rslts){
  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  }
  
  Xd<-rslts$plot.del
  Xa<-rslts$plot.amp
  pdf(get.file("Figures/Fig5E.pdf"))
  violin.split(Xd$TIL,Xd$Del,conditions = Xd$Gene,
               ylab = "TIL levels",xlab = "mTIL UP Genes",
               main = "TCGA, mTIL CNA vs. TIL levels",cex.axis = 0.5)
  violin.split(Xa$TIL,Xa$Amp,conditions = Xa$Gene,col1 = "darkred",
               ylab = "TIL levels",xlab = "mTIL DOWN Genes",
               main = "TCGA, mTIL CNA vs. TIL levels")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

