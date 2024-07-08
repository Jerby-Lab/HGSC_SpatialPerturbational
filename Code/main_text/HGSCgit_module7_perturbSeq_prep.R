HGSC_Fig7e_prep<-function(rslts,sig){
  # select out the significant data 
  sig1<-get.top.elements(rslts$sum$sum[,1:2],7,min.ci = 0.01)
  p.hot<-c(paste0(sig1$Repressor,".up"),paste0(sig1$Activators,".down"))
  genes<-get.abundant(unlist(intersect.list1(rslts$sig[p.hot],sig$Mtil.up)),5)
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
  X.lab<-cbind.data.frame(Mtil = ifelse(b1,"Repressor","Activator"),
                          Condition = get.strsplit(prtA,".",1))
  row.names(X.lab) <- colnames(X1)
  rslts<-list(DEGs = X1,labels = X.lab)
  saveRDS(rslts,file = get.file("Results/SourceData_fig7e.rds"))
  return()
  
}