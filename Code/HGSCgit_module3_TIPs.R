#### Results Section 3 ###
# Figure 3. T/NK cell states reflect their tumor infiltration status.
# Table S5. Immune tumor infiltration programs (TIPs) derived from SMI data. Associated with Figure 3. 

# Figure 3A. CD8 T cell umaps
# Figure 3B. CD8 T cell TIP genes' association with tumor infiltration in different cell subtypes.
# Figure 3C. TIPs overlaps (barplots)
# Figure 3D. Xenium spatial maps with CD8 TIP
# Figure 3E. Ligand-receptor interactions (generated via CytoScape)

HGSC_Figure3_TIP<-function(r1.smi,r.xenium,r1.xenium,R){
  if(missing(r1.smi)){
    r1.smi<-HGSC_SMI.process.CD8()
    r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
    r1.xenium<-HGSC_Xenium.process.CD8.NK(r = r.xenium)
    R<-TIP_find_all()
    HGSC_Figure3.regenerate_TIP(r1.smi,r.xenium,r1.xenium,R)
    return()
  }
  
  
  #1 Regenerate Figure 3A: UMAPs of CD8 T cells with TIP scores.
  TIP_Fig3A_umaps(r1 = r1.smi)
  #2 Regenerate Figure 3B: Dot-plot showing CD8 T cell TIP genes association with infiltration in different immune subsets.
  TIP_Fig3B_dotplot(R = R)
  #3 Regenerate Figure 3C: Barplot showing the overlap across the different TIPs 
  TIP_Fig3C_barplot(R)
  #4 Regenerate Figure 3D: Spatial maps showing the CD8 TIP scores in Xenium data 
  TIP_Fig3D_Xenium(r = r.xenium,r1 = r1.xenium)
  #5 Test for statistical significance in the Xenium data 
  P<-TIP_CD8.Xenium.test(r1 = r1.xenium)
  
  #6 Optional: Regenerate an extended version of Figure 3D
  TIP_Fig3D_Xenium.extended.version(r = r.xenium,r1 = r1.xenium)
  
  #7 Write TIPs to Table S5.
  write.xlsx(x = as.data.frame(list.2.mat(R$sig)),file = get.file("Tables/TableS5.xlsx"),asTable = T)
  
  return()
  
}

TIP_Fig3A_umaps<-function(r1){
  if(missing(r1)){r1<-HGSC_SMI.CD8.process()}
  l1<-umap.ggplot(r1$umap,labels = r1$plot[,c("Malignant","Fibroblast","TIP")],main = "CD8 T cells (all genes)")
  l2<-umap.ggplot(r1$cts$umap,labels = r1$plot[,c("Malignant","TIP","clusters2")],main = "CD8 T cells (T cell genes)")
  
  pdf(get.file("Figures/Fig3A.pdf"))
  call.multiplot(c(l1,l2,"NULL"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

TIP_Fig3B_dotplot<-function(R){
  if(missing(R)){R<-TIP_find_all()}
  # Load immune checkpoint gene set
  ICs<-readRDS(get.file("/Data/ICs.rds"))
  # Load GCPR gene set
  GPCRs<-readRDS(get.file("/Data/GPCRs.rds"))
  
  X<-lapply(names(R$cell.type.specific), function(x){
    X1<-R$cell.type.specific[[x]]$HLM[,c("HLM1.Z","HLM1.Estimate","HLM2.Z","HLM2.Estimate","SpecificB")]
    X1$Gene<-rownames(X1)
    X1$cell.type<-x
    return(X1)
  })
  names(X)<-names(R$cell.type.specific)
  X1<-X$CD8.T.cell
  
  sig<-R$sig[c("CD8.T.cell.TUMOR.up","CD8.T.cell.TUMOR.down")]
  g<-unlist(sig)
  sig<-get.top.cor(X1[g,c("HLM1.Z","HLM2.Z")],q = 25,min.ci = 3)
  g<-unique(unlist(sig))
  g<-g[order(X1[g,"HLM1.Z"])]
  ICs<-ICs[order(X1[ICs,"HLM1.Z"])]
  GPCRs<-intersect(g,GPCRs)
  GPCRs<-GPCRs[order(X1[GPCRs,"HLM1.Z"])]
  g<-c(ICs,GPCRs,setdiff(g,c(ICs,GPCRs)))
  
  X0<-rbind(X$CD8.T.cell,X$CD4.T.cell,X$NK.cell,
            X$B.cell,X$Treg,X$Monocyte,X$Mast.cell)
  X0<-X0[is.element(X0$Gene,g),]
  X0$Gene <- factor(X0$Gene, levels = g)
  X0$cell.type <- factor(X0$cell.type, 
                         levels = c("CD8.T.cell","CD4.T.cell","NK.cell",
                                    "Treg","Monocyte"))
  X0$Z<-abs(X0$HLM1.Z)
  X0$Z[X0$Z>40]<-40
  X0$Estimate<-X0$HLM1.Estimate
  p<-call.dotPlot(X0,cex = 8)
  
  pdf(get.file("Figures/Fig3B.pdf"))
  print(p)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

TIP_Fig3C_barplot<-function(R){
  sig<-lapply(R$cell.type.specific,function(X) return(intersect.lists.by.idx(X$sig[c(1,3)],X$sig[c(2,4)])))
  sig<-unlist(sig,recursive = F)
  sig<-sig[laply(sig,length)>2]
  names(sig)<-gsub("HLM1.Z","TIP",names(sig))
  
  pdf(get.file("Figures/Fig3C.pdf"))
  print(upset(fromList(sig[grepl("up",names(sig))]), order.by = "freq"))
  print(upset(fromList(sig[grepl("down",names(sig))]), order.by = "freq"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

TIP_Fig3D_Xenium<-function(r,r1){
  fun_color_range <- colorRampPalette(c("blue","grey","red"))
  f<-function(x){
    n1<-min(100,floor(length(x)/3))
    v<-call.discretize(x,n1)
    return(fun_color_range(n1)[v])}
  
  X<-apply(cap.mat(r1$scores[,1:2],cap = 0.01,MARGIN = 2),2,f)
  rownames(X)<-rownames(r1$scores)
  colnames(X)<-colnames(r1$scores[,1:2])
  
  # Expand to all cells, including non-malignant ones
  idx<-match(r$cells,r1$cells)
  X1<-X[idx,] 
  X1[is.na(idx)&r$cell.types=="Malignant",]<-"grey" # T cells are marked in black
  # Non-malignant cells other than CD8 T cells and NK cells are not shown
  X1[is.na(idx)&r$cell.types!="Malignant",]<-"white" 
  rownames(X1)<-r$cells
  r$cd8.b<-is.element(r$cells,r1$cells)
  r$b.plot<-r$cd8.b|r$cell.types=="Malignant"
  
  samples<-paste0("XEN_T10_",c("H4","C4","F2","F3","E5","G3","H1","I2"))
  
  pdf(get.file("Figures/Fig3D.pdf"))
  par(mfrow=c(2,2),oma = c(0, 0, 0, 1),xpd = T)
  for(x in samples){
    b1<-r$samples==x&r$b.plot;
    if(x=="XEN_T10_F3"){b1<-b1&r$coor[,1]>4000}
    if(x=="XEN_T10_G3"){b1<-b1&r$coor[,"y"]>3000}
    call.plot.plus(r$coor[b1,],labels = r$scores[b1,3],my.col = X1[b1,1],
                 b.top = r$cd8.b[b1],cex = ifelse(r$cd8.b[b1],0.4,0.2),
                 legend.flag = F,main = x,cex.axis = 1e-10)
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return(X1)
}

TIP_Fig3D_Xenium.extended.version<-function(r,r1){
  l<-NULL
  idx<-paste0("XEN_T10_",c("H4","C4","F2","F3","E5","G3","H1","I2"))
  for(x in idx){
    b<-r$samples==x
    b1<-r1$samples==x
    if(sum(b1)<40){next()}
    if(x=="XEN_T10_F3"){b<-b&r$coor[,"x"]>4000;b1<-b1&r1$coor[,"x"]>4000}
    if(x=="XEN_T10_G3"){b<-b&r$coor[,"y"]>3000;b1<-b1&r1$coor[,"y"]>3000}
    X<-cap.mat(r1$scores[b1,],MARGIN = 2,cap = 0.05)
    l[[x]]<-umap.ggplot(r1$coor[b1,],labels = X[,"R"],main = x,size = 1,remove.legend = F)
    l[[paste0(x,"_boxplot")]]<-call.boxplot(r1$scores[b1,"R"],main = x,ylab = "TIP score",xlab = "Malignant Env",
                                          add.n.of.samples(discretize.3.labels(r1$tme[b1,"Malignant"],q = 0.1)))
    l[[paste0(x,"_all")]]<-umap.ggplot(r$coor[b,],labels = r$cell.types[b],size = 0.1,main = x,remove.legend = F)
  }
  pdf(get.file("Figures/Fig3D_split.version.pdf"))
  call.multiplot(l[!grepl("boxplot",names(l))],nplots = 2,cols = 1)
  call.multiplot(l[grepl("boxplot",names(l))],nplots = 6,cols = 3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

HGSC_Figure3.data.regenerate<-function(){
  #1. Download the SMI and Xenium datasets 
  r.smi<-readRDS(get.file("Data/SMI_data.rds"))
  r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
  #2 Compute the tumor infiltration programs (TIPs) for different immune cell types.
  rslts<-TIP_find_all(r.smi)
  #3 Process CD8 T cells and add TIP scores 
  r1.smi<-HGSC_SMI.process.CD8(r = r.smi,rslts = R$cell.type.specific$CD8.T.cell)
  #4 Process CD8 T cells and add TIP scores 
  r1.xenium<-HGSC_Xenium.process.CD8.NK(r = r.xenium,rslts = R$cell.type.specific$CD8.T.cell)
   
  return()
}

TIP_find_all<-function(r){
  file1<-get.file("Results/TIP_all.rds")
  if(file.exists(file1)){
    return(readRDS(file1))
  }
  motile.cell<-c("CD8.T.cell","CD4.T.cell","NK.cell","Treg","Monocyte")
  R<-list()
  for(x in setdiff(motile.cell,names(R))){
    R[[x]]<-TIP_find(r,motile.cell = x,overwrite = F)
  }
  X<-lapply(names(R), function(x){
    X<-R[[x]]$HLM[,c("HLM1.Z","HLM1.Estimate","HLM2.Z","HLM2.Estimate","SpecificB")]
    X$Gene<-rownames(X)
    X$cell.type<-x
    return(X)
  })
  names(X)<-names(R)
  summary(X)
  
  X1<-union.multiple.mats(X)
  Xb<-X1[,grepl("SpecificB",colnames(X1),fixed = T)]
  X1<-X1[,grepl(".Z",colnames(X1),fixed = T)]
  X1$n.up<-rowSums(X1>3,na.rm = T)
  X1$n.down<-rowSums(X1<(-3),na.rm = T)
  Xb$n<-rowSums(Xb,na.rm = T)
  
  sig<-lapply(R,function(X) return(intersect.lists.by.idx(X$sig[c(1,3)],X$sig[c(2,4)])))
  sig<-unlist(sig,recursive = F)
  summary(sig)
  sig<-sig[laply(sig,length)>2]
  names(sig)<-gsub("HLM1.Z","TUMOR",names(sig))
  rslts<-list(cell.type.specific = R,sig = sig,sumZ = X1,sum.Specific = Xb)
  saveRDS(rslts,file = file1)
  return(rslts)
}

TIP_find<-function(r,motile.cell = "CD8.T.cell",overwrite = F){
  file1<-paste0(get.file("Results/TIP_"),motile.cell,".rds")
  if(!overwrite && file.exists(file1)){return(readRDS(file1))}
  
  r1<-set.list(r,r$cell.subtypes==motile.cell,name = motile.cell)
  r1$gene.dr<-rowSums(r1$tpm>0)
  r1<-set.list(r1,r1$gene.dr>10)
  r1$tme<-r$frames.tme[r1$frames,]
  r1$mal<-r1$tme[,"Malignant"]
  HLM1<-apply.formula.HLM(r1,X = r1$mal,Y = r1$tpm,formula = "y ~ (1 | frames) + x")
  HLM2<-apply.formula.HLM(r1,X = r1$mal,Y = r1$tpm,formula = "y ~ (1 | samples) + x")
  
  f<-function(x){
    b<-is.element(r$cell.subtypes,c(motile.cell,x))
    if(x=="TNK.cell"){
      b<-is.element(r$cell.subtypes,motile.cell)|is.element(r$cell.types,x)
    }
    Z<-t.test.mat(r$tpm[,b],r$cell.subtypes[b]==motile.cell)[,3]
    return(Z)
  }
  cell.types<-setdiff(c("Fibroblast","Malignant","B.cell","TNK.cell","Monocyte"),unique(r1$cell.types))
  X1<-t(laply(cell.types,f))
  colnames(X1)<-cell.types
  assertthat::are_equal(rownames(HLM1),rownames(HLM2))
  X2<-cbind.data.frame(HLM1 = HLM1,HLM2 = HLM2,X1[rownames(HLM1),],
                       Specific = rowMin(X1[rownames(HLM1),]),
                       SpecificB = rowMin(X1[rownames(HLM1),])>3)
  Xp<-X2[order(-X2$HLM1.Z),]
  Xp<-Xp[Xp$SpecificB&abs(Xp$HLM1.Z)>2,]
  z<-Xp$HLM1.Z;names(z)<-rownames(Xp)
  sig<-get.top.cor(Xp[,c("HLM1.Z","HLM2.Z")],q = 100,min.ci = 2)
  rslts<-list(motile.cell = motile.cell,HLM = X2,sig = sig,z = z)
  saveRDS(rslts,file = file1)
  return(rslts)
}

call.dotPlot<-function(X,cex = 12){
  p<-ggplot(data = X,aes(x = cell.type, y = Gene, color = Estimate, size = Z)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = cex),
          axis.text.y = element_text(size=cex)) +
    scale_colour_gradient2(low = "blue",high = "red",midpoint = 0,
                           oob = scales::squish, name = 'Effect size')+
    geom_point(shape = 1,colour = "black")
  return(p)
}

HGSC_SMI.process.CD8<-function(r,rslts,recompute = F){
  datafile<-get.file("/Data/SMI_data_CD8.T.cells.rds")
  
  if(file.exists(datafile)&!recompute){return(readRDS(datafile))}
  
  if(missing(r)){r<-readRDS(get.file("Data/SMI_data.rds"))}
  if(missing(rslts)){rslts<-readRDS(get.file("Results/TIP_CD8.T.cell.rds"))}
  
  motile.cell<-"CD8.T.cell"
  r1<-set.list(r,r$cell.subtypes==motile.cell,name = motile.cell)
  assertthat::are_equal(unique(r1$cell.subtypes),rslts$motile.cell)
  
  r1<-prep4OE(r1)
  r1$scores<-get.OE(r1,rslts$sig)
  r1$scoresC<-cap.mat(r1$scores,cap = 0.05,MARGIN = 2)
  r1$tme<-r$frames.tme[r1$frames,]
  r1$tmeC<-cap.mat(r1$tme,MARGIN = 2,cap = 0.05)
  
  r1<-seuratW_get.embedding(r1,no.genes = 950,n.pcs = 15,cd.flag = F,umap.flag = T,norm.flag = F)
  r2<-set.list(r1,is.element(r1$genes,rownames(rslts$HLM)[rslts$HLM$Specific>(-1)]))
  r2<-seuratW_get.embedding(r2,no.genes = nrow(r2$cd),n.pcs = 15,cd.flag = F,umap.flag = T,norm.flag = F,resolution = 1)
  r1$plot<-cbind.data.frame(r1$tmeC[,c("Malignant","Fibroblast")],
                            TIP = r1$scoresC[,"HLM1.Z"],clusters1 = r1$clusters,clusters2 = r2$clusters)
  r1$cts<-r2[c("umap","pca","pca.load","clusters")]
  saveRDS(r1,file = get.file("/Data/SMI_data_CD8.T.cells.rds"))
  return(r1)
}

HGSC_Xenium.process.CD8.NK<-function(r,rslts){
  file1<-get.file("Data/Xenium_data_CD8.T.NK.cells.rds")
  if(file.exists(file1)){return(readRDS(file1))}
  
  if(missing(r)){
    r<-readRDS(get.file("Data/Xenium_data.rds"))
  }
  
  if(missing(rslts)){
    rslts<-TIP_find(motile.cell = "CD8.T.cell")
    rslts$scRNA<-readRDS(get.file("Results/TIP_CD8.T.cell_scRNA.rds"))
  }
  
  b<-is.element(rownames(rslts$scRNA),r$genes)
  sig<-get.top.cor(rslts$scRNA[b&rslts$scRNA[,"P"]<1e-6,],q = 100,min.ci = 0.1,idx = "R")
  sig1<-rslts$sig[c("CD8.T.cell.TUMOR.up","CD8.T.cell.TUMOR.down")]
  sig1<-intersect.list1(sig1,r$genes)
  setdiff.lists.by.idx(sig,sig1)
  setdiff.lists.by.idx(sig1,sig)
  
  b<-is.element(r$cell.subtypes,c("CD8.T.cell","NK.cell"))
  r1<-set.list(r,b)
  r1<-prep4OE(r1,n.cat = 20)
  r1$scores<-get.OE(r1,sig)
  saveRDS(r1,file = file1)
  
  return(r1)
  
  
}

TIP_CD8.Xenium.test<-function(r1){
  r1$tumor<-r1$tme[,"Malignant"]>median(r1$tme[,"Malignant"])
  samples<-c("All",unique(r1$samples))
  P<-laply(samples,function(x){
    if(x=="All"){
      p1<-t.test.mat(t(r1$scores),r1$tumor)
      return(c(ttest = p1[,"zscores"],n = length(r1$cells)))
    }
    b<-r1$samples==x
    p1<-t.test.mat(t(r1$scores[b,]),r1$tme[b,"Malignant"]>median(r1$tme[b,"Malignant"]))
    return(c(ttest = p1[,"zscores"],n = sum(b)))
  })
  rownames(P)<-samples
  return(P)
}

