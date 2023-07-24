#### Results Section 4 ###
# Figure 4. mTIL
# Table S6. mTIL program and its regulators.

# Figure 4A. CD8 T cell umaps
# Figure 4B. CD8 T cell TIP genes' association with tumor infiltration in different cell subtypes.
# Figure 4C. TIPs overlaps (barplots)
# Figure 4D. Xenium spatial maps with CD8 TIP
# Figure 4E. Ligand-receptor interactions (generated via CytoScape)

HGSC_Figure4_mTIL<-function(r,r1,rslts){
  if(missing(r)){
    r<-readRDS(get.file("Data/SMI_data.rds"))
    r1<-readRDS(get.file("Data/SMI_data_malignant.rds"))
    r1<-mTIL_Fig4_prepData(r1)
    rslts<-readRDS(get.file("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds"))
  }
  
  #1 Regenerate Figure 4A: mTIL heatmap
  mTIL_Fig4A(r1 = r1,rslts = rslts)
  #2 Regenerate Figure 4C: mTIL spatial maps
  mTIL_Fig4C(r = r,r1 = r1,rslts = rslts)
  #3 Regenerate Figure 4D: mTIL as a function of TIL proximity and abundance (boxplot)
  mTIL_Fig4D(r1 = r1,rslts = rslts)
  #4 Regenerate Figure 4E: ROCs
  mTIL_Fig4E(r1 = r1,rslts = rslts)
  return()
}

mTIL_Fig4A<-function(r1,rslts){
  g<-intersect(c("CXCL10","CXCL2","CXCL5","CXCL9","LGALS1","NR3C1",
                 "BCL2","FGFR1","HDAC1","ITGB5","RELA","HDAC5"),unlist(rslts$sig))
  idx1<-order(r1$scoresAv[,"hot100"])
  idx2<-c(rslts$sig$hot100.up,rslts$sig$hot100.down)
  gene.cor<-spearman.cor(r1$tpmAv[,idx2],r1$scoresAv[,"hot"])
  sig1<-get.top.cor(gene.cor,q = 50,idx = "R")
  sig1<-union.lists(sig1,get.top.cor(gene.cor[g,],q = 50,idx = "R"),unique.flag = T)
  idx2<-c(sig1$R.up,sig1$R.down)
  X1<-cap.mat(t(r1$tpmAv[idx1,idx2]),cap = 0.05,MARGIN = 1)
  col.labels <-cbind.data.frame(TILs = r1$tmeAv[,"TNK.cell"]>median(r1$tmeAv[,"TNK.cell"]),
                                CD4 = r1$tmeAv[,"CD4.T.cell"]>0,
                                CD8 = r1$tmeAv[,"CD8.T.cell"]>0,
                                NK = r1$tmeAv[,"NK.cell"]>0,
                                samples = get.strsplit(rownames(r1$tmeAv),"_X",1))[idx1,]
  col.labels$sites<-r1$sites_binary[match(col.labels$samples,r1$samples)]
  row.labels <- ifelse(is.element(idx2,sig1$R.up),"Hot","Cold")
  row.labels <- cbind.data.frame(row.labels,row.labels)
  # pdf("/Volumes/Resource2/HGSC/Figures/HGSC.plot2_mTIL.heatmap.pdf")
  pdf(get.file("Figures/Fig4A.pdf"))
  call.heatmap(X1,scale = "row",col.labels = col.labels,row.labels = row.labels,cexRow = 0.1,legend.flag = F)
  barplot(r1$scoresAv[idx1,"hot100"],cex.names = 1e-10,xlab = "Frames",ylab = "Immune-rich signature: Overall Expression")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return(sig1)
}

mTIL_Fig4C<-function(r,r1,rslts){
  #1 Use the mTIL scores to color the malignant cells.
  X<-get.mat(r1$cells,c("hot","hot100"),data = NA)
  for(x in unique(r1$samples)){
    b<-r1$samples==x
    X[b,]<-apply(cap.mat(r1$scores[b,c("hot","hot100")],cap = 0.01,MARGIN = 2),2,function(x) scores.2.colors(x))
  }
  
  #2 Color TNK cells in black and other cells in grey.
  idx<-match(r$cells,r1$cells)
  X1<-X[idx,] 
  X1[is.na(idx)&r$cell.types=="TNK.cell",]<-"black" # T cells are marked in black
  X1[is.na(idx)&r$cell.types!="TNK.cell",]<-"grey" # Non-malignant cells other than T cells are shown in grey
  rownames(X1)<-r$cells
  
  #3 Plot the spatial maps per sample to regenerate Figure 4D.
  r$ids<-paste(r$patients,r$sites,sep = ", ")
  rslts$ttest.per.sample$ids<-r$ids[match(rownames(rslts$ttest.per.sample),r$samples)]
  samples<-c('SMI_T10_F001','SMI_T10_F015','SMI_T12_F001',
             'SMI_T12_F009','SMI_T12_F016','SMI_T13_F001')
  pdf(get.file("Figures/Fig4C.pdf"))
  par(mfrow=c(2,2),oma = c(0, 1, 0, 1),xpd = T)
  for(x in samples){
    b1<-r$samples==x
    call.plot.plus(-r$coor[b1,c("y","x")],labels = r$scores[b1,3],my.col = X1[b1,1],
                   b.top = r$cell.types[b1]=="TNK.cell",cex = 0.5,legend.flag = F,
                   main = paste(x,"\nZ =",round(rslts$ttest.per.sample[x,"zero.hot"],1),
                                "\nZ.median =",round(rslts$ttest.per.sample[x,"median.hot"],1),
                                rslts$ttest.per.sample[x,"ids"]),
                   pch = ifelse(r$cell.types[b1]=="TNK.cell",15,16))
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

mTIL_Fig4D<-function(r1,rslts,q1 = 0.75){
  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_data_malignant.rds"))
    r1<-mTIL_Fig4_prepData(r1)
    rslts<-readRDS("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds")
  }
  
  P<-list(samples = r1$scoresSamples[,"hot100"],
          patches = r1$scoresAv[,"hot100"],
          cells = r1$scores[,"hot100"])
  L3<-list(samples = discretize.3.labels(r1$tmeSamples[,"TNK.cell"],1-q1),
           patches = discretize.3.labels(r1$tmeAv[,"TNK.cell"],1-q1),
           cells = discretize.3.labels(r1$tme[,"TNK.cell"],1-q1))
  titles<-c(samples = "Sample-level",patches = "Frame-level",cells = "Cell-level")
  
  f<-function(x){
    l<-list(c("Low","Moderate"),c("Moderate","High"))
    p1<-call.boxplot(P[[x]],L3[[x]],ylab = "mTIL program OE",
                     xlab = "TIL levels",main = titles[x])
    p1<-p1+stat_compare_means(comparisons = l,label = "p.signif",method = "t.test")+theme(legend.position="none")
    return(p1)
  }
  
  p<-lapply(names(L3),f)
  
  pdf(get.file("Figures/Fig4D.pdf"))
  p1<-call.boxplot(r1$scores[,"hot100"],paste0("L",r1$envTIL1),main = "Fig3D",
                   ylab = "mTIL program OE",xlab = "TIL levels (radius)",blank.flag = T,labels = NULL)
  p[[4]]<-p1+stat_compare_means(comparisons = list(c("L0","L1"),c("L1","L2"),c("L2","L3")),label = "p.signif",method = "t.test")+theme(legend.position="none")
  print(call.multiplot(p[c(4,1,3,2)],cols = 3,nplots = 6))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

mTIL_Fig4E<-function(r1,rslts,q1 = 0.75){
  if(missing(r1)){
    r1<-mTIL_Fig4_prepData()
    rslts<-readRDS("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds")
  }
  
  P<-list(samples = r1$scoresSamples[,"hot100"],
          patches = r1$scoresAv[,"hot100"],
          cells = r1$scores[,"hot100"])
  Y<-list(samples = r1$tmeSamples[,"TNK.cell"]>quantile(r1$tmeSamples[,"TNK.cell"],q1),
          patches = r1$tmeAv[,"TNK.cell"]>quantile(r1$tmeAv[,"TNK.cell"],q1),
          cells = r1$tme[,"TNK.cell"]>quantile(r1$tme[,"TNK.cell"],q1))
  
  pdf(get.file("Figures/Fig4E.pdf"))
  plot.multi.ROCs(P,Y)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

HGSC_SMI.process.mal<-function(r1){
  
  if(missing(r1)){r1<-readRDS(get.file("Data/SMI_data_malignant.rds"))}
  
  envTIL0<-paste(r1$tmeCell[[1]][,"TNK.cell"]>0,
                 r1$tmeCell[[2]][,"TNK.cell"]>0,
                 r1$tmeCell[[3]][,"TNK.cell"]>0)
  
  r1$envTIL1<-rowSums(cbind(r1$tmeCell[[1]][,"TNK.cell"]>0,
                            r1$tmeCell[[2]][,"TNK.cell"]>0,
                            r1$tmeCell[[3]][,"TNK.cell"]>0))
  
  r1$envTIL2<-rowSums(cbind(r1$tmeCell[[2]][,"TNK.cell"]>0,
                            r1$tmeCell[[2]][,"CD4.T.cell"]>0,
                            r1$tmeCell[[2]][,"CD8.T.cell"]>0,
                            r1$tmeCell[[2]][,"NK.cell"]>0))
  
  r1$b.cold<-r1$tmeCell[[2]][,"TNK.cell"]==0&r1$tmeCell[[2]][,"TNK.cell_LC"]==0
  r1$cd8.pos<-r1$tmeCell[[2]][,"CD8.T.cell"]>0
  r1$cd4.pos<-r1$tmeCell[[2]][,"CD4.T.cell"]>0
  r1$nk.pos<-r1$tmeCell[[2]][,"NK.cell"]>0
  return(r1)
  
}

plot.multi.ROCs<-function(P,Y,b,main = ""){
  p1<-plot.auc(P[[1]],Y[[1]],plot.flag = F,subplotF = F)
  p2<-plot.auc(P[[2]],Y[[2]],plot.flag = F,subplotF = F)
  p3<-plot.auc(P[[3]],Y[[3]],plot.flag = F,subplotF = F)
  a<-c(get.auc(P[[1]],Y[[1]]),
       get.auc(P[[2]],Y[[2]]),
       get.auc(P[[3]],Y[[3]]))
  a<-round(a,2)
  plot(p1,ylim = c(0,1),main = main)
  plot(p2,ylim = c(0,1),col = "red",add = T)
  plot(p3,ylim = c(0,1),col = "blue",add = T)
  abline(a = 0,b = 1,col = "gray30")
  legend(x = 0.4,y = 0.3,col = c("black","red","blue"),
         legend = paste0(names(P)," (AUC = ",a,")"),lty = 1,lwd = 3,
         cex = 0.8)
}

scores.2.colors<-function(x.class){
  palette("default")
  call_col<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = "hsv")
  return(call_col)
}

#### ********** Extra ********** ####


mTIL_Fig_supp<-function(r1){
  r1<-mTIL_Fig4_prepData()
  
  pdf("/Volumes/ljerby/ljerby/HGSC/Figures/HGSC_Fig3_supp.mTIL.per.site.pdf")
  print("Supplemantary Figure showing the mTIL-TIL connections in different sites.")
  print(call.boxplot(r1$scores[,"hot100"],r1$sites,labels = paste0(r1$envTIL1),
                     ylab = "mTIL program OE",xlab = "Site",legend.name = "TIL level"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  pdf("/Volumes/Resource2/HGSC/Figures/HGSC.plot4_mTIL.boxplot2.pdf")
  
  print("Different levels of TILs, based on distance and cell types")
  p1<-call.boxplot(r1$scores[,"hot100"],paste0("L",r1$envTIL1),ylab = "mTIL program OE",xlab = "TIL levels (radius)",blank.flag = T,labels = NULL)
  p1<-p1+stat_compare_means(comparisons = list(c("L0","L1"),c("L1","L2"),c("L2","L3")),label = "p.signif",method = "t.test")+theme(legend.position="none")
  p2<-call.boxplot(r1$scores[,"hot100"],paste0("L",r1$envTIL2),ylab = "mTIL program OE",xlab = "No. types of TILs",labels = NULL)
  p2<-p2+stat_compare_means(comparisons = list(c("L0","L1"),c("L1","L2"),c("L2","L3"),c("L3","L4")),label = "p.signif",method = "t.test")+theme(legend.position="none")
  call.multiplot(list(p1,NULL,p2,NULL,NULL),nplots = 6,cols = 3)
  
  
  pdf("/Volumes/Resource2/HGSC/Figures/HGSC.plot4_mTIL.boxplot2.pdf")
  print("mTIL as a function of CD8/CD4 T cells and NK cells.")
  par(mfrow=c(2,3),oma = c(1, 1, 0, 1),xpd = T)
  b<-r1$b.cold|r1$cd4.pos
  boxplot.test(r1$scores[b,"hot100"],ifelse(r1$cd4.pos[b],"CD4+","CD4-"),xlab = "Cell environment",
               main = "CD4 T cell",ylab = "mTIL expression",legend.flag = F)
  b<-r1$b.cold|r1$cd8.pos
  boxplot.test(r1$scores[b,"hot100"],ifelse(r1$cd8.pos[b],"CD8+","CD8-"),xlab = "Cell environment",
               main = "CD8 T cell",ylab = "mTIL expression",legend.flag = F)
  b<-r1$b.cold|r1$nk.pos
  boxplot.test(r1$scores[b,"hot100"],ifelse(r1$nk.pos[b],"NK+","NK-"),xlab = "Cell environment",
               main = "NK cell",ylab = "mTIL expression",legend.flag = F)
  
  # library(ggpubr)
  
  
  print("Supplementary Figure showing the baseline variation across patients")
  p3<-call.boxplot(r1$scores[r1$b.cold,"hot100"],r1$patients[r1$b.cold],ylab = "Program Expression",
                   xlab = "Patients",add.anova = T,main = "TIL-high Program, cold niches")
  print(p3)
  
  # call.boxplot(r1$scores[r1$b.cold,"hot100"],r1$samples[r1$b.cold],ylab = "Program Expression",
  #            xlab = "Patients",labels = r1$sites[r1$b.cold],add.anova = T,main = "TIL-high Program, cold niches")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

SMI_mTIL.detect<-function(){
  r<-readRDS(get.file("Data/SMI_data.rds"))
  r1<-cell2env_prep.data(r,cell.type = "Malignant")
  assertthat::are_equal(unique(r1$cell.types),"Malignant")
  cell2env1<-cell2env_main(r1,tme.cell.type = "TNK.cell",cell.type = "Malignant",name = "mTIL")
}

cell2env_main<-function(r1,tme.cell.type = "TNK.cell",cell.type = "Malignant",name = "mTIL"){
  
  file1<-paste0(get.file("Results/HGSC_"),name,"_",paste(cell.type,collapse = "."),
                "2env_",paste(tme.cell.type,collapse = "."),".rds")
  if(file.exists(file1)){return(rslts<-readRDS(file1))}
  
  print(paste(cell.type,"<-> TME",tme.cell.type))
  if(length(tme.cell.type)>1){
    r1$tme.v<-rowSums(r1$tme[,tme.cell.type])
  }else{
    r1$tme.v<-r1$tme[,tme.cell.type]
  }
  
  if(any(grepl("_LC",unique(r1$cell.types)))){
    r0<-set.list(r1,r1$b.subsample)
    print(ggbarplot.frq(r1$b.subsample,r1$samples,freq.flag = F,l2.name = "Sample",l1.name = "Include"))
  }else{
    r0<-set.list(r1,r1$b.subsample.hc)
    # print(ggbarplot.frq(r1$b.subsample.hc,r1$samples,freq.flag = F,l2.name = "Sample",l1.name = "Include"))
  }
  
  rslts<-list(cell.type = unique(r1$cell.types),tme.cell.type = tme.cell.type,
              subsample2.flag = subsample2.flag, cells = r0$cells)
  
  Z1<-apply.formula.all.HLM(r = r0,Y = r0$tpm,X = r0$tme.v>0,
                            MARGIN = 1,formula = "y ~ (1 | frames) + x + comp",ttest.flag = F)
  
  rslts$hlm1<-Z1
  rslts$sum<-Z1$HLM[-1]
  
  # Add signatures of the top markers identified 
  rslts$sig<-get.top.cor(rslts$sum[,1:2],min.ci = -4,q = 100)[c(1,3)]
  names(rslts$sig)<-paste(name,c("up","down"),sep = ".")
  saveRDS(rslts,file = file1)
  # SMI_heatmap_hot.cold(tme.cell.type)
  return(rslts)
}

cell2env_prep.data<-function(r,cell.type){
  r1<-set.list(r,is.element(r$cell.types,cell.type))
  r1<-log.comp(r1)
  r1$tme<-r$frames.tme[r1$frames,]
  r1$tpmAv<-average.mat.rows(t(r1$tpm),r1$frames)
  r1$tmeAv<-r$frames.tme[rownames(r1$tpmAv),]
  return(r1)
}

cell2env_SVM<-function(r1,cv.type,tme.cell.type){
  if(cv.type=="all"){
    rslts<-list(svm.pts.cv = cell2env_SVM(r1,cv.type = "patient.level",tme.cell.type),
                svm.cells.cv = cell2env_SVM(r1,cv.type = "cell.level",tme.cell.type))
    return(rslts)
  }
  
  if(length(tme.cell.type)>1){
    r1$tme.v<-rowSums(r1$tme[,tme.cell.type])
    r1$tmeAv.v<-rowSums(r1$tmeAv[,tme.cell.type])
  }else{
    r1$tme.v<-r1$tme[,tme.cell.type]
    r1$tmeAv.v<-r1$tmeAv[,tme.cell.type]
  }
  
  if(cv.type=="patient.level"){
    train.patients<-sample(unique(r1$patients),40)
    train.frames<-unique(r1$frames[is.element(r1$patients,train.patients)])
    b<-is.element(r1$patients,train.patients)
  }else{
    train.frames<-sample(unique(r1$frames),round(length(unique(r1$frames))/2))
    b<-is.element(r1$frames,train.frames)
  }
  
  bAv<-is.element(rownames(r1$tmeAv),train.frames)
  rslts<-list(cell.type = unique(r1$cell.types),tme.cell.type = tme.cell.type,
              cv.type = cv.type,
              train = r1$cells[b],test = r1$cells[!b])
  rslts$tme.cor.train<-spearman.cor(r1$tpmAv[bAv,],r1$tmeAv.v[bAv])
  sig1<-get.top.cor(rslts$tme.cor.train,idx = "R",q = 50)
  sig2<-get.top.cor(rslts$tme.cor.train,idx = "R",q = 100)
  
  # library(e1071)
  y<-r1$tmeAv.v>0
  rslts$svm1<-SMI_call.svm(y = y,X = r1$tpmAv[,unlist(sig1)],bAv)
  rslts$svm2<-SMI_call.svm(y = y,X = r1$tpmAv[,unlist(sig2)],bAv)
  rslts$svm3<-SMI_call.svm(y = y,X = r1$tpmAv,bAv)
  print(c(rslts$svm1$auc.test,rslts$svm2$auc.test,rslts$svm3$auc.test))
  return(rslts)
}

SMI_call.svm<-function(y,X,b,cost = 1,main = "",X3){
  model <- svm(X[b,], as.factor(y[b]),probability = T,cost = cost)
  pred1 <- predict(model, X[b,],probability = T)
  pred1p <- attr(pred1, "probabilities")[,"TRUE"]
  predAll <- predict(model, X,probability = T)
  predAllp <- attr(predAll, "probabilities")[,"TRUE"]
  
  # d<-ifelse(mean(pred1p[pred1])>mean(pred1p[pred1==FALSE]),1,-1)
  # pred1p<-pred1p*d
  table(pred1,y[b])
  # boxplot(pred1p~pred1)
  pred2 <- predict(model, X[!b,],probability = T)
  pred2p <- attr(pred2, "probabilities")[,"TRUE"]
  # boxplot(pred2p~pred2)
  
  rslts1<-list(train = mean(y[b]==pred1),
               test = mean(y[!b]==pred2),
               pred.train = cbind(prob = pred1p,y = y[b]),
               pred.test = cbind(prob = pred2p,y = y[!b]),
               predAll = cbind(prob = predAllp, y = y),
               auc.train = plot.auc(pred1p,y[b],main = "Train"),
               auc.test = plot.auc(pred2p,y[!b],main = "Test"))
  par(mfrow=c(1,2),oma = c(1, 1, 0, 1),xpd = T)
  boxplot.test(pred1p,y[b],main = main)
  if(sum(y[!b]==1)>2|sum(y[!b]==0)>2){
    boxplot.test(pred2p,y[!b],ref.label = TRUE,alternative = "greater")
  }
  
  if(!missing(X3)){
    pred3<-predict(model, X3,probability = T)
    pred3p <- attr(pred3, "probabilities")[,"TRUE"]
    rslts1$pred3<-pred3p
  }
  
  return(rslts1)
}



