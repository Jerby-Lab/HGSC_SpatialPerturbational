
#' Discovery Dataset Malignant Cells
#' Processes Discovery dataset for malignant cells.
#' @param r1 Data frame with SMI data (optional).
#' @return Data frame with processed SMI data.
HGSC_SMI.process.mal<-function(r1){
  
  if(missing(r1)){r1<-readRDS(get.file("/Data/ST_Discovery_malignant.rds"))}
  
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

#' Extract M_{TIL} and T/NK cell levels
#' Detects M_{TIL} and T/NK cell levels in Discovery dataset using overall
#' expression calcualtions. 
#' @return List of results from mTIL detection.
SMI_mTIL.detect<-function(){
  r<-readRDS(get.file("Data/SMI_data.rds"))
  r1<-cell2env_prep.data(r,cell.type = "Malignant")
  assertthat::are_equal(unique(r1$cell.types),"Malignant")
  cell2env1<-cell2env_main(r1,tme.cell.type = "TNK.cell",cell.type = "Malignant",name = "mTIL")
}

#' Cell <> Environment Main Analysis
#' Main function for cell to environment analysis (uses mixed effects modeling)
#' to discover gene sets associated with abundance of a certain cell type in a 
#' cells environment. 
#' @param r1 List object with Discovery data (malignant cells only)
#' @param tme.cell.type Type of TME cells (default is "TNK.cell").
#' @param cell.type Type of cells to analyze (default is "Malignant").
#' @param name Name for the analysis (default is "mTIL").
#' @return List of results from the cell to environment analysis.
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

#' Prep Data for Cell <> Environment analysis 
#' Prepares data for cell to environment analysis.
#' @param r List object with Discovery data (malignant cells only)
#' @param cell.type Type of cells to analyze.
#' @return Data frame with prepared data.
cell2env_prep.data<-function(r,cell.type){
  r1<-set.list(r,is.element(r$cell.types,cell.type))
  r1<-log.comp(r1)
  r1$tme<-r$frames.tme[r1$frames,]
  r1$tpmAv<-average.mat.rows(t(r1$tpm),r1$frames)
  r1$tmeAv<-r$frames.tme[rownames(r1$tpmAv),]
  return(r1)
}

#' Cell to Environment SVM
#' Performs SVM analysis for cell to environment data.
#' @param r1 List object with Discovery data (malignant cells only)
#' @param cv.type Type of cross-validation ("all", "patient.level", "cell.level").
#' @param tme.cell.type Type of TME cells.
#' @return List of SVM results.
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

#' Wrapper around SVM analysis for Discovery dataset
#' Calls SVM analysis for Discovery dataset
#' @param y Response variable.
#' @param X Predictor matrix.
#' @param b Boolean vector for training data.
#' @param cost Cost parameter for SVM (default is 1).
#' @param main Title for the plot (default is "").
#' @param X3 Additional predictor matrix for validation (optional).
#' @return List of SVM results.
SMI_call.svm<-function(y,X,b,cost = 1,main = "",X3){
  model <- svm(X[b,], as.factor(y[b]),probability = T,cost = cost)
  pred1 <- predict(model, X[b,],probability = T)
  pred1p <- attr(pred1, "probabilities")[,"TRUE"]
  predAll <- predict(model, X,probability = T)
  predAllp <- attr(predAll, "probabilities")[,"TRUE"]
  
  table(pred1,y[b])
  pred2 <- predict(model, X[!b,],probability = T)
  pred2p <- attr(pred2, "probabilities")[,"TRUE"]
  
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