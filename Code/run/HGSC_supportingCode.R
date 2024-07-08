#### Supporting Code for Statistics, Plotting, and Visualization

#' Plot UMAP (via ggplot)
#'
#' Wrapper around ggplot to visualize the UMAP. 
#' 
#' @param umapX matrix of UMAP coordinates 
#' @param labels vector of labels of each point 
#' @param labels.name (default = "") string of name of point levels 
#' @param main (default = "") string of the name of the plot
#' @param size (default = 0.2) size of text
#' @param xlim1 bounds on x-axis to plot 
#' @param ylim1 bounds on y-axis to plot 
#' @param reorder.flag (default = F) boolean to shuffle the order the points are plotted 
#' @param remove.legend (default = F) boolean to remove legend.
#' @return ggplot object of finalized plot. 
#' @export
umap.ggplot<-function(umapX,labels,labels.name = "",main = "",size = 0.2,
                      xlim1,ylim1,reorder.flag = F,remove.legend = F){
  if((is.matrix(labels)|is.data.frame(labels))&&ncol(labels)>1){
    p<-lapply(colnames(labels),function(x){
      p1<-umap.ggplot(umapX,labels = labels[,x],labels.name = labels.name,main = x,size = size)
      return(p1)
    })
    names(p)<-colnames(labels)
    return(p)
  }
  if(reorder.flag){
    print("reordering")
    idx<-order(labels)
    labels<-labels[idx]
    umapX<-umapX[idx,]
  }
  
  xylabs <- colnames(umapX)
  colnames(umapX)<-c("UMAP1","UMAP2")
  X <- cbind.data.frame(umapX,col = labels)
  p <- ggplot(X, aes(x = UMAP1, y = UMAP2, color = col)) + geom_point(size = size) + labs(color = labels.name)
  p <- p + ggtitle(main)
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab(xylabs[1]) + ylab(xylabs[2])
  if(!missing(xlim1)){
    xlim(xlim1)
    ylim(ylim1)
  }
  
  if(!is.numeric(labels)){
    if(remove.legend){p<-p+theme(legend.position = "none")}
    return(p)
  }
  p<-p+scale_color_gradient2(midpoint=mean(labels),
                             low="blue", mid="gray",high="red", space ="Lab")
  
  if(remove.legend){
    p<-p+theme(legend.position = "none")
  }
  return(p)
}

#' Call Multiplot
#'
#' Calls the multiplot function to print multiple plots in a grid layout.
#'
#' @param plotlist List of ggplot objects to be plotted.
#' @param nplots Number of plots to display per call (default is 4).
#' @param cols Number of columns in the layout (default is 2).
#' @return None. This function is used for its side effect of printing plots.
#' @export
call.multiplot<-function(plotlist,nplots = 4,cols = 2){
  flag<-F
  while(!is.null(plotlist)&!flag){
    print(multiplot(plotlist = plotlist[1:min(nplots,length(plotlist))],cols = cols))
    flag<-(min(nplots,length(plotlist))+1)>length(plotlist)
    plotlist<-plotlist[(min(nplots,length(plotlist))+1):length(plotlist)]
  }
}

#' Multiplot
#'
#' Arranges multiple ggplot objects in a grid layout.
#'
#' @param ... ggplot objects.
#' @param plotlist List of ggplot objects (optional).
#' @param file File path to save the plot (optional).
#' @param cols Number of columns in the layout (default is 1).
#' @param layout Custom layout matrix (optional).
#' @return None. This function is used for its side effect of printing plots.
#' @export
multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL){
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Kaplan-Meier Plot
#'
#' Creates a Kaplan-Meier plot with options for customization.
#'
#' @param r Data frame containing survival information.
#' @param v Vector of variable values for stratification.
#' @param main Title of the plot (default is '').
#' @param X Additional covariates (optional).
#' @param qua Quantile threshold for defining high/low groups (default is 0.2).
#' @param xlim X-axis limits (optional).
#' @param direction Direction for one-sided p-values (default is 0).
#' @param legend.flag Boolean to show legend (default is TRUE).
#' @param ylab Y-axis label (default is "Survival probability").
#' @param four.levels Boolean to use four levels for stratification (default is FALSE).
#' @return Vector of group labels.
#' @export
km.plot3 <- function(r,v,main = '',X = NULL,qua = 0.2,xlim = NULL,direction = 0,
                     legend.flag = T,ylab = "Survival probability",four.levels = F){
  M1<-summary(coxph(r$survival ~ v))$coefficients
  coxD<-M1[1,"coef"]
  cox.p<-M1[1,"Pr(>|z|)"]
  if(!is.null(X)){
    Mc<-summary(coxph(r$survival ~ cbind(v,X)))$coefficients
    cox.p.c<-Mc[1,"Pr(>|z|)"]
    coxD.c<-Mc[1,"coef"]
  }
  # q1<-mean(v,na.rm=T)+(sd(v,na.rm = T)*0.75);q2<-mean(v,na.rm=T)-(sd(v,na.rm = T)*0.75);b1<-v>=q1;b2<-v<=q2
  b1<-v>=quantile(v,1-qua,na.rm = T);b2<-v<=quantile(v,qua,na.rm = T)
  G<-ifelse(b1,"High","Moderate")
  col<-c("red","blue","darkgreen")
  if(four.levels){
    b1m<-!b1&!b2&v>=median(v,na.rm = T)
    b2m<-!b1&!b2&v<median(v,na.rm = T)
    G[b1m]<-"Moderate.high";G[b2m]<-"Moderate.low"
    col<-c("red","blue","darkgreen","black")
  }
  G[b2]<-"Low"
  km2<-npsurv(r$survival ~ G)
  sdf2<-survdiff(r$survival ~ G)
  sdf2<-(1 - pchisq(sdf2$chisq, length(sdf2$n) - 1))/3
  if(four.levels){
    l<-paste0(c("High","Low","Moderate.high","Moderate.low")," (",km2$n,")")
  }else{
    l<-paste0(c("High","Low","Moderate")," (",km2$n,")")
  }
  
  if(is.null(xlim)){
    survplot(km2,col = col,lty = c(1,1), xlab = 'Years',label.curves = F,
             ylab = ylab,n.risk = T)
  }else{
    survplot(km2,col = col,lty = c(1,1), xlab = 'Years',label.curves = F,xlim = c(0,xlim),
             ylab = ylab,n.risk = T)
  }
  if(legend.flag){
    legend("topright",fill = col[c(setdiff(1:length(col),2),2)],cex = 0.8,
           legend = l[c(setdiff(1:length(col),2),2)])
  }
  
  if(!is.null(X)){
    if(direction==0){
      P<-c(cox.p,cox.p.c,sdf2)
    }else{
      P<-get.onesided.p.value(direction*c(coxD,coxD.c,coxD),c(cox.p,cox.p.c,sdf2))
    }
    P<-format(P,scientific = T,digits = 2)
    main<-paste0(main,"\nP=",P[1],", Pc=",P[2],"\nlogrank=",P[3])
  }else{
    if(direction==0){
      P<-c(cox.p,sdf2)
    }else{
      P<-get.onesided.p.value(direction*c(coxD,coxD),c(cox.p,sdf2))
    }
    P<-format(P,scientific = T,digits = 2)
    main<-paste0(main,"\nP=",P[1],", logrank=",P[2])
  }
  title(main,cex.main =1)
  return(G)
}

#' Call Plot Plus
#'
#' Creates a scatter plot with additional features like highlighting top points.
#'
#' @param x Numeric vector or matrix of x-coordinates.
#' @param y Numeric vector of y-coordinates (optional).
#' @param labels Vector of labels for points.
#' @param b.top Logical vector indicating top points.
#' @param red.top Boolean to color top points in red (default is FALSE).
#' @param regression.flag Boolean to add regression line (default is FALSE).
#' @param my.col Custom colors for points (optional).
#' @param set.flag Boolean to set plot parameters (default is FALSE).
#' @param cor.flag Boolean to show correlation (default is FALSE).
#' @param pch Plotting character (default is 16).
#' @param cex Size of points (default is 0.3).
#' @param main Title of the plot (default is "").
#' @param ylab Y-axis label (default is "tSNE2").
#' @param xlab X-axis label (default is "tSNE1").
#' @param cex.axis Size of axis labels (default is 0.6).
#' @param add.N Boolean to add sample size to labels (default is FALSE).
#' @param grey.zeros Boolean to grey out zero values (default is FALSE).
#' @param legend.flag Boolean to show legend (default is TRUE).
#' @return Regression line values.
#' @export
call.plot.plus<-function(x, 
                         y = NULL,
                         labels,
                         b.top
                         ,red.top = F,
                         regression.flag = F,
                         my.col = NULL,
                         set.flag = F,
                         cor.flag = F,
                         pch=16,
                         cex=0.3,
                         main="",
                         ylab = "tSNE2",
                         xlab = "tSNE1", 
                         cex.axis = 0.6,
                         add.N = F,
                         grey.zeros = F,
                         legend.flag = T){
  
  regl<-call.plot(x = x,y = y,labels,regression.flag,my.col = my.col,
                  set.flag = set.flag,cor.flag = cor.flag,
                  pch = pch,cex = cex,main = main,ylab = ylab,xlab = xlab,
                  cex.axis = cex.axis,
                  add.N = add.N,legend.flag = legend.flag)
  if(is.null(y)){
    v<-colnames(x)
    if(xlab==""){xlab<-v[1]}
    if(ylab==""){ylab<-v[2]}
    y<-x[,2];x<-x[,1]
  }
  if(red.top){
    points(x[b.top],y[b.top],cex = cex,col = "red",pch = 1)
  }else{
    if(is.null(my.col)){
      if(grey.zeros){
        my.col<-rep("grey",length(labels))
        my.col[labels>0]<-labels.2.colors(labels[labels>0])
      }else{
        my.col <- labels.2.colors(labels)
      }
    }
    points(x[b.top],y[b.top],cex = cex,col = my.col[b.top],pch = 16)
  }
  return(regl)
  
}

#' Call Plot
#'
#' Creates a scatter plot with options for customization.
#'
#' @param x Numeric vector or matrix of x-coordinates.
#' @param y Numeric vector of y-coordinates (optional).
#' @param labels Vector of labels for points.
#' @param regression.flag Boolean to add regression line (default is FALSE).
#' @param my.col Custom colors for points (optional).
#' @param set.flag Boolean to set plot parameters (default is FALSE).
#' @param cor.flag Boolean to show correlation (default is FALSE).
#' @param legend.flag Boolean to show legend (default is TRUE).
#' @param pch Plotting character (default is 16).
#' @param cex Size of points (default is 0.5).
#' @param main Title of the plot (default is "").
#' @param ylab Y-axis label (default is "UMAP2").
#' @param xlab X-axis label (default is "UMAP1").
#' @param cex.axis Size of axis labels (default is 0.6).
#' @param add.N Boolean to add sample size to labels (default is FALSE).
#' @param cex.main Size of main title (default is 1).
#' @param color.spec Color specification for points (default is "rgb").
#' @return Optional regression line values.
#' @export
call.plot<-function(x, y = NULL,labels,regression.flag = F,my.col = NULL,
                    set.flag = F,cor.flag = F,legend.flag = T,
                    pch=16,cex=0.5,main="",ylab = "UMAP2",xlab = "UMAP1", 
                    cex.axis = 0.6,add.N = F,cex.main = 1,
                    color.spec = "rgb"){
  main<-capitalize(main)
  if(add.N&length(unique(labels))<30){
    labels<-add.n.of.samples(labels)
  }
  if(set.flag){
    par(mar=c(8, 7, 4.1, 12.1), xpd=TRUE)
  }
  if(is.null(my.col)){
    my.col<-labels.2.colors(labels,color.spec = color.spec)
  }
  if(is.null(y)){
    if(missing(xlab)){xlab<-colnames(x)[1]}
    if(missing(ylab)){ylab<-colnames(x)[2]}
    y<-x[,2];x<-x[,1]
  }
  
  if(cor.flag){
    xy.cor<-spearman.cor(y,x)
    main <- paste(main, "\nR =",format(xy.cor[1],digits = 2),"P =",format(xy.cor[2],scientific = T,digits = 2))
  }
  plot(x,y,col=my.col,pch=pch,cex=cex,main=main,ylab=ylab,xlab = xlab,cex.axis = cex.axis,cex.main = cex.main)
  
  labels<-gsub(" ","_",labels)
  l<-(max(x,na.rm = T)-min(x,na.rm = T))/20
  if(length(unique(labels))<30&legend.flag){
    if(length(pch)==length(labels)){
      map<-unique(paste(labels,my.col,pch))
      labels.n<-as.matrix(table(labels))
      idx<-match(get.strsplit(map,' ',1),names(labels.n))
      map[,1]<-paste0(map[,1]," (N = ",m[idx],")")
      print(as.integer(get.strsplit(map,' ',3)))
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),
             legend = get.strsplit(map,' ',1),
             col = get.strsplit(map,' ',2),
             inset=c(-0.5,0),
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }else{
      map<-unique(paste(labels,my.col,pch))
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),inset = c(-0.5,0),
             legend = gsub("_"," ",get.strsplit(map,' ',1)),
             col = get.strsplit(map,' ',2),xpd = T,
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }
    
  }
  if(regression.flag ==1){
    b<-!is.na(x)&!is.na(y)
    v<-lowess(x[b],y[b])
    lines(v)
    return(v)
  }
  if(regression.flag ==2){
    b<-!is.na(x)&!is.na(y)
    ulabels<-unique(labels)
    for(i in ulabels){
      bi<-b&labels==i
      v<-lowess(x[bi],y[bi])
      lines(v)
    }
    
  }
  
  
}

#' Call Boxplot
#'
#' Creates a boxplot with additional customization options.
#'
#' @param y Numeric vector of y-values.
#' @param x Factor vector of x-values.
#' @param unique.x Unique levels of x (optional).
#' @param f Function to apply to each group (default is median).
#' @param ylab Y-axis label (default is '').
#' @param xlab X-axis label (default is '').
#' @param main Title of the plot (default is '').
#' @param labels Vector of labels for points (optional).
#' @param legend.name Name for the legend (default is "").
#' @param add.anova Boolean to add ANOVA p-value to title (default is FALSE).
#' @param blank.flag Boolean to remove grid lines (default is TRUE).
#' @param cex Size of axis labels (default is 0.7).
#' @param order.flag Boolean to order x levels by function values (default is TRUE).
#' @param b.ref Reference group for comparison (optional).
#' @param p.val.show Boolean to show p-values (default is NULL).
#' @return ggplot object of the boxplot.
#' @export
call.boxplot<-function (y,x,unique.x, f = median,
                        ylab = '',xlab = '',main = '',labels=NULL,
                        legend.name = "",add.anova = F,blank.flag = T,
                        cex = 0.7,order.flag = T,b.ref = NULL,
                        p.val.show = is.null(b.ref)){
  b<-is.infinite(y)|is.na(y)
  if(length(labels)==length(x)){
    labels<-labels[!b]
  }
  y<-y[!b];x<-x[!b]
  if(p.val.show){
    a<-aov(y ~ as.factor(x))
    kt<-kruskal.test(y~as.factor(x))
    p.kt = format(kt$p.value,digits=2)
    p.anova=format(unlist(summary(a))['Pr(>F)1'],digits=2)
    if(add.anova){
      main = paste(main,'\n(ANOVA p-value = ',p.anova,'\nKruskal p-value =',p.kt,')',sep = '')
    }
  }
  if(!is.null(b.ref)){
    main <- paste0(main,"\n(",my.format.pval(t.test.labels(y,b.ref,alternative = "greater")),", ",
                   "AUC = ",round(get.auc(y,b.ref),2),")")
  }
  
  if(missing(unique.x)){
    unique.x<-unique(x)
    x.med<-apply(as.matrix(unique.x),1,function(xi) f(y[is.element(x,xi)]))
    names(x.med)<-unique.x
    if(order.flag){
      unique.x<-unique.x[order(x.med)]
    }
    #labels<-labels[unique.x]
  }
  
  #idx<-order(y.med[y])
  x <- data.frame(name = x, val = y)
  x$name <- factor(x$name, levels = unique.x)
  p <- ggplot(x, aes(x = name, y = val,fill = labels)) + #geom_point(stat="identity") +
    labs(y = ylab, title = main, x = xlab) +
    geom_boxplot() + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = rel(cex)),
          axis.text.y = element_text(size = rel(cex)),
          axis.title.y = element_text(size = rel(cex)),
          plot.title = element_text(size = rel(1)))
  p <- p + labs(color = legend.name) + scale_fill_discrete(name = legend.name)
  if(blank.flag){
    p <- p+ theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())
  }
  return(p)
  #labs(title = paste(screen.name,'ANOVA p-value = ',p.anova)) +
}

#' Call Heatmap
#'
#' Creates a heatmap with various customization options.
#'
#' @param m Matrix of values to be plotted.
#' @param main Title of the plot (default is '').
#' @param col.labels Column labels (optional).
#' @param row.labels Row labels (optional).
#' @param k Number of clusters for hierarchical clustering (default is 3).
#' @param filter.na Boolean to filter out NA values (default is TRUE).
#' @param cexCol Size of column labels (default depends on number of columns).
#' @param cexRow Size of row labels (default depends on number of rows).
#' @param m.value Value to display in the key (default is '').
#' @param scale Scaling method (default is "none").
#' @param cluster.flag Clustering method (default is "none").
#' @param sym.flag Boolean for symmetric clustering (default is FALSE).
#' @param xlab X-axis label (default is "").
#' @param ylab Y-axis label (default is "").
#' @param row.col Row colors (optional).
#' @param legend.flag Boolean to show legend (default is TRUE).
#' @param method Distance method for clustering (default is 'euclidean').
#' @param symm Boolean for symmetric distance matrix (default is FALSE).
#' @param palette Color palette for heatmap (default is "redblue").
#' @return Hierarchical clustering object of rows.
#' @export
call.heatmap<-function(m,main = '',col.labels = NULL,
                       row.labels = NULL,k = 3,filter.na = T,
                       cexCol = ifelse(ncol(m)>70,0.00001,1),
                       cexRow = ifelse(nrow(m)>70,0.00001,1),
                       m.value = '',scale = "none",
                       cluster.flag = "none",sym.flag = F,
                       xlab = "",ylab = "",row.col = NULL,legend.flag = T,
                       method = 'euclidean',symm = F, palette = "redblue"){
  if(filter.na){
    b.row<-!is.na(rowSums(m))
    m<-m[b.row,]
    if(!is.null(row.labels)){row.labels<-subset(row.labels,b.row)}
    b.col<-!is.na(colSums(m))
    m<-m[,b.col]
    if(!is.null(col.labels)){col.labels<-subset(col.labels,b.col)}
  }
  
  if(cluster.flag!="none"){
    if(method =="cor"){
      hc <- hclust(as.dist(2-cor(m)), method="complete");
      hr <- hclust(as.dist(2-cor(t(m))), method="complete");
    }else{
      hc <- hclust(dist(t(m),method = method), method="complete");
      hr <- hclust(dist(m,method =method), method="complete");
    }
    if(sym.flag){
      hc<-hr
    }
    Rowv <- as.dendrogram(hr)
    Colv <- as.dendrogram(hc)
    if(!is.null(col.labels)){
      col.labels<-cbind.data.frame(col.labels,
                                   clusters = paste0("C",cutree(hc, k = k)))
    }else{
      col.labels<-cbind.data.frame(clusters = paste0("C",cutree(hc, k = k)))
    }
    
  }else{
    Rowv <- NA;Colv <- NA
    hc<-T;hr<-T
  }
  
  col.col <- NULL
  if(!is.null(col.labels)){
    col.labels<-as.data.frame(col.labels)
    col.col<-(apply(col.labels,2,labels.2.colors))
    colnames(col.col)<-colnames(col.labels)
  }
  if(!is.null(row.labels)&is.null(row.col)){
    row.labels<-as.data.frame(row.labels)
    row.col<-t(
      as.matrix(t(laply(1:ncol(row.labels),
                        function(i) return(labels.2.colors(row.labels[,i]))))))
    if(ncol(row.labels)==1){row.col<-t(row.col)}
    row.labels<-t(as.matrix(row.labels))
    rownames(row.col)<-rownames(row.labels)
  }
  
  myheatcol <- redblue(50)
  myheatcol <- c(rep(myheatcol[1],10),myheatcol,rep(myheatcol[length(myheatcol)],10))
  myheatcol<-myheatcol[seq(length(myheatcol),1,-1)]
  plot.heatmap(m,main,Rowv=Rowv, Colv=Colv,
               m.value,cexRow,cexCol,myheatcol,scale,col.col,row.col,
               cluster.flag = cluster.flag,
               xlab = xlab,ylab = ylab,symm = symm)
  if(!legend.flag||(is.null(col.labels)&&is.null(row.labels))){return(col.labels)}
  legend.place<-rep(c("topright","right","bottomright","left",
                      "bottomleft","topleft"),5)
  if(!is.null(row.labels)&&identical(col.labels,t(row.labels))){
    row.labels<-NULL
  }
  
  if(!is.null(col.labels)){
    coltitles<-colnames(col.labels)[1]
    col.col<-as.matrix(col.col[,!duplicated(t(col.labels))])
    col.labels<-as.matrix(col.labels[,!duplicated(t(col.labels))])
    n1<-ifelse(is.matrix(col.labels),ncol(col.labels),1)
    
    for(i in 1:n1){
      v<-unique(paste(col.col[,i],col.labels[,i],sep = "?"))
      legend(legend.place[i],legend = get.strsplit(v,"?",2),col=get.strsplit(v,"?",1),
             pch = 15,cex = 0.4,title = coltitles[i],title.adj = 0.1)
    }
    legend.place<-legend.place[(i+1):length(legend.place)]
  }
  if(!is.null(row.labels)){
    rowtitles<-rownames(row.labels)
    n1<-ifelse(is.matrix(row.labels),nrow(row.labels),1)
    for(i in 1:n1){
      v<-unique(paste(row.col[i,],row.labels[i,],sep = "?"))
      legend(legend.place[i],legend = get.strsplit(v,"?",2),
             col=get.strsplit(v,"?",1),pch = 15,cex = 0.4,
             title = rowtitles[i],title.adj = 0.1)
    }
  }
  return(hr)
}

#' Plot Heatmap
#'
#' Helper function to create a heatmap with specified parameters.
#'
#' @param m Matrix of values to be plotted.
#' @param main Title of the plot.
#' @param Rowv Row dendrogram.
#' @param Colv Column dendrogram.
#' @param m.value Value to display in the key.
#' @param cexRow Size of row labels.
#' @param cexCol Size of column labels.
#' @param myheatcol Color palette for heatmap.
#' @param scale Scaling method.
#' @param col.col Column colors (optional).
#' @param row.col Row colors (optional).
#' @param cluster.flag Clustering method.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param symm Boolean for symmetric distance matrix.
#' @return None. This function is used for its side effect of printing a heatmap.
#' @export
plot.heatmap<-function(m,main,Rowv,Colv,m.value,cexRow,cexCol,
                       myheatcol,scale,col.col = NULL,row.col = NULL,
                       cluster.flag = "none",
                       xlab = "",ylab = "",symm = F){
  ColSideColorsSize = 1
  if(cluster.flag=="none"){Rowv = NA;Colv = NA}
  if(cluster.flag=="row"){Colv = NA}
  if(cluster.flag=="col"){Rowv = NA}
  if(is.null(col.col)&is.null(row.col)){
    heatmap.3(m, main = main,  Rowv = Rowv,Colv = Colv,
              col = myheatcol, density.info="none",margins=c(10,10),
              ColSideColorsSize = ColSideColorsSize,
              RowSideColorsSize = ColSideColorsSize,scale = scale,#RowAxisColors = 1,
              key=TRUE,
              KeyValueName = m.value,
              symm = symm,cexRow=cexRow,cexCol=cexCol,
              dendrogram = cluster.flag,
              xlab = xlab,ylab = ylab,symbreaks = T)
    return()
  }
  if(is.null(col.col)){col.col<-as.matrix(rep("black",ncol(m)));ColSideColorsSize = 0.5}
  if(is.null(row.col)){row.col<-t(as.matrix(rep("black",nrow(m))));RowSideColorsSize = 0.5}
  heatmap.3(m, main = main,  Rowv = Rowv,Colv = Colv,
            col=myheatcol, density.info="none",margins=c(10,10),
            ColSideColorsSize = ColSideColorsSize,
            RowSideColorsSize = ColSideColorsSize,scale = scale,#RowAxisColors = 1,
            key=TRUE , KeyValueName = m.value,symm = symm,cexRow=cexRow,cexCol=cexCol,
            ColSideColors=col.col,RowSideColors = row.col,dendrogram = cluster.flag,
            xlab = xlab,ylab = ylab)
}

#' Plot AUC
#'
#' Plots the ROC curve and calculates the AUC.
#'
#' @param p1 Predicted values.
#' @param y1 True values.
#' @param main Title of the plot (default is "").
#' @param subplotF Boolean to create subplot (default is TRUE).
#' @param precF Boolean to plot precision-recall curve (default is FALSE).
#' @param add Boolean to add to existing plot (default is FALSE).
#' @param col Color of the plot (default is "black").
#' @param plot.flag Boolean to plot the curve (default is TRUE).
#' @return AUC value.
#' @export
plot.auc<-function(p1,y1,main = "",subplotF = T,precF = F,add = F,
                   col = "black",plot.flag = T){
  if(subplotF){
    par(mfrow=c(1,2),oma = c(0, 0, 3, 0))
  }
  pr <- ROCR::prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  if(!plot.flag){return(prf)}
  plot(prf,main = paste0(main,'\nAUC =',round(auc,digits = 2)),cex.main=1,add = add,col = col)
  abline(a = 0,b = 1,col = "gray30")
  if(precF){
    prf <- performance(pr, measure = "prec", x.measure = "rec")
    plot(prf,ylim = c(0,1))
    abline(h = mean(y1))
  }
  return(auc)
}

#' Plot Multiple ROCs
#'
#' Plots multiple ROC curves with AUC values.
#'
#' @param P List of predicted values.
#' @param Y List of true values.
#' @param b Boolean vector indicating group membership.
#' @param main Title of the plot (default is "").
#' @return None. This function is used for its side effect of printing a plot.
#' @export
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

#' Violin Split
#'
#' Creates a split violin plot.
#'
#' @param scores Vector of scores.
#' @param treatment Factor vector of treatment groups.
#' @param conditions Factor vector of conditions.
#' @param main Title of the plot (default is "").
#' @param xlab X-axis label (default is "Sample").
#' @param ylab Y-axis label (default is "Scores").
#' @param legend.flag Boolean to show legend (default is TRUE).
#' @param show.pval Boolean to show p-value (default is FALSE).
#' @param col1 Color for the first group (default is "lightblue").
#' @param cex.axis Size of axis labels (default is 1).
#' @return None. This function is used for its side effect of printing a plot.
#' @export
violin.split<-function(scores, treatment, conditions, main = "",
                       xlab = "Sample",ylab = "Scores",legend.flag = T,
                       show.pval = F,col1 = "lightblue",cex.axis = 1){
  if(length(unique(conditions))==1){
    p<-t.test.mat(m = rbind(scores,scores),
                  b = treatment == sort(treatment,decreasing = T)[1])[1,1]
  }else{
    p<-t.test.groups(x = rbind(scores,scores),
                     b = treatment == sort(treatment,decreasing = T)[1],
                     g = conditions)[1,]
    p<-p[sort(names(p))]
  }
  if(show.pval){
    p<-my.format.pval(2*(10^-abs(p)))
    main<-paste(main,p,sep = "\n")
  }
  treatment<-as.factor(treatment)
  beanplot(scores ~ treatment*conditions, ll = 0.0,las = 2,
           main = main, side = "both", xlab=xlab,ylab = ylab,
           col = list(c(col1, "black"),"gray"),
           axes=T,cex.main = 1,cex.axis = cex.axis)
  if(legend.flag){
    legend("bottomright", fill = c(col1,"gray"),
           legend = levels(treatment), box.lty=0)
  }
}

#' Call Dot Plot
#'
#' Generates a dot plot.
#'
#' @param X Data frame for plotting.
#' @param cex Size of the points (default is 12).
#' @return ggplot object for the dot plot.
#' @export
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

#' Cap Matrix
#'
#' Caps the values in a matrix to specified quantiles.
#'
#' @param M Matrix of values.
#' @param cap Quantile cap value (default is 0.01).
#' @param MARGIN Margin to apply the function (1 for rows, 2 for columns).
#' @return Matrix with capped values.
#' @export
cap.mat<-function(M,cap = 0.01,MARGIN = 1){
  Z<-apply(M,MARGIN = MARGIN,function(x){
    q9<-quantile(x,1-cap)
    q1<-quantile(x,cap)
    x[x>q9]<-q9;x[x<q1]<-q1
    return(x)
  })
  if(MARGIN==1){Z<-t(Z)}
  return(Z)
}

#' Call Discretize
#'
#' Discretizes a vector into specified number of categories.
#'
#' @param v Numeric vector to be discretized.
#' @param n.cat Number of categories.
#' @param q1 Quantile values for discretization.
#' @return Discretized vector.
#' @export
call.discretize<-function(v,n.cat,q1){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)),na.rm = T)
  u<-matrix(data = 1,nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<=q1[i])]<-i
  }
  return(u)
}

#' Discretize to 3 Labels
#'
#' Discretizes a vector or matrix into three labels: High, Moderate, Low.
#'
#' @param X Numeric vector or matrix.
#' @param q Quantile threshold for high/low (default is 0.1).
#' @param verbose Boolean to print messages (default is FALSE).
#' @return Discretized vector or matrix with three labels.
#' @export
discretize.3.labels<-function(X,q = 0.1,verbose = F){
  if(q>0.5){q<-(1-q)}
  if(verbose){
    print(paste("Low quantile<=",q))
    print(paste("High quantile>=",1-q))
  }
  
  f<-function(v){
    b.low<-v<=quantile(v,q,na.rm = T)
    b.high<-v>=quantile(v,1-q,na.rm = T)
    labels<-ifelse(b.high,"High",ifelse(b.low,"Low","Moderate"))
    if(!any(is.na(labels))){
      labels<-factor(labels,levels = c("High","Moderate","Low"))
    }
    return(labels)
  }
  if(!is.matrix(X)){return(f(X))}
  B<-apply(X,2,f)
  return(B)
}

#' Get Top Correlations
#'
#' Retrieves the top correlations from a matrix.
#'
#' @param m Matrix of correlation values.
#' @param q Number of top correlations to retrieve (default is 100).
#' @param min.ci Minimum correlation value (default is 0).
#' @param idx Indices to retrieve specific correlations (optional).
#' @param add.prefix Prefix to add to the correlation names (default is "").
#' @param sort.flag Boolean to sort the correlations (default is TRUE).
#' @return List of top correlations.
#' @export
get.top.cor<-function(m,q = 100,min.ci = 0,idx = NULL,
                      add.prefix ="",sort.flag = T){
  m<-as.matrix(m)
  if(is.null(colnames(m))){colnames(m)<-1:ncol(m)}
  m.pos<-(-m);m.neg<-m
  colnames(m.pos)<-paste0(colnames(m.pos),".up")
  colnames(m.neg)<-paste0(colnames(m.neg),".down")
  v<-get.top.elements(cbind(m.pos,m.neg),q,min.ci = (-abs(min.ci)),sort.flag = sort.flag)
  names(v)<-c(colnames(m.pos),colnames(m.neg))
  if(!is.null(idx)){
    v<-v[paste(idx,c("up","down"),sep = ".")]
  }
  names(v)<-paste0(add.prefix,names(v))
  return(v)
}

#' Get Top Elements
#'
#' Retrieves the top elements from a matrix based on specified criteria.
#'
#' @param m Matrix of values.
#' @param q Number of top elements to retrieve (default is 100).
#' @param min.ci Minimum value to consider (optional).
#' @param main Prefix for the output list (optional).
#' @param sort.flag Boolean to sort the elements (default is TRUE).
#' @return List of top elements.
#' @export
get.top.elements<-function(m,q = 100,min.ci = NULL,
                           main = "",sort.flag = T){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    if(sort.flag){
      top.l[[i]]<-sort(v[b])
    }else{
      top.l[[i]]<-v[b][order(m[b,i])]
    }
    
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}

#' Intersect Lists by Index
#'
#' Intersects two lists element-wise by index.
#'
#' @param l1 First list.
#' @param l2 Second list.
#' @param remove.empty Boolean to remove empty intersections (default is FALSE).
#' @return List of intersected elements.
#' @export
intersect.lists.by.idx<-function(l1,l2,remove.empty = F){
  L<-lapply(1:length(l1), function(x) intersect(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  if(!remove.empty){return(L)}
  L<-L[laply(L,length)>0]
  return(L)
}

#' Intersect List with Elements
#'
#' Intersects a list with a vector of elements.
#'
#' @param l List of elements.
#' @param g Vector of elements to intersect with.
#' @param n1 Minimum number of elements in the intersection (default is 0).
#' @return List of intersected elements.
#' @export
intersect.list1<-function(l,g,n1=0){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[laply(l1,length)>n1]
  return(l1)
}

#' Union of Lists
#'
#' Combines two lists element-wise by index.
#'
#' @param l1 First list.
#' @param l2 Second list.
#' @param unique.flag Boolean to keep unique elements only (default is TRUE).
#' @param disregard.names Boolean to disregard names of the lists (default is FALSE).
#' @return List of combined elements.
#' @export
union.lists<-function(l1,l2,unique.flag = T,disregard.names = F){
  if(disregard.names){
    names(l2)<-names(l1)
  }else{
    idx<-robust.match(names(l1),names(l2))
  }
  if(unique.flag){
    L<-lapply(names(l1), function(x) unique(sort(c(l1[[x]],l2[[x]]))))
  }else{
    L<-lapply(names(l1), function(x) c(l1[[x]],l2[[x]]))
  }
  
  names(L)<-names(l1)
  return(L)
}

#' Robust Match
#'
#' Matches elements of two vectors after processing them to a common format.
#'
#' @param v1 First vector.
#' @param v2 Second vector.
#' @return Index of matches in the second vector.
#' @export
robust.match<-function(v1,v2){
  rmv<-c('_',"-",'.'," ",":","+")
  v1<-casefold(multi.gsub(pattern = rmv,replacement = '_',x = v1))
  v2<-casefold(multi.gsub(pattern = rmv,replacement = '_',x = v2))
  idx<-match(v1,v2)
  return(idx)
}

#' Multiple Global Substitution
#'
#' Performs multiple global substitutions on a character vector.
#'
#' @param pattern Vector of patterns to replace.
#' @param replacement Replacement string (default is '').
#' @param x Character vector to process.
#' @return Character vector with substitutions made.
#' @export
multi.gsub<-function(pattern,replacement = '',x){
  for(i in 1:length(pattern)){
    x<-gsub(pattern = pattern[i],replacement = replacement ,x = x,fixed = T)
  }
  return(x)
}

#' Get Substring by Index
#'
#' Splits strings by a separator and returns a specific part.
#'
#' @param v Character vector to split.
#' @param sep Separator string.
#' @param idx Index of the substring to return.
#' @return Character vector of substrings.
#' @export
get.strsplit<-function(v,sep,idx){
  v<-as.character(v)
  vi<-laply(strsplit(v,split = sep,fixed = T),function(x) x[idx])
  return(vi)
}

#' t-test Matrix
#'
#' Performs t-tests on each row of a matrix.
#'
#' @param m Matrix of values.
#' @param b Boolean vector indicating group membership.
#' @param two.sided Boolean for two-sided t-test (default is FALSE).
#' @param rankf Boolean to return ranks (default is FALSE).
#' @param fold.changeF Boolean to return fold change (default is FALSE).
#' @param BH.flag Boolean to adjust p-values using Benjamini-Hochberg (default is FALSE).
#' @return Data frame of t-test results.
#' @export
t.test.mat<-function(m,b,two.sided=F,rankf = F,fold.changeF = F,BH.flag = F){
  if(length(b)!=ncol(m)){
    print("Error. Inconsistent no. of samples.")
    return()
  }
  if(sum(b)<2||sum(!b)<2){
    return(get.mat(rownames(m),c('more','less',"zscores")))
  }
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) t.test(x[b],x[!b])$p.value))
  }else{
    p<-t(apply(m,1,function(x) c(t.test(x[b],x[!b],alternative = 'greater')$p.value,
                                 t.test(x[b],x[!b],alternative = 'less')$p.value)))
    colnames(p)<-c('more','less')
    p<-cbind(p,get.p.zscores(p))
    colnames(p)[3]<-"zscores"
  }
  if(rankf){
    p<-cbind(p,rank(p[,1]),rank(p[,2]))
    colnames(p)[4:5]<-c("rank.more","rank.less")
  }
  if(fold.changeF){
    p<-cbind.data.frame(p,pos.mean = rowMeans(m[,b]),neg.mean = rowMeans(m[,!b]))
    p$logFC<-log2(p$pos.mean/p$neg.mean)
  }
  if(BH.flag){
    p<-as.data.frame(p)
    p$BH.more<-p.adjust(p$more,method = "BH")
    p$BH.less<-p.adjust(p$less,method = "BH")
  }
  
  return(p)
}

#' t-test Groups
#'
#' Performs t-tests on groups defined by a factor.
#'
#' @param x Matrix of values.
#' @param b Boolean vector indicating group membership.
#' @param g Factor vector defining groups.
#' @param min.n Minimum number of samples per group (default is 1).
#' @param combine.n Threshold for combining p-values (optional).
#' @return Data frame of t-test results.
#' @export
t.test.groups<-function(x,b,g,min.n = 1,combine.n){
  x<-as.matrix(x)
  gu<-intersect(get.abundant(g[!b],min.n),get.abundant(g[b],min.n))
  if(is.null(rownames(x))){
    rownames(x)<-1:nrow(x)
  }
  v<-get.mat(rownames(x),gu)
  for (i in 1:length(gu)){
    b.g<-is.element(g,gu[i]);
    v[,i]<-t.test.mat(x[,b.g],b[b.g])[,3]
  }
  if(!missing(combine.n)){
    v<-cbind.data.frame(Z.up = rowSums(v>combine.n),
                        Z.down = rowSums(v<(-combine.n)),
                        P.up = fisher.combine(get.pval.from.zscores(v)),
                        P.down = fisher.combine(get.pval.from.zscores(v)),v)
  }
  return(v)
}

#' Get Z-scores from p-values
#'
#' Converts p-values to z-scores. 
#'
#' @param p Data frame of p-values.
#' @return Vector of z-scores.
#' @export
get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

#' Get Matrix
#'
#' Creates a matrix with specified row and column names and initial values.
#'
#' @param m.rows Vector of row names.
#' @param m.cols Vector of column names.
#' @param data Initial value for the matrix elements (default is NA).
#' @return Matrix with specified dimensions and initial values.
#' @export
get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

#' Add Number of Samples to Labels
#'
#' Adds the number of samples to labels.
#'
#' @param l Vector of labels.
#' @param n.flag Boolean to include "n=" prefix (default is TRUE).
#' @param sep Separator string (default is " ").
#' @return Vector of labels with sample numbers added.
#' @export
add.n.of.samples<-function(l,n.flag = T,sep = " "){
  num.samples<-table(l)
  idx<-match(l,names(num.samples))
  if(n.flag){
    l<-paste0(l,sep,"(n = ",num.samples[idx],")")
  }else{
    l<-paste0(l,sep,"(",num.samples[idx],")")
  }
  return(l)
}

#' List to Matrix
#'
#' Converts a list of vectors to a matrix.
#'
#' @param l List of vectors.
#' @return Matrix with elements of the list.
#' @export
list.2.mat<-function(l){
  n1<-max(laply(l,length))
  m<-t(laply(l,function(x) c(x,matrix(data = "",nrow = n1-length(x)+1))))
  m<-m[1:n1,]
  colnames(m)<-names(l)
  return(m)
}

#' Spearman Correlation
#'
#' Calculates the Spearman correlation between two vectors or matrices.
#'
#' @param v1 First vector or matrix.
#' @param v2 Second vector or matrix (optional).
#' @param method Correlation method (default is 'spearman').
#' @param use Method for handling missing values (default is "pairwise.complete.obs").
#' @param match.flag Boolean to match columns by name (default is FALSE).
#' @param alternative Alternative hypothesis for correlation test (default is "two.sided").
#' @param upper.tri.flag Boolean to return upper triangle of the correlation matrix (default is FALSE).
#' @return List of correlation results.
#' @export
spearman.cor<-function(v1,v2 = NULL,method = 'spearman',
                       use = "pairwise.complete.obs",match.flag = F,
                       alternative = "two.sided",upper.tri.flag = F){
  if(is.null(v2)){
    v2<-v1
  }
  if(!is.matrix(v1)){v1<-as.matrix(v1)}
  if(!is.matrix(v2)){v2<-as.matrix(v2)}
  if(match.flag){
    n=ncol(v1)
    if(is.null(colnames(v1))){colnames(v1)<-1:ncol(v1)}
    results<-get.mat(m.cols = c("R","P"),m.rows = colnames(v1))
    for(i in 1:ncol(v1)){
      c.i <- cor.test(v1[,i],v2[,i],method = method,
                      use = use, alternative = alternative)
      results[i,1] <- c.i$estimate
      results[i,2] <- c.i$p.value
    }
  }else{
    n1=ncol(v1)
    m<-matrix(nrow = n1,ncol = ncol(v2))
    rownames(m)<-colnames(v1)
    colnames(m)<-colnames(v2)
    results<-list(cor = m, p = m)
    for(i in 1:n1){
      f<-function(x){
        c.i<-cor.test(v1[,i],x,method = method,
                      use = use, alternative = alternative);
        c(c.i$estimate,c.i$p.value)}
      c.i <- apply(v2,2,f)
      results$cor[i,] <- c.i[1,]
      results$p[i,] <- c.i[2,]
    }
    if(ncol(v2)==1){
      results<-cbind(results$cor,results$p)
      colnames(results)<-c('R','P')
    }
  }
  if(upper.tri.flag){
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}

#' Labels to Colors
#'
#' Converts labels to colors.
#'
#' @param x.class Vector of labels.
#' @param x Vector of values (optional).
#' @param number.flag Boolean to return numeric values (default is FALSE).
#' @param color.spec Color specification for the output (default is "hsv").
#' @return Vector of colors corresponding to the labels.
#' @export
labels.2.colors<-function(x.class,x = NULL,number.flag = F,color.spec = "hsv"){
  palette("default")
  call_col<-c("black","red","cadetblue","gray",
              "darkgreen","darkorange","darkviolet","gold3",
              "lightpink","deeppink2","deepskyblue",
              palette(),rainbow(20))
  
  no.classes<-length(unique(x.class))
  if(number.flag){
    call_col<-match(x.class,sort(unique(x.class)))
  }else{
    if(is.numeric(x.class[1])){
      call_col<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = color.spec)
    }else{
      call_col<-call_col[match(x.class,sort(unique(x.class)))]
    }
  }
  if(!is.null(x)){names(call_col)<-x  }
  return(call_col)
}

#' Get AUC
#'
#' Calculates the area under the ROC curve.
#'
#' @param p1 Predicted values.
#' @param y1 True values.
#' @return AUC value.
#' @export
get.auc<-function(p1,y1){
  pr <- prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  return(auc)
}

#' Call Format P-Value
#'
#' Formats p-values for display.
#'
#' @param p Vector of p-values.
#' @param prnt.flag Boolean to print p-values (default is FALSE).
#' @param d Separator string (default is "=").
#' @return Formatted p-value string.
#' @export
call.format.pval<-function(p,prnt.flag = F,d = "="){
  if(length(p)>1){
    P<-laply(p,my.format.pval)
    P<-gsub("P = ","",P)
    P<-paste("P =",paste(P,collapse = ", "))
    return(P)
  }
  if(abs(p)>1){
    p<-10^(-abs(p))
  }
  if(p>0.05){p<-paste("P",d,round(p,3));return(p)}
  p<-gsub("e","*10",paste("P",d,format(p,scientific = T,digits= 3)))
  p<-gsub("-0","-",p)
  if(prnt.flag){
    p<-paste0("(",p,")")
  }
  return(p)
}

#' Get Abundant Elements
#'
#' Retrieves abundant elements from a vector.
#'
#' @param v Vector of values.
#' @param abn.c Minimum abundance count (default is 2).
#' @param boolean.flag Boolean to return a logical vector (default is FALSE).
#' @param top Index of top elements to retrieve (optional).
#' @param decreasing Boolean to sort in decreasing order (default is TRUE).
#' @return Vector of abundant elements or logical vector.
#' @export
get.abundant<-function(v,abn.c = 2,boolean.flag = F,top,decreasing = T){
  m<-as.matrix(table(v))
  m<-as.matrix(m[order(m,decreasing = decreasing),])
  if(!missing(top)){
    abn.c<-m[top]
  }
  m<-m[m>=abn.c,]
  abn.names<-names(m)
  if(boolean.flag){
    b<-is.element(v,abn.names)
    return(b)
  }
  return(abn.names)
}

#' "Not in" Magic
#'
#' Negation of the `%in%` Magic function
#'
#' @export
`%!in%` <- Negate(`%in%`)

#' Spatial Co-Occurrence
#'
#' Hypergeometric Tests for spatial co-occurrence in single-cell data.
#'
#' @param r Data frame with single-cell data.
#' @param type Type of co-occurrence to test ("Att" for attraction, "Rej" for rejection).
#' @return Matrix of p-values for co-occurrence tests.
#' @export
spatial.co.occur<-function(r,type = "Att"){
  if(!identical(rownames(r$frames.metadata),rownames(r$frames.tme))){
    print("Error. The frames in r$frames.metadata do not match those in r$frames.tme.")
    return()
  }
  p<-t(combn(colnames(r$frames.tme),2))
  f1<-function(x,y){
    return(get.hyper.p.value(x,y)[1])
  }
  f<-function(x,p){
    b<-r$frames.metadata$samples==x
    Y<-r$frames.tme[b,]
    if(type == "Att"){
      J<-apply(p,1,function(p1) f1(Y[,p1[1]]>0,Y[,p1[2]]>0))
      # J <- p.adjust(J, "BH")
    }else{
      J<-apply(p,1,function(p1) f1(Y[,p1[1]]>0,Y[,p1[2]]==0))
      # J <- p.adjust(J, "BH")
    }
    return(J)
  }
  J<-t(plyr::laply(unique(r$samples),function(x) f(x,p)))
  colnames(J)<-unique(r$samples)
  rownames(J)<-paste(p[,1],p[,2],sep = "_")
  J[is.na(J)] <- 1
  return(J)
}

#' Transfer Data to Seurat Object
#'
#' Transfers metadata from a list to a Seurat object.
#'
#' @param r List of metadata.
#' @param so Seurat object.
#' @param transfer_list List of fields to transfer (default includes several common fields).
#' @return Seurat object with updated metadata.
#' @export
transfer_data_list_to_so <- function(r,
                                     so,
                                     transfer_list = c("coor",
                                                       "samples",
                                                       "TMAs",
                                                       "patients",
                                                       "samples",
                                                       "sites",
                                                       "treatment",
                                                       "cell.types",
                                                       "cell.types2"))
{
  for (field in transfer_list) {
    if (is.null(dim(r[[field]]))) {
      so@meta.data[[field]] <- r[[field]]
    }
    else {
      so@meta.data <- cbind(so@meta.data, r[[field]])
    }
  }
  return(so)
}

#' Seuratify
#'
#' Creates a Seurat object from count data.
#'
#' @param counts Matrix of count data.
#' @param prj_name Project name (default is "").
#' @return Seurat object.
#' @export
seuratify <- function(counts, prj_name = "") {
  so <- Seurat::CreateSeuratObject(counts, project = prj_name)
  so@meta.data$cellid <- row.names(so@meta.data)
  return(so)
}

#' Cap Object
#'
#' Caps the values in an object to specified quantiles.
#'
#' @param X Numeric vector or matrix.
#' @param q Quantile cap value.
#' @return Vector or matrix with capped values.
#' @export
cap_object <- function (X, q) {
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}

#' RGB to Hex
#'
#' Converts RGB values to hex color code.
#'
#' @param x Vector of RGB values.
#' @return Hex color code.
#' @export
rgb2hex <- function (x) {
  hex = rgb(x[1], x[2], x[3], maxColorValue = 255)
  return(hex)
}

#' Get Overall Expression Scores Wrapper
#'
#' Wrapper that calculates overall expression scores for gene sets.
#'
#' @param r Data frame with single-cell data.
#' @param sig List of gene sets.
#' @return Matrix of scores.
#' @export
get.OE <- function (r, sig){
  scores <- get.OE1(r, sig)
  names(sig) <- gsub(" ", ".", names(sig))
  two.sided <- unique(gsub(".up", "", gsub(".down", "", names(sig))))
  b <- is.element(paste0(two.sided, ".up"), names(sig)) & is.element(paste0(two.sided,
                                                                            ".down"), names(sig))
  if (any(b)) {
    two.sided <- two.sided[b]
    scores2 <- as.matrix(scores[, paste0(two.sided, ".up")] -
                           scores[, paste0(two.sided, ".down")])
    colnames(scores2) <- two.sided
    scores <- cbind(scores2, scores)
  }
  if (!is.null(r$cells)) {
    rownames(scores) <- r$cells
  }
  else {
    if (!is.null(r$samples)) {
      rownames(scores) <- r$samples
    }
  }
  return(scores)
}

#' Get Overall Expression Scores
#'
#' Calculates overall expression scores for gene sets.
#'
#' @param r Data frame with single-cell data.
#' @param sig List of gene sets.
#' @return Matrix of scores.
#' @export
get.OE1 <- function (r, sig){
  if (is.list(sig)) {
    scores <- t(plyr::laply(sig, function(g) get.OE1(r, g)))
    rownames(scores) <- r$cells
    colnames(scores) <- names(sig)
    return(scores)
  }
  g <- sig
  b <- is.element(r$genes, g)
  assertthat::is.string(rownames(r$binZ)[1])
  n1 <- plyr::laply(rownames(r$binZ), function(x) sum(b[r$genes.dist.q ==
                                                          x]))
  rand.scores <- t(r$binZ) %*% n1
  if (sum(b) == 1) {
    raw.scores <- r$zscores[b, ]
  }
  else {
    raw.scores <- colSums(r$zscores[b, ])
  }
  scores <- (raw.scores - rand.scores)/sum(b)
  return(scores)
}

#' Scale and Center
#'
#' Scales and centers a matrix.
#'
#' @param X Matrix of values.
#' @param MARGIN Margin to apply the function (1 for rows, 2 for columns).
#' @return Scaled and centered matrix.
#' @export
scale_and_center <- function (X, MARGIN){
  X_norm = t(apply(X, MARGIN, function(x) {
    loc = mean(x, na.rm = T)
    center = x - loc
    scale = center/sd(x)
    return(scale)
  }))
  return(X_norm)
}

#' Subset List
#'
#' Subsets a list based on specified cells.
#'
#' @param r List of data.
#' @param subcells Vector of cell IDs to subset.
#' @return Subsetted list.
#' @export
subset_list <- function (r, subcells) {
  n_cells <- length(r$cells)
  n_genes <- length(r$genes)
  q <- lapply(r, function(x) {
    if (is.null(dim(x))) {
      if (length(x) == n_cells) {
        return(x[r$cells %in% subcells])
      }
      else {
        return(x)
      }
    }
    else if (dim(x)[2] == n_cells) {
      return(x[, r$cells %in% subcells])
    }
    else if (dim(x)[1] == n_cells) {
      return(x[r$cells %in% subcells, ])
    }
    else {
      return(x)
    }
  })
  return(q)
}

#' Spearman Correlation with Detailed Output
#'
#' Calculates the Spearman correlation between two vectors or matrices with detailed output.
#'
#' @param v1 First vector or matrix.
#' @param v2 Second vector or matrix (optional).
#' @param method Correlation method (default is 'spearman').
#' @param use Method for handling missing values (default is "pairwise.complete.obs").
#' @param match.flag Boolean to match columns by name (default is FALSE).
#' @param alternative Alternative hypothesis for correlation test (default is "two.sided").
#' @param upper.tri.flag Boolean to return upper triangle of the correlation matrix (default is FALSE).
#' @return List of correlation results.
#' @export
spearman.cor <- function (v1, v2 = NULL,
                          method = "spearman",
                          use = "pairwise.complete.obs",
                          match.flag = F,
                          alternative = "two.sided",
                          upper.tri.flag = F)
{
  if (is.null(v2)) {
    v2 <- v1
  }
  if (!is.matrix(v1)) {
    v1 <- as.matrix(v1)
  }
  if (!is.matrix(v2)) {
    v2 <- as.matrix(v2)
  }
  if (match.flag) {
    n = ncol(v1)
    if (is.null(colnames(v1))) {
      colnames(v1) <- 1:ncol(v1)
    }
    results <- get.mat(m.cols = c("R", "P"), m.rows = colnames(v1))
    for (i in 1:ncol(v1)) {
      c.i <- cor.test(v1[, i], v2[, i], method = method,
                      use = use, alternative = alternative)
      results[i, 1] <- c.i$estimate
      results[i, 2] <- c.i$p.value
    }
  }
  else {
    n1 = ncol(v1)
    m <- matrix(nrow = n1, ncol = ncol(v2))
    rownames(m) <- colnames(v1)
    colnames(m) <- colnames(v2)
    results <- list(cor = m, p = m)
    for (i in 1:n1) {
      f <- function(x) {
        c.i <- cor.test(v1[, i], x, method = method,
                        use = use, alternative = alternative)
        c(c.i$estimate, c.i$p.value)
      }
      c.i <- apply(v2, 2, f)
      results$cor[i, ] <- c.i[1, ]
      results$p[i, ] <- c.i[2, ]
    }
    if (ncol(v2) == 1) {
      results <- cbind(results$cor, results$p)
      colnames(results) <- c("R", "P")
    }
  }
  if (upper.tri.flag) {
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}

#' Get Top Correlations with Detailed Output
#'
#' Retrieves the top correlations from a matrix with detailed output.
#'
#' @param m Matrix of correlation values.
#' @param q Number of top correlations to retrieve (default is 100).
#' @param min.ci Minimum correlation value (default is 0).
#' @param idx Indices to retrieve specific correlations (optional).
#' @param add.prefix Prefix to add to the correlation names (default is "").
#' @param sort.flag Boolean to sort the correlations (default is TRUE).
#' @return List of top correlations.
#' @export
get.top.cor <- function (m, q = 100, min.ci = 0, idx = NULL, add.prefix = "") {
  m <- as.matrix(m)
  if (is.null(colnames(m))) {
    colnames(m) <- 1:ncol(m)
  }
  m.pos <- (-m)
  m.neg <- m
  colnames(m.pos) <- paste0(colnames(m.pos), ".up")
  colnames(m.neg) <- paste0(colnames(m.neg), ".down")
  v <- get.top.elements(cbind(m.pos, m.neg), q, min.ci = (-abs(min.ci)))
  names(v) <- c(colnames(m.pos), colnames(m.neg))
  if (!is.null(idx)) {
    v <- v[paste(idx, c("up", "down"), sep = ".")]
  }
  names(v) <- paste0(add.prefix, names(v))
  return(v)
}

#' Spatial Sample Visualization
#'
#' Creates a visualization of spatial samples with cell types.
#'
#' @param seg_path Path to the segmentation file.
#' @param celltypes Vector of cell types.
#' @param cell2rgb Mapping of cell types to RGB colors.
#' @param samplename Name of the sample.
#' @param background Background color (default is "black").
#' @param outpath Output path for the visualization (default is "~/").
#' @param outfile Output file name (default is "out.jpg").
#' @param cont_field Continuous field for coloring (optional).
#' @param low_qc_color Color for low-quality cells (default is 0).
#' @param contvals Continuous values for coloring (optional).
#' @return None. This function is used for its side effect of creating a visualization.
#' @export
spatial_sample_visualization <- function (seg_path,
                                          celltypes,
                                          cell2rgb,
                                          samplename,
                                          background = "black",
                                          outpath = "~/",
                                          outfile = "out.jpg",
                                          cont_field = "",
                                          low_qc_color = 0,
                                          contvals = NULL)
{
  cellseg = read.csv(seg_path)
  colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {
    strsplit(x, split = "X0.")[[1]][2]
  }))
  colnames(cellseg)[1] <- "0"
  cellmask <- cellseg
  cellmask[cellmask != 0] <- low_qc_color
  if (background == "white") {
    cellmask[cellmask == 0] <- 1
  }
  else {
    cellmask[cellmask == 0] <- 0
  }
  cellmask <- EBImage::Image(as.matrix(cellmask))
  cellmask <- EBImage::channel(cellmask, "rgb")
  levels <- setdiff(unique(celltypes), cont_field)
  for (type in levels) {
    cellids = names(celltypes[celltypes == type])
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x,
                                                                     split = "c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x) {
      x %in% cellidx
    }))
    celltype_rgb_1 <- cellmask@.Data[, , 1]
    celltype_rgb_1[celltype_mask] <- as.numeric(cell2rgb[[type]][1])/255
    cellmask[, , 1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[, , 2]
    celltype_rgb_2[celltype_mask] <- as.numeric(cell2rgb[[type]][2])/255
    cellmask[, , 2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[, , 3]
    celltype_rgb_3[celltype_mask] <- as.numeric(cell2rgb[[type]][3])/255
    cellmask[, , 3] <- celltype_rgb_3
  }
  if (cont_field != "") {
    cellids = names(celltypes[celltypes == cont_field])
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x,
                                                                     split = "c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x) {
      x %in% cellidx
    }))
    cellids <- paste0(samplename, "_c", cellseg[celltype_mask])
    colormap <- t(sapply(unique(contvals), simplify = T,
                         grDevices::col2rgb))/255
    celltype_rgb_1 <- cellmask@.Data[, , 1]
    celltype_rgb_1[celltype_mask] <- colormap[contvals[cellids],
                                              1]
    cellmask[, , 1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[, , 2]
    celltype_rgb_2[celltype_mask] <- colormap[contvals[cellids],
                                              2]
    cellmask[, , 2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[, , 3]
    celltype_rgb_3[celltype_mask] <- colormap[contvals[cellids],
                                              3]
    cellmask[, , 3] <- celltype_rgb_3
  }
  
  EBImage::writeImage(cellmask, files = outfile)
}

#' Intersect Lists
#'
#' Intersects two lists element-wise.
#'
#' @param sig1 First list of elements.
#' @param sig2 Second list of elements.
#' @return List of intersected elements.
#' @export
intersect.lists <- function(sig1, sig2){
  out <- lapply(1:length(sig1), function(x){
    intersect(sig1[[x]], sig2[[x]])
  })
  names(out) <- names(sig1)
  return(out)
}

#' De Novo Cell Type Markers for scRNA
#'
#' Identifies de novo cell type markers for single-cell RNA data.
#'
#' @param r Data frame with single-cell data.
#' @param n.non.mal Minimum number of non-malignant cells (default is 50).
#' @param q.dr Quantile threshold for dropout rate (default is 0.2).
#' @return List of identified cell type markers.
#' @export
scRNA_denovo.cell.type.markers<-function(r,n.non.mal,q.dr = 0.2){
  get.FC1<-function(x1,x2){
    Z<-log2(gene.av[sig[[x1]],x1]/gene.av[sig[[x1]],x2])
    return(Z)
  }
  get.ttest1<-function(x1,x2){
    b1<-is.element(r$cell.types,x1)
    b2<-is.element(r$cell.types,x2)
    b<-b1|b2
    Z<-t.test.mat(r$tpm[sig[[x1]],b],b1[b])
    return(Z)
  }
  get.ttest<-function(x1){
    if(length(sig[[x1]])<1){
      return(NA)
    }
    Z<-t(laply(cell.types,function(x2){
      get.ttest1(x1,x2)[,3]
    }))
    colnames(Z)<-cell.types
    b.mal<-startsWith(cell.types,"Malignant")
    Z<-cbind.data.frame(n = rowSums(Z>10,na.rm = T),
                        n.mal = rowSums(Z[,b.mal]>10,na.rm = T),
                        n.non.mal = rowSums(Z[,!b.mal]>10,na.rm = T),
                        n.non.T = rowSums(Z[,!b.tcell]>10,na.rm = T),Z)
    return(Z)
  }
  get.FC<-function(x1){
    Z<-t(laply(cell.types,function(x2) get.FC1(x1,x2)))
    colnames(Z)<-cell.types
    b.mal<-startsWith(cell.types,"Malignant")
    Z<-cbind.data.frame(n = rowSums(Z>0.2,na.rm = T),
                        n.mal = rowSums(Z[,b.mal]>0.2,na.rm = T),
                        n.non.mal = rowSums(Z[,!b.mal]>0.2,na.rm = T),
                        n.non.T = rowSums(Z[,!b.tcell]>0.2,na.rm = T),Z)
    return(Z)
  }
  
  r$b.mal<-r$cell.types=="Malignant"
  b<-!r$b.mal|is.element(r$samples,get.abundant(r$samples[r$b.mal],abn.c = 50))
  b<-b&get.abundant(r$cell.types,abn.c = 50,boolean.flag = T)
  # print(paste("Removing",sum(!b),"cells."))
  r<-set.list(r,b)
  
  r$cell.types[r$b.mal]<-paste(r$cell.types[r$b.mal],r$patients[r$b.mal],sep = "_")
  cell.types<-unique(r$cell.types)
  gene.av <- t(laply(cell.types,function(x) return(rowMeans(r$tpm[,r$cell.types==x]))))
  gene.dr <- t(laply(cell.types,function(x) return(rowMeans(r$tpm[,r$cell.types==x]>0))))
  
  colnames(gene.av)<-cell.types
  colnames(gene.dr)<-cell.types
  genes<-r$genes
  sigDR<-apply(gene.dr,2,function(x) genes[x>q.dr])
  n.mal<-min(length(unique(r$patients[r$b.mal])),40)
  sigDR$Malignant<-get.abundant(unlist(sigDR[grepl("Malignant",names(sigDR))]),n.mal)
  # print(summary(sigDR))
  sig<-sigDR
  
  b.mal<-grepl("Malignant",cell.types)
  b.tcell<-is.element(cell.types,c("T.cell","CD4.T","CD8.T"))
  if(missing(n.non.mal)){n.non.mal<-(sum(!b.mal)-1)}
  n.non.T<-sum(!b.tcell)-1
  
  Z1<-lapply(cell.types,get.FC)
  names(Z1)<-cell.types;summary(Z1)
  # Identify genes which are up-regulated in a non-malignant cell type compared
  # to all other (or at least n.non.mal) non-malignant cells types, and compared to malignant cells in at least n.mal patients
  sig<-lapply(Z1, function(X) sort(rownames(X)[X$n.mal>=n.mal&X$n.non.mal>=n.non.mal]));summary(sig)
  sig[b.tcell]<-lapply(Z1[b.tcell], function(X) sort(rownames(X)[X$n.non.T>=n.non.T&X$n.mal>=n.mal]))
  sig[b.mal]<-lapply(Z1[b.mal], function(X) sort(rownames(X)[X$n.non.mal>=n.non.mal]))
  sig$Malignant<-get.abundant(unlist(sig[b.mal]),n.mal)
  # print(summary(sig[c("Malignant",cell.types[!b.mal])]))
  sigFC<-sig
  
  Z2<-lapply(cell.types,get.ttest)
  names(Z2)<-cell.types
  sig<-lapply(Z2, function(X){
    if(identical(X,NA)){return(NA)}
    return(sort(rownames(X)[X$n.mal>=n.mal&X$n.non.mal>=n.non.mal]))});summary(sig)
  sig[b.tcell]<-lapply(Z2[b.tcell], function(X) sort(rownames(X)[X$n.non.T>=n.non.T&X$n.mal>=n.mal]))
  sig[b.mal]<-lapply(Z2[b.mal], function(X) sort(rownames(X)[X$n.non.mal>=n.non.mal]))
  sig$Malignant<-sort(get.abundant(v = unlist(sig[b.mal]),
                                   abn.c = min(sum(b.mal),40),boolean.flag = F))
  # print(summary(sig[c("Malignant",cell.types[!b.mal])]))
  
  sig.strict<-intersect.lists(sigFC,sigDR)
  sig.strict<-intersect.lists(sig,sig.strict)
  sig.strict<-sig.strict[laply(sig.strict,length)>0]
  rslts<-list(cohort = r$cohortName,
              gene.dr = gene.dr,
              gene.av = gene.av,
              FC = Z1,ttest = Z2,
              sigDR = sigDR,
              sigFC = sigFC,
              sig.strict = sig.strict,
              sig.mal = sig[grepl("Malignant_",names(sig))],
              sig = sig[!grepl("Malignant_",names(sig))])
  if(sum(is.element(cell.types,c("Myeloid","DC")))==2){
    rslts<-scRNA_denovo.cell.type.markers.similar.cell.types(rslts,"Myeloid","DC")
  }
  max.nonmal<-rowMax(rslts$gene.dr[,!grepl("Malignant",colnames(rslts$gene.dr))])
  rslts$sig$Malignant.strict<-rslts$sig$Malignant[max.nonmal[rslts$sig$Malignant]<0.2]
  # print(summary(rslts$sig))
  return(rslts)
}

#' Set List
#'
#' Sets elements of a list based on a boolean vector.
#'
#' @param r List of data.
#' @param b Boolean vector indicating elements to keep.
#' @param name Name for the subsetted list (optional).
#' @return Subsetted list.
#' @export
set.list<-function (r,b,name){
  set.field<-function (v,b){
    d <- dim(v)
    d.b<-length(b)
    if(!is.null(d)){
      if(d[1]==d.b){v <- subset(v,subset = b)}
      if(d[2]==d.b){v <- v[,b]}
    }else{if(length(v)==d.b){v <- v[b]}}
    return(v)
  }
  rn<-lapply(r, set.field, b = b)
  if(!missing(name)){rn$name<-name}
  return(rn)
}

#' Row Maximum
#'
#' Calculates the maximum value for each row of a matrix.
#'
#' @param X Matrix of values.
#' @return Vector of row maximum values.
#' @export
rowMax <- function (X) {
  y <- apply(X, 1, function(x) max(x, na.rm = T))
  return(y)
}

#' Row Minimum
#'
#' Calculates the minimum value for each row of a matrix.
#'
#' @param m Matrix of values.
#' @return Vector of row minimum values.
#' @export
rowMin <- function (m) {
  return(-rowMax(-m))
}

#' Prepare Data for OE Calculation
#'
#' Prepares data for overall expression score calculation.
#'
#' @param r Data frame with single-cell data.
#' @param n.cat Number of categories for discretization (default is 50).
#' @return Data frame with prepared data.
#' @export
prep4OE <- function (r, n.cat = 50) {
  r$zscores <- center.matrix(r$tpm, dim = 1, sd.flag = T)
  X <- 10 * ((2^r$tpm) - 1)
  r$genes.dist <- log2(rowMeans(X, na.rm = T) + 1)
  r$genes.dist.q <- discretize.prvt(r$genes.dist, n.cat = n.cat)
  b <- rowSums(is.na(r$zscores)) == 0
  if (any(!b)) {
    r <- set.list(r, b)
  }
  r$binZ <- average.mat.rows(r$zscores, r$genes.dist.q, f = colMeans)
  return(r)
}

#' Discretize Data
#'
#' Discretizes a vector into specified number of categories
#'
#' @param v Numeric vector to be discretized.
#' @param n.cat Number of categories.
#' @param q1 Quantile values for discretization.
#' @return Discretized vector.
#' @export
discretize.prvt <- function (v, n.cat, q1) {
  q1 <- quantile(v, seq(from = (1/n.cat), to = 1, by = (1/n.cat)),
                 na.rm = T)
  u <- matrix(data = 1, nrow = length(v))
  for (i in 2:n.cat) {
    u[(v >= q1[i - 1]) & (v <= q1[i])] <- i
  }
  u <- paste0("Q", u)
  return(u)
}

#' Center Matrix
#'
#' Centers a matrix by subtracting the mean.
#'
#' @param m Matrix of values.
#' @param dim Dimension to center (1 for rows, 2 for columns).
#' @param sd.flag Boolean to standardize by standard deviation (default is FALSE).
#' @return Centered matrix.
#' @export
center.matrix <- function (m, dim = 1, sd.flag = F) {
  if (dim == 1) {
    zscores <- sweep(m, 1, rowMeans(m, na.rm = T), FUN = "-")
  }
  else {
    zscores <- sweep(m, 2, colMeans(m, na.rm = T), FUN = "-")
  }
  if (sd.flag) {
    zscores <- sweep(zscores, dim, apply(m, dim, function(x) (sd(x,
                                                                 na.rm = T))), FUN = "/")
  }
  return(zscores)
}

#' Average Matrix Rows
#'
#' Averages rows of a matrix based on specified IDs.
#'
#' @param m Matrix of values.
#' @param ids Vector of IDs.
#' @param f Function to apply to each group (default is colMeans).
#' @return Averaged matrix.
#' @export
average.mat.rows <- function (m, ids, f = colMeans)
{
  ids.u <- sort(unique(ids))
  m1 <- get.mat(ids.u, colnames(m))
  for (x in ids.u) {
    b <- is.element(ids, x)
    if (sum(b) == 1) {
      m1[x, ] <- m[b, ]
    }
    else {
      m1[x, ] <- f(m[b, ])
    }
  }
  return(m1)
}

#' Cast Sites
#'
#' Casts a data frame by anatomical sites and specific columns.
#'
#' @param df Data frame to cast.
#' @param column_idx Indices of columns to cast.
#' @return Cast data frame.
#' @export
cast_sites <- function(df, column_idx){
  out <- Reduce(function(x, y) merge(x, y, by = c("patients", "treatment")),
                lapply(column_idx, function(x){
                  col = colnames(df)[x]
                  out <- df[,c(1:3, x)] %>%
                    spread(sites_binary, .data[[col]])
                  colnames(out)[3:4] <- paste0(col, "_", colnames(out)[3:4])
                  return(out)
                }))
  return(out)
}

#' Get Frames
#'
#' Retrieves frames for a single-cell spatial data list object
#'
#' @param r Data frame with single-cell spatial data.
#' @param n1 Number of frames in total
#' @return Data frame with updated frames.
#' @export
get.frames<-function(r,n1){
  r$x<-ceil(r$coor[,1]/n1)
  r$y<-ceil(r$coor[,2]/n1)
  b<-r$samples==r$samples[1]
  # my.plot(r$coor[b,],labels = paste(r$x,r$y)[b])
  r$frames<-paste0(r$samples,"_X",r$x,"_Y",r$y)
  return(r)
}

#' Get Hypergeometric P-Value
#'
#' Calculates the hypergeometric p-value for two sets.
#'
#' @param b1 Boolean vector for the first set. 
#' @param b2 Boolean vector for the second set.
#' @param full.flag Boolean to return full output (default is TRUE).
#' @return Vector of p-values and expected values.
#' @export
get.hyper.p.value<-function(b1,b2,full.flag = T){
  p1<-NA;p2<-NA;e<-0;
  if(any(b1)&&any(b2)){
    p1<-max(1-phyper(sum(b1&b2)-1, sum(b1), sum(!b1), sum(b2)),1e-17)
    e<-sum(b2)*(sum(b1)/length(b2))
    p2<-sum(b1&b2)/e
  }
  if (full.flag){
    p<-c(p1,p2,sum(b1&b2),e)
    names(p)<-c('hyper.p.value','ob.vs.exp','ob','exp')
    return(p)  
  }
  return(p1)
}

#' Get One-Sided P-Value
#'
#' Calculates one-sided p-values from z-scores.
#'
#' @param c Vector of z-scores.
#' @param p Vector of p-values.
#' @return Vector of one-sided p-values.
#' @export
get.onesided.p.value <- function(c,p){
  p[p==0] <-1e-17
  p.one.side <- p
  p.one.side[] <- NA
  b<-c>0&!is.na(c)
  p.one.side[b]=p[b]/2
  b<-c<=0&!is.na(c)
  p.one.side[b]=1-(p[b]/2)
  return(p.one.side)
}

#' Scores to Colors
#'
#' Converts scores to colors.
#'
#' @param x.class Vector of scores.
#' @return Vector of colors corresponding to the scores.
#' @export
scores.2.colors<-function(x.class){
  palette("default")
  call_col<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = "hsv")
  return(call_col)
}

#' Create Volcano Plot
#'
#' Creates a volcano plot with specified parameters.
#'
#' @param x Vector of x-values.
#' @param y Vector of y-values.
#' @param dot_names Vector of dot names.
#' @param zcat Threshold for categorization (default is 1.3).
#' @param x_label Label for the x-axis (default is "X-axis").
#' @param y_label Label for the y-axis (default is "Y-axis").
#' @param quadrant_labels Labels for the quadrants (default is c("Quadrant 1", "Quadrant 2", "Quadrant 3", "Quadrant 4")).
#' @param cex Size of the text labels (default is 1).
#' @return ggplot object of the volcano plot.
#' @export
create_volcano_plot <- function(x, y,dot_names, zcat = 1.3,
                                x_label = "X-axis", y_label = "Y-axis",
                                quadrant_labels = c("Quadrant 1", "Quadrant 2", 
                                                    "Quadrant 3", "Quadrant 4"),
                                cex = 1) {
  # create a data frame
  df <- data.frame(x, y, dot_names)
  
  df$significant <- "No"
  df$significant[(x > 0 & y > zcat)|(x > zcat & y > 0)] <- "Resistance"
  df$significant[(x < 0 & y < (-zcat))|x < (-zcat) & y < 0] <- "Response"
  df$significant[(x < (-zcat) & y > 0)|(x < 0 & y > zcat)] <- "Opposing1"
  df$significant[(x > 0 & y < (-zcat))|x > zcat & y < 0] <- "Opposing2"
  
  
  # create the plot with four quadrants
  p <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 3) +
    scale_color_manual(values = c("No" = "grey", "Resistance" = "blue", "Response" = "green", "Opposing1" = "orange","Opposing2" = "red"), guide = FALSE) +
    geom_vline(xintercept = zcat, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = zcat, color = "gray50", linetype = "dashed") +
    geom_vline(xintercept = -zcat, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = -zcat, color = "gray50", linetype = "dashed") +
    theme_minimal() +
    labs(x = x_label, y = y_label) +
    annotate("text", x = max(df$x) - 0.1, y = max(df$y) - 0.1, label = quadrant_labels[1], hjust = 1, vjust = 1, color = "blue") +
    annotate("text", x = min(df$x) + 0.1, y = max(df$y) - 0.1, label = quadrant_labels[2], hjust = 0, vjust = 1, color = "orange") +
    annotate("text", x = min(df$x) + 0.1, y = min(df$y) + 0.1, label = quadrant_labels[3], hjust = 0, vjust = 0, color = "green") +
    annotate("text", x = max(df$x) - 0.1, y = min(df$y) + 0.1, label = quadrant_labels[4], hjust = 1, vjust = 0, color = "red") +
    coord_cartesian(xlim = range(df$x), ylim = range(df$y)) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    geom_text_repel(data = subset(df, significant != "No"), aes(label = dot_names), size = cex, box.padding = 0.5)+
    guides(color = guide_legend(title = "Significance"))
  return(p)
  
}

#' Seurat Wrapper for Embedding
#'
#' Wrapper function for creating embeddings using Seurat.
#'
#' @param r Data frame with single-cell data.
#' @param no.genes Number of highly variable genes to use (default is 2000).
#' @param n.pcs Number of principal components to use (default is 10).
#' @param cd.flag Boolean to use CD data (default is TRUE).
#' @param PCs Precomputed principal components (optional).
#' @param resolution Resolution for clustering (default is 0.4).
#' @param norm.flag Boolean to normalize data (default is TRUE).
#' @param tsne.flag Boolean to create t-SNE embedding (default is FALSE).
#' @param umap.flag Boolean to create UMAP embedding (default is TRUE).
#' @param full.flag Boolean to return full Seurat object (default is FALSE).
#' @return Data frame with updated embeddings.
#' @export
seuratW_get.embedding<-function(r,no.genes = 2000,n.pcs = 10,cd.flag = T,PCs = NULL,
                                resolution = 0.4,norm.flag = T,tsne.flag = F,
                                umap.flag = T,full.flag = F){
  
  D1<-seuratW_process(r = r,
                      no.genes = no.genes,
                      scores = NULL,
                      plot.flag = T,
                      n.pcs = n.pcs,
                      cd.flag = cd.flag,
                      PCA.approx = T,resolution = resolution,
                      norm.flag = norm.flag,
                      tsne.flag = tsne.flag,
                      umap.flag = umap.flag,
                      PCs = PCs)
  
  idx<-my.match(r$cells,rownames(D1@meta.data));table(is.na(idx))
  r$clusters<-paste0("C",FetchData(D1,vars = "seurat_clusters")[idx,])
  r$pca<-Embeddings(object = D1, reduction = "pca")[idx,]
  r$pca.load<-Loadings(object = D1,reduction = "pca")
  if(tsne.flag){r$tsne<-Embeddings(object = D1, reduction = "tsne")[idx,]}
  if(umap.flag){r$umap<-Embeddings(object = D1, reduction = "umap")[idx,]}
  if(full.flag){r$seurat<-D1}
  return(r)
}

#' Seurat Wrapper for Processing
#'
#' Wrapper function for processing data using Seurat.
#'
#' @param r Data frame with single-cell data.
#' @param no.genes Number of highly variable genes to use (default is 1000).
#' @param n.pcs Number of principal components to use (default is 10).
#' @param tsne.flag Boolean to create t-SNE embedding (default is TRUE).
#' @param umap.flag Boolean to create UMAP embedding (default is TRUE).
#' @param cluster.iter Boolean to iterate clustering (default is FALSE).
#' @param cd.flag Boolean to use CD data (default is FALSE).
#' @param scores Additional scores to include (optional).
#' @param plot.flag Boolean to create plots (default is FALSE).
#' @param fileName File name for saving plots (optional).
#' @param D1 Precomputed Seurat object (optional).
#' @param resolution Resolution for clustering (default is 0.4).
#' @param norm.flag Boolean to normalize data (default is TRUE).
#' @param tsne.method Method for t-SNE (default is "Rtsne").
#' @param PCA.approx Boolean to use approximate PCA (default is TRUE).
#' @param PCs Precomputed principal components (optional).
#' @return Processed Seurat object.
#' @export
seuratW_process<-function(r,no.genes = 1000, n.pcs = 10,tsne.flag = T,umap.flag = T,
                          cluster.iter = F,cd.flag = F,scores = NULL,plot.flag = F,fileName = NULL,
                          D1,resolution = 0.4,norm.flag = T,tsne.method = "Rtsne",PCA.approx = T,PCs){
  print(paste(n.pcs,"PCs used."))
  if(missing(D1)){
    if(cd.flag){D1.data<-round(r$cd,2)}else{D1.data<-round(r$tpm,2)}
    D1<-seuratW_createObject(r$genes,D1.data,r$sampleName,norm.flag = norm.flag)
  }
  
  m = D1@assays$RNA@meta.data
  row.names(m) <- r$genes
  hvg.D1 <- get.top.elements(-m[,1:4], q = no.genes)$vf_vst_counts_variance.standardized
  # D1 = FindVariableFeatures(D1, nfeatures = no.genes, method = "vst")
  # hvg.D1 = VariableFeatures(D1)
  D1 <- RunPCA(object = D1,features = hvg.D1,npcs = n.pcs,approx = PCA.approx)
  X<-Embeddings(object = D1, reduction = "pca")[,1:n.pcs]
  # v<-Loadings(object = D1, reduction = "pca")
  
  if(!missing(PCs)&!is.null(PCs)){
    print("Using previous PCs!")
    idx.pca<-match(rownames(D1@reductions$pca@cell.embeddings),
                   rownames(PCs$cell.embeddings))
    D1@reductions$pca@feature.loadings<-PCs$feature.loadings
    D1@reductions$pca@cell.embeddings<-PCs$cell.embeddings[idx.pca,]
  }
  
  b<-duplicated(X)
  if(any(b)){
    print(paste("Removing",sum(b),"cells."))
    D1<-subset(D1,cells = rownames(D1@meta.data)[!b])
  }
  
  if(tsne.flag){D1 <- RunTSNE(object = D1, reduction.use = "pca",dims = 1:n.pcs, do.fast = TRUE)}
  if(umap.flag){D1<-RunUMAP(D1,dims = 1:n.pcs)}
  
  D1 <- FindNeighbors(object = D1,reduction = "pca",dims = 1:n.pcs)
  D1 <- FindClusters(D1,resolution = resolution)
  # reduction.type = "pca",
  # dims.use = 1:n.pcs, save.SNN = T, resolution = resolution)
  if(cluster.iter){
    while(length(unique(D1@ident))==1){
      resolution<-resolution+0.05
      print(paste("Clustering with resolution",resolution))
      D1 <- FindClusters(D1,resolution = resolution)
    }
    D1@resolution<-resolution
  }
  
  idx<-my.match(rownames(D1@meta.data),colnames(D1.data))
  if(!is.null(r$samples)){D1@meta.data<-cbind.data.frame(D1@meta.data,samples = r$samples[idx])}
  if(!is.null(scores)){D1@meta.data<-cbind.data.frame(D1@meta.data,scores[idx,])}
  return(D1)
}

#' Seurat Wrapper for Creating Object
#'
#' Wrapper function for creating a Seurat object.
#'
#' @param genes Vector of gene names.
#' @param D1.data Matrix of data.
#' @param D1.name Name for the Seurat object.
#' @param find.var Boolean to find variable features (default is TRUE).
#' @param norm.flag Boolean to normalize data (default is TRUE).
#' @return Seurat object.
#' @export
seuratW_createObject<-function(genes,D1.data,D1.name,find.var = T,norm.flag = T){
  print(paste("Processing",D1.name,"..."))
  D1.data<-round(D1.data,2)
  D1.data<-D1.data[genes,!duplicated(colnames(D1.data))];dim(D1.data)
  D1 <- CreateSeuratObject(counts = D1.data)
  if(norm.flag){
    D1 <- NormalizeData(object = D1)
  }
  D1 <- ScaleData(object = D1)
  if(find.var){D1 <- FindVariableFeatures(object = D1, do.plot = FALSE)}
  D1@meta.data$source <- D1.name
  return(D1)
}

#' Set Difference of Lists by Index
#'
#' Computes the set difference of two lists element-wise by index.
#'
#' @param l1 First list.
#' @param l2 Second list.
#' @return List of set differences.
#' @export
setdiff.lists.by.idx<-function(l1,l2){
  L<-lapply(1:length(l1), function(x) setdiff(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  print(summary(L))
  return(L)
}

#' Apply Formula for Hierarchical Linear Model (Mixed Effects Model)
#'
#' Applies a hierarchical linear model (mixed effects model) formula to each 
#' row of a matrix.
#'
#' @param r Data frame with single-cell data.
#' @param X Matrix of predictor variables.
#' @param Y Matrix of response variables.
#' @param MARGIN Margin to apply the function (1 for rows, 2 for columns).
#' @param formula Formula for the hierarchical linear model.
#' @param ttest.flag Boolean to perform t-tests (default is FALSE).
#' @return Data frame of model results.
#' @export
apply.formula.HLM<-function(r,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x",ttest.flag = F){
  if(is.matrix(Y)){
    if(ttest.flag){
      m1<-t.test.mat(Y,X)
      b<-rowSums(p.adjust.mat(m1[,1:2])<0.1,na.rm = T)>0
      m1<-m1[b,];Y<-Y[b,]
      print(paste("Testing",sum(b),"genes that show a signal."))
    }
    m<-t(apply(Y,MARGIN = MARGIN,function(y){formula.HLM(y,X,r,formula = formula)}))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,function(x){formula.HLM(Y,x,r,formula = formula)}))
  }
  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  if(ttest.flag){
    m<-cbind.data.frame(m,ttest = m1)
  }
  return(m)
}

#' Hierarchical Linear Model Formula (Mixed Effects Model)
#'
#' Applies a hierarchical linear model (mixed effects model) formula to data.
#'
#' @param y Response variable.
#' @param x Predictor variable.
#' @param r0 Data frame with additional variables.
#' @param formula Formula for the hierarchical linear model.
#' @param val Value for the predictor variable (default is determined automatically).
#' @param return.all Boolean to return all coefficients (default is FALSE).
#' @return Vector of model coefficients and p-values.
#' @export
formula.HLM<-function(y,x,r0, formula = "y ~ (1 | samples) + x",
                      val = ifelse(is.numeric(x),"","TRUE"),return.all = F){
  r0$x<-x;r0$y<-y
  f<-function(r0){
    M1 <- with(r0, lmer (formula = formula))
    if(return.all){
      c1<-summary(M1)$coef[,c("Estimate","Pr(>|t|)")]
    }else{
      c1<-summary(M1)$coef[paste0("x",val),]
      idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
      c1<-c1[idx]
    }
    return(cbind(c1,singular = isSingular(M1)))
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA,NA))})
  # print(c1)
  return(c1)
}
