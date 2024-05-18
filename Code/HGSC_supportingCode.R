umap.ggplot<-function(umapX,labels,labels.name = "",main = "",size = 0.2,xlim1,
                      ylim1,reorder.flag = F,remove.legend = F){
  if((is.matrix(labels)|is.data.frame(labels))&&ncol(labels)>1){
    p<-lapply(colnames(labels),function(x){
      print(x)
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

call.multiplot<-function(plotlist,nplots = 4,cols = 2){
  flag<-F
  while(!is.null(plotlist)&!flag){
    print(multiplot(plotlist = plotlist[1:min(nplots,length(plotlist))],cols = cols))
    flag<-(min(nplots,length(plotlist))+1)>length(plotlist)
    plotlist<-plotlist[(min(nplots,length(plotlist))+1):length(plotlist)]
  }
}

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

call.discretize<-function(v,n.cat,q1){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)),na.rm = T)
  u<-matrix(data = 1,nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<=q1[i])]<-i
  }
  return(u)
}

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

call.plot.plus<-function(x, y = NULL,labels,b.top,red.top = F,regression.flag = F,my.col = NULL,set.flag = F,cor.flag = F,
                         pch=16,cex=0.3,main="",ylab = "tSNE2",xlab = "tSNE1", cex.axis = 0.6,
                         add.N = F,grey.zeros = F,legend.flag = T){

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

call.plot<-function(x, y = NULL,labels,regression.flag = F,my.col = NULL,set.flag = F,cor.flag = F,legend.flag = T,
                    pch=16,cex=0.5,main="",ylab = "UMAP2",xlab = "UMAP1", cex.axis = 0.6,add.N = F,cex.main = 1,
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

call.boxplot<-function (y,x,unique.x, f = median,
                        ylab = '',xlab = '',main = '',labels=NULL,
                        legend.name = "",add.anova = F,blank.flag = T,
                        cex = 0.7,order.flag = T,b.ref = NULL,p.val.show = is.null(b.ref)){
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

get.top.cor<-function(m,q = 100,min.ci = 0,idx = NULL, add.prefix ="",sort.flag = T){
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

get.top.elements<-function(m,q = 100,min.ci = NULL,main = "",sort.flag = T){
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

intersect.lists.by.idx<-function(l1,l2,remove.empty = F){
  L<-lapply(1:length(l1), function(x) intersect(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  if(!remove.empty){return(L)}
  L<-L[laply(L,length)>0]
  return(L)
}

intersect.list1<-function(l,g,n1=0){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[laply(l1,length)>n1]
  return(l1)
}

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

robust.match<-function(v1,v2){
  rmv<-c('_',"-",'.'," ",":","+")
  v1<-casefold(multi.gsub(pattern = rmv,replacement = '_',x = v1))
  v2<-casefold(multi.gsub(pattern = rmv,replacement = '_',x = v2))
  idx<-match(v1,v2)
  return(idx)
}

multi.gsub<-function(pattern,replacement = '',x){
  for(i in 1:length(pattern)){
    x<-gsub(pattern = pattern[i],replacement = replacement ,x = x,fixed = T)
  }
  return(x)
}

get.strsplit<-function(v,sep,idx){
  v<-as.character(v)
  vi<-laply(strsplit(v,split = sep,fixed = T),function(x) x[idx])
  return(vi)
}

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
    v<-cbind.data.frame(Z.up = rowSums(v>combine.n),Z.down = rowSums(v<(-combine.n)),
                        P.up = fisher.combine(get.pval.from.zscores(v)),
                        P.down = fisher.combine(get.pval.from.zscores(v)),v)
  }
  return(v)
}

get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

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

list.2.mat<-function(l){
  n1<-max(laply(l,length))
  m<-t(laply(l,function(x) c(x,matrix(data = "",nrow = n1-length(x)+1))))
  m<-m[1:n1,]
  colnames(m)<-names(l)
  return(m)
}

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
      c.i <- cor.test(v1[,i],v2[,i],method = method,use = use, alternative = alternative)
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
        c.i<-cor.test(v1[,i],x,method = method,use = use, alternative = alternative);
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

call.heatmap<-function(m,main = '',col.labels = NULL,row.labels = NULL,k = 3,filter.na = T,
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
      col.labels<-cbind.data.frame(col.labels,clusters = paste0("C",cutree(hc, k = k)))
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
    row.col<-t(as.matrix(t(laply(1:ncol(row.labels),function(i) return(labels.2.colors(row.labels[,i]))))))
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

labels.2.colors<-function(x.class,x = NULL,number.flag = F,color.spec = "hsv"){
  palette("default")
  call_col<-c("black","red","cadetblue","gray","darkgreen","darkorange","darkviolet","gold3",
              "lightpink","deeppink2","deepskyblue",palette(),rainbow(20))

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

plot.auc<-function(p1,y1,main = "",subplotF = T,precF = F,add = F,col = "black",plot.flag = T){
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

get.auc<-function(p1,y1){
  pr <- prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  return(auc)
}

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

violin.split<-function(scores, treatment, conditions, main = "",xlab = "Sample",ylab = "Scores",legend.flag = T,
                       show.pval = F,col1 = "lightblue",cex.axis = 1){
  # require(beanplot)
  if(length(unique(conditions))==1){
    p<-t.test.mat(m = rbind(scores,scores),b = treatment == sort(treatment,decreasing = T)[1])[1,1]
  }else{
    p<-t.test.groups(x = rbind(scores,scores),b = treatment == sort(treatment,decreasing = T)[1],g = conditions)[1,]
    p<-p[sort(names(p))]
  }
  print(p)
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
  return(p)
}

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

`%!in%` <- Negate(`%in%`)

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

seuratify <- function(counts, prj_name = "") {
  so <- Seurat::CreateSeuratObject(counts, project = prj_name)
  so@meta.data$cellid <- row.names(so@meta.data)
  return(so)
}

cap_object <- function (X, q) {
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}

rgb2hex <- function (x) {
  hex = rgb(x[1], x[2], x[3], maxColorValue = 255)
  return(hex)
}

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

scale_and_center <- function (X, MARGIN){
  X_norm = t(apply(X, MARGIN, function(x) {
    loc = mean(x, na.rm = T)
    center = x - loc
    scale = center/sd(x)
    return(scale)
  }))
  return(X_norm)
}

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

get.top.elements <- function (m, q = 100, min.ci = NULL, main = "") {
  top.l <- list()
  v <- rownames(m)
  for (i in 1:ncol(m)) {
    mi <- m[, i]
    mi <- mi[!is.na(mi)]
    idx <- order(mi, decreasing = F)
    ci <- mi[idx[min(q, length(mi))]]
    ci <- min(ci, min.ci)
    b <- m[, i] <= ci
    b[is.na(m[, i])] <- F
    top.l[[i]] <- sort(v[b])
  }
  if (main != "") {
    main <- paste0(main, ".")
  }
  names(top.l) <- paste0(main, colnames(m))
  return(top.l)
}

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
    # print(type)
    cellids = names(celltypes[celltypes == type])
    # print(length(cellids))
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
    print(cont_field)
    print(length(cellids))
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

intersect.lists <- function(sig1, sig2){
  out <- lapply(1:length(sig1), function(x){
    intersect(sig1[[x]], sig2[[x]])
  })
  names(out) <- names(sig1)
  return(out)
}

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

rowMax <- function (X) {
  y <- apply(X, 1, function(x) max(x, na.rm = T))
  return(y)
}

rowMin <- function (m) {
  return(-rowMax(-m))
}

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

spatial.co.occur<-function(r,type = "Att"){
  # r - single cell data structure with the following frame information:
  #     1) frames.tme (k x m1) - the relative abundance of m1 different cell types across k frames
  #     2) frames.metadata (k x m2) - k frames metadata including "samples" information that denotes which sample each frame is from.
  # type -  whether the co-occurrence should be tested for "attraction" (i.e., higher than expected)
  #         or rejections ("i.e., lower than expected). Should be either "Att" or "Rej", respectively.
  
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

get.frames<-function(r,n1){
  r$x<-ceil(r$coor[,1]/n1)
  r$y<-ceil(r$coor[,2]/n1)
  b<-r$samples==r$samples[1]
  # my.plot(r$coor[b,],labels = paste(r$x,r$y)[b])
  r$frames<-paste0(r$samples,"_X",r$x,"_Y",r$y)
  return(r)
}

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

cap_object <- function (X, q) 
{
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}

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
