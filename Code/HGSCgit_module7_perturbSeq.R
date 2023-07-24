#### Results Section 5 ###
# Figure 7. Perturb-Seq meta-analyses

# Figure 7B. CD8 T cell umaps
# Figure 7C. CD8 T cell TIP genes' association with tumor infiltration in different cell subtypes.

# Figure 7. High content CRISPR screen in ovarian cancer cells in monoculture and co-culture with targeting NK and T cells.
# 7A. Overview of experimental design.
# 7B. CRISPR-based differential fitness.
# 7C. Perturb-Seq UMAP.
# 7D. Perturb-Seq-based mTIL regulators.
# 7E. mTIL DEGs under different perturbations.
# 7F. KOs altering the response to NK cells.
# 7G. Differential expression of KO signatures in monoculture and co-culture
# 7H. Perturb-Seq UMAPs with cell colored based on gene KO signatures.

HGSC_Figure7_perturbOC<-function(r,rslts,sig,fitnessR){
  if(missing(r)){
    r<-readRDS("Data/PerturbSeq_TYKnuNK.rds")
    rslts<-readRDS("Results/PerturbSeq_TYKnuNK_DEGs.rds")
    sig<-readRDS("Results/HGSC_mTIL.rds")
    fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  }
  
  r$scores1c<-cap.mat(r$scores1,cap = 0.05,MARGIN = 2)
  r$condL<-r$cond
  r$condL[r$cond=="Mono"]<-"Monoculture"
  r$condL[r$cond=="1to1"]<-"Co-culture, 1:1"
  r$condL[r$cond=="25to1"]<-"Co-culture, 2.5:1"
  
  HGSC_Fig7B(fitnessR)
  HGSC_Fig7C(r = r)
  HGSC_Fig7D(rslts = rslts)
  HGSC_Fig7E(rslts = rslts,sig = sig)
  HGSC_Fig7F(rslts = rslts)
  HGSC_Fig7G(r = r,rslts = rslts)
  HGSC_Fig7H(r = r)
  
  
  return()
}

HGSC_Fig7B<-function(fitnessR){
  if(missing(fitnessR)){
    fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  }
  pdf(get.file("Figures/Fig7B.pdf"))
  p<-create_volcano_plot(x = fitnessR$Tcell.z, y = fitnessR$NK.z,
                         dot_names = rownames(fitnessR), x_label = "T cell selection", y_label = "NK selection",
                         quadrant_labels = c("Resistance", "T cell response\nNK resistance",
                                             "Response", "NK response\nT cell resistance"),zcat = 2)
  print(p)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig7C<-function(r){
  pdf(get.file("Figures/Fig7C.pdf"))
  print(umap.ggplot(r$umap,labels = add.n.of.samples(r$condL)))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig7D<-function(rslts){
  P<-rslts$sum$sum
  P<-P[P$N<2,]
  z<-P$Z;names(z)<-rownames(P)
  z[abs(z)>8]<-8*sign(z[abs(z)>8])
  
  pdf(get.file("Figures/Fig7D.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "mTIL differential expression",
          xlab = "Gene KO",main = "Perturbation impact on mTIL program",
          col = ifelse(z>2,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Activator","Repressor","NS"), pch = 15,
         col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig7E<-function(rslts,sig){
  sig1<-get.top.elements(rslts$sum$sum[,1:2],7,min.ci = 0.01)
  p.hot<-c(paste0(sig1$Repressor,".up"),paste0(sig1$Activators,".down"))
  genes<-get.abundant(unlist(intersect.list1(rslts$sig[p.hot],sig$mTIL.up)),5)
  ids<-c(unique(get.strsplit(p.hot[grepl("up",p.hot)],".",1)),"NTC",unique(get.strsplit(p.hot[grepl("down",p.hot)],".",1)))
  
  prt<-unlist(sig1)
  X<-cbind.data.frame(Mono = rslts$Mono$deg.ttest[match(genes,rownames(rslts$Mono$deg.ttest)),prt],
                      Co1 = rslts$Co1$deg.ttest[genes,prt],
                      Co2 = rslts$Co2$deg.ttest[genes,prt])
  rownames(X)<-genes
  X[is.na(X)]<-0
  X[abs(X)>4]<-4*sign(X[abs(X)>4])
  b<-is.element(get.strsplit(colnames(X),".",2),sig1$Repressor)
  
  hc1 <- hclust(dist(t(X[,b]),method = 'euclidean'), method="complete")
  hc2 <- hclust(dist(t(X[,!b]),method = 'euclidean'), method="complete")
  prt1<-hc1$labels[hc1$order]
  prt1<-c(prt1[grepl("Co",prt1)],prt1[!grepl("Co",prt1)])
  prt2<-hc2$labels[hc2$order]
  prt2<-c(prt2[grepl("Co",prt2)],prt2[!grepl("Co",prt2)])
  prtA<-c(prt1,prt2)
  
  X1<-X[,prtA]
  b1<-is.element(get.strsplit(prtA,".",2),sig1$Repressor)
  X.lab<-cbind.data.frame(mTIL = ifelse(b1,"Repressor","Activator"),
                          Condition = get.strsplit(prtA,".",1))
  
  pdf(get.file("Figures/Fig7E.pdf"))
  call.heatmap(X1,cluster.flag = "row",cexRow = 0.5,col.labels = X.lab,cexCol = 0.5,
               xlab = "Pertrubations",ylab = "Genes",m.value = "Z score",legend.flag = F)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
  
}

HGSC_Fig7F<-function(rslts){
  z<-rslts$prt.deNK1[,"zscores"]
  names(z)<-rownames(rslts$prt.deNK1)
  z<-sort(z)
  z[abs(z)>30]<-30*sign(z[abs(z)>30])
  
  pdf(get.file("Figures/Fig7F.pdf"))
  barplot(abs(z),las=2,cex.names = 0.3,ylab = "Differential expression,\nmono- vs. co-culture",
          xlab = "Perturbation signature",main = "",cex.axis = 0.1,
          col = ifelse(z>4,"blue",ifelse(z<(-2),"darkred","grey")))
  legend("top",legend = c("Higher","Lower","NS"), pch = 15,col = c("blue","darkred","grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig7G<-function(r,rslts){
  sig1<-get.top.elements(rslts$prt.deNK1[,c("BH.more","BH.less")],7,min.ci = 0.01)
  par(mfrow=c(2,3),oma = c(2, 1, 0, 5),xpd = T)
  X<-NULL
  for(x in unlist(sig1)){
    X1<-cbind.data.frame(scores = r$scores1[r$b.ctrl,x],cond = r$cond1[r$b.ctrl],
                         prt = x)
    X<-rbind(X,X1)
  }
  X1<-X[is.element(X$prt,sig1$BH.more),]
  X2<-X[!is.element(X$prt,sig1$BH.more),]
  
  pdf(get.file("Figures/Fig7G.pdf"))
  violin.split(scores = X1$scores, treatment = ifelse(X1$cond=="Mono","Monoculture","Co-culture"),
               conditions = X1$prt,ylab = "Pertrbation signature OE",xlab = "Perturbation",
               main = "Genetic pertrubation mimicking NK effects")
  violin.split(scores = X2$scores, treatment = ifelse(X2$cond=="Mono","Monoculture","Co-culture"),
               conditions = X2$prt,ylab = "Pertrbation signature OE",xlab = "Perturbation",
               main = "Genetic pertrubation counteracting NK effects")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

HGSC_Fig7H<-function(r){
  f1<-function(x){
    call.boxplot(r$scores1[,x],paste0(r$targets,r$cond1),cex = 0.5,
                 labels = paste0(r$cond1,
                                 ifelse(r$targets==x,paste(", ",x,"KO"),"")),xlab = "Perturbations",
                 ylab = paste0(x," signature OE"),
                 main = x)}
  
  X<-cbind.data.frame(Conditions = r$condL,
                      r$scores1c[,c("ACTR8.up","IRF1","MED12.up","STAT1")])
  colnames(X)<-gsub(".up","",colnames(X))
  
  pdf(get.file("Figures/Fig7H.pdf"))
  print(umap.ggplot(r$umap,labels = X[,-1],size = 0.1))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

create_volcano_plot <- function(x, y,dot_names, zcat = 1.3,x_label = "X-axis", y_label = "Y-axis",
                                quadrant_labels = c("Quadrant 1", "Quadrant 2", "Quadrant 3", "Quadrant 4"),cex = 1) {
  # Create a data frame
  df <- data.frame(x, y, dot_names)
  
  df$significant <- "No"
  df$significant[(x > 0 & y > zcat)|(x > zcat & y > 0)] <- "Resistance"
  df$significant[(x < 0 & y < (-zcat))|x < (-zcat) & y < 0] <- "Response"
  df$significant[(x < (-zcat) & y > 0)|(x < 0 & y > zcat)] <- "Opposing1"
  df$significant[(x > 0 & y < (-zcat))|x > zcat & y < 0] <- "Opposing2"
  
  
  # Create the plot with four quadrants
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
  # geom_text(data = subset(df, significant != "No"), aes(label = dot_names), vjust = -1, hjust = 0, size = 3)
  
  # p <- p + guides(color = guide_legend(title = "Significance", 
  #                                      override.aes = list(shape = c(16,17, 18, 19, 20))))
  return(p)
  
}

