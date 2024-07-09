#### Results Section 3 ###
# Figure 3. MTIL (Malignant Transcriptional program that robustly marks the 
# presence of Infiltrating Lymphocytes)

# Figure 3a. MTIL heatmap
# Figure 3b. MTIL Gene Set Enrichment Analysis
# Figure 3c. MTIL spatial maps (composite made in Adobe Illustrator)
# Figure 3d: MTIL as a function of TIL proximity and abundance (boxplot)
# Figure 3e: MTIL ROCs
# Figure 3f: MTIL in Representative Whole Tissue (Test 2)
# Figure 3f: MTIL in Test 2 Dataset

#' Figure 3 Wrapper Function
#' Reproduce main text Figures 3a-f. 
#' @param r Discovery data.
#' @param r1 Discovery data, malignant cells.
#' @param Mtil.sig Mtil signature. 
#' @return Regenerates Figure 3 from (Yeh et al., 2024) in the designated Figures folder. 
HGSC_Figure3_mTIL<-function(r,r1,Mtil.sig){
  #1 Figure 3a: MTIL heatmap
  mTIL_Fig3a(r1 = r1,Mtil.sig = Mtil.sig)
  #2 Figure 3b: Gene Set Enrichment Analysis
  mTIL_Fig3b(r = r,Mtil.sig = Mtil.sig)
  #4 Figure 3c: MTIL spatial maps
  mTIL_Fig3c(r = r,r1 = r1, Mtil.sig = Mtil.sig)
  #5 Figure 3d: MTIL as a function of TIL proximity and abundance (boxplot)
  mTIL_Fig3d(r1 = r1)
  #6 Figure 3e: ROCs
  mTIL_Fig3e(r1 = r1)
  #7 Figure 4f: MTIL in Whole Tissue (Test 2)
  mTIL_Fig3f(r = r,mal = r1)
  
  return()
}

#' Figure 3a. MTIL gene expression heatmap in the Discovery dataset. 
#' @param r1 Discovery data, malignant cells.
#' @param Mtil.sig Mtil signature. 
#' @return Fig. 3a (Figures folder). 
mTIL_Fig3a<-function(r1,Mtil.sig){
  # get genes
  g<-intersect(c("CXCL10","CXCL2","CXCL5","CXCL9","LGALS1","NR3C1",
                 "BCL2","FGFR1","HDAC1","ITGB5","RELA","HDAC5"),
               unlist(Mtil.sig))
  
  # get expression matrix
  idx1<-order(r1$scoresAv[,"hot100"])
  idx2<-c(Mtil.sig$Mtil.up,Mtil.sig$Mtil.down)
  gene.cor<-spearman.cor(r1$tpmAv[,idx2],r1$scoresAv[,"hot"])
  sig1<-get.top.cor(gene.cor,q = 50,idx = "R")
  sig1<-union.lists(sig1,get.top.cor(gene.cor[g,],q = 50,idx = "R"),
                    unique.flag = T)
  idx2<-c(sig1$R.up,sig1$R.down)
  X1<-cap.mat(t(r1$tpmAv[idx1,idx2]),cap = 0.05,MARGIN = 1)
  
  # get column and row annotation labels
  col.labels <-cbind.data.frame(TILs = r1$tmeAv[,"TNK.cell"]>
                                  median(r1$tmeAv[,"TNK.cell"]),
                                CD4 = r1$tmeAv[,"CD4.T.cell"]>0,
                                CD8 = r1$tmeAv[,"CD8.T.cell"]>0,
                                NK = r1$tmeAv[,"NK.cell"]>0,
                                samples = get.strsplit(rownames(r1$tmeAv),
                                                       "_X",1))[idx1,]
  col.labels$sites<-r1$sites_binary[match(col.labels$samples,r1$samples)]
  row.labels <- ifelse(is.element(idx2,sig1$R.up),"Hot","Cold")
  row.labels <- cbind.data.frame(row.labels,row.labels)
  
  # print to disk
  pdf(get.file("Figures/Fig3a.pdf"))
  call.heatmap(X1,scale = "row",
               col.labels = col.labels,
               row.labels = row.labels,
               cexRow = 0.1,
               legend.flag = F)
  barplot(r1$scoresAv[idx1,"hot100"],
          cex.names = 1e-10,
          xlab = "Frames",
          ylab = "Immune-rich signature: Overall Expression")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 3b. MTIL gene ontology enrichment analysis 
#' @param r Discovery data, malignant cells.
#' #' @param Mtil.sig Mtil signature. 
#' @return Fig. 3b (Figures folder). 
mTIL_Fig3b <- function(r,Mtil.sig){
  # run mtil up GO enrichment
  results1 <- (gprofiler2::gost(Mtil.sig$Mtil.up,
                                custom_bg = r$genes[!grepl("NegPrb", r$genes)],
                                user_threshold = 0.1,
                                correction_method = "fdr",
                                sources = c("GO:BP", "GO:MF", "GO:CC")))
  
  # compile results into dataframe
  dfup <- results1$result %>%
    filter(term_id %!in% unlist(results1$result$parents))
  dfup <- dfup[c(6, 4, 3, 17, 28, 11, 42, 46),] %>%
    mutate(neglog10 = -log10(p_value)) %>%
    mutate(term = paste0(term_name, "\n(",term_id, ")")) %>%
    arrange(neglog10)
  dfup$term <- factor(dfup$term, levels = dfup$term)
  
  # run mtil down GO enrichment
  results2 <- gprofiler2::gost(Mtil.sig$Mtil.down,
                               custom_bg = r$genes[!grepl("NegPrb", r$genes)],
                               user_threshold = 0.1,
                               correction_method = "fdr",
                               sources = c("GO:BP", "GO:MF", "GO:CC"))
  
  # compile results into dataframe
  dfdown <- results2$result %>%
    filter(term_id %!in% unlist(results2$result$parents))
  dfdown <- dfdown[c(9,11,14, 18, 20, 21, 22, 19),]
  dfdown <- dfdown %>% mutate(neglog10 = -log10(p_value)) %>%
    mutate(term = paste0(term_name, "\n(",term_id, ")")) %>%
    arrange(neglog10)
  dfdown$term <- factor(dfdown$term, levels = dfdown$term)
  
  # format data frame for plotting
  dfup$direction = "Up"
  dfdown$direction = "Down"
  out <- rbind(dfup, dfdown)
  out$direction <- factor(out$direction, levels = c("Up", "Down"))
  out$intersection_size <- factor(unlist(lapply(out$intersection_size, function(x){
    if(x < 5) {
      return("<5")
    } else if ( x < 10) {
      return("6-10")
    } else if (x < 20) {
      return("11-20")
    } else {
      return(">21")
    }
  })), levels = c("<5", "6-10", "11-20", ">21"))
  
  # make plot
  p2 <- ggplot(out, aes(x = neglog10, y = term, col = source,
                        size = intersection_size)) +
    facet_wrap(~direction, scales = "free") +
    geom_point() +
    theme_classic() +
    guides(colour = guide_legend("GO Type" ,
                                 override.aes = list(fill="white",
                                                     size = 4,
                                                     shape = 15)),
           size = guide_legend("No. Genes")) +
    xlab("-log10(q_value)") +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_text(size = 12),
          strip.background = element_rect(colour = "white"))
  # plot to disk
  ggsave(get.file("Figures/Fig3b.pdf"), width = 15, height = 5, units = "in",plot = p2)
}

#' Figure 3c. MTIL gene ontology enrichment analysis plot.
#' @param r Discovery data.
#' @param r1 Discovery data, malignant cells.
#' #' @param Mtil.sig Mtil signature. 
#' @return Fig. 3c (Figures folder). 
mTIL_Fig3c <-function(r,r1, Mtil.sig = Mtil.sig){
  #1 Use the mTIL scores to color the malignant cells.
  X.ttest<-readRDS(get.file("Results/ST_Mtil_Discovery_ttest.per.sample.rds"))
  X<-get.mat(r1$cells,c("hot","hot100"),data = NA)
  
  for(x in unique(r1$samples)){
    b<-r1$samples==x
    X[b,]<-apply(cap.mat(r1$scores[b,c("hot","hot100")],
                         cap = 0.05,MARGIN = 2),2,
                 function(x){
                   p <- ggplot(data.frame(s = x, s = x),
                               aes(x = s,
                                   y = s,
                                   col = s)) +
                     geom_point() +
                     scale_color_gradient2(low = '#009392', 
                                           mid = '#f6edbd', 
                                           high = '#A01C00',
                                           midpoint = median(
                                             cap_object(x, 0.05))) +
                     coord_fixed()
                   colors <- data.frame(ggplot_build(p)$data[[1]])[,1]
                   return(colors)
                 })
  }
  
  #2 Color TNK cells in black and other cells in grey.
  idx<-match(r$cells,r1$cells)
  X1<-X[idx,]
  X1[is.na(idx)&r$cell.types=="TNK.cell",]<-"black" # T cells are marked in black
  X1[is.na(idx)&r$cell.types!="TNK.cell",]<-"grey" # Non-malignant cells other than T cells are shown in grey
  rownames(X1)<-r$cells
  
  #3 Plot the spatial maps per sample to regenerate Figure 4D. (centroids in situ)
  r$ids<-paste(r$patients,r$sites,sep = ", ")
  X.ttest$ids<-r$ids[match(rownames(X.ttest),r$samples)]
  samples<-c('SMI_T10_F001','SMI_T10_F015','SMI_T12_F001',
             'SMI_T12_F009','SMI_T12_F016','SMI_T13_F001')
  pdf(get.file("Figures/Fig3c.pdf"))
  par(mfrow=c(2,2),oma = c(0, 1, 0, 1),xpd = T)
  for(x in samples){
    b1<-r$samples==x
    call.plot.plus(-r$coor[b1,c("y","x")],labels = r$scores[b1,3],my.col = X1[b1,1],
                   b.top = r$cell.types[b1]=="TNK.cell",cex = 0.5,legend.flag = F,
                   main = paste(x,"\nZ =",round(X.ttest[x,"ttest_zero"],1),
                                "\nZ.median =",round(X.ttest[x,"ttest_median"],1),
                                X.ttest[x,"ids"]),
                   pch = ifelse(r$cell.types[b1]=="TNK.cell",15,16))
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  #4 Also generate space-filling spatial maps to regenerate Figure 3d.
  # a. load information
  cell2rgb <- list("Other" = c(200, 200, 200), # grey
                   "Malignant" = c(7, 224, 0), #green before color bar
                   "TNK.cell" = c(0, 0, 0)) # black
  
  # b. plot for each sample
  out <- lapply(samples, function(sample){
    # segpath
    seg_path = Sys.glob(paste0(get.file("Results/Segmentation/"), sample, "*"))[1]
    
    # subset out cells
    cells <- r$cells[r$samples == sample]
    q <- subset_list(r, cells)
    q$cell.types <- gsub("_LC", "", q$cell.types)
    celltypes <- q$cell.types
    names(celltypes) <- q$cells
    celltypes[celltypes != "TNK.cell" & celltypes != "Malignant"] <- "Other"
    
    # get scores for sample (new)
    mal <- subset_list(q, subcells = q$cells[q$cell.types == "Malignant"])
    mal$scores <- get.OE(mal, Mtil.sig)
    
    # set up the colors of the continuous values
    p <- ggplot(data.frame(s = cap_object(mal$scores[,"Mtil"], 0.05), 
                           mal$coor),
                aes(x = x,
                    y =y,
                    col = s)) +
      geom_point() +
      scale_color_gradient2(low = '#009392', 
                            mid = '#f6edbd', 
                            high = '#A01C00',
                            midpoint = median(
                              cap_object(mal$scores[,"Mtil"], 0.05))) +
      coord_fixed()
    colors <- data.frame(ggplot_build(p)$data[[1]])[,1]
    names(colors) <- names(mal$scores[,"Mtil"])
    leg <- ggpubr::as_ggplot(ggpubr::get_legend(p))
    contvals = colors
    
    if (sample == "SMI_T10_F001") {
      pdf(get.file("Figures/Fig3c_legend.pdf"),
          height = 4,
          width = 3)
      print(leg)
      dev.off()
    }
    
    # run visualization
    spatial_sample_visualization(seg_path,
                                 celltypes,
                                 cell2rgb,
                                 sample,
                                 background = "white",
                                 outpath = outpath,
                                 low_qc_color = 1,
                                 cont_field = "Malignant",
                                 contvals = contvals,
                                 outfile = get.file(paste0("Figures/Fig3c_",
                                                           sample,
                                                           ".png")))
    
    return("")
  })
  
  return()
}

#' Figure 3d. MTIL expression as a function of TIL proximity and abundance (boxplots)
#' @param r1 Discovery data, malignant cells.
#' @param q1 quartile cutoff for discretizing T/NK cell levels. 
#' @return Fig. 3d (Figures folder). 
mTIL_Fig3d<-function(r1,q1 = 0.75){
  P<-list(samples = r1$scoresSamples[,"hot100"],
          patches = r1$scoresAv[,"hot100"],
          cells = r1$scores[,"hot100"])
  L3<-list(samples = discretize.3.labels(r1$tmeSamples[,"TNK.cell"],1-q1),
           patches = discretize.3.labels(r1$tmeAv[,"TNK.cell"],1-q1),
           cells = discretize.3.labels(r1$tme[,"TNK.cell"],1-q1))
  titles<-c(samples = "Sample-level",
            patches = "Frame-level",
            cells = "Cell-level")
  
  # plot boxplots to compare TIL levels
  f<-function(x){
    l<-list(c("Low","Moderate"),c("Moderate","High"))
    p1<-call.boxplot(P[[x]],L3[[x]],ylab = "mTIL program OE",
                     xlab = "TIL levels",main = titles[x])
    p1<-p1+stat_compare_means(comparisons = l,
                              label = "p.signif",
                              method = "t.test") + 
      theme(legend.position="none")
    return(p1)
  }
  p<-lapply(names(L3),f)
  
  # plot to disk
  pdf(get.file("Figures/Fig3d.pdf"))
  p1<-call.boxplot(r1$scores[,"hot100"],
                   paste0("L",r1$envTIL1),main = "Fig3d",
                   ylab = "mTIL program OE",
                   xlab = "TIL levels (radius)",blank.flag = T,labels = NULL)
  p[[4]]<-p1+stat_compare_means(comparisons = list(c("L0","L1"),
                                                   c("L1","L2"),
                                                   c("L2","L3")),
                                label = "p.signif",method = "t.test")+
    theme(legend.position="none")
  print(call.multiplot(p[c(4,1,3,2)],cols = 3,nplots = 6))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

#' Figure 3e. MTIL as a predictor of TIL levels (ROC curves)
#' @param r1 Discovery data, malignant cells.
#' @param q1 quartile cutoff for discretizing T/NK cell levels. 
#' @return Fig. 3e (Figures folder). 
mTIL_Fig3e<-function(r1,q1 = 0.75){
  # prep scores and discretize T/NK cell levels
  P<-list(samples = r1$scoresSamples[,"hot100"],
          patches = r1$scoresAv[,"hot100"],
          cells = r1$scores[,"hot100"])
  Y<-list(samples = r1$tmeSamples[,"TNK.cell"]>quantile(r1$tmeSamples[,"TNK.cell"],q1),
          patches = r1$tmeAv[,"TNK.cell"]>quantile(r1$tmeAv[,"TNK.cell"],q1),
          cells = r1$tme[,"TNK.cell"]>quantile(r1$tme[,"TNK.cell"],q1))
  
  # plot ROC to disk
  pdf(get.file("Figures/Fig3e.pdf"))
  plot.multi.ROCs(P,Y)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

#' Figure 3f. MTIL in Representative Whole Tissue (Test 2)
#' Visualize MTIL in situ in whole tissue data from the Test 2 dataset. 
#' @param r Test 2 data, a single whole tissue section
#' @param mal Test 2 data, malignant cells
#' @return Fig. 3f (Figures folder). 
mTIL_Fig3f<- function(r,mal){
  # map MTIL onto whole tissue data.
  plt_mal <- data.frame(mal$coor, cap_object(mal$scores, q = 0.1))
  nonmal <- subset_list(r, subcells = r$cells[r$cell.types != "Malignant"])
  remove(r)
  plt_df <- data.frame(nonmal$coor, cell.types = nonmal$cell.types) %>% 
    mutate(cell.types = c("Other", "TNK.cell")[(cell.types == "TNK.cell") + 1]) 
  
  # make whole tissue spatial map based on cell centroids. 
  whole = ggplot() + 
    geom_point(data = filter(plt_df, cell.types == "Other"), 
               mapping = aes(x = x, y = y),
               color = "lightgrey", 
               size = 0.1) + 
    geom_point(data = plt_mal, 
               mapping = aes(x = x, y = y, col = hot100), 
               size = 0.1) + 
    scale_color_gradient2(low = '#009392', 
                          mid = '#f6edbd', 
                          high = '#A01C00', midpoint = 0) + 
    geom_point(data = filter(plt_df, cell.types == "TNK.cell"), 
               mapping = aes(x = x, y = y),
               color = "black", 
               size = 0.1) +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank(), 
          legend.position = "none") + 
    coord_fixed() 
  
  # write to disk
  png(get.file("Figures/Fig3f.png"),
      height = 16, width = 16, units = "in", res = 500)
  print(whole)
  dev.off()
}

#' Figure 3g. MTIL Overall Expression in Test 2 Dataset
#' 
#' This function generates boxplots of MTIL Overall Exrpression in whole 
#' tissue as a function of T/NK cell levels on the FOV-, Frame-, and Cell-level. 
#' 
#' @param i which whole tissue data to use. 
#' @return Fig. 3f (Figures folder). 
mTIL_Fig3g <- function(){
  out <- lapply(c("113", "8", "2", "1"), 
                function(i){readRDS(
                  get.file(paste0(
                    "Data/ST_Test2_HGSC", i , "_malignant.rds")))})
  
  q1 <- 0.75
  P<- list(samples = do.call("c", 
                             lapply(out, 
                                    function(r1){r1$scoresSamples[,"hot100"]})),
           patches = do.call("c", 
                             lapply(out, 
                                    function(r1){r1$scoresAv[,"hot100"]})), 
           cells = do.call("c", 
                           lapply(out, 
                                  function(r1){r1$scores[,"hot100"]})))
  L3<-list(
    samples = discretize.3.labels(
      do.call("c", lapply(out, function(r1){r1$tmeSamples[,"TNK.cell"]})),1-q1),
    patches = discretize.3.labels(
      do.call("c", lapply(out, function(r1){r1$tmeAv[,"TNK.cell"]})),1-q1),
    cells = discretize.3.labels(
      do.call("c", lapply(out, function(r1){r1$tme[,"TNK.cell"]})),1-q1))
  titles<-c(samples = "FOV-level",patches = "Frame-level",cells = "Cell-level")
  
  f<-function(x){
    print(x)
    l<-list(c("Low","Moderate"),c("Moderate","High"))
    p1<-call.boxplot(P[[x]],L3[[x]],ylab = "MTIL program OE",
                     xlab = "TIL levels",main = titles[x], cex = 1)
    p1<-p1+stat_compare_means(comparisons = l,
                              label = "p.format",
                              method = "t.test", 
                              method.args = 
                              list(alternative = "less"), 
                              p.format = "%.5e")+
      theme(legend.position="none")
    return(p1)
  }
  
  p<-lapply(names(L3),f)
  pdf(get.file("Figures/Fig3g.pdf"), width = 2, height = 4)
  print(p)
  dev.off()
}

