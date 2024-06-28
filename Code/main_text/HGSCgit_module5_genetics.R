#### Results Section 5 ###
# Figure 5. Copy Number Alterations (CNAs) mapping to mTIL and TIL levels 
# in spatial data and TCGA.
#
# Figure 5a: MTIL at baseline
# Figure 5b: Top MTIl <> CNA correlation genes (boxplots)
# Figure 5c: MTIL vs. CNAs Mixed Effects
# Figure 5d: Cartoon diagram generated outside of R (Illustrator)
# Figure 5e: TCGA: MTIL vs. CNAs
# Figure 5f: Experimental data (BMP-7 assay) figure generated in PRISM 
# Figure 5g: Figure generated in Cytoscape
# Figure 5h: Experimental data (Migration assay), figure generated in PRISM 

#' Figure 5 Wrapper Function
#'
#' This function calls code to reproduce main text Figures 5a, b, c, e. 
#' All other figures were produced outside of R. Source Data for 5f and 5h are
#' available with the published manuscript.  
#'
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure5_CNAs<-function(r,r1,rslts1,rslts2){
  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_mTIL_CNA.rds"))
    rslts1<-readRDS(get.file("Results/HGSC_CNAs.vs.mTIL_SMI.rds"))
    rslts2<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  }
  
  #1. Regenerate Figure 5a: mTIL at baseline
  HGSC_Fig5a()
  #2. Regenerate Figure 5b: Top mTIl <> CNA correlation genes (boxplots)
  HGSC_Fig5b(r1 = r1,rslts = rslts1)
  #3. Regenerate Figure 5c: mTIL at baseline
  HGSC_Fig5c(rslts = rslts1)
  #5. Regenerate Figure 5e: TCGA: mTIL vs. CNAs
  HGSC_Fig5e(rslts = rslts2)
}

#' Figure 5a. MTIL at baseline
#'
#' Here we plot the variation of MTIL expression at baseline (defined as MTIL
#' expression in cells that are in TIL-deprived niches) per patient. Data 
#' is from the Discovery and Test 2 datasets. 
#' 
#' @return this function returns nothing, but writes a figure in .pdf format 
#' in the Figures/ folder. 
HGSC_Fig5a <-function(){
  # load variation data. 
  df = readRDS(get.file("Results/HGSC_PatientVariation.rds"))
  
  # sort and label patients from different datasets 
  df$pat <- df$patients
  df$pat[grepl("HGSC", df$dataset)] <- paste0(
    df$patients[grepl("HGSC", 
                      df$dataset)],"_WT")
  df$pat <- gsub("HGSC", "", df$pat)
  order <- (df %>% group_by(pat) %>%
              summarize(mtil_m = median(hot100)) %>%
              arrange(mtil_m))$pat
  df$pat <- factor(df$pat, order)
  df$wt <- df$dataset
  df$wt[grepl("HGSC", df$dataset)] <- "Test 2 Dataset"
  pats =unique(filter(df, wt != "Discovery Dataset")$patients)
  df$wt[df$wt == "Discovery Dataset" & 
          df$patients %in% pats] <- "Discovery Dataset (Matched)"
  
  # plot
  p = ggplot(df, aes(x = pat, y = hot100, fill = wt)) +
    geom_boxplot(outlier.size = 0.2, size = 0.5) +  
    scale_fill_manual(values = c("white", 
                                 "#8ED2C6",
                                 "#C47AB3"), name = "Dataset") +
    theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          axis.line = element_line(size = 0.5)) +
    xlab("Tissue Profile") +
    ylab("MTIL Expression") + 
    ggtitle("MTIL expression per Tissue Profile, T/NK cell negative Niches")
  
  # print to disk. 
  pdf(get.file("Figures/Fig5a.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}

#' Figure 5b. Top MTIl <> CNA correlation genes (boxplots)
#'
#' Here we visualize specific examples of MTIL expression as a function of 
#' copy number alterations. 
#' 
#' @param r1 list object containing malignant cells from the Discovery Dataset
#' @param rslts list object containing scores from mixed effects modeling 
#' between gene expression and copy number alterations in the Discovery dataset
#' @return this function returns nothing, but writes a figure in .pdf format 
#' in the Figures/ folder. 
HGSC_Fig5b<-function(r1,rslts){
  # load data if missing 
  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_mTIL_CNA.rds"))
    rslts<-readRDS(get.file("Results/HGSC_mTIL_vsCNAs.rds"))
  }
  
  # function to make boxplots
  f1<-function(x){
    idx<-order(r1$cnv[,x])
    call.boxplot(r1$scores[idx,"mTIL.up"],
                 r1$cnv[idx,x],order.flag = F,
                 xlab = paste(x,"Copy Number"),ylab = "mTIL program (OE)",
                 main = paste(
                   x,"(",call.format.pval(rslts$HLM[x,"mTIL.up.P"]),")"))
  }
  
  # make boxplots for select genes
  l1<-lapply(c("TCF7L2","IFNGR2","AXL","IFNAR1","ACTA2","RUNX1"),f1)
  
  # plot to disk
  pdf(get.file("Figures/Fig5b.pdf"))
  call.multiplot(l1,nplots = 6,cols = 3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  return()
}

#' Figure 5c. MTIL vs. CNAs Mixed Effects
#'
#' Visualization of full list of MTIL genes with statistically significant 
#' association with CNA. 
#' 
#' @param rslts list object containing scores from mixed effects modeling 
#' between gene expression and copy number alterations in the Discovery dataset. 
#' @return this function returns nothing, but writes a figure in .pdf format 
#' in the Figures/ folder. 
HGSC_Fig5c<-function(rslts){
  # load data if missing
  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_mTIL_vsCNAs.rds"))
  }
  
  # adjust the Z-score for MHT 
  z<-rslts$Z[,"mTIL.up.Z"]
  names(z)<-rownames(rslts$Z)
  z<-z[p.adjust(10^-abs(z),method = "BH")<0.1]
  z<-sort(z)
  
  # plot to disk 
  pdf(get.file("Figures/Fig5c.pdf"))
  barplot(abs(z),las=2,cex.names = 0.5,xlab = "Genes CNA",
          col = ifelse(z<0,"lightblue","darkred"),
          ylab = "Z-score",main = "CNAs correlated with the mTIL program")
  legend(legend = c("Negative","Positive"), pch = 15,
         "topright",col = c("lightblue","darkred"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 5e. MTIL vs. CNAs Mixed Effects
#'
#' Visualization of full list of MTIL genes with statistically significant 
#' association with CNA. 
#' 
#' @param rslts list object containing scores from mixed effects modeling 
#' between gene expression and copy number alterations in the TCGA dataset. 
#' @return this function returns nothing, but writes a figure in .pdf format 
#' in the Figures/ folder. 
HGSC_Fig5e<-function(rslts){
  # load results if missiong
  if(missing(rslts)){
    rslts<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  }
  
  # split results into deletion and amplification
  Xd<-rslts$plot.del
  Xa<-rslts$plot.amp
  
  # make beanplots. 
  pdf(get.file("Figures/Fig5e.pdf"))
  violin.split(Xd$TIL,Xd$Del,conditions = Xd$Gene,
               ylab = "TIL levels",xlab = "mTIL UP Genes",
               main = "TCGA, mTIL CNA vs. TIL levels",cex.axis = 0.5)
  violin.split(Xa$TIL,Xa$Amp,conditions = Xa$Gene,
               ylab = "TIL levels",xlab = "mTIL DOWN Genes",
               main = "TCGA, mTIL CNA vs. TIL levels")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

