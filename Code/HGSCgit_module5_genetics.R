#### Results Section 5 ###
# Figure 5. Copy Number Alterations (CNAs) mapping to mTIL and TIL levels in SMI spatial data and TCGA.
#
# Figure 5A: mTIL at baseline
# Figure 5B: Top mTIl <> CNA correlation genes (boxplots)
# Figure 5C: mTIL at baseline
# Figure 5D: mTIL survival
# Figure 5E: TCGA: mTIL vs. CNAs

HGSC_Figure5_CNAs<-function(r,r1,rslts1,rslts2){

  if(missing(r1)){
    r1<-readRDS(get.file("Data/SMI_mTIL_CNA.rds"))
    rslts1<-readRDS(get.file("Results/HGSC_CNAs.vs.mTIL_SMI.rds"))
    rslts2<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  }

  #1. Regenerate Figure 5A: mTIL at baseline
  # HGSC_Fig5A(r1 = r1)
  #2. Regenerate Figure 5B: Top mTIl <> CNA correlation genes (boxplots)
  HGSC_Fig5B(r1 = r1,rslts = rslts1)
  #3. Regenerate Figure 5C: mTIL at baseline
  HGSC_Fig5C(rslts = rslts1)
  #4. Regenerate Figure 5D: mTIL survival
  HGSC_Fig5D(r=r)
  #5. Regenerate Figure 5E: TCGA: mTIL vs. CNAs
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

HGSC_Fig5D <- function(r){
  # Fx Get mTIL Overall Expression
  get_mtil_expression <- function(r, mtil){
    r1 <- subset_list(r, subcells = r$cells[r$cell.types == "Malignant"])
    r1 <- prep4OE(r1)
    r1$mTIL <- get.OE(r1, mtil)
    mtil_df_mean <- data.frame(
      data.frame(r1[c("patients",
                      "treatment",
                      "samples",
                      "sites_binary")])[row.names(r1$mTIL),],
      r1$mTIL[,grepl("100", colnames(r1$mTIL))]) %>%
      group_by(patients, treatment, samples, sites_binary) %>%
      summarize_all(mean) %>%
      group_by(patients, treatment, sites_binary) %>%
      select(-samples) %>%
      summarize_all(mean) %>%
      rename(mean_hot100 = hot100)

    mtil_df_med <- data.frame(
      data.frame(r1[c("patients",
                      "treatment",
                      "samples",
                      "sites_binary")])[row.names(r1$mTIL),],
      r1$mTIL[,grepl("100", colnames(r1$mTIL))]) %>%
      group_by(patients, treatment, samples, sites_binary) %>%
      summarize_all(median) %>%
      group_by(patients, treatment, sites_binary) %>%
      select(-samples) %>%
      summarize_all(mean) %>%
      rename(median_hot100 = hot100)


    mtil <- Reduce(function(x, y) merge(x, y, by = c("patients",
                                                     "treatment",
                                                     "sites_binary"), all=TRUE),
                   list(mtil_df_mean, mtil_df_med))

    mtil <- mtil[,c(1:4,7)]
    mtil <- cast_sites(mtil, 4:5)
    return(mtil)
  }

  # Fx for Format Survival Data
  format_survival_data <- function(features){
    clin <- readRDS(get.file("Data/Clinical.rds"))
    surv <- select(clin, patients, specimen, fu_time1, fu_time2, outcome) %>% unique()
    surv$outcomes_binary <- unlist(lapply(surv$outcome, function(x){
      if (is.na(x)) {
        return(NA)
      } else if (x == "Alive") {
        return(x)
      } else {
        return ("Dead")
      }
    }))
    surv$event <- (surv$outcomes_binary == "Dead") + 0
    surv <- filter(surv, !is.na(event))
    surv <- merge(surv, features, by = "patients")
    surv$fu_time1 <- as.numeric(surv$fu_time1)
    return(surv)
  }

  # Fx for Univariate Model
  univariate_modeling <- function(features, surv){
    univariate <- do.call("rbind", lapply(colnames(features)[3:6], function(x){
      tmp <- surv
      colnames(tmp)[colnames(tmp) == x] <- "x"
      coxm <- coxph(formula = Surv(fu_time1, event) ~ x, data = tmp)
      out <- summary(coxm)$coefficients
      return(out)
    }))
    row.names(univariate) <- colnames(features)[3:6]
    colnames(univariate) <- c("coef", "exp", "se", "z", "p")
    univariate <- data.frame(univariate) %>% arrange(p)
    univariate[grepl("Adnexa", row.names(univariate)),]
    return(univariate)
  }

  # Fx for Kaplan Meier Curve
  km_plots <- function(univariate, surv){
    fu_features <- row.names(
      univariate[grepl("Adnexa", row.names(univariate)),]
    )[c(1:2)]

    pdf(get.file("Figures/Fig5D.pdf"),
        width = 5, height = 5)
    for (feat in fu_features) {
      tmp <- surv[,c(c(1:8), which(colnames(surv) == feat))]
      colnames(tmp)[9] <- "feat"

      cut <- surv_cutpoint(tmp,
                           time = "fu_time1",
                           event = "event",
                           "feat",
                           minprop = 0.1)

      tmpcat <- surv_categorize(cut)

      p <- (survdiff(Surv(fu_time1, event) ~ feat, tmpcat))$pvalue
      fit <- survfit(Surv(fu_time1, event) ~ feat, tmpcat)
      leg = paste0(paste0(c("high ", "low ")),
                   " (n = ", unlist(table(filter(data.frame(tmpcat),
                                                 !is.na(feat))$feat)), ")")
      plot(fit,
           main = paste0(feat, "\n(p = ", format(round(p, digits = 5),
                                                 scientific = T), ")"),
           ylab = "Survival Probability",
           xlab = "Follow-up Time (Days)",
           col = c("darkred", "#9ee7ff"),
           lwd = 2)
      legend("bottomleft", legend=leg,fill=c("darkred", "#9ee7ff"), bty="n")
    }
    dev.off()
  }

  # load mTIL
  mtil <- readRDS(get.file("Results/mTIL_sig.rds"))

  # get mTIL expression
  mtil_df <- get_mtil_expression(r, mtil)

  # format survival
  surv <- format_survival_data(mtil_df)

  # univariate modeling
  univariate <- univariate_modeling(mtil_df, surv)

  # make KM plots
  km_plots(univariate, surv)
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

