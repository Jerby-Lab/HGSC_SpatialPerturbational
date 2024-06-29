#### Results Section 4 ###
# Figure 4. MTIL as prognostic marker and predictor of Immune Checkpoint Blockade (ICB). 

# Figure 4a. Forest Plot of Overall Survival in HGSC Spatial Cohort
# Figure 4b. Kaplan Meier Curves of Overall Survival in HGSC Spatial Cohort
# Figure 4c. MTIL as predictor of time to event clinical data, public data 
# Figure 4d: MTIL as predictor of response, I-SPY2 trial 
# Figure 4e: T/NK and TMB predictors of ICB

#' Figure 4 Wrapper Function
#'
#' This function calls code to reproduce main text Figures 4a-e
#'
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure4_ICB<-function(){
  #1 Regenerate Figure 4a: Forest Plot of Overall Survival in HGSC Spatial Cohort
  ICB_Fig4a()
  #2 Regenerate Figure 4b: Kaplan Meier Curves of Overall Survival in HGSC Spatial Cohort
  ICB_Fig4b()
  # #3 Regenerate Figure 4c: MTIL as predictor of time to event clinical data, public data 
  ICB_Fig4c()
  # #4 Regenerate Figure 4d: MTIL as predictor of response, I-SPY2 trial
  ICB_Fig4d()
  # #5 Regenerate Figure 4e: T/NK and TMB predictors of ICB
  ICB_Fig4e()
  return()
}

#' Figure 4a. Forest Plot of Overall Survival in HGSC Spatial Cohort
#'
#' This function takes fitted multivariate cox proportional hazards models and
#' visualizes hazards ratios in a forest plot. 
#'
#' @return null, but writes figures in .pdf format 
#' in the Figures/ folder. 
ICB_Fig4a <- function(){
  # sub-function for finding the confidence interval of hazards ratio
  get_CI <- function(cox_model){
    summary_cox <- summary(cox_model)
    # extract coefficients and standard errors
    coefficients <- summary_cox$coefficients
    se <- coefficients[, "se(coef)"]
    # compute 90% CIs
    z <- qnorm(0.95)  # 1.645 for 90% CI
    hr <- exp(coefficients[, "coef"])
    lower_90 <- exp(coefficients[, "coef"] - z * se)
    upper_90 <- exp(coefficients[, "coef"] + z * se)
    # combine into a data frame
    ci_data <- data.frame(
      Variable = rownames(coefficients),
      HR = hr,
      Lower_90 = lower_90,
      Upper_90 = upper_90
    )
    return(ci_data)
  }
  
  # construct the input data
  CI <- get_CI(readRDS(get.file("Results/HGSC_MultivariateCox1.rds")))
  CI2 <- get_CI(readRDS(get.file("Results/HGSC_MultivariateCox2.rds")))
  forest_data <- rbind(rbind(CI[1:6, ], CI2[c("BRCAP", "TMB"),]), CI[7,])
  
  # set capping thresholds
  max_hr_threshold = 100  # for example, cap at 15
  min_hr_threshold = 0.0001 # cap lower limit to avoid zero or negative
  
  # apply capping
  forest_data$HR[forest_data$HR > max_hr_threshold] <- max_hr_threshold
  forest_data$Lower[forest_data$Lower < min_hr_threshold] <- min_hr_threshold
  forest_data$Lower[forest_data$Lower > max_hr_threshold] <- max_hr_threshold
  forest_data$Upper[forest_data$Upper > max_hr_threshold] <- max_hr_threshold
  
  # format the base data 
  base_data <- tibble::tibble(
    mean  = forest_data$HR, 
    lower = forest_data$Lower,
    upper =  forest_data$Upper,
    variable = c("Age at Diagnosis (> 65 vs. < 65)",
                 "Disease Stage at Diagnosis (III vs. IV)", 
                 "Neoadjuvant Chemotherapy (+/-)",
                 "Immune Checkpoint Blockade (+/-)", 
                 "Maintenance Bevacizumab (+/-)",
                 "PARP Inhibitor Use (+/-)", 
                 "BRCA1/2 Mutations (P or NP/WT)", 
                 "Tumor Mutational Burden (m/kB)", 
                 "Average MTIL Expression"))
  # make the plot 
  p = base_data |>
    forestplot(labeltext = c(variable),
               clip = c(0.01, 10),
               xlog = TRUE, 
               title = "Overall Survival, HGSC") |>
    fp_set_style(box = "royalblue",
                 line = "darkblue",
                 summary = "royalblue") |> 
    fp_add_header(variable = c("", "Variable")) |> 
    fp_set_zebra_style("#EFEFEF")
  
  # print to disk 
  pdf(get.file("Figures/Fig4a.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
}

#' Figure 4b. Kaplan Meier Curves of Overall Survival in HGSC Spatial Cohort
#'
#' This function discretizes MTIL into bottom quartile, top quartile of patients
#' based on average MTIL expression and estimated T/NK cell levels. This function
#' then computes the log-rank p-value for predictability of survival and 
#' visualizes the results on a kaplan meier curve. 
#'
#' @return null, but writes figures in .pdf format 
#' in the Figures/ folder. 
ICB_Fig4b <- function() {
  # Load in the data. 
  surv <- readRDS(get.file("Data/HGSC_SurvivalData.rds"))
  surv$fu_time3 <- surv$fu_time1/365
  
  # Specifying the features to analyze
  fu_features <- c("MTIL", "TNK")
  leg.name = c("MTIL Expression", "T/NK Cell Density")
  names(leg.name) = fu_features
  title.name = c("Overall Survival, by MTIL Expression",
                 "Overall Survival, by T/NK Cell Density")
  names(title.name) = fu_features
  
  # Opening a PDF device to save plots
  pdf(get.file("Figures/Fig4b.pdf"), width = 7, height = 7)
  
  for (feat in fu_features) {
    # Identifying necessary columns
    relevant_columns <- surv[, c(which(colnames(surv) %in% c("patients",
                                                             "fu_time3", 
                                                             "event")), 
                                 which(colnames(surv) == feat))]
    colnames(relevant_columns)[4] <- "feat"
    
    # Calculating the 25th, 75th, and additional percentiles for the feature
    quantiles <- quantile(relevant_columns[,"feat"], probs = c(0.25, 0.75), na.rm = TRUE)
    low_quantile <- quantiles[1]
    high_quantile <- quantiles[2]
    
    # Categorizing data into quartiles including the moderate levels (26-74%)
    relevant_columns$categorized_feat <- ifelse(relevant_columns$feat > high_quantile, "Top 25%",
                                                ifelse(relevant_columns$feat < low_quantile, "Bottom 25%", "Middle 50%"))
    
    # Now all categories, including "Middle 50%", will be analyzed
    fit <- survfit(Surv(fu_time3, event) ~ categorized_feat, data = relevant_columns, conf.int = 0.9)
    p <- survdiff(Surv(fu_time3, event) ~ categorized_feat, data = relevant_columns)$pvalue
    
    count = table(relevant_columns$categorized_feat)
    custom_labels <- c("Bottom 25%" = paste0("Low (n = ", count[1], ")"),
                       "Middle 50%" = paste0("Moderate (n = ", count[2], ")"), 
                       "Top 25%" = paste0("High (n = ", count[3], ")") )
    
    # Plotting with ggsurvplot for flexible plotting and confidence intervals
    p1 <- ggsurvplot(fit, data = relevant_columns,
                     pval = T, # Automatically add p-value
                     risk.table = T, # Add a risk table
                     pval.method = TRUE, # Show method of p-value calculation
                     palette = c("#3A53A4", "#196533", "#ED2224"),
                     conf.int = TRUE, #Show 95% confidence intervall
                     conf.int.style = "ribbon", # Display confidence interval as a ribbon
                     conf.int.color = "lightgrey", # Set color of the confidence interval# 
                     xlab = "Follow-up Time (Days)",
                     ylab = "Survival Probability",
                     title = paste0(title.name[feat]), 
                     subtitle = paste0("(p = use Cox Reg. from Table 8"),
                     risk.table.y.text.col = F, # Color text in risk table
                     risk.table.y.text = F, censor = F, 
                     legend = "right", # Position legend on the left
                     legend.labs = custom_labels, 
                     legend.title = leg.name[feat]) # Custom labels for the legend)
    print(p1)
  }
  dev.off()
}

#' Figure 4c. MTIL as predictor of time to event clinical data, public data 
#'
#' This function plots kaplan meier curves of MTIL estimated from bulk 
#' transcriptomics data as a predictor of survival in publicly available 
#' Melanoma, Non-small cell lung cancer, and HGSC datasetts. 
#'
#' @return null, but writes figures in .pdf format 
#' in the Figures/ folder. 
ICB_Fig4c <- function() {
  # load bulk transcriptomics data 
  mel = readRDS(get.file("Data/ICB_melanoma_Liu2019.rds"))
  lc = readRDS(get.file("Data/ICB_NSCLC_Ravi2023.rds"))
  lc$survival <- Surv(lc$metadata$Harmonized_PFS_Days/365,
                      lc$metadata$Harmonized_OS_Event)
  hgsc = readRDS(get.file("Data/ICGC_HGSC_AU.rds"))
  hgsc<-set.list(hgsc,hgsc$metadata$donor_tumour_stage_at_diagnosis=="III")
  
  # quantize mtil expression 
  q9m <-quantile(mel$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(mel$survival[,1]))))
  q9l <-quantile(lc$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(lc$survival[,1]))))
  q9h <-quantile(hgsc$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(hgsc$survival[,1]))))
  
  # plot kaplan meier curves to disk
  pdf(get.file("Figures/Fig4c.pdf"), width = 8, height = 2.5)
  par(mfrow=c(1,3),oma = c(0, 1, 0, 1),xpd = T)
  out = km.plot3(mel,mel$scores[,"mTIL"],qua = 0.2,xlim = q9m,direction = -1,main = "Melanoma mTIL")
  out = km.plot3(lc,lc$scores[,"mTIL"],qua = 0.2,xlim = q9l,direction = -1,main = "NSCLC mTIL")
  out = km.plot3(hgsc,hgsc$scores[,"mTIL"],qua = 0.2,xlim = q9h,direction = -1,main = "HGSC mTIL")
  dev.off()
}

# Figure 4d: MTIL as predictor of response, I-SPY2 trial 

#' This function plots boxplots of MTIL expression (from bulk transcriptomics)
#' as a function of binary response (pCR vs. no pCR) in two arms of the I-SPY2 
#' clinical trials 
#'
#' @return null, but writes figures in .pdf format 
#' in the Figures/ folder. 
ICB_Fig4d <- function() {
  # load Pusztai et al dataset
  r <- readRDS("Data/Breast_Durv_Pusztai.rds")
  
  # make plotting data frame for the Pusztai et al dataset
  plt <- data.frame(r[c("pCR", "arm", "HR", "HER2")], r$MTIL, 
                    TNK.cell = r$cell.sig[,"TNK.cell"]) 
  plt <- filter(plt, arm != "control") # HER2- 71 patients 
  plt <- plt %>% 
    mutate(pCR = factor(c("No (n=42)", "Yes (n=29)")[(pCR + 1)]))
  
  # make boxplot for Pusztai et al dataset
  a = ggboxplot(plt, x = "pCR", y = "mTIL") + 
    geom_signif(comparisons = list(c("No (n=42)", "Yes (n=29)")), 
                map_signif_level=TRUE, test.args = list(alternative = "less")) + 
    ggtitle("I-SPY2 Trial\nDurvalumab+Olaparib") + 
    ylab("MTIL Overall Expression")
  
  # load Wolf et al dataset
  r <-readRDS("Data/Breast_Pembro_Wolf.rds")
  
  # make plotting data frame for the Wolf et al dataset
  plt <- data.frame(r[c("pCR", "arm", "HR", "HER2")], r$MTIL,
                    TNK.cell = r$cell.sigs[,"TNK.cell"]) 
  plt <- filter(plt, arm == "Pembro") # HER2- 71 patients 
  plt <- plt %>% 
    mutate(pCR = factor(c("No (n=38)", "Yes (n=31)")[(pCR + 1)]))
  
  # make boxplot for Wolf et al dataset 
  b = ggboxplot(plt, x = "pCR", y = "mTIL") + 
    geom_signif(comparisons = list(c("No (n=38)", "Yes (n=31)")), 
                map_signif_level=TRUE, test.args = list(alternative = "less")) + 
    ggtitle("I-SPY2 Trial\nPembrolizumab+Paclitaxel") + 
    ylab("MTIL Overall Expression")
  
  # plot to disk 
  pdf(get.file("Figures/Fig4d.pdf"), width = 7, height = 8)
  print(a + b)
  dev.off()
}

# Figure 4e: T/NK and TMB predictors of ICB

#' This function plots kaplan meier curves of other biomarkers (estimated T/NK 
#' cell levels from bulk transcriptomics data, and tumor mutational burden) to 
#' test their predictability of survival in publicly available 
#' Melanoma, Non-small cell lung cancer, and HGSC datasets. 
#'
#'
#' @return null, but writes figures in .pdf format 
#' in the Figures/ folder. 
ICB_Fig4e <- function(){
  # load bulk transcriptomics data 
  mel = readRDS(get.file("Data/ICB_melanoma_Liu2019.rds"))
  lc = readRDS(get.file("Data/ICB_NSCLC_Ravi2023.rds"))
  lc$survival <- Surv(lc$metadata$Harmonized_PFS_Days/365,
                      lc$metadata$Harmonized_OS_Event)
  hgsc = readRDS(get.file("Data/ICGC_HGSC_AU.rds"))
  hgsc<-set.list(hgsc,hgsc$metadata$donor_tumour_stage_at_diagnosis=="III")
  
  # quantize mtil expression 
  q9m <-quantile(mel$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(mel$survival[,1]))))
  q9l <-quantile(lc$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(lc$survival[,1]))))
  q9h <-quantile(hgsc$survival[,1],
                 na.rm = T,
                 probs = 1-(10/sum(!is.na(hgsc$survival[,1]))))
  
  # plot kaplan meier curves to disk
  pdf(get.file("Figures/Fig4e.pdf"), width = 10.5, height = 2.5)
  par(mfrow=c(1,4),oma = c(0, 1, 0, 1),xpd = T)
  out = km.plot3(mel,mel$tme[,"T.cell"],
                 qua = 0.2,xlim = q9m,direction = -1,main = "Melanoma TNK")
  out = km.plot3(mel,mel$conf[,"TMB"],
                 qua = 0.2,xlim = q9m,direction = -1,main = "Melanoma TMB")
  out = km.plot3(lc,lc$tme[,"T.cell"],
                 qua = 0.2,xlim = q9l,direction = -1,main = "NSCLC TNK")
  out = km.plot3(hgsc,hgsc$tme[,"T.cell"],
                 qua = 0.2,xlim = q9h,direction = -1,main = "HGSC TNK")
  dev.off()
}
