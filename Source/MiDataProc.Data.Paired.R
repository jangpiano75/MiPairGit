REFERENCE_CHECK_P <- function(data_transform = "", method_name = "", FDR = ""){
  reference_lists <- NULL
  if(data_transform == "CLR"){
    reference_lists = c(reference_lists, "Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-60.")
  }else if(data_transform == "Proportion"){
    reference_lists = c(reference_lists, "Martin-Fernandez J-A, Hron K, Templ M, Filzmoser P, Palarea-Albaladejo J. Bayesian-multiplicative treatment of count zeros in compositional data sets. Stat Modelling. 2015;15(2):134-58.")
  }else if(data_transform == "Count (Rarefied)"){
    "Sanders HL. Marine benthic diversity: a comparative study. Am Nat. 1968;102(925):243-82."
  }
  
  if(method_name == "Wilcoxon signed-rank test" |method_name ==  "Wilcoxon signed-rank test (default)"){
    reference_lists = c(reference_lists, "Wilcoxon F. Individual comparisons by ranking methods. Biometrics Bulletin. 1945;1(6):80.")}
  else if(method_name == "Multivariate Hotelling's t-squared test"){
    reference_lists = c(reference_lists, 
                        "Hotelling H. The generalization of Student's ratio. The Annals of Mathematical Statistics. 1931;2(3):360-378.")}
  else if(method_name == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
    reference_lists = c(reference_lists, "Tukey, John. Comparing Individual Means in the Analysis of Variance. Biometrics. 1949;5(2):99-114.", "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  else if(method_name == "Friedman's test (global) with Conover's test (pairwise)" | method_name == "Friedman's test (global) with Conover's test (pairwise)(default)"){
    reference_lists = c(reference_lists, "Friedman M. A comparison of alternative tests of significance for the problem of m rankings. The Annals of Mathematical Statistics. 1940;11(1):86-92.", "Conover WJ. (1999). Pratical Nonparametric Statistics. 3rd ed. New York: John Wiley & SonsAnderson MJ. A new method for non-parametric multivariate analysis of variance. Austral Ecol. 2001;26(1):32-46.",
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  
  else if(method_name == "Durbin's test (global) with Conover's test (pairwise)(default)" | method_name == "Durbin's test (global) with Conover's test (pairwise)"){
    reference_lists = c(reference_lists, 
                        "Conover WJ. (1999). Pratical Nonparametric Statistics. 3rd ed. New York: John Wiley & Sons Anderson MJ. A new method for non-parametric multivariate analysis of variance. Austral Ecol. 2001;26(1):32-46.",
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  
  else if(method_name == "LMM (LRT (global) with t-test (pairwise))"){
    reference_lists = c(reference_lists, "Laird NM, Ware JH. Random-effects models for longitudinal data. Biometrics. 1982;38(4):963-74.", 
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  else if(method_name == "PERMANOVA"){
    reference_lists = c(reference_lists, "Anderson MA. new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.", 
                        "McArdle BH, Anderson MJ. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.")
  }else if(method_name == "PERMANOVA_BASE"){
    reference_lists = c(reference_lists, "Anderson MA. new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.", 
                        "McArdle BH, Anderson MJ. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.",
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }else if(method_name == "PERMANOVA_ACROSS"){
    reference_lists = c(reference_lists, "Anderson MA. new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.", 
                        "McArdle BH, Anderson MJ. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.",
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  
  else if(method_name == "PERMANOVA-FL"){
    reference_lists = c(reference_lists, "Zhu Z, Satten GA, Mitchell C, Hu Y-J. Constraining permanova and LDM to within-set comparisons by projection improves the efficiency of analyses of matched sets of Microbiome Data. Microbiome, 2021:9(1).")
  }else if(method_name == "LDM"){
    reference_lists = c(reference_lists, "Zhu Z, Satten GA, Mitchell C, Hu Y-J. Constraining permanova and LDM to within-set comparisons by projection improves the efficiency of analyses of matched sets of Microbiome Data. Microbiome. 2021:9(1).",
                        "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }

  else if(method_name == "GLMM (Binomial)" | method_name == "GLMM (Negative Binomial)" | method_name == "GLMM (Beta)"){
    reference_lists = c(reference_lists, "Breslow NE, Clayton DG. Approximate inference in generalized linear mixed models. J Am Stat Assoc. 1993;88(421):9-25.")
  }else if(method_name == "GEE (Binomial)"){
    reference_lists = c(reference_lists, "Liang K-Y, Zeger SL. Longitudinal data analysis using generalized linear models. Biometrika. 1986;73(1):13-22.")
  }
  
  if(FDR == "Yes"){
    reference_lists = c(reference_lists, "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  
  if(length(reference_lists) == 0){
    return(NULL)
  }
  else{
    for (i in seq(1:length(reference_lists))){
      reference_lists[i] = paste(i, ". ", reference_lists[i], sep = "")
    }
    return(reference_lists)
  }
}