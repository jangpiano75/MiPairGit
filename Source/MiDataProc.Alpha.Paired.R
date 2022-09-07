library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
#library(quantreg)
library(fossil)
library(picante)
library(entropart)
library(lme4)
library(lmerTest)
#library(broom.mixed)
#library(gee)
#library(geepack)
library(ICSNP)

#####################
# Data Manipulation #
#####################

rem.tax.d <- c("", "metagenome", "gut metagenome", "mouse gut metagenome")
rem.tax.str.d <- c("uncultured", "incertae", "Incertae", "unidentified", "unclassified", "unknown")

tax.tab.clean <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- na.code
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  ind <- which(tax.tab.h[,2] != na.code)
  ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  ind <- which(tax.tab.h[,3] != na.code)
  ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  ind <- which(tax.tab.h[,4] != na.code)
  ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  ind <- which(tax.tab.h[,5] != na.code)
  ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  ind <- which(tax.tab.h[,6] != na.code)
  ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  ind <- which(tax.tab.h[,7] != na.code)
  ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

otu.tab.clean <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}

biom.clean <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  if (tax.tab.c) {
    tax.tab <- tax.tab.clean(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

num.tax.rank <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.cleaned <- tax.tab.clean(tax.tab, rem.tax, rem.tax.str, na.code = na.code)
  num.taxa <- c()
  for (i in 1:6) {
    taxa <- unique(tax.tab.cleaned[,i+1])
    uni.taxa <- sum(taxa == na.code)
    num.taxa[i] <- nrow(taxa) - uni.taxa
  }
  return(num.taxa)
}

lib.size.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  lib.size.sum <- c(mean(lib.size), quantile(lib.size))
  names(lib.size.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(lib.size = lib.size, lib.size.sum = lib.size.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

mean.prop.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  mean.prop.sum <- c(mean(mean.prop), quantile(mean.prop))
  names(mean.prop.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(mean.prop = mean.prop, mean.prop.sum = mean.prop.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

rarefy.func <- function(biom, cut.off, multi.rarefy = FALSE) {
  
  if (!multi.rarefy | multi.rarefy == 1) {
    biom <- rarefy_even_depth(biom, cut.off, rngseed = 487)
  } else {
    otu.tab <- otu_table(biom)
    tax.tab <- tax_table(biom)
    tree <- phy_tree(biom)
    sam.dat <- sample_data(biom)
    
    otu.tab.list <- list()
    for (i in 1:multi.rarefy) {
      otu.tab.list[[i]] <- otu_table(rarefy_even_depth(biom, cut.off, rngseed = i), taxa_are_rows = TRUE)
    }
    
    sum.otu.tab <- otu.tab.list[[1]]
    for (i in 2:multi.rarefy) {
      sum.otu.tab <- sum.otu.tab + otu.tab.list[[i]]
    }
    otu.tab <- otu_table(round(sum.otu.tab/multi.rarefy), taxa_are_rows = TRUE)
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat) 
  }
  
  return(biom)
}

alpha.pe.pqe.func <- function(x, tree, norm = TRUE) {
  ind <- which(x != 0)
  s.tree <- prune_taxa(names(x[ind]), tree)
  pe <- AllenH(x[ind], 1, s.tree, Normalize = norm, CheckArguments = FALSE)
  pqe <- AllenH(x[ind], 2, s.tree, Normalize = norm, CheckArguments = FALSE)
  return(c(pe, pqe))
}

alpha.v1.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice, alpha.pd))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  
  return(alpha.div)
}

cov.func <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind.pri <- colnames(sam.dat) == sel.pri.var
  ind.mon.sin.rev <- mon.sin.rev.bin.con$is.mon | mon.sin.rev.bin.con$is.rev | mon.sin.rev.bin.con$is.sin
  return(colnames(sam.dat)[!(ind.pri | ind.mon.sin.rev)])
}

alpha.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  
  return(bin.cat)
}

alpha.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

alpha.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

alpha.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, alpha.div) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  return(list(bin.var = bin.var, alpha.div = alpha.div))
}

alpha.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x, na.rm = TRUE), quantile(x, na.rm = TRUE))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

alpha.bin.sum.func <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  ref.sum <- matrix(NA, n.alpha, 7)
  com.sum <- matrix(NA, n.alpha, 7)
  
  for (i in 1:n.alpha) {
    ind.alpha <- alpha.div[,i]
    sum.out <- tapply(ind.alpha, bin.var, alpha.ind.sum.func)
    ref.sum[i,] <- sum.out[[1]]
    com.sum[i,] <- sum.out[[2]]
  }
  rownames(ref.sum) <- colnames(alpha.div)
  colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  rownames(com.sum) <- colnames(alpha.div)
  colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  
  for(i in 1:nrow(ref.sum)){
    ref.sum[i,] <- decimal_adjust(ref.sum[i,])
  }
  
  for(i in 1:nrow(com.sum)){
    ref.sum[i,] <- decimal_adjust(com.sum[i,])
  }
  
  out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
  names(out) <- levels(bin.var)
  return(out)
}

prim_vars_pair <- function(sam_dat){
  name_vec <- c() 
  for (name in names(sam_dat)){
    if (length(table(sam_dat[,name])) != 1){
      name_vec <- c(name_vec, name)
    }
  }
  return(name_vec)
}

data_mani_paired <- function(alpha_div, sam_dat, prim_id, block_id){
  dat = cbind(sam_dat[,block_id], sam_dat[,prim_id], alpha_div)
  dat = dat[order(dat[,block_id], dat[,prim_id]),]
  colnames(dat)[1] <- "block_var"
  colnames(dat)[2] <- "prime_var"
  return(dat)
}

alpha.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

is.mon.sin.rev.combined <- function(sam.dat){
  is.symmetric <- logical()
  nvar <- ncol(sam.dat)
  
  for (i in 1:nvar){
    table <- table(sam.dat[,i])
    if (length(unique(table)) == 1 & length(table) > 1){
      is.symmetric[i] = TRUE
    }
    else{
      is.symmetric[i] = FALSE
    }
  }
  return(list(is.symmetric = is.symmetric))
}


####################################################
# Comparative analysis for paired alpha diversity  #
####################################################

alpha.t.paired <- function(alpha_div, alpha_mani){

  if (length(unique(as.numeric(table(alpha_mani$prime_var)))) != 1){
    ind <- as.numeric(table(alpha_mani$block_var)) != length(table(alpha_mani$prime_var)) 
    alpha_mani <- alpha_mani[!(alpha_mani$block_var %in% names(table(alpha_mani$block_var))[ind]),]   
  }

  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  ref.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[1])
  com.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[2])
  
  out <- matrix(NA, ncol(alpha_div), 6)
  
  for (i in 1:ncol(alpha_div)){
    fit <- t.test(alpha_mani[,i+2][ref.ind], alpha_mani[,i+2][com.ind], paired = TRUE)
    out[i,] <- c(round(fit$statistic, 3), round(fit$stderr, 3), round(fit$parameter, 3), round(fit$conf.int,3) , round(fit$p.value, 3))
  }

  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha_div)
  colnames(out) <- c("t", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
}

alpha.wilcox.paired <- function(alpha_div, alpha_mani){
  
  if (length(unique(as.numeric(table(alpha_mani$prime_var)))) != 1){
    ind <- as.numeric(table(alpha_mani$block_var)) != length(table(alpha_mani$prime_var)) 
    alpha_mani <- alpha_mani[!(alpha_mani$block_var %in% names(table(alpha_mani$block_var))[ind]),]
  }
  
  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  ref.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[1])
  com.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[2])
  
  out <- matrix(NA, ncol(alpha_div), 2)
  
  for (i in 1:ncol(alpha_div)){
    fit <- wilcox.test(alpha_mani[,i+2][ref.ind], alpha_mani[,i+2][com.ind], paired = TRUE, correct = FALSE)
    out[i,] <- c(round(fit$statistic,3), round(fit$p.value, 3))
  }
  
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha_div)
  colnames(out) <- c("Wilcoxon", "P.value")
  return(out)
}

hotel.data.standard <- function (alpha_mani){
  
  if (length(unique(as.numeric(table(alpha_mani$prime_var)))) != 1){
    ind <- as.numeric(table(alpha_mani$block_var)) != length(table(alpha_mani$prime_var))     
    alpha_mani <- alpha_mani[!(alpha_mani$block_var %in% names(table(alpha_mani$block_var))[ind]),] 
  }
  
  alpha_mani_standard <- alpha_mani
  
  for (i in 3:ncol(alpha_mani)){
    alpha_mani_standard[,i] <- (alpha_mani[,i] - mean(alpha_mani[,i]))/sd(alpha_mani[,i])
  }
  
  return (alpha_mani_standard)
}


alpha.hotel.paired <- function(alpha_mani){
  time.p.cat <- names(table(alpha_mani[, "prime_var"]))
  
  ref.ind <- which(alpha_mani[, "prime_var"] == time.p.cat[1])
  com.ind <- which(alpha_mani[, "prime_var"] == time.p.cat[2])
  
  mat_1 <- alpha_mani[ref.ind, 3:ncol(alpha_mani)]
  mat_2 <- alpha_mani[com.ind, 3:ncol(alpha_mani)]
  mat <- mat_1 - mat_2 
  
  fit = HotellingsT2(mat)
  
  out <- c(round(fit$statistic,3), round(fit$p.value, 3) , fit$parameter) 
  names(out) <- c("T.2", "P.value", "DF1", "DF2")
  return(out)
}

hotel.data.mani.for.CI <- function(alpha_mani, par = "ref"){  #parameter can be "ref" or "com" 
  time.p.cat <- names(table(alpha_mani[, "prime_var"]))
  
  ref.ind <- which(alpha_mani[, "prime_var"] == time.p.cat[1])
  com.ind <- which(alpha_mani[, "prime_var"] == time.p.cat[2])
  
  mat_1 <- alpha_mani[ref.ind, 3:ncol(alpha_mani)]
  mat_2 <- alpha_mani[com.ind, 3:ncol(alpha_mani)]
  
  if (par == "ref"){
    return (mat_1)
  } else{
    return (mat_2)
  }
}

simul_CI <- function (dat_ref, dat_com, level = 0.95, type = "Simul") {
  p <- ncol(dat_ref)
  n <- nrow(dat_com)
  d <- dat_ref - dat_com   
  
  dbar <- apply(d, 2, mean)  
  
  ref_mean <- apply(dat_ref, 2, mean)
  com_mean <- apply(dat_com, 2, mean)
  
  scit <- matrix(rep(0, p * 5), nrow = p)
  csq <- (((n - 1) * p)/(n - p)) * qf(level, p, n - p)
  
  cov <- (1/(n-1)) * as.matrix((d - dbar)) %*% as.matrix(t(d - dbar))
  
  for (i in 1:p){
    scit[i, 1] <- dbar[i]
    scit[i, 2] <- dbar[i] - sqrt(((cov[i, i])^2/n) * csq)
    scit[i, 3] <- dbar[i] + sqrt(((cov[i, i])^2/n) * csq)
    
  }
  
  sum_dat_ci <- data.frame(Estimate = scit[, 1], Case_mean = ref_mean, Control_mean = com_mean, LowerCI = scit[, 2], UpperCI = scit[, 3])
  rownames(sum_dat_ci) <- colnames(dat_ref)
  
  return(sum_dat_ci)
}

global.alpha.lmm.p <- function(alpha_mani){
  alpha_relevel_dat <- alpha_mani 
  alpha_relevel_dat[,"prime_var"] <- as.factor(alpha_relevel_dat[,"prime_var"])
  alpha_relevel_dat <<- alpha_relevel_dat 
  
  p_val <- c() 
  for (i in 3:length(alpha_relevel_dat)){
    fit_full <- lmer(alpha_relevel_dat[,i] ~ alpha_relevel_dat[,2]+ (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    fit_nested <- lmer(alpha_relevel_dat[,i] ~ 1 + (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    p_val <- c(p_val, as.numeric(unlist(anova(fit_full, fit_nested))["Pr(>Chisq)2"]))
  }
  
  return(p_val)
}

paired.alpha.lmm.dat <- function(alpha_mani){
  alpha_relevel_dat <- alpha_mani 
  alpha_relevel_dat[,"prime_var"] <- as.factor(alpha_relevel_dat[,"prime_var"])
  alpha_relevel_dat <<- alpha_relevel_dat 
  
  mat <- matrix(NA, nrow = 9, ncol = 5)
  
  for (i in 3:length(alpha_relevel_dat)){
    fit <- lmer(alpha_relevel_dat[,i] ~ alpha_relevel_dat[,2]+ (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    mat[i-2,] <- coef(summary(fit))[2,] 
  }
  rownames(mat) <- names(alpha_mani)[3:11]
  mat <- apply(mat, 2, decimal_adjust)
  mat <- data.frame(mat)
  colnames(mat) <- c("Estimate", "Std Err", "DF", "t", "P.value")
  
  return(mat)
}

###########################
# Visualization (Boxplot) #
###########################

alpha.bin.paired <- function(alpha_div, alpha_mani, test_out){
  
  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  
  ref.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[1])
  com.ind <- which(alpha_mani[,"prime_var"] == time.p.cat[2])
  
  alpha_ref <- alpha_mani[ref.ind, 3:ncol(alpha_mani)]
  alpha_com <- alpha_mani[com.ind, 3:ncol(alpha_mani)] 
  
  n.alpha <- ncol(alpha_div)
  alpha.div.ind <- colnames(alpha_div)
  pvs <- as.numeric(test_out[,"P.value"]) 
  ind.p.sig <- which(pvs < 0.05)
  
  par(mfrow = c(3, 3))
  for (i in 1:n.alpha){
    if (is.element(i, ind.p.sig)){
      xlab.v <- paste("*p:", p.value.0.1(pvs[i]), sep = "")
    } else {
      xlab.v <- paste("p:", p.value.0.1(pvs[i]), sep = "")
    }
    boxplot(alpha_ref[,i], alpha_com[,i], at = c(1, 2), xlab = xlab.v, ylab = alpha.div.ind[i], names = c(time.p.cat[1], time.p.cat[2]), col=c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)), notch = TRUE, horizontal = FALSE) 
  }
}

Hotell_CI_plot_3 <- function(hotel_test_result, CI, type = "Simul", level = 0.95){
  plot_mat <- as.matrix(rbind(c("Alpha Diversity", "Mean Ref", "Mean Conf" ,"Mean Diff"), cbind(c(rownames(CI), "P-value"), c(round(CI$Case_mean, 3), NA), c(round(CI$Control_mean, 3), NA), c(round(CI$Estimate, 3), round(hotel_test_result[["P.value"]], 3)))))
  
  
  plot_mat_list <- list(list(), list(), list(), list())
  plot_mat_list[[1]] <- plot_mat[, 1]
  plot_mat_list[[2]] <- plot_mat[, 2]
  plot_mat_list[[3]] <- plot_mat[, 3]
  plot_mat_list[[4]] <- plot_mat[, 4]
  
  Plot <- forestplot(labeltext = plot_mat_list, mean = c(NA, CI$Estimate, NA), lower = c(NA, CI$LowerCI, NA), upper = c(NA, CI$UpperCI, NA), 
                     hrzl_lines = list("2" = gpar(lwd = 2.5, col = "#000044"), "11" = gpar(lwd = 2.5, col = "#000044"), "12" = gpar(lwd = 2.5, col = "#000044")), new_page = TRUE, boxsize = 0.25, line.margin = 0.3, grid = 0, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),is.summary = c(rep(TRUE, 1), rep(FALSE, 9), TRUE),
                     col = fpColors(box = rgb(1, 0, 0, 0.5), line = "black", summary = "red3"), xlab = paste(paste0(level*100, "%"), "Simultaneous Confidence Interval"),
                     txt_gp = fpTxtGp(label = list(gpar(fontfamilly = "", cex = 0.7), gpar(fontfamily = "", cex = 0.7)), ticks = gpar(fontfamily = "", cex = 0.7), xlab = gpar(fontfamily= "", cex = 0.7))) 
  
  return(Plot)
}
