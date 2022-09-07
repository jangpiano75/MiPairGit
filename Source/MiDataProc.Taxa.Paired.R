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
library(gridGraphics)
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggthemes)
#library(erer)
library(stringr)
library(devtools)
#library(betareg)
library(lme4)
library(lmerTest)
#library(MiRKAT)
library(DiagrammeR)
library(reticulate)
#library(NBZIMM)
library(nlme)
#library(gee)
#library(geepack)
#library(gridExtra)

########
# Taxa #
########

tax.trans <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "NANANA") {
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
  
}

taxa.names.rank <- function(taxa.out){ 
  taxon.names <- list()
  taxon.names <- lapply(taxa.out, function(x) str_split(names(x), ";"))
  dup.list <- list(NA,NA,NA,NA,NA,NA)
  
  ranks <- c("K_", "P_", "C_", "O_", "F_", "G_", "S_")
  
  taxon.names.rank <- list()
  for(rank in 1:6){
    taxon <- lapply(taxon.names[[rank]], function(x) str_sub(x,start = 3))
    taxon.names.rank[[rank]] <- sapply(taxon, tail, 1)
    
    if(length(taxon.names.rank[[rank]]) != length(unique(taxon.names.rank[[rank]]))){
      duplicated.taxons <- unique(taxon.names.rank[[rank]][duplicated(taxon.names.rank[[rank]])])
      
      for(i in 1:length(duplicated.taxons)){
        duplicated.taxon <- duplicated.taxons[i]
        ind.dup <- which(taxon.names.rank[[rank]] %in% duplicated.taxon)
        
        for(j in 1:length(ind.dup)){
          duplicated.taxon <- paste(duplicated.taxon,"*",collapse = "")
          taxon.names.rank[[rank]][ind.dup[j]] <- duplicated.taxon 
          dup.list[[rank]][j] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[j]]), sep = ""), collapse = " | "), sep = "")
        }
      }
    }
  }
  names(taxon.names.rank) <- names(taxa.out)
  return(list(names = taxon.names.rank, duplicates = dup.list))
}

#####################
# Data manipulation #
#####################

taxa.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x, na.rm = TRUE), quantile(x, na.rm = TRUE))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

taxa.bin.sum.func <- function(bin.var, taxa) {
  
  if(is.null(ncol(taxa))){
    out <- NULL
    
  }else{
    n.taxa <- ncol(taxa)
    ref.sum <- matrix(NA, n.taxa, 7)
    com.sum <- matrix(NA, n.taxa, 7)
    for (i in 1:n.taxa) {
      ind.taxa <- taxa[,i]
      sum.out <- tapply(ind.taxa, bin.var, taxa.ind.sum.func)
      ref.sum[i,] <- sum.out[[1]]
      com.sum[i,] <- sum.out[[2]]
    }
    rownames(ref.sum) <- colnames(taxa)
    colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
    rownames(com.sum) <- colnames(taxa)
    colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
    out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
    names(out) <- levels(as.factor(bin.var))
  }
  return(out)
}

taxa.bin.sum.united.func <- function(bin.var, taxa.out) {
  taxa.bin.sum <-list()
  for(i in 1:6) {
    taxa.bin.sum[[i]] <- taxa.bin.sum.func(bin.var, taxa.out[[i]])
  }
  names(taxa.bin.sum) <- names(taxa.out)
  return(taxa.bin.sum)
}

taxa.bin.cat.ref.united.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  return(list(bin.var = bin.var, taxa = taxa.out))
}

bin.q.func <- function(out, method = c("BH", "BY")) {
  if(is.null(out)){
    return(NULL)
  }else{
    Q.value <- p.adjust(out$P.value, method = method)
    return(cbind(out, Q.value))
  }
}

data_mani_p_taxa <- function(taxa, sam_dat, prim_id, block_id){
  dat = cbind(sam_dat[,block_id], sam_dat[,prim_id], taxa)
  dat = dat[order(dat[,block_id], dat[,prim_id]),]
  colnames(dat)[1] <- "block_var"
  colnames(dat)[2] <- "prime_var"
  return(dat)
}

data_mani_p_sample <- function(sam_dat, prim_id, block_id){
  dat = cbind(sam_dat[,block_id], sam_dat[,prim_id])
  dat = dat[order(dat[,block_id], dat[,prim_id]),]
  colnames(dat)[1] <- "block_var"
  colnames(dat)[2] <- "prime_var"
  return(dat)
}

bin.q.united.func <- function(taxa.out, method = "BH") {
  q.out <- list()
  for(i in 1:6) {
    q.out[[i]] <- bin.q.func(taxa.out[[i]], method)
  }
  names(q.out) <- names(taxa.out)
  return(q.out)
}

taxa.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  if (length(bin.cat) != 2) {
    stop(paste(sel.bin.var, " is not binary", sep = ""))
  }
  return(bin.cat)
}

taxa.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

taxa.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

taxa.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  return(list(bin.var = bin.var, taxa = taxa))
}

taxa.bin.cat.ref.united.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, taxa = taxa.out))
}

###################################################
# Comparative analysis for paired taxa diversity  #
###################################################

taxa.t.paired <- function(taxa, data_mani_taxa){
  
  if (length(unique(as.numeric(table(data_mani_taxa$prime_var)))) != 1){
    ind <- as.numeric(table(data_mani_taxa$block_var)) != length(table(data_mani_taxa$prime_var))
    data_mani_taxa <- data_mani_taxa[!(data_mani_taxa$block_var %in% names(table(data_mani_taxa$block_var))[ind]),]
  }
  
  p.cat <- names(table(data_mani_taxa[,"prime_var"]))
  ref.ind <- which(data_mani_taxa[, "prime_var"] == p.cat[1])
  com.ind <- which(data_mani_taxa[, "prime_var"] == p.cat[2])
  
  out <- matrix(NA, ncol(taxa), 7)
  
  for (i in 1:ncol(taxa)){
    fit <- t.test(data_mani_taxa[,i+2][com.ind], data_mani_taxa[,i+2][ref.ind], paired = TRUE)
    out[i,] <- c(fit$statistic, fit$estimate, fit$stderr, fit$parameter, fit$conf.int, fit$p.value)
  }
  out <- data.frame(out)
  rownames(out) <- colnames(taxa)
  colnames(out) <- c("t", "Est", "SE", "DF", "Lower", "Upper", "P.value")
  return(out)
}

taxa.bin.t.paired.united <- function(taxa, sam_dat, prim_id, block_id){
  
  t.test.p <- list() 
  for (i in 1:6){
    taxon <- taxa[[i]]
    data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
    t.test.p[[i]] <- taxa.t.paired(taxon, data.mani.taxa)
  }
  names(t.test.p) <- names(taxa)
  return(t.test.p)
}

bin.q.united.func <- function(taxa.out, method = "BH") {
  q.out <- list()
  for(i in 1:6) {
    q.out[[i]] <- bin.q.func(taxa.out[[i]], method)
  }
  names(q.out) <- names(taxa.out)
  return(q.out)
}


taxa.wilcox.sign.paired= function(taxa, data_mani_taxa){
  
  if (length(unique(as.numeric(table(data_mani_taxa$prime_var)))) != 1){
    ind <- as.numeric(table(data_mani_taxa$block_var)) != length(table(data_mani_taxa$prime_var))
    data_mani_taxa <- data_mani_taxa[!(data_mani_taxa$block_var %in% names(table(data_mani_taxa$block_var))[ind]),]
  }

  p.cat <- names(table(data_mani_taxa[,"prime_var"]))
  ref.ind <- which(data_mani_taxa[, "prime_var"] == p.cat[1])
  com.ind <- which(data_mani_taxa[, "prime_var"] == p.cat[2])
  
  out <- matrix(NA, ncol(taxa), 2)
  
  for (i in 1:ncol(taxa)){
    fit <- wilcox.test(data_mani_taxa[,i+2][com.ind], data_mani_taxa[,i+2][ref.ind], paired = TRUE, correct = FALSE, conf.int = TRUE)
    out[i,] <- c(fit$statistic, fit$p.value)
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(taxa)
  colnames(out) <- c("W", "P.value")
  return(out)
}


taxa.wilcox.test.est.func <- function(bin.var, taxa, sel.ref, sel.com, q.out) {
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  
  taxa.added <- list()
  for(i in 1:6) {
    if(!is.null(taxa[[i]])){
      n.tax <- length(taxa[[i]])
      med.diff <- numeric()
      for(j in 1:n.tax) {  
        taxon <- taxa[[i]][,j]
        
        med.ref <- median(taxon[ind.ref], na.rm = TRUE)
        med.com <- median(taxon[ind.com], na.rm = TRUE)
        
        med.diff[j] <- med.com - med.ref
      }
      q.out[[i]] <- cbind(Est = med.diff, q.out[[i]])
    }
  }
  return(q.out)
}


taxa.bin.wilcox.paired.united <- function(taxa, sam_dat, prim_id, block_id){
  wilcox.test.p <- list() 
  for (i in 1:6){
    taxon <- taxa[[i]]
    data.mani.taxa <- data_mani_p_taxa(taxa[[i]], sam_dat, prim_id, block_id)
    wilcox.test.p[[i]] <- taxa.wilcox.sign.paired(taxon, data.mani.taxa)
  }
  names(wilcox.test.p) <- names(taxa)
  return(wilcox.test.p)
}

no_constant_ldm <- function(taxa){
  result <- c() 
  for (i in 1:ncol(taxa)){
    result_ind <- unique(taxa[,i] == 0)
    if (length(result_ind) == 1){
      if(result_ind == TRUE){
        result <- c(result, i)
      }
    }
  }
  if(length(result) == 0){
    taxa.out <- taxa
  }else {
    taxa.out <- taxa[,-result]
  }
  return(taxa.out)
}

taxa.ldm.paired.united_trans <- function(sam_dat, taxa, prim, block, n.perm = 3000){
  
  ldm.test.p <- list() 
  data_mani <- data.frame(prim.var = sam_dat[[prim]], block.id = sam_dat[[block]], row.names = rownames(sam_dat))
  ldm_1 <- ldm(formula = taxa_1 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_1 <- ldm_1$F.otu.tran
  p_val_1 <- ldm_1$p.otu.tran
  q_val_1 <- p.adjust(ldm_1$p.otu.tran, "BH")
  ldm.test.p[[1]] <- data.frame(matrix(c(unlist(f_val_1), unlist(p_val_1), unlist(q_val_1)), nrow=length(q_val_1), dimnames= list(colnames(q_val_1), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[1]]) <- names(taxa[[1]])
  
  ldm_2 <- ldm(formula = taxa_2 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_2 <- ldm_2$F.otu.tran
  p_val_2 <- ldm_2$p.otu.tran
  q_val_2 <- p.adjust(ldm_2$p.otu.tran, "BH")
  ldm.test.p[[2]] <- data.frame(matrix(c(unlist(f_val_2), unlist(p_val_2), unlist(q_val_2)), nrow=length(q_val_2), dimnames= list(colnames(q_val_2), c("F", "P.value","Q.value"))))
  rownames(ldm.test.p[[2]]) <- names(taxa[[2]])
  
  ldm_3 <- ldm(formula = taxa_3 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_3 <- ldm_3$F.otu.tran
  p_val_3 <- ldm_3$p.otu.tran
  q_val_3 <- p.adjust(ldm_3$p.otu.tran, "BH")
  ldm.test.p[[3]] <- data.frame(matrix(c(unlist(f_val_3),unlist(p_val_3),  unlist(q_val_3)), nrow=length(q_val_3), dimnames= list(colnames(q_val_3), c("F","P.value", "Q.value"))))
  rownames(ldm.test.p[[3]]) <- names(taxa[[3]])
  
  ldm_4 <- ldm(formula = taxa_4 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_4 <- ldm_4$F.otu.tran
  p_val_4 <- ldm_4$p.otu.tran
  q_val_4 <- p.adjust(ldm_4$p.otu.tran, "BH")
  ldm.test.p[[4]] <- data.frame(matrix(c(unlist(f_val_4), unlist(p_val_4), unlist(q_val_4)), nrow=length(q_val_4), dimnames= list(colnames(q_val_4), c("F","P.value", "Q.value"))))
  rownames(ldm.test.p[[4]]) <- names(taxa[[4]])
  
  ldm_5 <- ldm(formula = taxa_5 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_5 <- ldm_5$F.otu.tran
  p_val_5 <- ldm_5$p.otu.tran
  q_val_5 <- p.adjust(ldm_5$p.otu.tran, "BH")
  ldm.test.p[[5]] <- data.frame(matrix(c(unlist(f_val_5), unlist(p_val_5), unlist(q_val_5)), nrow=length(q_val_5), dimnames= list(colnames(q_val_5), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[5]]) <- names(taxa[[5]])
  
  ldm_6 <- ldm(formula = taxa_6 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_6 <- ldm_6$F.otu.tran
  p_val_6 <- ldm_6$p.otu.tran
  q_val_6 <- p.adjust(ldm_6$p.otu.tran, "BH")
  ldm.test.p[[6]] <- data.frame(matrix(c(unlist(f_val_6), unlist(p_val_6),  unlist(q_val_6)), nrow=length(q_val_6), dimnames= list(colnames(q_val_6), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[6]]) <- names(taxa[[6]])
  
  names(ldm.test.p) <- names(taxa)
  return(ldm.test.p)
}

taxa.ldm.paired.united_freq <- function(sam_dat, taxa, prim, block, n.perm = 3000){
  
  ldm.test.p <- list() 
  data_mani <- data.frame(prim.var = sam_dat[[prim]], block.id = sam_dat[[block]], row.names = rownames(sam_dat))
  
  ldm_1 <- ldm(formula = taxa_1 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_1 <- ldm_1$F.otu.freq
  p_val_1 <- ldm_1$p.otu.freq
  q_val_1 <- p.adjust(ldm_1$p.otu.freq, "BH")
  ldm.test.p[[1]] <- data.frame(matrix(c(unlist(f_val_1), unlist(p_val_1), unlist(q_val_1)), nrow=length(q_val_1), dimnames= list(colnames(q_val_1), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[1]]) <- names(taxa[[1]])
  
  ldm_2 <- ldm(formula = taxa_2 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_2 <- ldm_2$F.otu.freq
  p_val_2 <- ldm_2$p.otu.freq
  q_val_2 <- p.adjust(ldm_2$p.otu.freq, "BH")
  ldm.test.p[[2]] <- data.frame(matrix(c(unlist(f_val_2), unlist(p_val_2), unlist(q_val_2)), nrow=length(q_val_2), dimnames= list(colnames(q_val_2), c("F", "P.value","Q.value"))))
  rownames(ldm.test.p[[2]]) <- names(taxa[[2]])
  
  ldm_3 <- ldm(formula = taxa_3 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_3 <- ldm_3$F.otu.freq
  p_val_3 <- ldm_3$p.otu.freq
  q_val_3 <- p.adjust(ldm_3$p.otu.freq, "BH")
  ldm.test.p[[3]] <- data.frame(matrix(c(unlist(f_val_3),unlist(p_val_3),  unlist(q_val_3)), nrow=length(q_val_3), dimnames= list(colnames(q_val_3), c("F","P.value", "Q.value"))))
  rownames(ldm.test.p[[3]]) <- names(taxa[[3]])
  
  ldm_4 <- ldm(formula = taxa_4 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_4 <- ldm_4$F.otu.freq
  p_val_4 <- ldm_4$p.otu.freq
  q_val_4 <- p.adjust(ldm_4$p.otu.freq, "BH")
  ldm.test.p[[4]] <- data.frame(matrix(c(unlist(f_val_4), unlist(p_val_4), unlist(q_val_4)), nrow=length(q_val_4), dimnames= list(colnames(q_val_4), c("F","P.value", "Q.value"))))
  rownames(ldm.test.p[[4]]) <- names(taxa[[4]])
  
  ldm_5 <- ldm(formula = taxa_5 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_5 <- ldm_5$F.otu.freq
  p_val_5 <- ldm_5$p.otu.freq
  q_val_5 <- p.adjust(ldm_5$p.otu.freq, "BH")
  ldm.test.p[[5]] <- data.frame(matrix(c(unlist(f_val_5), unlist(p_val_5), unlist(q_val_5)), nrow=length(q_val_5), dimnames= list(colnames(q_val_5), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[5]]) <- names(taxa[[5]])
  
  ldm_6 <- ldm(formula = taxa_6 | as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
  f_val_6 <- ldm_6$F.otu.freq
  p_val_6 <- ldm_6$p.otu.freq
  q_val_6 <- p.adjust(ldm_6$p.otu.freq, "BH")
  ldm.test.p[[6]] <- data.frame(matrix(c(unlist(f_val_6), unlist(p_val_6),  unlist(q_val_6)), nrow=length(q_val_6), dimnames= list(colnames(q_val_6), c("F","P.value","Q.value"))))
  rownames(ldm.test.p[[6]]) <- names(taxa[[6]])
  
  names(ldm.test.p) <- names(taxa)
  return(ldm.test.p)
}

taxa_pair_lmm <- function(taxon, sam_dat, prim_id, block_id){
  taxa.prop.ind <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
  taxa.prop.ind <<- taxa.prop.ind
  
  result <- NA
  for (i in 3:ncol(taxa.prop.ind)){
    fit<- lmer(taxa.prop.ind[[i]] ~ taxa.prop.ind[,2]  + (1|taxa.prop.ind[,1]), data = taxa.prop.ind)  
    result_ind <- coef(summary(fit))[2,]
    result <- rbind(result, result_ind)
  }
  result <- result[2:nrow(result),]
  rownames(result) <- names(taxon)
  
  result <- cbind(result, p.adjust(result[,5], "BH"))
  colnames(result) <- c("Estimate", "Std Err", "DF", "t", "P.value", "Q.value")
  
  result <- data.frame(result)
  return(result)
}

taxa.pair.lmm.united <- function(taxa, sam_dat, prim_id, block_id){
  list <- list() 
  for (i in 1:length(taxa)){
    list[[i]] <- taxa_pair_lmm(taxa[[i]], sam_dat, prim_id, block_id)
  }
  names(list) <- names(taxa)
  return(list)
}

taxa.pair.bin.boxplot <- function(bin.var, taxa.out, t.test.q.out, taxa.names.out, page, mult.test.cor = TRUE, sam_dat, prim_id, block_id) {  ## the taxa.out is actually taxa.out$clr
  
  sig.list <- list()
  
  if(mult.test.cor) {
    for(i in 1:6) {
      sig.list[[i]] <- which(t.test.q.out[[i]]$Q.value < 0.05)
    }
  } else {
    for(i in 1:6) {
      sig.list[[i]] <- which(t.test.q.out[[i]]$P.value < 0.05)
    }
  }
  
  j <- page
  
  nrow <- ceiling(length(sig.list[[j]])/4)
  
  if(nrow > 0){
    nrow <- ceiling(length(sig.list[[j]])/4)
    par(mfrow = c(nrow,4))
    id <- 0
    if(page > 1) {
      for(m in 1:(page-1)){
        id <- id+length(sig.list[[m]])
      }
    }
    
    for(k in 1:length(sig.list[[j]])) {
      if (grepl(":",taxa.names.out$names[[j]][sig.list[[j]][k]])) {
        name<- gsub(":","_", taxa.names.out$names[[j]][sig.list[[j]][k]])
      } else {
        name <- taxa.names.out$names[[j]][sig.list[[j]][k]]
      }
      if(mult.test.cor){
        qval.clr <- t.test.q.out[[j]][sig.list[[j]][k],"Q.value"]
        round.qval.clr <- format(round(qval.clr, digits = 3), nsmall = 3)
        xlab.v = paste("*q:", round.qval.clr, sep="")
      } else {
        if(!mult.test.cor){
          pval.clr <- t.test.q.out[[j]][sig.list[[j]][k],"P.value"]
          round.pval.clr <- format(round(pval.clr, digits = 3), nsmall = 3)
          xlab.v = paste("*p:", round.pval.clr, sep="")
        }
      }
      
      id <- id+1 
      
      taxon <- as.matrix(taxa.out[[j]][sig.list[[j]][k]]) 
      
      data.mani.taxon <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id) 
      
      time.p.cat <- names(table(data.mani.taxon[,"prime_var"]))
      
      ref.ind <- which(data.mani.taxon[,"prime_var"] == time.p.cat[1])
      com.ind <- which(data.mani.taxon[,"prime_var"] == time.p.cat[2])
      
      taxa_ref <- data.mani.taxon[ref.ind, 3:ncol(data.mani.taxon)]
      taxa_com <- data.mani.taxon[com.ind, 3:ncol(data.mani.taxon)]
      boxplot(taxa_ref, taxa_com, main = id, xlab = xlab.v, ylab = name, names = c(time.p.cat[1], time.p.cat[2]), notch = TRUE, col=c(rgb(0,1,0,0.5), rgb(1,0,0,0.5))) 
      #names = substr(levels(bin.var), 1, 8)
    }
  } else {
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot.new()
    text(x = 0.5, y = 0.5, paste("No significant taxa are found in ", ranks[page], sep = ""), 
         cex = 1.2, col = "black")
  }
}

taxa.sig.dend <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(NA, 1, 1)
  for (i in 1:6) {
    ind.sig <- which(out[[i]]$Q.value < 0.05)
    if (length(ind.sig) >= 1) {
      sig.out <- out[[i]][ind.sig,]
      taxa <- rownames(sig.out)
      print(sig.out)
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.value"], 3), nsmall = 3))
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    } 
  }
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(!(tax.tab.num[,tax.rank] %in% sig.taxa))
      non.sig.taxa <- names(table(tax.tab.num[ind.non.sig,tax.rank]))
      
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0) 
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!is.na(Family.desc)]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        
        connection.with.upper[[5]] <- connection 
      }
    }
  }else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
    
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  sig.Family <- Family[!grepl("Family", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  #ind <- text.tab.all[,1] %in% Genus[!grepl("Genus", Genus)]
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  sig.Species <- Species[!grepl("Species", Species)]
  #ind <- text.tab.all[,1] %in% Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Species <- sig.Species[ind.pos.sig] 
  neg.sig.Species <- sig.Species[ind.neg.sig]
  
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa)] 
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa)] 
  }
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 5
  
  neg <- 1
  non <- 1
  
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa)){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa)){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  for(i in 1:length(dummy)){
    flow.text.dum[[5]][i] <- paste(dummy[i], "[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " ")
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:5){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  print(grViz(flow.text))
  
}

taxa.sig.dend_2 <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE", is.LDM = FALSE) {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  if (is.LDM == FALSE){
    text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
    ci.tab.all <- matrix(NA, 1, 1)
    for (i in 1:6) {
      ind.sig <- which(out[[i]]$Q.value < 0.05)
      if (length(ind.sig) >= 1) {
        sig.out <- out[[i]][ind.sig,]
        taxa <- rownames(sig.out)
        text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                          format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                          format(round(sig.out[,"Q.value"], 3), nsmall = 3))
        ci.tab <- cbind(sig.out[,1])
        text.tab.all <- rbind(text.tab.all, text.tab)
        ci.tab.all <- rbind(ci.tab.all, ci.tab)
      }
    }
  }else{
    text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "Q.value"), 1, 4)
    ci.tab.all <- matrix(NA, 1, 1)
    for (i in 1:6){
      ind.sig <- which(out[[i]]$Q.value < 0.05)
      if (length(ind.sig) >= 1){
        sig.out <- out[[i]][ind.sig, ]
        taxa <- rownames(sig.out)
        text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3),
                          format(round(sig.out[,"Q.value"], 3), nsmall = 3))
        ci.tab <- cbind(sig.out[,1])
        text.tab.all <- rbind(text.tab.all, text.tab) 
        ci.tab.all <- rbind(ci.tab.all, ci.tab) 
      }
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    } 
  }
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(!(tax.tab.num[,tax.rank] %in% sig.taxa))
      non.sig.taxa <- names(table(tax.tab.num[ind.non.sig,tax.rank]))
      
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0) 
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!is.na(Family.desc)]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        
        connection.with.upper[[5]] <- connection 
      }
    }
  }else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
    
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  sig.Family <- Family[!grepl("Family", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  #ind <- text.tab.all[,1] %in% Genus[!grepl("Genus", Genus)]
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  sig.Species <- Species[!grepl("Species", Species)]
  #ind <- text.tab.all[,1] %in% Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Species <- sig.Species[ind.pos.sig] 
  neg.sig.Species <- sig.Species[ind.neg.sig]
  
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa)] 
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa)] 
  }
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 5
  
  neg <- 1
  non <- 1
  
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa)){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa)){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  for(i in 1:length(dummy)){
    flow.text.dum[[5]][i] <- paste(dummy[i], "[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " ")
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:5){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  print(grViz(flow.text))
}


add_NA <- function(taxa.out, tax.tab) {
  taxa.out.ori <- taxa.out
  tax.tab <- as.data.frame(tax.tab)
  for(type in 1:length(taxa.out)) {
    na.num <- 1
    for (rank in 1:6) {
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))) {
        if(substring(colnames(taxa.out[[type]][[rank]])[i], nchar(colnames(taxa.out[[type]][[rank]])[i])-2) == "NA ") {
          colnames(taxa.out[[type]][[rank]])[i] <- paste0(colnames(taxa.out[[type]][[rank]])[i], na.num)
          na.num <- na.num + 1
        }
      }
    }
    for (rank in 1:5) {
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))){
        colnames(taxa.out[[type]][[rank+1]]) <- str_replace(colnames(taxa.out[[type]][[rank+1]]), colnames(taxa.out.ori[[type]][[rank]])[i], colnames(taxa.out[[type]][[rank]])[i])
      }
    }
    
    for(rank in 1:6) {
      for(i in 1:length(unique(colnames(taxa.out[[type]][[rank]])))) {
        ind <- tax.tab[[rank+1]] == unique(colnames(taxa.out.ori[[type]][[rank]]))[i]
        
        tax.tab[[rank+1]][ind] <- unique(colnames(taxa.out[[type]][[rank]]))[i]
      }
    }
  }
  return(list(taxa.out = taxa.out, tax.tab = tax.tab))
}



taxa.sig.dend <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(NA, 1, 1)
  for (i in 1:6) {
    ind.sig <- which(out[[i]]$Q.value < 0.05)
    if (length(ind.sig) >= 1) {
      sig.out <- out[[i]][ind.sig,]
      taxa <- rownames(sig.out)
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.value"], 3), nsmall = 3))
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}



tax.tab.cleanNA <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  # print(head(tax.tab.c))
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- "NA "
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
  
  ranks <- c("p_", "c_", "o_", "f_", "g_", "s_")
  for (i in 1:6) {
    na.num <- 1
    ind <- which(tax.tab.h[,i+1] != na.code)
    ind_omit <- which(tax.tab.h[,i+1] != na.code & tax.tab.h[,i] == na.code)
    tax.tab.h[ind ,i+1] <- paste(tax.tab.h[ind,i],paste(ranks[i],tax.tab.h[ind,i+1], sep = ""), sep = ";")
    # tax.tab.h[ind ,i+1] <- apply(as.matrix(tax.tab.h[ind ,i+1]), 1, function(x) 
    #   if(substring(x, nchar(x)-2) == "_NA") {
    #     x <- paste0(x, na.num); na.num <- na.num +1
    #     print(na.num)
    #   } else { 
    #     x <- x
    #   })
    # for(k in 1:length(tax.tab.h[ind, i+1])){
    #   if(substring(tax.tab.h[ind, i+1][k], nchar(tax.tab.h[ind, i+1][k])-2) == "_NA") {
    #     print(tax.tab.h[ind, i+1][k])
    #     tax.tab.h[ind, i+1] <- paste0(substring(tax.tab.h[ind, i+1], na.num))
    #     na.num <- na.num + 1          
    #   }
    # }
    # if(substring(tax.tab.h[ind ,i+1], nchar(tax.tab.h[ind ,i+1])-2))
    
    if(length(ind_omit)!=0) tax.tab.h[ind_omit,c((i+1):7)] = na.code
  }
  
  # ind <- which(tax.tab.h[,2] != na.code)
  # ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  # tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,3] != na.code)
  # ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  # tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,4] != na.code)
  # ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  # tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,5] != na.code)
  # ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  # tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,6] != na.code)
  # ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  # tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,7] != na.code)
  # ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  # tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}



biom.cleanSNA <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, surv.time, censor, follow, subgroup.var = NULL, subgroup = NULL) {
  #surv.time, censor, follow, subgroup = NULL
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

  if (!is.null(subgroup.var)){
    ind <- sam.dat [[subgroup.var]] == subgroup
    sam.dat  <- sam.dat [ind, ]
  }
  
  otu.tab <- otu.tab[,colnames(otu.tab) %in% rownames(sam.dat)]
  
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
    tax.tab <- tax.tab.cleanNA(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
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



tax.trans.na <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "") {
  
  # tax.tab <- matrix(nrow = nrow(tax.tab.ori), ncol = ncol(tax.tab.ori))
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    
    na.num <- 1
    for(i in 1:length(tax)) {
      na.taxa <- substring(tax[i], nchar(tax[i])-2)
      if(na.taxa == "_NA") {
        new.tax <- paste0(tax[i], na.num)
        ind.na.tab <- which(tax.tab[,j+1] == tax[i])
        tax.tab[ind.na.tab, j+1] <<- new.tax
        rare.tax.tab[ind.na.tab, j+1] <<- new.tax
        tax[i] <- new.tax
        na.num <- na.num + 1
      }
      
    }
    
    
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
  
  
}
