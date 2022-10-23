#####################
# Data manipulation #
#####################
taxa.bin.cat.ref.united.func.mult.pairwise <- function(sel.bin.var, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  taxa.out <- taxa
  return(list(bin.var = bin.var, taxa = taxa.out))
}

reduced_data <- function(sam_dat, prim_id, levels){
  sam.dat <- sam_dat[sam_dat[[prim_id]] %in% levels, ]
  return(sam.dat)
}

reduced_taxa <- function(taxa, sam_dat, prim_id, block_id, levels){
  taxa.reduced <- list() 
  for (i in 1:length(taxa)){
    taxa.ind <- taxa[[i]]
    taxa.tran.ind <- cbind(prim.var = sam_dat[[prim_id]], block.id = sam_dat[[block_id]], taxa.ind)
    taxa.ind <- taxa.tran.ind[taxa.tran.ind[,"prim.var"] %in% levels,]
    taxa.reduced[[i]] <- taxa.ind[, 3:ncol(taxa.ind)]
  }
  names(taxa.reduced) <- names(taxa)
  return(taxa.reduced)
}

reduced_tax_ad <- function(re_tax){
  tax <- list()
  for (i in 1:length(re_tax)){
    tax[[i]] <- re_tax[[i]][, 3:length(re_tax[[i]])]
  }
  return(tax[[i]])
}

taxa_mani <- function(taxa, sam_dat, prim_var, species = TRUE){
  taxa.list <- list() 
  if (species){
    for (i in 1:length(taxa)){
      taxa.list[[i]] <- cbind(taxa[[i]], sam_dat)
      taxa.list[[i]][, prim_var] <- as.factor(taxa.list[[i]][, prim_var])
    }
    names(taxa.list) <- names(taxa)
  }else{
    for (i in 1:length(taxa)){
      taxa.list[[i]] <- cbind(taxa[[i]], sam_dat)
      taxa.list[[i]][, prim_var] <- as.factor(taxa.list[[i]][, prim_var])
    }
    names(taxa.list) <- names(taxa)[1:5]
  }
  
  return(taxa.list)
}

taxa_summary<- function(summary, ref_level, prim_var, taxa_prop_ind){
  summary<- summary[2:nrow(summary),]
  result <- matrix(NA, nrow = nrow(summary) , ncol = 5)
  
  levels <- levels(as.factor(taxa_prop_ind[, prim_var]))

  for (i in 1:nrow(summary)){
    ref <- rep(ref_level, length(levels)-1)
    com <- levels[2:length(levels)]
  }
  
  for (i in 1:ncol(summary)){
    result[,i] <- noquote(decimal_adjust(summary[,i]))
  }
  
  result <- cbind(ref, com, result)
  
  colnames(result) <- c("Ref", "Com", "Estimate", "StdErr", "DF", "t", "P.value")
  result <- data.frame(result)
  return(result)
}


taxa.bin.mult.paired <- function(taxa, sam_dat, prim_id, block_id, page, q_val_result){
  
  taxon <- taxa[[page]]
  data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
  
  time.p.cat <- names(table(data.mani.taxa[,"prime_var"]))
  
  q_val_result[[page]][!complete.cases(as.numeric(q_val_result[[page]]))] <- 1 
  sig.num <- length(q_val_result[[page]][q_val_result[[page]] < 0.05])
  
  nrow <- ceiling(sig.num/4)
  
  if (nrow > 0){
    par(mfrow = c(nrow, 4))
    
    ind <- list()
    for (i in 1:length(time.p.cat)){
      ind[[i]] <- which(data.mani.taxa[,"prime_var"] == time.p.cat[i])
    }
    
    each_in_taxon <- list() 
    meaningful = c() 
    taxa_dat_list <- list() 
    id <- 0 
    
    
    for (j in 1:length(taxon))
    {
      if (j %in% which(q_val_result[[page]] < 0.05)){ 
        meaningful = c(meaningful, j)
        id <- id+1
        taxa_dat <- c() 
        for (i in 1:length(time.p.cat)){
          taxa_dat <- cbind(taxa_dat, data.mani.taxa[, 2+j][ind[[i]]])
        }
        colnames(taxa_dat) <- time.p.cat
        taxa_dat_list[[id]] <- taxa_dat
        
        for (k in 1:length(meaningful)){
          each_in_taxon[[k]] <- taxa_dat
        }
      }
    }
    
    xlab.v <- c() 
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      xlab.v[i] <- paste0("*q:", p.value.0.1_char(decimal_adjust(q_val_result[[page]][q_val_result[[page]]<0.05][i])), sep = "")
    }
    
    ylab.v <- c() 
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      ylab.v[i] <- rev(strsplit(names(taxon)[q_val_result[[page]] < 0.05][[i]], ";")[[1]])[1]
    }
    
    for (i in 1:length(q_val_result[[page]][q_val_result[[page]]<0.05])){
      boxplot(taxa_dat_list[[i]], main = i, xlab = xlab.v[i], ylab = ylab.v[i], names = time.p.cat, col = "Paleturquoise1",  notch = TRUE , boxwex=0.9, horizontal = FALSE)
    }
    
  }else{
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot.new()
    text(x = 0.5, y = 0.5, paste("No significant taxa are found in ", ranks[page], sep = ""), cex = 1.2, col = "black")
  }
}

#########################################################
# Comparative analysis for multi-paired taxa diversity  #
#########################################################

taxa.t.pair.mult.ind <- function(taxa_ind, data_mani_taxa, method_adj_ind = "BH"){
  ind <- list()
  time.p.cat <- names(table(data_mani_taxa[,"prime_var"]))
  
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(data_mani_taxa[,"prime_var"] == time.p.cat[i])
  }
  comb <- combn(length(time.p.cat), 2)
  
  each_rank <- list() 
  for (i in 1:ncol(taxa_ind)){
    out <- matrix(NA, ncol(comb), 8)
    for (j in 1:ncol(comb)){
      combination <- comb[,j]
      fit <- t.test(data_mani_taxa[,i+2][ind[[combination[1]]]], data_mani_taxa[,i+2][ind[[combination[2]]]], paired = TRUE)
      out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], -fit$statistic, fit$stderr, fit$parameter, -fit$conf.int, fit$p.value)
    }
    q_val <- p.adjust(out[,8], method = method_adj_ind, n = length(out[,8]))
    
    out <- as.data.frame(cbind(out, q_val))
    colnames(out) <- c("Ref", "Com", "t", "Std Err", "DF", "Lower", "Upper", "P.value", "Adj. P.value")
    each_rank[[i]] <- out 
  }
  names(each_rank) <- colnames(taxa_ind)  
  return(each_rank)
}

taxa.t.pair.mult.united <- function(taxa, prim_id, block_id, sam_dat, method_adj = "BH" ){
  t.test.p <- list() 
  
  for (i in 1:length(taxa)){
    taxon <- taxa[[i]]
    data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)  
    t.test.p[[i]] <- taxa.t.pair.mult.ind(taxon, data.mani.taxa, method_adj_ind = method_adj )
  }
  
  names(t.test.p) <- names(taxa) 
  return(t.test.p)
}

taxa.wilcox.pair.mult.ind <- function(taxa_ind, data_mani_taxa, method_adj_ind = "BH"){
  ind <- list()
  time.p.cat <- names(table(data_mani_taxa[,"prime_var"]))
  
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(data_mani_taxa[,"prime_var"] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  each_rank <- list() 
  
  for (i in 1:ncol(taxa_ind)){
    out <- matrix(NA, ncol(comb), 4)
    for (j in 1:ncol(comb)){
      combination <- comb[,j]
      fit <- wilcox.test(data_mani_taxa[,i+2][ind[[combination[1]]]], data_mani_taxa[,i+2][ind[[combination[2]]]], paired = TRUE, correct = FALSE)
      out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], -fit$statistic, fit$p.value)
    }
    
    q_val <- p.adjust(out[,4], method = method_adj_ind, n = length(out[,4]))
    out <- as.data.frame(cbind(out, q_val))
    colnames(out) <- c("Ref", "Com", "t", "P.value", "Adj.P.value")
    each_rank[[i]] <- out 
  }
  names(each_rank) <- colnames(taxa_ind)
  return(each_rank)
}

taxa.wilcox.pair.mult.united <- function(taxa, prim_id, block_id, sam_dat, method_adj = "BH" ){
  wilcox.test.p <- list() 
  
  for (i in 1:length(taxa)){
    taxon <- taxa[[i]]
    data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
    wilcox.test.p[[i]] <- taxa.wilcox.pair.mult.ind(taxon, data.mani.taxa, method_adj_ind = method_adj)
  }
  
  names(wilcox.test.p) <- names(taxa) 
  return(wilcox.test.p)
}

taxa.ldm.paired.mult.ind <- function(sam_dat, taxon, prim_id, block_id, n.perm){
  
  data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
  time.p.cat <- names(table(data.mani.taxa[,"prime_var"]))
  
  ind <- list()
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(data.mani.taxa[,"prime_var"] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  ref_conf_list <- c()
  each_rank <- list() 
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    taxa_data <- data.mani.taxa[c(ind[[combination[1]]], ind[[combination[2]]]),]
    taxa_data_in <<- taxa_data[, 3:ncol(taxa_data)]
    
    sam_dat_mani <- data_mani_p_sample(sam_dat, prim_id, block_id)
    data_ind <- sam_dat_mani[c(ind[[combination[1]]], ind[[combination[2]]]),]
    data_mani <- data.frame(prim.var = data_ind[, "prime_var"], block.id = data_ind[, "block_var"], row.names = rownames(data_ind))
    
    ldm <- ldm(formula = taxa_data_in| as.factor(block.id) ~ prim.var, data= data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
    
    dat <- cbind(rep(time.p.cat[combination[1]], ncol(taxa_data_in)), rep(time.p.cat[combination[2]], ncol(taxa_data_in)))
    rownames(dat) <- colnames(taxa_data_in)
    dat_2 <- cbind(ldm$F.otu.freq, ldm$p.otu.freq, p.adjust(as.vector(ldm$p.otu.freq), "BH"))

    
    merged_dat <- merge(dat, dat_2, by= 0, all = TRUE)
    rownames(merged_dat) <- merged_dat$Row.names
    merged_dat <- merged_dat[, 2:ncol(merged_dat)]
    colnames(merged_dat) <- c("Ref", "Com", "F", "P.value", "Adj.P.value")
    
    merged_dat[,3] <- decimal_adjust(merged_dat[,3])
    merged_dat[,4] <- decimal_adjust(merged_dat[,4])
    merged_dat[,5] <- decimal_adjust(merged_dat[,5])
    each_rank[[j]] <- merged_dat
  } 
  return(each_rank)
}

taxa.ldm.paired.mult.united <- function(sam_dat, taxa, prim_id, block_id, n.perm){
  united_list <- list() 
  
  for (i in 1:length(taxa)){
    taxa_result <- taxa.ldm.paired.mult.ind(sam_dat, taxa[[i]], prim_id, block_id, n.perm)
    taxa_list <- list() 
    for (k in 1:nrow(taxa_result[[1]])){
      taxa_result_k <- NULL
      for (j in 1:length(taxa_result)){
        taxa_result_k <- rbind(taxa_result_k, taxa_result[[j]][k,])
        rownames(taxa_result_k) <- NULL
      }
      taxa_list[[k]] <- taxa_result_k
    }
    names(taxa_list) <- rownames(taxa_result[[1]])
    united_list[[i]] <- taxa_list
  }
  names(united_list) <- names(taxa)
  return(united_list)
}

taxa.ldm.paired.mult.ind_trans <- function(sam_dat, taxon, prim_id, block_id, n.perm){
  
  data.mani.taxa <- data_mani_p_taxa(taxon, sam_dat, prim_id, block_id)
  time.p.cat <- names(table(data.mani.taxa[,"prime_var"]))
  
  ind <- list()
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(data.mani.taxa[,"prime_var"] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  ref_conf_list <- c()
  each_rank <- list() 
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    taxa_data <- data.mani.taxa[c(ind[[combination[1]]], ind[[combination[2]]]),]
    taxa_data_in <<- taxa_data[, 3:ncol(taxa_data)]
    
    sam_dat_mani <- data_mani_p_sample(sam_dat, prim_id, block_id)
    data_ind <- sam_dat_mani[c(ind[[combination[1]]], ind[[combination[2]]]),]
    data_mani <- data.frame(prim.var = data_ind[, "prime_var"], block.id = data_ind[, "block_var"], row.names = rownames(data_ind))
    
    ldm <- ldm(formula = taxa_data_in| as.factor(block.id) ~ prim.var, data = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
    
    dat <- cbind(rep(time.p.cat[combination[1]], ncol(taxa_data_in)), rep(time.p.cat[combination[2]], ncol(taxa_data_in)))
    rownames(dat) <- colnames(taxa_data_in)
    dat_2 <- cbind(ldm$F.otu.tran, ldm$p.otu.tran, p.adjust(as.vector(ldm$p.otu.tran)))
    
    merged_dat <- merge(dat, dat_2, by= 0, all = TRUE)
    rownames(merged_dat) <- merged_dat$Row.names
    merged_dat <- merged_dat[, 2:ncol(merged_dat)]
    colnames(merged_dat) <- c("Ref", "Com", "F", "P.value", "Adj.P.value")
    
    merged_dat[,3] <- decimal_adjust(merged_dat[,3])
    merged_dat[,4] <- decimal_adjust(merged_dat[,4])
    merged_dat[,5] <- decimal_adjust(merged_dat[,5])
    each_rank[[j]] <- merged_dat
  } 
  return(each_rank)
}

taxa.ldm.paired.mult.united_trans <- function(sam_dat, taxa, prim_id, block_id, n.perm){
  united_list <- list() 
  
  for (i in 1:length(taxa)){
    taxa_result <- taxa.ldm.paired.mult.ind_trans(sam_dat, taxa[[i]], prim_id, block_id, n.perm)
    taxa_list <- list() 
    for (k in 1:nrow(taxa_result[[1]])){
      taxa_result_k <- NULL
      for (j in 1:length(taxa_result)){
        taxa_result_k <- rbind(taxa_result_k, taxa_result[[j]][k,])
        rownames(taxa_result_k) <- NULL
      }
      taxa_list[[k]] <- taxa_result_k
    }
    names(taxa_list) <- rownames(taxa_result[[1]])
    united_list[[i]] <- taxa_list
  }
  names(united_list) <- names(taxa)
  return(united_list)
}

remove_zero_col <- function(data){
  num_zero <- sapply(data, function(x)sum(x != 0))
  col_zero <- c() 
  for (i in 1:length(num_zero)){
    if (num_zero[[i]] == 0){
      col_zero <- c(col_zero, i)
    }
  }
  return(col_zero)
}

taxa.ldm.mult.glob <- function(sam_dat, taxa, prim_id, block_id, n.perm = 3000, trans = FALSE) {
  output_result <- list()
  result <- list()  
  q_val_result <- list()
  for (i in 1:length(names(taxa))){
    
    data_mani <- data.frame(prim.var = sam_dat[[prim_id]], block.id = sam_dat[[block_id]], row.names = rownames(sam_dat))
    taxa.ind <<- taxa[[i]]
    out <- matrix(NA, length(taxa.ind), 3)
    
    fit <- ldm(taxa.ind | as.factor(block.id) ~ prim.var, dat = data_mani, cluster.id = block.id, perm.within.type = "free", perm.between.type = "none", n.perm.max = n.perm)
    if (trans){
      out <- cbind(as.vector(fit[["F.otu.tran"]]), as.vector(fit[["p.otu.tran"]]), as.vector(p.adjust(fit[["p.otu.tran"]], "BH")))  #as.vector(fit[["q.otu.tran"]]))
      each_p_result <- as.vector(fit$p.otu.tran)
      each_q_result <- decimal_adjust(p.adjust(each_p_result, "BH"))
    }else{
      out <-cbind(as.vector(fit[["F.otu.freq"]]), as.vector(fit[["p.otu.freq"]]), as.vector(p.adjust(fit[["p.otu.freq"]], "BH"))) #as.vector(fit[["q.otu.freq"]]))
      each_p_result <- as.vector(fit$p.otu.freq)
      each_q_result <- decimal_adjust(p.adjust(each_p_result, "BH"))
    }
    
    out <- data.frame(out)
    colnames(out) <- c("F", "P.value", "Adj.P.value")
    
    if (length(remove_zero_col(taxa.ind)) != 0){
      rownames(out) <- names(taxa.ind)[-remove_zero_col(taxa.ind)]
    }else{
      rownames(out) <- names(taxa.ind)
    }
    result[[i]] <- out 
    q_val_result[[i]] <- each_q_result 
  }
  names(result) <- names(taxa)
  names(q_val_result) <- names(taxa)
  
  output_result$global_test <- result 
  output_result$pval <- q_val_result
  return(output_result)
}  

taxa.data.mani <- function(sam_dat, taxa, prim_id, block_id, species = TRUE){
  taxa.transformed <- list() 
  if (species){
    for (i in 1:length(taxa)){
      taxa.ind <- taxa[[i]]
      taxa.tran.ind <- cbind(prim.var = sam_dat[[prim_id]], block.id = sam_dat[[block_id]], taxa.ind)
      taxa.transformed[[i]] <- taxa.tran.ind
    }
    names(taxa.transformed) <- names(taxa)
  }else{
    for (i in 1:5){
      taxa.ind <- taxa[[i]]
      taxa.tran.ind <- cbind(prim.var = sam_dat[[prim_id]], block.id = sam_dat[[block_id]], taxa.ind)
      taxa.transformed[[i]] <- taxa.tran.ind
    }
    names(taxa.transformed) <- names(taxa)[1:5]
  }
  
  return(taxa.transformed)
}

taxa.f.pair.mult.overall <- function(taxa_mani){
  output_list <- list() 
  summary_list <- list() 
  pval_list <- list() 
  
  for (i in 1:length(taxa_mani)){
    taxa.ind <- taxa_mani[[i]]
    pval_vec <- c() 
    
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    out <- matrix(NA, nrow = length(taxa.ind)-2, ncol = 2)
    for (j in 1:(length(taxa.ind)-2)){
      summary_aov <- unlist(summary(aov(taxa.ind[,j+2] ~ prim.var + block.id, data = taxa.ind)))
      p_val <- as.vector(summary_aov[["Pr(>F)2"]])
      pval_vec <- c(pval_vec, p_val)
      summary_ind <- c(as.vector(summary_aov[["F value2"]]), as.vector(summary_aov[["Pr(>F)2"]]))
      out[j,] <- summary_ind
    }
    
    pval_list[[i]] <- decimal_adjust(p.adjust(pval_vec, "BH"))
    rownames(out) <- colnames(taxa.ind[3:length(taxa.ind)])
    colnames(out) <- c("F", "P.value")
    summary_list[[i]] <- out
  }
  
  names(summary_list) <- names(taxa_mani)
  names(pval_list) <- names(taxa_mani)
  output_list$global_test <- summary_list 
  output_list$pval <- pval_list 
  
  return(output_list) 
}

taxa.f.mult.glob.p_val <- function(taxa_mani){
  pval_list <- list()
  for (i in 1:length(taxa_mani)){
    pval_vec <- c()
    taxa.ind <- taxa_mani[[i]]
    for (j in 1:(length(taxa.ind)-2)){
      summary_aov <- unlist(summary(aov(taxa.ind[,j+2] ~ prim.var + block.id, data = taxa.ind)))
      p_val <-  as.vector(summary_aov[["Pr(>F)2"]])
      pval_vec <- c(pval_vec, p_val)
    }
    pval_list[[i]] <- decimal_adjust(p.adjust(pval_vec, "BH"))
  }
  names(pval_list) <- names(taxa_mani)
  return(pval_list)
}

q_val_combined_table.anova <- function(result_table, q_val_list){
  new_list_with_q <- list() 
  for (i in 1:length(result_table)){
    new_list_with_q[[i]] <- cbind(result_table[[i]], as.vector(q_val_list[[i]]))
    colnames(new_list_with_q[[i]]) <- c("F", "P.value", "Adj.P.value")
  }
  names(new_list_with_q) <- names(result_table)
  return(new_list_with_q)
}

q_val_combined_table.friedman <- function(result_table, q_val_list){
  new_list_with_q <- list() 
  for (i in 1:length(result_table)){
    new_list_with_q[[i]] <- cbind(result_table[[i]], as.vector(q_val_list[[i]]))
    colnames(new_list_with_q[[i]]) <- c("Chisq", "DF", "P.value", "Adj.P.value")
  }
  names(new_list_with_q) <- names(result_table)
  return(new_list_with_q)
}

sig.taxa <- function(p_val){
  
  sig_list <- list() 
  for (i in 1:length(p_val)){
    sig_list[[i]] <- as.numeric(p_val[[i]]) < 0.05 
    sig_list[[i]][!complete.cases(as.numeric(p_val[[i]]))] <- FALSE
  }
  return(sig_list)
}

rank_sig <- function(pairwise_result, sig.taxa.list){
  sig.global.list <- list() 
  for (i in 1:length(pairwise_result)){
    sig.global.list[[i]] <- pairwise_result[[i]][sig.taxa.list[[i]]]
  }
  names(sig.global.list) <- names(pairwise_result)
  return(sig.global.list)
}

list_to_dat <- function(rank_sig_list){
  dat_list <- list() 
  
  for (i in 1:length(rank_sig_list)){
    names(rank_sig_list[[i]]) <- NULL
    dat_list[[i]] <- do.call(rbind.data.frame, rank_sig_list[[i]])
  }
  names(dat_list) <- names(rank_sig_list)
  return(dat_list)
}

taxa.f.pair.mult.tukey <- function(taxa_mani){
  summary_list <- list() 
  for (i in 1:length(taxa_mani)){
    taxa.ind <- taxa_mani[[i]]
    
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    each_list <- list() 
    for (j in 1:(length(taxa.ind) -2)){
      rmanova_result <- aov(taxa.ind[,j+2] ~ prim.var + block.id, data = taxa.ind)
      tukey_result <- TukeyHSD(rmanova_result)$prim.var
      ind_result <- c() 
      for (k in 1:nrow(tukey_result)){
        ref_conf_list <- strsplit(rownames(TukeyHSD(rmanova_result)$prim.var), split = "-")[[k]]
        ref <- ref_conf_list[2]
        conf <- ref_conf_list[1]
        
        ref_mean <- as.vector(mean(taxa.ind[,j+2][taxa.ind$prim.var == ref]))
        conf_mean <- as.vector(mean(taxa.ind[,j+2][taxa.ind$prim.var == conf]))
        
        tukey_result_ind <- c(ref, conf, decimal_adjust(c(ref_mean, conf_mean, tukey_result[k,])))
        ind_result <- rbind(ind_result, tukey_result_ind)
      }
      rownames(ind_result) <- 1:nrow(ind_result)
      colnames(ind_result) <- c("Ref", "Com", "Mean (Ref)", "Mean (Com)", "Diff", "Lower", "Upper", "Adj.P.value")
      each_list[[j]] <- as.data.frame(ind_result)
    }
    names(each_list) <-  colnames(taxa.ind[3:length(taxa.ind)])
    summary_list[[i]] <- each_list 
  }
  names(summary_list) <- names(taxa_mani)
  return(summary_list)
}

taxa.friedman.pair.mult.overall <- function(taxa_mani){
  output_list <- list() 
  summary_list <- list() 
  pval_list <- list() 
  for (i in 1:length(taxa_mani)){
    pval_vec <- c() 
    taxa.ind <- taxa_mani[[i]]
    
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    out <- matrix(NA, nrow = length(taxa.ind)-2, ncol = 3)
    for (j in 1:(length(taxa.ind)-2)){
      friedman_result <- friedmanTest(y = taxa.ind[, j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id)
      p_val <- as.vector(friedman_result$p.value)
      pval_vec <- c(pval_vec, p_val) 
      summary_ind <- c(decimal_adjust(as.vector(friedman_result$statistic)), friedman_result$parameter ,as.vector(friedman_result$p.value))
      out[j,] <- summary_ind
    }
  
    rownames(out) <- colnames(taxa.ind[3:length(taxa.ind)])
    colnames(out) <- c("Chisq", "DF", "P.value")
    summary_list[[i]] <- out
    pval_list[[i]] <- decimal_adjust(p.adjust(pval_vec, "BH"))
  }
  names(summary_list) <- names(taxa_mani)
  names(pval_list) <- names(taxa_mani)
  
  output_list$global_test <- summary_list 
  output_list$pval <- pval_list 
  return(output_list) 
}

fill_missing_na_taxa <- function(data.mani){
  
  data.mani$block.id <- as.character(data.mani$block.id)
  data.mani$prim.var <- as.character(data.mani$prim.var)
  
  tab <- table(data.mani$block.id, data.mani$prim.var)
  miss_rc <- which(tab == 0, arr.ind = TRUE) 
  miss_rc <- data.frame(miss_rc)
  
  miss_prim <- colnames(tab)[miss_rc$col]
  miss_block <- rownames(tab)[miss_rc$row]
  
  length <- nrow(data.mani)
  for (i in 1:length(miss_prim)){
    data.mani[length+i, ] <- c(miss_prim[i], miss_block[i], rep(NA, ncol(data.mani)-2))
  }
  
  data.mani <- data.mani[order(data.mani$block.id, data.mani$prim.var),]
  return(data.mani)
}

taxa.durbin.pair.mult.overall <- function(taxa_mani){
  output_list <- list() 
  summary_list <- list() 
  pval_list <- list() 
  for (i in 1:length(taxa_mani)){
    pval_vec <- c() 
    taxa.mani.ind <- taxa_mani[[i]]
    
    taxa.ind <- fill_missing_na_taxa(taxa.mani.ind)
    
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    out <- matrix(NA, nrow = length(taxa.ind)-2, ncol = 3)
    
    for (j in 1:(length(taxa.ind)-2)){
      durbin_result <- PMCMRplus::durbinTest(y = taxa.ind[, j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id)
      p_val <- as.vector(durbin_result$p.value)
      pval_vec <- c(pval_vec, p_val)
      summary_ind <- c(as.vector(durbin_result$statistic), durbin_result$parameter ,as.vector(durbin_result$p.value))
      out[j,] <- summary_ind
    }
    pval_list[[i]] <- pval_vec 
    rownames(out) <- colnames(taxa.ind[3:length(taxa.ind)])
    colnames(out) <- c("Chisq", "DF", "P.value")
    summary_list[[i]] <- out
  }
  names(summary_list) <- names(taxa_mani)
  names(pval_list) <- names(taxa_mani)
  
  output_list$global_test <- summary_list 
  output_list$pval <- pval_list 
  
  return(output_list) 
}

taxa.pair.mult.conover <- function(taxa_mani, p_adjustment = "BH"){
  summary_list <- list() 
  for (i in 1:length(taxa_mani)){
    taxa.ind <- taxa_mani[[i]]
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    each_list <- list() 
    for (j in 1:(length(taxa.ind)-2)){
      
      conover_statistic <- frdAllPairsConoverTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = p_adjustment)$statistic 
      
      conover_qval <- frdAllPairsConoverTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = p_adjustment)$p.value 
      
      conover_qval[is.nan(conover_qval)] <- ""
      
      conover_pval <- frdAllPairsConoverTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = "none")$p.value 
      
      ref <- c() 
      conf <- c() 
      t <- c() 
      Q.val <- c() 
      P.val <- c() 
      
      for (k in 1:length(rownames(conover_qval))){
        ref_ind <- rep(colnames(conover_qval)[k], length(na.omit(conover_pval[,k])))
        ref <- c(ref, ref_ind)
        
        conf_ind <- rownames(conover_qval)[k:length(rownames(conover_pval))]
        conf <- c(conf, conf_ind)
        
        t_ind <- as.vector(na.omit(conover_statistic[,k]))
        t <- c(t, t_ind)
        
        Q_ind <- na.omit(conover_qval[,k])
        P_ind <- decimal_adjust(na.omit(conover_pval[,k]))
        Q.val <- c(Q.val, Q_ind)
        P.val <- c(P.val, P_ind)
      }
      conover_result <- cbind(Ref = ref, Com = conf, t = decimal_adjust(t), P.value = P.val, Adj.P.value =  Q.val) 
      each_list[[j]] <- data.frame(conover_result)
      rownames(each_list[[j]]) <- 1:nrow(conover_result)
    }
    names(each_list) <- names(taxa.ind[3:length(taxa.ind)])
    summary_list[[i]] <- each_list 
  }
  names(summary_list) <- names(taxa_mani)
  return(summary_list)
}

taxa.pairwise.mult.durbin <- function(taxa_mani, p_adjustment = "BH"){
  summary_list <- list() 
  for (i in 1:length(taxa_mani)){
    
    taxa.mani.ind <- taxa_mani[[i]]
    taxa.ind <- fill_missing_na_taxa(taxa.mani.ind)
    
    taxa.ind$prim.var <- as.factor(taxa.ind$prim.var)
    taxa.ind$block.id <- as.factor(taxa.ind$block.id)
    
    each_list <- list() 
    for (j in 1:(length(taxa.ind)-2)){
      
      durbin_statistic <- PMCMRplus::durbinAllPairsTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = p_adjustment)$statistic 
      
      durbin_qval <- PMCMRplus::durbinAllPairsTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = p_adjustment)$p.value 
      
      durbin_pval <- PMCMRplus::durbinAllPairsTest(y = taxa.ind[,j+2], groups = taxa.ind$prim.var, blocks = taxa.ind$block.id, p.adjust.method = "none")$p.value 
      
      ref <- c() 
      conf <- c() 
      t <- c() 
      Q.val <- c() 
      P.val <- c() 
      
      for (k in 1:length(rownames(durbin_qval))){
        ref_ind <- rep(colnames(durbin_qval)[k], length(na.omit(durbin_pval[,k])))
        ref <- c(ref, ref_ind)
        
        conf_ind <- rownames(durbin_qval)[k:length(rownames(durbin_pval))]
        conf <- c(conf, conf_ind)
        
        t_ind <- as.vector(na.omit(durbin_statistic[,k]))
        t <- c(t, t_ind)
        
        Q_ind <- decimal_adjust(na.omit(durbin_qval[,k]))
        P_ind <- decimal_adjust(na.omit(durbin_pval[,k]))
        Q.val <- c(Q.val, Q_ind)
        P.val <- c(P.val, P_ind)
      }
      durbin_result <- cbind(Ref = ref, Com = conf, t = decimal_adjust(t), P.value = P.val, Adj.P.value =  Q.val) 
      each_list[[j]] <- data.frame(durbin_result)
      rownames(each_list[[j]]) <- 1:nrow(durbin_result)
    }
    names(each_list) <- names(taxa.ind[3:length(taxa.ind)])
    summary_list[[i]] <- each_list 
  }
  names(summary_list) <- names(taxa_mani)
  return(summary_list)
}

taxa.lmm <- function(taxa_mani, sam.dat, div_num, prim_var, block_var, ref_level, method_adj = "BH" ){
  taxa.prop.ind <- taxa_mani[[div_num]]
  taxa.prop.ind[, prim_var] <- relevel(taxa.prop.ind[,prim_var], ref_level)
  taxa.prop.ind <<- taxa.prop.ind
  
  out_p <- c() 
  out <- list() 
  
  prim_var_col <<- match(prim_var, names(taxa.prop.ind))
  block_var_col <<- match(block_var, names(taxa.prop.ind))
  
  n_sam <- ncol(sam.dat)
  for (i in 1:(ncol(taxa.prop.ind)- n_sam)){
    result <- lmer(taxa.prop.ind[[i]] ~ taxa.prop.ind[,prim_var_col]  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  #taxa.prop.ind[,prim_var_col] 
    prf_sum <- coef(summary(result))[2:nrow(coef(summary(result))), 5]
    out_p <- c(out_p, prf_sum)
    
    result <- taxa_summary(coef(summary(result)), ref_level, prim_var, taxa.prop.ind)  
    out[[i]] <- result
  }
  
  ind <- list() 
  num_level <- length(table(taxa.prop.ind[, prim_var]))
  
  for (j in 1:(num_level - 1)){
    name <- paste0("taxa.prop.ind[, prim_var_col]", names(table(taxa.prop.ind[, prim_var]))[j+1])
    ind[[j]] <- which(names(out_p) == name) 
  }
  
  p_val <- list() 
  for (k in 1:(num_level - 1)){
    p_val[[k]] <- out_p[ind[[k]]]
  }
  
  q_val <- list() 
  for (f in 1:(num_level - 1)){
    q_val[[f]] <- p.adjust(p_val[[f]], method = method_adj, n = length(p_val[[f]]))
  }
  
  names(out) <- names(taxa.prop.ind[1:(length(names(taxa.prop.ind)) - n_sam)])
  
  mat <- matrix(unlist(q_val), nrow = ncol(taxa.prop.ind)- n_sam) 
  
  for (i in 1:(ncol(taxa.prop.ind) - n_sam)){
    q.val <- mat[i,]
    out[[i]] <- cbind(out[[i]], q.val)
    colnames(out[[i]]) <- c("Ref", "Com", "Est", "SE", "DF", "t", "P.value", "Adj.P.value")
  }
  return(out)
}

taxa_total <- function(taxa_mani, sam_dat, prim_var, block_var, ref_level, method_adj = "BH" ){
  total_list <- list() 
  for (i in 1:length(taxa_mani)){
    total_list[[i]] <- taxa.lmm(taxa_mani, sam_dat,  i, prim_var, block_var, ref_level, method_adj)
  }
  names(total_list) <- names(taxa_mani)
  return(total_list)
}

no_constant_acr <- function(taxa_arc, sam_dat){
  for (i in 1:length(taxa_arc)){
    taxa.arc.ind <- taxa_arc[[i]]
    ind <- c()
    
    n_sam <- ncol(sam_dat)
    for (j in 1:(ncol(taxa.arc.ind)-n_sam)){
      if (length(taxa.arc.ind[[j]][taxa.arc.ind[[j]] == 0]) == (length(taxa.arc.ind[[j]]))){
        ind <- c(ind, j)
      }
    }
    if(length(ind) !=0){
      taxa_arc[[i]] <- taxa_arc[[i]][,-ind]
    }
  }
  return(taxa_arc)
}

global_list <- function(taxa_origin, taxa_mani, div_num, prim_var, block_var, ref_level){
  taxa.prop.ind <- taxa_mani[[div_num]]
  taxa.prop.ind[, prim_var] <- relevel(taxa.prop.ind[,prim_var], ref_level)
  taxa.prop.ind <<- taxa.prop.ind
  
  taxa.ori.ind <- taxa_origin[[div_num]]
  
  prim_var_col <<- match(prim_var, names(taxa.prop.ind))
  block_var_col <<- match(block_var, names(taxa.prop.ind))
  
  mat <- matrix(NA, nrow = ncol(taxa.ori.ind), ncol = 5)
  for (i in 1:ncol(taxa.ori.ind)){
    fit_full<- lmer(taxa.prop.ind[[i]] ~ taxa.prop.ind[,prim_var_col]  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    fit_nested <- lmer(taxa.prop.ind[[i]] ~ 1  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    mat[i,] <- unlist(anova(fit_full, fit_nested))[c("logLik1", "logLik2", "Chisq2", "Df2", "Pr(>Chisq)2")]
  }
  mat <- cbind(mat, p.adjust(mat[,5], "BH"))
  mat <- data.frame(mat)
  rownames(mat) <- names(taxa.ori.ind)
  colnames(mat) <- c("logLik : Nested", "logLik : Complex", "Chisq", "DF", "P.value", "Adj.P.value")
  
  return(mat)
}

global_total_list <- function(taxa_origin, taxa_mani, prim_var, block_var, ref_level, arc = FALSE){
  total_list <- list() 
  
  if(arc){
    for (i in 1:length(taxa_mani)){
      total_list[[i]] <- global_list_arc(taxa_origin, taxa_mani, i, prim_var, block_var, ref_level)
    }
    names(total_list) <- names(taxa_mani)
  }
  else{
    for (i in 1:length(taxa_mani)){
      total_list[[i]] <- global_list(taxa_origin, taxa_mani, i, prim_var, block_var, ref_level)
    }
    names(total_list) <- names(taxa_mani)
  }
  
  return(total_list)
}

global_list_arc <- function(taxa_origin, sam_dat, taxa_mani, div_num, prim_var, block_var, ref_level){
  taxa.prop.ind <- taxa_mani[[div_num]]
  taxa.prop.ind[, prim_var] <- relevel(taxa.prop.ind[,prim_var], ref_level)
  taxa.prop.ind <<- taxa.prop.ind
  
  taxa.ori.ind <- taxa_origin[[div_num]]
  
  prim_var_col <<- match(prim_var, names(taxa.prop.ind))
  block_var_col <<- match(block_var, names(taxa.prop.ind))
  
  mat <- matrix(NA, nrow = ncol(taxa.ori.ind), ncol = 5)
  for (i in 1:(ncol(taxa.prop.ind[[i]]) - ncol(sam_dat))){
    fit_full<- lmer(taxa.prop.ind[[i]] ~ taxa.prop.ind[,prim_var_col]  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    fit_nested <- lmer(taxa.prop.ind[[i]] ~ 1  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    mat[i,] <- unlist(anova(fit_full, fit_nested))[c("logLik1", "logLik2", "Chisq2", "Df2", "Pr(>Chisq)2")]
  }
  
  mat <- cbind(mat, p.adjust(mat[,5], "BH"))
  mat <- data.frame(mat)
  rownames(mat) <- names(taxa.ori.ind)
  colnames(mat) <- c("logLik : Nested", "logLik : Complex", "Chisq", "DF", "P.value", "Adj.P.value")
  return(mat)
}

global_p_lmm <- function(taxa_mani, sam_dat, div_num, prim_var, block_var, ref_level){
  
  taxa.prop.ind <- taxa_mani[[div_num]]
  taxa.prop.ind[, prim_var] <- relevel(taxa.prop.ind[,prim_var], ref_level)
  taxa.prop.ind <<- taxa.prop.ind
  
  prim_var_col <<- match(prim_var, names(taxa.prop.ind))
  block_var_col <<- match(block_var, names(taxa.prop.ind))
  
  p_val <- c() 
  for (i in 1:(ncol(taxa.prop.ind)-ncol(sam_dat))){
    fit_full <- lmer(taxa.prop.ind[[i]] ~ taxa.prop.ind[,prim_var_col]  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    fit_nested <- lmer(taxa.prop.ind[[i]] ~ 1  + (1|taxa.prop.ind[,block_var_col]), data = taxa.prop.ind)  
    p_val[i] <- unlist(anova(fit_full, fit_nested))["Pr(>Chisq)2"]
  }
  
  return(p_val)
}

global_total_p <- function(taxa_mani, sam_dat, prim_var, block_var, ref_level){
  p_val <- list() 
  for (i in 1:length(taxa_mani)){
    p_val[[i]] <- decimal_adjust(p.adjust(global_p_lmm(taxa_mani, sam_dat, i, prim_var, block_var, ref_level), "BH"))
  }
  names(p_val) <- names(taxa_mani)
  return(p_val)
}

make_p_val_list <- function(result_table){
  list_total <- list() 
  for (i in 1:length(result_table)){
    list_total_ind <- list()
    for (j in 1:nrow(result_table[[i]][[1]])){
      vec <- NULL
      for (k in 1:length(result_table[[i]])){
        vec <- c(vec, result_table[[i]][[k]]$Adj.P.value[j])
      }
      list_total_ind[[j]] <- vec 
    }
    list_total[[i]] <- list_total_ind
  }
  return(list_total)
}

make_q_val_list <- function(p_val_list){
  q_val_result <- list() 
  for (p in 1:length(p_val_list)){
    q_val_result_ind <- list()
    for (k in 1:length(p_val_list[[p]])){
      q_val_result_ind[[k]] <- p.adjust(p_val_list[[p]][[k]], "BH")
    }
    q_val_result[[p]] <- q_val_result_ind
  }
  return(q_val_result)
}

q_val_dat <- function(q_val_list, result_table, p.val = TRUE){
  q_val_result_dat <- list()
  
  for (i in 1:length(q_val_list)){
    mat <- matrix(unlist(q_val_list[[i]]), ncol = length(q_val_list[[i]]), nrow = length(q_val_list[[i]][[1]]))
    q_val_result_dat[[i]] <- mat 
  }
  
  if(p.val){
    for (k in 1:length(result_table)){
      for (p in 1:length(result_table[[k]])){
        
        result_table[[k]][[p]]$P.value <- p.value.0.1_char(result_table[[k]][[p]]$P.value)
        result_table[[k]][[p]]$Adj.P.value <- p.value.0.1_char(decimal_adjust(q_val_result_dat[[k]][p,]))  
      }
    }
  }else{
    for (k in 1:length(result_table)){
      for (p in 1:length(result_table[[k]])){
        
        result_table[[k]][[p]]$Adj.P.value <- p.value.0.1_char(decimal_adjust(q_val_result_dat[[k]][p,]))  
      }
    }
  }
  
  return(result_table)
}

q_val_dat_fc <- function(q_val_list, result_table){
  q_val_result_dat <- list()
  
  for (i in 1:length(q_val_list)){
    mat <- matrix(unlist(q_val_list[[i]]), ncol = length(q_val_list[[i]]), nrow = length(q_val_list[[i]][[1]]))
    q_val_result_dat[[i]] <- mat 
  }
  
  for (k in 1:length(result_table)){
    for (p in 1:length(result_table[[k]])){
      result_table[[k]][[p]][, "Adj.P.value"] <- q_val_result_dat[[k]][p,] 
    }
  }
  return(result_table)
}


q_val_dat_original <- function(q_val_list, result_table){
  q_val_result_dat <- list()
  
  for (i in 1:length(q_val_list)){
    mat <- matrix(unlist(q_val_list[[i]]), ncol = length(q_val_list[[i]]), nrow = length(q_val_list[[i]][[1]]))
    q_val_result_dat[[i]] <- mat 
  }
  
  for (k in 1:length(result_table)){
    for (p in 1:length(result_table[[k]])){
      result_table[[k]][[p]]$Adj.P.value <- as.numeric(q_val_result_dat[[k]][p,])  #07111 수정 
    }
  }
  return(result_table)
}

q_val_convert <- function(pairwise_result){
  
  for (i in 1:length(pairwise_result)){
    q_val_ind <- list() 
    
    for (j in 1:length(pairwise_result[[i]])){
      pairwise_result[[i]][[j]]$P.value <- p.value.0.1_char(decimal_adjust(pairwise_result[[i]][[j]]$P.value))
      pairwise_result[[i]][[j]]$Adj.P.value <- p.value.0.1_char(decimal_adjust(p.adjust(pairwise_result[[i]][[j]]$Adj.P.value, "BH"))) #0718 수정
    }
  }
  return(pairwise_result)
}

q_val_convert_qc <- function(pairwise_result){
  
  for (i in 1:length(pairwise_result)){
    q_val_ind <- list() 
    
    for (j in 1:length(pairwise_result[[i]])){
      pairwise_result[[i]][[j]]$Adj.P.value <-p.adjust(pairwise_result[[i]][[j]]$Adj.P.value, "BH")
    }
  }
  return(pairwise_result)
}

#####################
# Volcano Pairwise  #
#####################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# library(DESeq2)
# 
# install.packages("GUniFrac")
# library(GUniFrac)

log2_fold_change <- function(taxa_mani){
  final_list <- list() 
  for (k in 1:length(taxa_mani)){
    
    dim <- k
    otu_tab <- t(taxa_mani[[dim]][,3:length(taxa_mani[[dim]])])
    time.p.cat <- names(table(taxa_mani[[dim]][,"prim.var"]))
    combination <- combn(time.p.cat, 2)
    
    total_list <- list() 
    
    for (i in 1:ncol(combination)){
      otu_tab <- t(taxa_mani[[dim]][,3:length(taxa_mani[[dim]])])
      otu_tab <- otu_tab[,taxa_mani[[dim]][,"prim.var"] %in% combination[,i]]
      meta <- data.frame(group = factor(taxa_mani[[dim]][,"prim.var"][taxa_mani[[dim]][,"prim.var"] %in% combination[,i]], levels = sort(combination[,i])), row.names = colnames(otu_tab))
      
      #mat <- matrix(NA, nrow = nrow(otu_tab), ncol = ncol(otu_tab))
      #for (j in 1:ncol(otu_tab)){
      #  mat[,j] <- (otu_tab[,j] / estimateSizeFactorsForMatrix(otu_tab)[j])
      #}
      
      normalized <- apply(otu_tab, 1, function(x){log2(mean(x[meta$group == combination[2, i]])/mean(x[meta$group == combination[1, i]]))}) #otu_tab <- mat 
      normalized <- replace(normalized, is.infinite(normalized), NA)
      total_list[[i]] <- normalized
    }
    names(total_list) <- apply(combination, 2, function(x)paste0(x[1], "-", x[2]))
    final_list[[k]] <- total_list
  }
  names(final_list) <- names(taxa_mani)
  return(final_list)
}
# log2_result <- log2_fold_change(taxa.mani)
# sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)

#taxa.mani.lmm 말고 그냥 taxa.mani 
log2_fold_change_lmm <- function(taxa_mani, ref_level){
  final_list <- list() 
  for (k in 1:length(taxa_mani)){
    
    dim <- k
    otu_tab <- t(taxa_mani[[dim]][,3:length(taxa_mani[[dim]])])
    time.p.cat <- names(table(taxa_mani[[dim]][,"prim.var"]))
    combination <- combn(time.p.cat, 2)
    ref_exist <- apply(combn(names(table(taxa_mani[[1]][,"prim.var"])), 2), 2,  function(x) ref_level %in% x)
    combination <- combination[,ref_exist]
    
    total_list <- list() 
    
    for (i in 1:ncol(combination)){
      otu_tab <- t(taxa_mani[[dim]][,3:length(taxa_mani[[dim]])])
      otu_tab <- otu_tab[,taxa_mani[[dim]][,"prim.var"] %in% combination[,i]]
      meta <- data.frame(group = factor(taxa_mani[[dim]][,"prim.var"][taxa_mani[[dim]][,"prim.var"] %in% combination[,i]], levels = sort(combination[,i])), row.names = colnames(otu_tab))
      
      normalized <- apply(otu_tab, 1, function(x){log2(mean(x[meta$group == combination[2, i]])/mean(x[meta$group == combination[1, i]]))}) #otu_tab <- mat 
      normalized <- replace(normalized, is.infinite(normalized), NA)
      total_list[[i]] <- normalized
    }
    names(total_list) <- apply(combination, 2, function(x)paste0(x[1], "-", x[2]))
    final_list[[k]] <- total_list
  }
  names(final_list) <- names(taxa_mani)
  return(final_list)
}


size_convert <- function(factor){
  if (factor == TRUE){
    size <- "glob: sig"
  }else{
    size <- "glob: non-sig"
  }
  return(size)
}

volcano_dat <- function(pairwise_result, log2fc_result, sig_global){
  total_list <- list() 
  
  for (i in 1:length(pairwise_result)){
    list <- list() 
    for (k in 1:length(pairwise_result[[i]])){
      mat <- matrix(NA, nrow(pairwise_result[[i]][[1]]), 9)
      colnames(mat) <- c("Ref", "Com", "Ref-Com", "Adj.P.value", "logFC", "global_sig", "pairwise_sig", "rank", "name")
      mat[,1] <- pairwise_result[[i]][[k]][,"Ref"]
      mat[,2] <- pairwise_result[[i]][[k]][,"Com"]
      mat[,3] <- paste0(pairwise_result[[i]][[k]][,"Ref"], "-", pairwise_result[[i]][[k]][,"Com"])
      mat[,4] <- as.numeric(pairwise_result[[i]][[k]][,"Adj.P.value"])   
      mat[,5] <- as.numeric(matrix(unlist(log2fc_result[[i]]), nrow = nrow(pairwise_result[[1]][[1]]), byrow = TRUE)[,k])
      mat[,6] <- size_convert(sig_global[[i]][k]) 
      mat[,7] <- ifelse(as.numeric(pairwise_result[[i]][[k]][,"Adj.P.value"]) <= 0.05, "significant", "not significant")
      mat[,8] <- str_to_title(names(pairwise_result)[i])
      mat[,9] <- names(pairwise_result[[i]])[k]
      list[[k]] <- data.frame(mat)
    }
    total_list[[i]] <- list
    names(total_list[[i]]) <- names(pairwise_result[[i]])
  }
  names(total_list) <- names(pairwise_result)
  return(total_list)
}



#volcanodat <- volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa)

volcano_dat_lmm <- function(pairwise_result, log2fc_result, sig_global){
  total_list <- list() 
  
  for (i in 1:length(pairwise_result)){
    list <- list() 
    for (k in 1:length(pairwise_result[[i]])){
      mat <- matrix(NA, nrow(pairwise_result[[i]][[1]]), 7)
      colnames(mat) <- c("Ref-Com", "Adj.P.value", "logFC", "global_sig", "pairwise_sig", "rank", "name")
      mat[,1] <- paste0(pairwise_result[[i]][[k]]$Ref, "-", pairwise_result[[i]][[k]]$Com)
      mat[,2] <- as.numeric(pairwise_result[[i]][[k]][,"Adj.P.value"])   
      mat[,3] <- as.numeric(matrix(unlist(log2fc_result[[i]]), nrow = nrow(pairwise_result[[1]][[1]]), byrow = TRUE)[,k])
      mat[,4] <- size_convert(sig_global[[i]][k]) 
      mat[,5] <- ifelse(as.numeric(pairwise_result[[i]][[k]][,"Adj.P.value"]) <= 0.05, "significant", "not significant")
      mat[,6] <- str_to_title(names(pairwise_result)[i])
      mat[,7] <- names(pairwise_result[[i]])[k]   
      list[[k]] <- data.frame(mat)
    }
    total_list[[i]] <- list
    names(total_list[[i]]) <- names(pairwise_result[[i]])
  }
  names(total_list) <- names(pairwise_result)
  return(total_list)
}


volcano_list_to_dat <- function(volcano_dat, lmm = FALSE, species = TRUE){
  
  list_pairwise <- list() 
  if(species){
    ncomb <- nrow(volcano_dat[[1]][[1]])
    length <- sum(as.numeric(sapply(volcano_dat, function(y)length(y))))
    for (k in 1:ncomb){
      phylum <- do.call(rbind, (lapply(volcano_dat[[1]], function(x) x[k,])))
      class <- do.call(rbind, (lapply(volcano_dat[[2]], function(x) x[k,])))
      order <- do.call(rbind, (lapply(volcano_dat[[3]], function(x) x[k,])))
      family <- do.call(rbind, (lapply(volcano_dat[[4]], function(x) x[k,])))
      genus <- do.call(rbind, (lapply(volcano_dat[[5]], function(x) x[k,])))
      species <- do.call(rbind, (lapply(volcano_dat[[6]], function(x) x[k,])))
      
      list_pairwise[[k]] <- rbind(phylum, class, order, family, genus, species)
    }
    if (lmm){
      names(list_pairwise) <- volcano_dat[[1]][[1]]$Ref.Com
    }else{
      names(list_pairwise) <-  paste0(volcano_dat[[1]][[1]]$Ref, "-", volcano_dat[[1]][[1]]$Com)
    }
  }else{
    
    ncomb <- nrow(volcano_dat[[1]][[1]])
    length <- sum(as.numeric(sapply(volcano_dat, function(y)length(y))))
    for (k in 1:ncomb){
      phylum <- do.call(rbind, (lapply(volcano_dat[[1]], function(x) x[k,])))
      class <- do.call(rbind, (lapply(volcano_dat[[2]], function(x) x[k,])))
      order <- do.call(rbind, (lapply(volcano_dat[[3]], function(x) x[k,])))
      family <- do.call(rbind, (lapply(volcano_dat[[4]], function(x) x[k,])))
      genus <- do.call(rbind, (lapply(volcano_dat[[5]], function(x) x[k,])))
      
      list_pairwise[[k]] <- rbind(phylum, class, order, family, genus)
    }
    if (lmm){
      names(list_pairwise) <- paste0(volcano_dat[[1]][[1]]$Ref.Com)
    }else{
      names(list_pairwise) <-  paste0(volcano_dat[[1]][[1]]$Ref, "-", volcano_dat[[1]][[1]]$Com)
    }
  }
  
  return(list_pairwise)
}

#showgrid 

###############
# Volcano 3D  #
###############

make_tripe_pvalues <- function(global_test, pairwise_result, taxa_mani, species = TRUE){
  
  q_list <- list() 
  if(species == TRUE){
    for (i in 1:length(pairwise_result)){
      mat <- matrix(data = NA, nrow = length(as.numeric(sapply(pairwise_result[[i]], `[`, 1, "Adj.P.value"))), ncol = nrow(pairwise_result[[i]][[1]])+1)
      
      mat[,1] <- as.numeric(global_test[[i]][,"Adj.P.value"])
      for (j in 1:nrow(pairwise_result[[i]][[1]])){
        mat[,j+1] <- as.numeric(sapply(pairwise_result[[i]],  `[`, j, "Adj.P.value"))
      }
      
      colnames(mat) <- c("global", "pairwise_1", "pairwise_2", "pairwise_3")
      rownames(mat) <- names(taxa_mani[[i]])[3:length(names(taxa_mani[[i]]))]
      q_list[[i]] <- mat 
    }
  }else{
    for (i in 1:5){
      mat <- matrix(data = NA, nrow = length(as.numeric(sapply(pairwise_result[[i]], `[`, 1, "Adj.P.value"))), ncol = nrow(pairwise_result[[i]][[1]])+1)
      
      mat[,1] <- as.numeric(global_test[[i]][,"Adj.P.value"])
      for (j in 1:nrow(pairwise_result[[i]][[1]])){
        mat[,j+1] <- as.numeric(sapply(pairwise_result[[i]],  `[`, j, "Adj.P.value"))
      }
      
      colnames(mat) <- c("global", "pairwise_1", "pairwise_2", "pairwise_3")
      rownames(mat) <- names(taxa_mani[[i]])[3:length(names(taxa_mani[[i]]))]
      q_list[[i]] <- mat 
    }
  }
  
  names(q_list) <- names(pairwise_result)
  
  return(q_list)
}

volcano_visualize_individual <- function(volcano_dat_fin, species){
  fig_list <- list()
  annotation_list <- list()
  title_vec <- names(volcano_dat_fin)

  q_val_list <- lapply(volcano_dat_fin, "[", , "Adj.P.value")
  ylim <- min(as.numeric(sapply(q_val_list, function(x)min(na.omit(x)))))
  ylim_val <- ceiling(-log10(ylim + 0.0001)) + 1

  volcano_dat_fin[[1]]$rank <- factor(volcano_dat_fin[[1]]$rank, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
  volcano_dat_fin[[1]]$global_sig <- factor(volcano_dat_fin[[1]]$global_sig, levels = c("glob: sig", "glob: non-sig"))

  x_lim_1 <- max(max(na.omit(as.numeric(volcano_dat_fin[[1]]$logFC))), abs(min(na.omit(as.numeric(volcano_dat_fin[[1]]$logFC)))))

  fig_list[[1]] <- plot_ly(data = volcano_dat_fin[[1]]) %>% add_trace(data = volcano_dat_fin[[1]][volcano_dat_fin[[1]]$global_sig == "glob: non-sig",], x = ~as.numeric(logFC), y = ~-log10(as.numeric(Adj.P.value) + 0.0001), symbols = 'cross', colors = "Set1", color = ~rank, type = "scatter",  mode = "markers", marker = list(size=7.5), text = ~paste("<b>logFC:</b>", as.numeric(logFC), "<b><br>Adj. P.value:</b>", as.numeric(Adj.P.value) , "<br>", rank, ":", name), hoverinfo = 'text') %>% add_trace(data = volcano_dat_fin[[1]][volcano_dat_fin[[1]]$global_sig == "glob: sig",], x = ~as.numeric(logFC), y = ~-log10(as.numeric(Adj.P.value) + 0.0001), symbol = factor("glob: sig", levels = c("glob: sig", "glob: non-sig")), symbols = 27, colors = "Set1", color = ~rank, type = "scatter",  mode = "markers", marker = list(size = 7.5), text = ~paste("<b>logFC:</b>", as.numeric(logFC), "<b><br>Adj. P.value:</b>", as.numeric(Adj.P.value) , "<br>", rank, ":", name), hoverinfo = 'text', showlegend = FALSE)  %>% layout(xaxis = list(title = paste0("<b>", title_vec[1], "</b>"), ticktext = list("<b>log<sub>2</sub>FC</b>"), tickvals = list(0)), yaxis = list(title =  list(text = "<b>-log<sub>10</sub>Q</b>", font = list(size = 13))), plot_bgcolor = "#FFFFFF", shapes = list(list(type = "rect", fillcolor = "Blue", line = list(color = "Blue"), opacity = 0, y0 = 1.301029996, y1 = ylim_val, x0 = - x_lim_1 - 1, x1= -2), list(type = "rect", fillcolor = "red", line = list(color = "red"), opacity = 0, y0 = 1.301029996, y1 = ylim_val, x0 = 2, x1 = x_lim_1 + 1)), showgrid = FALSE)

  for (i in 2:(length(volcano_dat_fin))){
    volcano_dat_fin[[i]]$rank <- factor(volcano_dat_fin[[i]]$rank, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
    volcano_dat_fin[[i]]$global_sig <- factor(volcano_dat_fin[[i]]$global_sig, levels = c("glob: sig", "glob: non-sig"))
    x_lim <- max(max(na.omit(as.numeric(volcano_dat_fin[[i]]$logFC))), abs(min(na.omit(as.numeric(volcano_dat_fin[[i]]$logFC)))))
    fig_list[[i]] <- plot_ly(volcano_dat_fin[[i]], x = ~as.numeric(logFC), y = ~-log10(as.numeric(Adj.P.value) + 0.0001), symbol = ~ global_sig, symbols = c('cross', 27), colors = "Set1", color = ~rank,  type = "scatter",  mode = "markers", marker = list(size=7.5), text = ~paste("<b>logFC:</b>", as.numeric(logFC), "<b><br>Adj. P.value:</b>", as.numeric(Adj.P.value) , "<br>", rank, ":", name), showlegend = FALSE, hoverinfo = 'text') %>% layout(xaxis = list(title = paste0("<b>", title_vec[i], "</b>"), ticktext = list("<b>log<sub>2</sub>FC</b>"), tickvals = list(0)), yaxis = list(title =  "-log<sub>10</sub>Q", font = list(size = 3)), plot_bgcolor = "#FFFFFF", shapes = list(list(type = "rect", fillcolor = "Blue",line = list(color = "Blue"), opacity = 0, y0 = 1.301029996, y1 = ylim_val,  x0 = -x_lim-1, x1= -2), list(type = "rect", fillcolor = "red", line = list(color = "red"), opacity = 0, y0 = 1.301029996, y1 = ylim_val, x0 = 2, x1 = x_lim + 1)))
  }

  number_row <- ceiling(length(volcano_dat_fin)/3)

  fig <- subplot(fig_list,titleX = TRUE, shareY = TRUE, nrows = number_row,  margin = c(0.02, 0.02, 0.02, 0.07)) %>% layout(showlegend = TRUE)
  fig
}

triple_list_to_dat <- function(p_value_result){
  dat <- p_value_result[[1]]
  
  for (i in 2:length(p_value_result)){
    dat <- rbind(dat, p_value_result[[i]])
  }
  dat <- data.frame(dat)
  return(dat)
}

rename_volcano_3d <- function(outcome){
  abbrev <- if (length(labs) == 3) labs else abbreviate(levels(outcome), 1)
  
  labs <- c("Not significant",
            paste0(abbrev[2], "+"),
            paste0(abbrev[2], "+", abbrev[3], "+"),
            paste0(abbrev[3], "+"),
            paste0(abbrev[1], "+", abbrev[3], "+"),
            paste0(abbrev[1], "+"),
            paste0(abbrev[1], "+", abbrev[2], "+"))
  return(labs)
}



volcano3D_RE <- function (polar, type = 1, label_rows = c(), label_size = 14,
                          arrow_length = 100, colour_code_labels = FALSE, label_colour = "black",
                          grid_colour = "grey80", grid_width = 2, grid_options = NULL,
                          axis_colour = "black", axis_width = 2, marker_size = 3, marker_outline_width = 0,
                          marker_outline_colour = "white", z_axis_title_offset = 1.2,
                          z_axis_title_size = 12, z_axis_angle = 0.5, radial_axis_title_size = 14,
                          radial_axis_title_offset = 1.2, xy_aspectratio = 1, z_aspectratio = 0.8,
                          camera_eye = list(x = 0.9, y = 0.9, z = 0.9), ...)
{

  options(scipen = 999)

  if (is(polar, "polar")) {
    args <- as.list(match.call())[-1]
    return(do.call(volcano3D_v1, args))
  }

  if (!is(polar, "volc3d"))
   stop("Not a 'volc3d' class object")

  args <- list(r_vector = polar@df[[type]]$r, z_vector = polar@df[[type]]$z)
  args <- append(args, grid_options)

  grid <- do.call(polar_grid, args)

  polar_grid <- grid@polar_grid
  axes <- grid@axes
  axis_labels <- grid@axis_labs
  h <- grid@z
  R <- grid@r
  max_offset <- max(c(z_axis_title_offset, radial_axis_title_offset))
  xyrange <- c(-1.05 * (max_offset) * R, 1.05 * max_offset *
                 R)
  axis_settings <- list(title = "", zeroline = FALSE, showline = FALSE,
                        showticklabels = FALSE, showgrid = FALSE, autotick = FALSE,
                        showspikes = FALSE)
  axis_settings_xy <- list(title = "", zeroline = FALSE, showline = FALSE,
                           showticklabels = FALSE, showgrid = FALSE, autotick = FALSE,
                           showspikes = FALSE, range = xyrange)
  df <- polar@df[[type]]
  if (length(label_rows) != 0) {
    if (!all(is.numeric(label_rows))) {
      if (!all(label_rows %in% rownames(df))) {
        stop("label_rows must be in rownames(polar_df)")
      }
    }
    if (all(is.numeric(label_rows))) {
      if (!all(label_rows < nrow(df))) {
        stop("label_rows not in 1:nrow(polar_df)")
      }
    }
    annot <- lapply(label_rows, function(i) {
      row <- df[i, ]
      if (colour_code_labels)
        ac <- row$col
      else ac <- label_colour
      annot <- list(x = row$x, y = row$y, z = row$z, text = rownames(row),
                    textangle = 0, ax = arrow_length, ay = 0, arrowcolor = ac,
                    font = list(color = ac), arrowwidth = 1, arrowhead = 0,
                    arrowsize = 1.5, yanchor = "middle")
    })
  }
  else {
    annot <- list()
  }
  plot_ly(polar@df[[type]], x = ~x, y = ~y, z = ~z, color = ~lab,
          colors = polar@scheme, hoverinfo = "text", text = ~paste0(rownames(polar@df[[type]]), "<br>Adj. P.value = ", format(round(pvalue, 3), scientific = FALSE)), key = rownames(polar@df[[type]]),
          marker = list(size = marker_size, line = list(color = marker_outline_colour,
                                                        width = marker_outline_width)), type = "scatter3d",
          mode = "markers", ...) %>% add_trace(x = polar_grid$x,
                                               y = polar_grid$y, z = polar_grid$z, color = I(grid_colour),
                                               line = list(width = grid_width), showlegend = FALSE,
                                               type = "scatter3d", mode = "lines", hoverinfo = "none",
                                               inherit = FALSE) %>% add_trace(x = axes$x, y = axes$y,
                                                                              z = axes$z, color = I(axis_colour), line = list(width = axis_width),
                                                                              showlegend = FALSE, type = "scatter3d", mode = "lines",
                                                                              hoverinfo = "none", inherit = FALSE) %>% add_text(x = radial_axis_title_offset *
                                                                                                                                  axis_labels$x, y = radial_axis_title_offset * axis_labels$y,
                                                                                                                                z = 0, text = levels(polar@outcome), color = I(axis_colour),
                                                                                                                                type = "scatter3d", mode = "text", textfont = list(size = radial_axis_title_size),
                                                                                                                                textposition = "middle center", hoverinfo = "none", showlegend = FALSE,
                                                                                                                                inherit = FALSE) %>% add_text(x = c(rep(R * sinpi(z_axis_angle),
                                                                                                                                                                        grid@n_z_breaks), R * z_axis_title_offset * sinpi(z_axis_angle)),
                                                                                                                                                              y = c(rep(R * cospi(z_axis_angle), grid@n_z_breaks),
                                                                                                                                                                    R * z_axis_title_offset * cospi(z_axis_angle)), z = c(grid@z_breaks,
                                                                                                                                                                                                                          h/2), text = c(grid@z_breaks, "-log<sub>10</sub>Q"),
                                                                                                                                                              textposition = "middle left", textfont = list(size = z_axis_title_size),
                                                                                                                                                              color = I(axis_colour), hoverinfo = "none", showlegend = FALSE,
                                                                                                                                                              inherit = FALSE) %>% add_text(x = grid@text_coords$x,
                                                                                                                                                                                            y = grid@text_coords$y, z = h * 0.03, text = grid@text_coords$text,
                                                                                                                                                                                            textposition = "middle center", textfont = list(size = 10),
                                                                                                                                                                                            color = I(axis_colour), hoverinfo = "none", showlegend = FALSE,
                                                                                                                                                                                            inherit = FALSE) %>% layout(margin = list(0, 0, 0, 0),
                                                                                                                                                                                                                        paper_bgcolor = "rgba(0, 0, 0, 0)", plot_bgcolor = "rgba(0, 0, 0, 0)",
                                                                                                                                                                                                                        scene = list(camera = list(eye = camera_eye), aspectratio = list(x = xy_aspectratio,
                                                                                                                                                                                                                                                                                         y = xy_aspectratio, z = z_aspectratio), dragmode = "turntable",
                                                                                                                                                                                                                                     xaxis = axis_settings_xy, yaxis = axis_settings_xy,
                                                                                                                                                                                                                                     zaxis = axis_settings, annotations = annot), xaxis = list(title = "x"),
                                                                                                                                                                                                                        yaxis = list(title = "y"))
}
