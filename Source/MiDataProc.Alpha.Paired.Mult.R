#####################
# Data Manipulation #
#####################

block_ids <- function(sam_dat){
  block_vec <- c() 
  for (i in 1:length(colnames(sam_dat))){
    length_element <- length(unique(as.numeric(table(sam_dat[,colnames(sam_dat)[i]])))) 
    if (length_element == 1){
      block_vec <- c(block_vec, colnames(sam_dat)[i])
    }
  }
  return(block_vec)
}

decimal_adjust <- function(x){
  x <- as.numeric(x)
  round_num <- round(x, 3)
  decimal <- formatC(round_num, 3, format = "f")
  return(decimal)
}

bin.cat.recode.func.mult <- function(sam.dat, sel.bin.var, ori.cat, rename){
  ind <- list()
  for (i in 1:length(table(sam.dat[,sel.bin.var]))){
    ind <- which(sam.dat[,sel.bin.var] == ori.cat[i])
    sam.dat[ind, sel.bin.var] <- rename[i]
  }
  return(sam.dat)
}

alpha.mult.pair.cat.ref.func <- function(sel.bin.var, rename.var, sam.dat, alpha.div){
  bin.var <- unlist(sam.dat[,sel.bin.var], use.names = FALSE)

  bin.var.vec <- c()
  alpha_div <- c() 
  ind <- list()
  for (i in 1:length(rename.var)){
    ind[[i]] <- which(bin.var == rename.var[i])

    bin.var.vec <- c(bin.var.vec, bin.var[ind[[i]]])
    alpha_div <- rbind(alpha.div, alpha.div[ind[[i]],])
  }
  
  bin.var.vec <- factor(bin.var.vec)
  return(list(bin.var = bin.var.vec, alpha.div = alpha_div))
}

alpha.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, alpha.div) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  return(list(bin.var = bin.var, alpha.div = alpha.div))
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}


p.value.0.1_char <- function(x){
  ind.0 <- which(x == "0.000")
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000")
  x[ind.1] <- ">.999"
  return(x)
}

##########################################################
# Comparative analysis for multi-paired alpha diversity  #
##########################################################
alpha.t.paired.mult <- function(alpha_mani, div_num, method_adj = "BH"){
  
  total_list <- list() 
  if (length(unique(as.numeric(table(alpha_mani$prime_var)))) != 1){
    ind <- as.numeric(table(alpha_mani$block_var)) != length(table(alpha_mani$prime_var))     
    alpha_mani <- alpha_mani[!(alpha_mani$block_var %in% names(table(alpha_mani$block_var))[ind]),] 
  }
  
  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  
  ind_list <- list()
  for (i in 1:length(time.p.cat)){
    ind_list[[i]] <- which(alpha_mani$prime_var == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, ncol(comb), 8)
  out_2 <- matrix(NA, ncol(comb), 8)
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    fit <- t.test(alpha_mani[,div_num+2][ind_list[[combination[1]]]], alpha_mani[,div_num+2][ind_list[[combination[2]]]], paired = TRUE)
    out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], decimal_adjust(fit$statistic), decimal_adjust(fit$stderr), decimal_adjust(fit$parameter), decimal_adjust(fit$conf.int), fit$p.value)
    out_2[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], fit$statistic, fit$stderr, fit$parameter, fit$conf.int, fit$p.value)
  }
  
  q_val <- p.adjust(out[,8], method = method_adj)
  q_val_2 <- p.adjust(out[,8], method = method_adj, n = length(out[,8]))
  
  out <- as.data.frame(cbind(out[,1:7], decimal_adjust(as.numeric(out[,8])), p.value.0.1_char(decimal_adjust(q_val)))) 
  colnames(out) <- c("Ref", "Com", "t", "Std Err", "Df", "Lower", "Upper", "P.value", "Q.value")
  
  out_2 <- as.data.frame(cbind(out_2[,1:7], as.numeric(out_2[,8]), q_val_2)) 
  colnames(out_2) <- c("Ref", "Com", "t", "Std Err", "Df", "Lower", "Upper", "P.value", "Q.value")
  
  total_list$download <- out_2 
  total_list$table <- out 
  
  return(total_list)
}

alpha.t.paired.mult.united <- function(alpha_div, alpha_mani, method_adj = "BH"){
  paired_t_list <- list() 
  
  table_list <- list() 
  download_list <- list() 
  
  for (i in 1:ncol(alpha_div)){
    div_list<- alpha.t.paired.mult(alpha_mani, i, method_adj)
    table_list[[i]] <- div_list$table
    download_list[[i]] <- div_list$download
  }
  
  names(table_list) <- names(alpha_div)
  names(download_list) <- names(alpha_div)
  
  paired_t_list$table <- table_list
  paired_t_list$download <- download_list
  return(paired_t_list)
}

alpha.wilcox.paired.mult <- function(alpha_mani, div_num, method_adj = "BH"){
  
  total_list <- list()
  if(length(unique(as.numeric(table(alpha_mani$prime_var)))) != 1){
    ind <- as.numeric(table(alpha_mani$block_var)) != length(table(alpha_mani$prime_var))
    alpha_mani <- alpha_mani[!(alpha_mani$block_var %in% names(table(alpha_mani$block_var))[ind]),]
  }
  
  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  ind_list <- list()
  for (i in 1:length(time.p.cat)){
    ind_list[[i]] <- which(alpha_mani[,"prime_var"] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, ncol(comb), 4)
  out_2 <- matrix(NA, ncol(comb), 4)
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    fit <- wilcox.test(alpha_mani[,div_num+2][ind_list[[combination[1]]]], alpha_mani[,div_num+2][ind_list[[combination[2]]]], paired = TRUE, correct = FALSE)
    out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], decimal_adjust(fit$statistic), as.numeric(fit$p.value))
    out_2[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], fit$statistic, as.numeric(fit$p.value))
  }
  
  q_val <- p.adjust(out[,4], method = method_adj, n = length(out[,4]))
  q_val_2 <- p.adjust(out[,4], method = method_adj, n = length(out[,4]))
  
  out <- as.data.frame(cbind(out[,1:3], decimal_adjust(out[,4]), p.value.0.1_char(decimal_adjust(q_val))))
  out_2 <- as.data.frame(cbind(out_2[,1:3], out_2[,4], q_val_2)) 
  
  colnames(out) <- c("Ref", "Com", "W", "P.value", "Q.value")
  colnames(out_2) <- c("Ref", "Com", "W", "P.value", "Q.value")
  
  total_list$download <- out_2 
  total_list$table <- out 
  
  return(total_list)
}

alpha.wilcox.paired.mult.united <- function(alpha_div, alpha_mani, method_adj = "BH"){
  wilcox_p_list <- list() 
  
  table_list <- list() 
  download_list <- list()
  
  for (i in 1:ncol(alpha_div)){
    div_list <- alpha.wilcox.paired.mult(alpha_mani, i, method_adj)
    table_list[[i]] <- div_list$table 
    download_list[[i]] <- div_list$download 
  }
  
  names(table_list) <- names(alpha_div)
  names(download_list) <- names(alpha_div)
  
  wilcox_p_list$table <- table_list 
  wilcox_p_list$download <- download_list 
  
  return(wilcox_p_list)
}

alpha.f.pair.mult.overall <- function(alpha_mani){
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  summary_vec <- c()
  
  for (i in 1:(length(alpha_mani)-2)){
    summary_aov <- unlist(summary(aov(alpha_mani[, i+2] ~ prime_var + block_var, data = alpha_mani)))
    summary_ind <- c(decimal_adjust(summary_aov[["F value2"]]), p.value.0.1_char(decimal_adjust(summary_aov[["Pr(>F)2"]])))  #revise 
    summary_vec <- rbind(summary_vec, summary_ind)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[3:11]
  
  colnames(summary_vec) = c("F", "P.value")
  return(summary_vec)
}

alpha.f.pair.mult.overall_download <- function(alpha_mani){
  summary_vec <- c()
  
  for (i in 1:(length(alpha_mani)-2)){
    summary_aov <- unlist(summary(aov(alpha_mani[, i+2] ~ prime_var + block_var, data = alpha_mani)))
    summary_ind <- c(summary_aov[["F value2"]], summary_aov[["Pr(>F)2"]]) #revise 
    summary_vec <- rbind(summary_vec, summary_ind)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[3:11]
  colnames(summary_vec) = c("F", "P.value")
  
  return(summary_vec)
}

alpha.f.pair.mult <- function(alpha_mani){
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  p.val <- c() 
  for (i in 1:(length(alpha_mani) - 2)){
    summary <- summary(aov(alpha_mani[, i + 2] ~ prime_var + block_var, data = alpha_mani))
    p.val.ind <- decimal_adjust(unlist(summary)[["Pr(>F)2"]])
    p.val <- c(p.val , p.val.ind)
  }
  return(p.val)
}

alpha.f.pair.mult.tukey <- function(alpha_mani){
  div_list <- list() 
  
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  
  for (i in 1:(length(alpha_mani)-2)){
    
    rmanova_result <- aov(alpha_mani[, i+2] ~ prime_var + block_var, data = alpha_mani)
    dat <- as.data.frame(TukeyHSD(rmanova_result)$prime_var)

    tukey_result <- data.frame(lapply(dat, decimal_adjust), row.names = rownames(dat)) 
    ind_result <- c() 

    for (j in 1:nrow(tukey_result)){
      ref_conf_list <- strsplit(rownames(TukeyHSD(rmanova_result)$prime_var), split = "-")[[j]]
      ref <- ref_conf_list[2]
      conf <- ref_conf_list[1]
      ref_mean <- decimal_adjust(mean(alpha_mani[,i+2][alpha_mani$prime_var == ref]))
      conf_mean <- decimal_adjust(mean(alpha_mani[,i+2][alpha_mani$prime_var == conf]))
      tukey_result_ind <- c(ref, conf, ref_mean, conf_mean, tukey_result[j,1:3], p.value.0.1_char(decimal_adjust(tukey_result[j, 4])))
      ind_result <- rbind(ind_result, tukey_result_ind)
    }
    
    rownames(ind_result) <- 1:nrow(ind_result)
    colnames(ind_result) <- c("Ref", "Com", "Mean (Ref) ", "Mean (Com)", "Diff", "Lower", "Upper", "Adj. P.value")
    div_list[[i]] <- as.data.frame(ind_result)
  }
  
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

alpha.f.pair.mult.tukey_download <- function(alpha_mani){
  div_list <- list() 
  for (i in 1:(length(alpha_mani)-2)){
    
    alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
    alpha_mani$block_var <- as.factor(alpha_mani$block_var)
    
    rmanova_result <- aov(alpha_mani[, i+2] ~ prime_var + block_var, data = alpha_mani)
    dat <- as.data.frame(TukeyHSD(rmanova_result)$prime_var)
    
    tukey_result <- data.frame(dat, row.names = rownames(dat)) 
    ind_result <- c() 
    
    for (j in 1:nrow(tukey_result)){
      ref_conf_list <- strsplit(rownames(TukeyHSD(rmanova_result)$prime_var), split = "-")[[j]]
      ref <- ref_conf_list[2]
      conf <- ref_conf_list[1]
      ref_mean <- mean(alpha_mani[,i+2][alpha_mani$prime_var == ref])
      conf_mean <- mean(alpha_mani[,i+2][alpha_mani$prime_var == conf])
      tukey_result_ind <- c(ref, conf, ref_mean, conf_mean, tukey_result[j,1:3], tukey_result[j, 4])
      ind_result <- rbind(ind_result, tukey_result_ind)
    }
    
    rownames(ind_result) <- 1:nrow(ind_result)
    colnames(ind_result) <- c("Ref", "Com", "Mean (Ref) ", "Mean (Com)", "Diff", "Lower", "Upper", "Adj. P.value")
    
    div_list[[i]] <- as.data.frame(ind_result)
  }
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

alpha.f.pair.mult.friedman <- function(alpha_mani){
  
  summary_vec <- c()
  
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  
  for (i in 1:(length(alpha_mani)-2)){
    friedman_result <- friedman.test(y = alpha_mani[, i+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var)
    friedman_ind_result <- c(decimal_adjust(friedman_result$statistic), friedman_result$parameter , p.value.0.1_char(decimal_adjust(friedman_result$p.value)))
    
    summary_vec <- rbind(summary_vec, friedman_ind_result)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[3:11]
  colnames(summary_vec) = c("Chisq", "Df", "P.value")
  
  return(summary_vec)
}

alpha.f.pair.mult.friedman_download <- function(alpha_mani){
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  summary_vec <- c()
  for (i in 1:(length(alpha_mani)-2)){
    friedman_result <- friedman.test(y = alpha_mani[, i+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var)
    friedman_ind_result <- c(friedman_result$statistic, friedman_result$parameter , friedman_result$p.value)
    
    summary_vec <- rbind(summary_vec, friedman_ind_result)
  }
  rownames(summary_vec) = colnames(alpha_mani)[3:11]
  colnames(summary_vec) = c("Chisq", "Df", "P.value")
  return(summary_vec)
}

alpha.f.pair.mult.friedman.p_val <- function(alpha_mani){
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  p_val <- c()
  for (i in 1:(length(alpha_mani)-2)){
    friedman_result <- friedman.test(y = alpha_mani[, i+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var)
    friedman_ind_p <- decimal_adjust(friedman_result$p.value)
    
    p_val <- c(p_val, friedman_ind_p)
  }
  return(p_val)
}

fill_missing_na <- function(data.mani){
  
  data.mani$block_var <- as.character(data.mani$block_var)
  data.mani$prime_var <- as.character(data.mani$prime_var)
  
  tab <- table(data.mani$block_var, data.mani$prime_var)
  miss_rc <- which(tab == 0, arr.ind = TRUE) 
  miss_rc <- data.frame(miss_rc)
  
  miss_prim <- colnames(tab)[miss_rc$col]
  miss_block <- rownames(tab)[miss_rc$row]
  
  length <- nrow(data.mani)
  for (i in 1:length(miss_prim)){
    data.mani[length+i, ] <- c(miss_block[i], miss_prim[i], rep(NA, 9))
  }
  data.mani <- data.mani[order(data.mani$block_var, data.mani$prime_var),]
  return(data.mani)
}

alpha.pair.mult.durbin <- function(alpha_mani){
  
  durbin_list <- list() 
  
  summary_vec <- c()
  summary_vec_2 <- c() 
  p_val <- c() 
  
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  
  for (i in 1:(length(alpha_mani)-2)){
    durbin_result <- PMCMRplus::durbinTest(y = alpha_mani[, i+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var)
    
    durbin_ind_result <- c(decimal_adjust(durbin_result$statistic), durbin_result$parameter , p.value.0.1_char(decimal_adjust(durbin_result$p.value)))
    durbin_ind_result_2 <- c(durbin_result$statistic, durbin_result$parameter , durbin_result$p.value)
    durbin_ind_p <- decimal_adjust(durbin_result$p.value)
    
    summary_vec <- rbind(summary_vec, durbin_ind_result)
    summary_vec_2 <- rbind(summary_vec_2, durbin_ind_result_2)
    
    p_val <- c(p_val, durbin_ind_p)
  }
  
  rownames(summary_vec) = colnames(alpha_mani)[3:11]
  colnames(summary_vec) = c("Chisq", "Df", "P.value")
  
  rownames(summary_vec_2) <- colnames(alpha_mani)[3:11]
  colnames(summary_vec_2) = c("Chisq", "Df", "P.value")
  
  durbin_list$table <- summary_vec
  durbin_list$download <- summary_vec_2 
  durbin_list$p.val <- p_val 
  
  return(durbin_list)
}

alpha.f.pair.mult.conover <- function(alpha_mani, p_adjustment = "BH"){
  div_list <- list() 
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  for (j in 1:(length(alpha_mani)-2)){
    conover_statistic <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$statistic
    conover_qval <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$p.value
    conover_pval <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = "none")$p.value
    
    ref<- c() 
    conf <- c()
    t <- c()
    Q.val <- c()
    P.val <- c()
    
    for (i in 1:length(rownames(conover_qval))){
      ref_ind <- rep(colnames(conover_qval)[i], length(na.omit(conover_pval[,i])))
      ref <- c(ref, ref_ind)
      
      conf_ind <- rownames(conover_qval)[i:length(rownames(conover_pval))]
      conf <- c(conf, conf_ind)
      
      t_ind<- decimal_adjust(na.omit(conover_statistic[,i]))
      t <- c(t, t_ind)
      
      Q_ind <- p.value.0.1_char(decimal_adjust(na.omit(conover_qval[,i])))
      P_ind <- p.value.0.1_char(decimal_adjust(na.omit(conover_pval[,i])))
      Q.val <- c(Q.val, Q_ind)
      P.val <- c(P.val, P_ind)
    }
    conover_result <- cbind(ref, conf, t,  Q.val)
    colnames(conover_result) <- c("Ref", "Com", "t", "Adj. P.value")
    rownames(conover_result) <- 1:nrow(conover_result)
    div_list[[j]] <- as.data.frame(conover_result)
  }
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

alpha.f.pair.mult.conover_download <- function(alpha_mani, p_adjustment = "BH"){
  div_list <- list() 
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  for (j in 1:(length(alpha_mani)-2)){
    conover_statistic <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$statistic
    conover_qval <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$p.value
    conover_pval <- frdAllPairsConoverTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = "none")$p.value
    
    ref<- c() 
    conf <- c()
    t <- c()
    Q.val <- c()
    P.val <- c()
    
    for (i in 1:length(rownames(conover_qval))){
      ref_ind <- rep(colnames(conover_qval)[i], length(na.omit(conover_pval[,i])))
      ref <- c(ref, ref_ind)
      
      conf_ind <- rownames(conover_qval)[i:length(rownames(conover_pval))]
      conf <- c(conf, conf_ind)
      
      t_ind<- decimal_adjust(na.omit(conover_statistic[,i]))
      t <- c(t, t_ind)
      
      Q_ind <- na.omit(conover_qval[,i])
      P_ind <- na.omit(conover_pval[,i])
      Q.val <- c(Q.val, Q_ind)
      P.val <- c(P.val, P_ind)
    }
    conover_result <- cbind(ref, conf, t, Q.val)
    colnames(conover_result) <- c("Ref", "Com", "t", "Adj. P.value")
    rownames(conover_result) <- 1:nrow(conover_result)
    div_list[[j]] <- as.data.frame(conover_result)
  }
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

alpha.pairwise.durbin <- function(alpha_mani, p_adjustment = "BH"){
  div_list <- list() 
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  for (j in 1:(length(alpha_mani)-2)){
    durbin_statistic <- PMCMRplus::durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$statistic
    durbin_qval <- durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$p.value
    durbin_pval <- durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = "none")$p.value
    
    ref<- c() 
    
    conf <- c()
    t <- c()
    Q.val <- c()
    P.val <- c()
    
    for (i in 1:length(rownames(durbin_qval))){
      ref_ind <- rep(colnames(durbin_qval)[i], length(na.omit(durbin_pval[,i])))
      ref <- c(ref, ref_ind)
      
      conf_ind <- rownames(durbin_qval)[i:length(rownames(durbin_pval))]
      conf <- c(conf, conf_ind)
      
      t_ind<- decimal_adjust(na.omit(durbin_statistic[,i]))
      t <- c(t, t_ind)
      
      Q_ind <- p.value.0.1_char(decimal_adjust(na.omit(durbin_qval[,i])))
      P_ind <- p.value.0.1_char(decimal_adjust(na.omit(durbin_pval[,i])))
      Q.val <- c(Q.val, Q_ind)
      P.val <- c(P.val, P_ind)
    }
    
    durbin_result <- cbind(ref, conf, t,  Q.val)
    colnames(durbin_result) <- c("Ref", "Com", "t", "Adj. P.value")
    rownames(durbin_result) <- 1:nrow(durbin_result)
    div_list[[j]] <- as.data.frame(durbin_result)
  }
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

alpha.pairwise.durbin_download <- function(alpha_mani, p_adjustment = "BH"){
  div_list <- list() 
  
  alpha_mani$block_var <- as.factor(alpha_mani$block_var)
  alpha_mani$prime_var <- as.factor(alpha_mani$prime_var)
  
  for (j in 1:(length(alpha_mani)-2)){
    durbin_statistic <- PMCMRplus::durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$statistic
    durbin_qval <- durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = p_adjustment)$p.value
    durbin_pval <- durbinAllPairsTest(y = alpha_mani[, j+2], groups = alpha_mani$prime_var, blocks = alpha_mani$block_var, p.adjust.method = "none")$p.value
    
    ref<- c() 
    conf <- c()
    t <- c()
    Q.val <- c()
    P.val <- c()
    
    for (i in 1:length(rownames(durbin_qval))){
      ref_ind <- rep(colnames(durbin_qval)[i], length(na.omit(durbin_pval[,i])))
      ref <- c(ref, ref_ind)
      
      conf_ind <- rownames(durbin_qval)[i:length(rownames(durbin_pval))]
      conf <- c(conf, conf_ind)
      
      t_ind<- decimal_adjust(na.omit(durbin_statistic[,i]))
      t <- c(t, t_ind)
      
      Q_ind <- durbin_qval[,i]
      P_ind <- durbin_pval[,i]
      Q.val <- c(Q.val, Q_ind)
      P.val <- c(P.val, P_ind)
    }
    
    durbin_result <- cbind(ref, conf, t, P.val, Q.val)
    colnames(durbin_result) <- c("Ref", "Com", "t", "P.value", "Q.value")
    rownames(durbin_result) <- 1:nrow(durbin_result)
    div_list[[j]] <- as.data.frame(durbin_result)
  }
  names(div_list) <- names(alpha_mani)[3:11]
  return(div_list)
}

##################
# Visualization  #
##################

alpha.bin.multi.paired <- function(alpha_div, alpha_mani){
  ind <- list()
  time.p.cat <- names(table(alpha_mani[,"prime_var"]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(alpha_mani[,"prime_var"] == time.p.cat[i])
  }
  
  each_alpha <- list()
  
  for (j in 1:length(alpha_div)){
    alpha_dat<- c()
    for (i in 1:length(time.p.cat)){
      alpha_dat <- cbind(alpha_dat, alpha_mani[, 2+j][ind[[i]]])
      }
    colnames(alpha_dat) <- time.p.cat
    each_alpha[[j]] <- alpha_dat
  }
  
  length <- length(time.p.cat)
  if (length <= 3){
    par(mfrow = c(3, 3), mar = c(4.1, 4.1, 1, 1.9))
    
    for (j in 1:length(names(alpha_div))){
      boxplot(each_alpha[[j]], ylab = names(alpha_div)[j], names = time.p.cat,   notch = TRUE, horizontal = FALSE, col = "Paleturquoise1",boxwex=0.7) #col = rainbow(length),
    }
    
  }else if (4 <= length & length <= 6){
    par(mfrow = c(5, 2), mar = c(4.1, 4.1, 1, 1.9))
    for (j in 1:length(alpha_div)){
      boxplot(each_alpha[[j]], ylab = names(alpha_div)[j], names = time.p.cat,  notch = TRUE, horizontal = FALSE, col = "Paleturquoise1", boxwex=0.7) # col = rainbow(length),
    }
  }else {
    par(mfrow = c(9, 1), mar = c(4.1, 4.1, 1, 1.9))
    for (j in 1:length(alpha_div)){
      boxplot(each_alpha[[j]], ylab = names(alpha_div)[j],  names = time.p.cat,  notch = TRUE, horizontal = FALSE, col = "Paleturquoise1", boxwex=0.7) # col = rainbow(length),
    }
  }
}

alpha.bin.f.mult.paired <- function(alpha_div, alpha_mani, test_result){
  ind <- list()
  time.p.cat <- names(table(alpha_mani[, "prime_var"]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(alpha_mani[, "prime_var"] == time.p.cat[i])
  }
  
  each_alpha <- list() 
  for (j in 1:length(alpha_div)){
    alpha_dat <- c() 
    for (i in 1:length(time.p.cat)){
      alpha_dat <- cbind(alpha_dat, alpha_mani[, 2+j][ind[[i]]])
    }
    colnames(alpha_dat) <- time.p.cat
    each_alpha[[j]] <- alpha_dat
  }
  
  ind.p.sig <- which(test_result < 0.05)
  
  length <- length(time.p.cat)
  if (length <= 3){
    par(mfrow = c(3, 3), mar = c(4.3, 4.1, 1, 1.9)) 
    
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      }

      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v, names = time.p.cat,  notch = TRUE, horizontal = FALSE, col = "Paleturquoise1", boxwex=0.7) #col = rainbow(length),
    }
    
  }else if (4 <= length & length <= 6){
    par(mfrow = c(5, 2), mar = c(4.1, 4.1, 1, 1.9)) 
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(test_result[i])), sep = "")
      }
      
      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v, names = time.p.cat, notch = TRUE, col = "Paleturquoise1", horizontal = FALSE, boxwex=0.7) #col = rainbow(length),
    }
  }else {
    par(mfrow = c(9, 1), mar = c(4.1, 4.1, 1, 1.9))
    for (i in 1:length(alpha_div)){
      if (is.element(i, ind.p.sig)){
        xlab.v <- paste("*p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      } else {
        xlab.v <- paste("p:", p.value.0.1_char(decimal_adjust(as.numeric(test_result[i]))), sep = "")
      }
      boxplot(each_alpha[[i]], ylab = names(alpha_div)[i], xlab = xlab.v,  names = time.p.cat,  notch = TRUE, col = "Paleturquoise1", horizontal = FALSE, boxwex=0.7) #col = rainbow(length), 
    }
  }
}

global.alpha.lmm.p <- function(alpha_mani, alpha_relevel){
  alpha_relevel_dat <- alpha_mani 
  alpha_relevel_dat[,"prime_var"] <- as.factor(alpha_relevel_dat[,"prime_var"])
  alpha_relevel_dat[,"prime_var"] <- relevel(alpha_relevel_dat[,"prime_var"], alpha_relevel)
  alpha_relevel_dat <<- alpha_relevel_dat 
  
  p_val <- c() 
  
  for (i in 3:length(alpha_relevel_dat)){
    fit_full <- lmer(alpha_relevel_dat[,i] ~ alpha_relevel_dat[,2]+ (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    fit_nested <- lmer(alpha_relevel_dat[,i] ~ 1 + (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    p_val <- c(p_val, as.numeric(unlist(anova(fit_full, fit_nested))["Pr(>Chisq)2"]))
  }
  return(p_val)
}

global.alpha.lmm.dat <- function(alpha_mani, alpha_relevel, download = FALSE){
  alpha_relevel_dat <- alpha_mani 
  alpha_relevel_dat[,"prime_var"] <- as.factor(alpha_relevel_dat[,"prime_var"])
  alpha_relevel_dat[,"prime_var"] <- relevel(alpha_relevel_dat[,"prime_var"], alpha_relevel)
  alpha_relevel_dat <<- alpha_relevel_dat 
  
  mat <- matrix(NA, nrow = 9, ncol = 5)
  
  for (i in 3:length(alpha_relevel_dat)){
    fit_full <- lmer(alpha_relevel_dat[,i] ~ alpha_relevel_dat[,2]+ (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    fit_nested <- lmer(alpha_relevel_dat[,i] ~ 1 + (1|alpha_relevel_dat[,1]), data = alpha_relevel_dat)
    if (download){
      mat[i-2,] <- unlist(anova(fit_full, fit_nested, test = "LRT"))[c("logLik1", "logLik2", "Chisq2", "Df2", "Pr(>Chisq)2")]
    }else{
      mat[i-2,] <- c(decimal_adjust(unlist(anova(fit_full, fit_nested, test = "LRT"))[c("logLik1", "logLik2", "Chisq2", "Df2")]), p.value.0.1_char(decimal_adjust(unlist(anova(fit_full, fit_nested, test = "LRT"))["Pr(>Chisq)2"])))
    }
  }
  
  rownames(mat) <- names(alpha_mani)[3:11]
  colnames(mat) <- c("logLik : Nested", "logLik : Complex", "Chisq", "Df", "P.value")
  return(mat)
}

alpha.pairwise.lmm <- function(alpha_mani, ref_level, method_adj = "BH", download = FALSE){
  alpha_relevel_dat <- alpha_mani
  alpha_relevel_dat[,"prime_var"] <- as.factor(alpha_relevel_dat[,"prime_var"])
  alpha_relevel_dat[,"prime_var"] <- relevel(alpha_relevel_dat[,"prime_var"], ref_level)
  alpha_relevel_dat <<- alpha_relevel_dat 
  
  n_row <- ncol(alpha_mani) -2
  result_list <- list() 
  
  for (i in 1:n_row){
    result <- lmer(alpha_relevel_dat[,i+2] ~ alpha_relevel_dat[,2] + (1|alpha_relevel_dat[,1]), dat = alpha_relevel_dat)
    mat <- coef(summary(result))[2:nrow(coef(summary(result))), -5]
    
    if(download){
      q.val <- p.adjust(coef(summary(result))[2:nrow(coef(summary(result))),5], method_adj)
      mat <- cbind(mat, q.val)
      colnames(mat) <- c("Est", "SE", "DF", "t", "Adj. P.value")
      rownames(mat) <- vapply(names(table(alpha_relevel_dat$prime_var))[2:length(names(table(alpha_relevel_dat$prime_var)))], function(x)paste0(ref_level, " - ", x), character(1))
      result_list[[i]] <- as.data.frame(mat)
    }else{
      q.val <- p.adjust(coef(summary(result))[2:nrow(coef(summary(result))),5], method_adj)
      mat <- cbind(mat, q.val)
      mat <- apply(mat, 2, function(x) decimal_adjust(x))
      mat <- cbind(ref_level, names(table(alpha_relevel_dat$prime_var))[2:length(names(table(alpha_relevel_dat$prime_var)))], mat)
      
      colnames(mat) <- c("Ref", "Com", "Est", "SE", "DF", "t",  "Adj. P.value")
      rownames(mat) <- NULL
      mat[, "Adj. P.value"] <- p.value.0.1_char(mat[, "Adj. P.value"])
      result_list[[i]] <- as.data.frame(mat)
    }
  }
  names(result_list) <- colnames(alpha_mani)[3:11]
  return(result_list)
}

reduced_alpha <- function(alpha.dat, sam.dat, prim_id, block_id, levels){
  
  alpha_united <- cbind(prim.var = sam.dat[[prim_id]], block.id = sam.dat[[block_id]], alpha.dat)
  alpha_ind <- alpha_united[which(alpha_united[,"prim.var"] %in% levels),]
  alpha_result <- alpha_ind[,3:ncol(alpha_ind)]

  return(alpha_result)
}




######################
#######Shiny##########
######################
prettyRadioButtons_new <- function (
  inputId, label, choices = NULL, selected = NULL, status = "primary", 
  shape = c("round", "square", "curve"), outline = FALSE, 
  fill = FALSE, thick = FALSE, animation = NULL, icon = NULL, 
  plain = FALSE, bigger = FALSE, inline = FALSE, width = NULL, 
  choiceNames = NULL, choiceValues = NULL) 
{
  status <- match.arg(status, c("default", "primary", "success", 
                                "info", "danger", "warning"))
  shape <- match.arg(shape)
  if (is.null(choices) && is.null(choiceNames) && is.null(choiceValues)) {
    choices <- character(0)
  }
  args <- shinyWidgets:::normalizeChoicesArgs(choices, choiceNames, choiceValues)
  selected <- shiny::restoreInput(id = inputId, default = selected)
  selected <- if (is.null(selected)) {
    args$choiceValues[[1]]
  }
  else {
    as.character(selected)
  }
  if (length(selected) > 1) 
    stop("The 'selected' argument must be of length 1")
  options1 <- shinyWidgets:::generatePretty(
    inputId = inputId, selected = selected, 
    inline = inline, type = "radio", choiceNames = args$choiceNames[1:2], 
    choiceValues = args$choiceValues[1:2], status = status, shape = shape, 
    outline = outline, fill = fill, thick = thick, animation = animation, 
    icon = icon, plain = plain, bigger = bigger
  )
  options2 <- shinyWidgets:::generatePretty(
    inputId = inputId, selected = selected, 
    inline = inline, type = "radio", choiceNames = args$choiceNames[3], 
    choiceValues = args$choiceValues[3], status = status, shape = shape, 
    outline = outline, fill = fill, thick = thick, animation = animation, 
    icon = icon, plain = plain, bigger = bigger
  )
  
  options <- tags$div(
    tags$div(
      tags$fieldset(
        tags$legend("More than two-group comparison (across groups)", style="font-size:11pt;font-weight:bold;"),
        options1
      )
    ),
    tags$div(
      tags$fieldset(
        tags$legend("More than two-group comparison (baseline to other groups)", style="font-size:11pt;font-weight:bold;"),
        tags$div(
          style = "display: inline-block;",
          options2
        ),
        
      )
    )
  )
  divClass <- "form-group shiny-input-radiogroup shiny-input-container"
  if (inline) 
    divClass <- paste(divClass, "shiny-input-container-inline")
  radioTag <- htmltools::tags$div(id = inputId, style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), class = divClass, 
    tags$label(class = "control-label", `for` = inputId, 
               class = if (is.null(label)) 
                 "shiny-label-null", label), options)
  shinyWidgets:::attachShinyWidgetsDep(radioTag, "pretty")
}

prettyRadioButtons_4 <- function (
  inputId, label, choices = NULL, selected = NULL, status = "primary", 
  shape = c("round", "square", "curve"), outline = FALSE, 
  fill = FALSE, thick = FALSE, animation = NULL, icon = NULL, 
  plain = FALSE, bigger = FALSE, inline = FALSE, width = NULL, 
  choiceNames = NULL, choiceValues = NULL) 
{
  status <- match.arg(status, c("default", "primary", "success", 
                                "info", "danger", "warning"))
  shape <- match.arg(shape)
  if (is.null(choices) && is.null(choiceNames) && is.null(choiceValues)) {
    choices <- character(0)
  }
  args <- shinyWidgets:::normalizeChoicesArgs(choices, choiceNames, choiceValues)
  selected <- shiny::restoreInput(id = inputId, default = selected)
  selected <- if (is.null(selected)) {
    args$choiceValues[[1]]
  }
  else {
    as.character(selected)
  }
  if (length(selected) > 1) 
    stop("The 'selected' argument must be of length 1")
  options1 <- shinyWidgets:::generatePretty(
    inputId = inputId, selected = selected, 
    inline = inline, type = "radio", choiceNames = args$choiceNames[1:3], 
    choiceValues = args$choiceValues[1:3], status = status, shape = shape, 
    outline = outline, fill = fill, thick = thick, animation = animation, 
    icon = icon, plain = plain, bigger = bigger
  )
  options2 <- shinyWidgets:::generatePretty(
    inputId = inputId, selected = selected, 
    inline = inline, type = "radio", choiceNames = args$choiceNames[4], 
    choiceValues = args$choiceValues[4], status = status, shape = shape, 
    outline = outline, fill = fill, thick = thick, animation = animation, 
    icon = icon, plain = plain, bigger = bigger
  )
  
  options <- tags$div(
    tags$div(
      tags$fieldset(
        tags$legend("More than two-group comparison (across groups)", style="font-size:11pt;font-weight:bold;"),
        options1
      )
    ),
    tags$div(
      tags$fieldset(
        tags$legend("More than two-group comparison (baseline to other groups)", style="font-size:11pt;font-weight:bold;"),
        tags$div(
          style = "display: inline-block;",
          options2
        ),
        
      )
    )
  )
  divClass <- "form-group shiny-input-radiogroup shiny-input-container"
  if (inline) 
    divClass <- paste(divClass, "shiny-input-container-inline")
  radioTag <- htmltools::tags$div(id = inputId, style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), class = divClass, 
    tags$label(class = "control-label", `for` = inputId, 
               class = if (is.null(label)) 
                 "shiny-label-null", label), options)
  shinyWidgets:::attachShinyWidgetsDep(radioTag, "pretty")
}

