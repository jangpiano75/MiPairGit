#########################################################
# Comparative analysis for multi-paired beta diversity  #
#########################################################

reduced_beta <- function(beta.dat, sam.dat, prim_id, block_id, levels){
  
  mat_vec <- match(rownames(beta.dat[[1]]), rownames(beta.dat[[1]]))
  dat <- list() 
  if(length(mat_vec[mat_vec == FALSE]) != 0){
    
    for (i in 1:length(beta.dat)){
      ind <- match(rownames(beta.dat[[i]], rownames(sam.dat)))
      ind_2 <- which(sam.dat[[prim_id]] %in% levels)
      dat[[i]] <- (beta.dat[[i]][ind, ind])[ind_2, ind_2]
    }
  }else{
    ind_2 <- which(sam.dat[[prim_id]] %in% levels)

    for (i in 1:length(beta.dat)){
      dat[[i]] <- beta.dat[[i]][ind_2, ind_2]
    }
  }
  names(dat) <- names(beta.dat)
  return(dat)
}

beta.bin.mult.cat.ref.func_sub <- function(sam.dat, prim_id, block_id, Ds.Ks, levels){
  bin.var <- factor(unlist(sam.dat[which(sam.dat[[prim_id]] %in% levels), prim_id]))
  
  Ds <- reduced_beta(Ds.Ks$Ds, sam.dat, prim_id, block_id, levels) 
  Ks <- reduced_beta(Ds.Ks$Ks, sam.dat, prim_id, block_id, levels) 
  
  return(list(bin.var = bin.var, Ds = Ds, Ks = Ks))
}

beta.permanova.paired.multi <- function(num_perm = 3000, div_num, beta_div, sam_dat, prim_id, block_id, method_adj = "BH"){
  
  beta_list <- list() 
  
  ind <- list() 
  time.p.cat <- names(table(sam_dat[,prim_id]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(sam_dat[,prim_id] == time.p.cat[i])
  }
  
  comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, ncol(comb), 5)
  out_2 <- matrix(NA, ncol(comb), 5)
  beta_div_ind <- beta_div[[div_num]]
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    combined_ind <- c(ind[[combination[[1]]]], ind[[combination[[2]]]])
    
    beta.dat.pair <- beta_div_ind[combined_ind, combined_ind]
    sam.dat.pair <- sam_dat[combined_ind,]
    
    perm <- how(nperm = num_perm)
    setBlocks(perm) <- sam.dat.pair[[block_id]]
    
    fit <- adonis2(beta.dat.pair ~ sam.dat.pair, data = data.frame(sam.dat.pair = sam.dat.pair[[prim_id]]), permutations = perm)
    
    out[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
    out_2[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], decimal_adjust(fit[1, "F"]), decimal_adjust(fit[1,"R2"]), fit[1,"Pr(>F)"])
  }
  
  q_val <- p.adjust(out[,5], method = method_adj)
  out <- as.data.frame(cbind(out[,1:3], out[,5], q_val))
  
  q_val_2 <- p.value.0.1_char(decimal_adjust(p.adjust(out_2[,5], method = method_adj, n = length(out_2[,5]))))  
  out_2 <- as.data.frame(cbind(out_2[,1:3],  q_val_2))
  
  colnames(out) <- c("Ref", "Com", "F", "P.value", "Adj. P.value")
  colnames(out_2) <- c("Ref", "Com", "F", "Adj. P.value")
  
  beta_list$download <- out
  beta_list$table <- out_2 
  
  return(beta_list)
}


beta.permanova.paired.base <- function(num_perm = 3000, level, div_num, beta_div, sam_dat, prim_id, block_id, method_adj = "BH"){
  
  beta_list <- list() 
  ind <- list()
  
  base_ind <- which(sam_dat[, prim_id] == level)
  time.p.cat <- names(table(sam_dat[-base_ind, prim_id]))
  
  
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(sam_dat[,prim_id] == time.p.cat[i])
  }
  
  #comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, length(ind), 5)
  out_2 <- matrix(NA, length(ind), 5)
  beta_div_ind <- beta_div[[div_num]]
  
  for (j in 1:length(ind)){
    
    combined_ind <- c(base_ind, ind[[j]])
    
    beta.dat.pair <- beta_div_ind[combined_ind, combined_ind]
    sam.dat.pair <- sam_dat[combined_ind,]
    
    perm <- how(nperm = num_perm)
    setBlocks(perm) <- sam.dat.pair[[block_id]]
    
    fit <- adonis2(beta.dat.pair ~ sam.dat.pair, data = data.frame(sam.dat.pair = sam.dat.pair[[prim_id]]), permutations = perm)
    
    out[j,] <- c(level, time.p.cat[j], fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
    out_2[j,] <- c(level, time.p.cat[j], decimal_adjust(fit[1, "F"]), decimal_adjust(fit[1,"R2"]), fit[1,"Pr(>F)"])
  }
  
  q_val <- p.adjust(out[,5], method = method_adj)
  out <- as.data.frame(cbind(out[,1:3], out[,5], q_val))
  
  q_val_2 <- p.value.0.1_char(decimal_adjust(p.adjust(out_2[,5], method = method_adj, n = length(out_2[,5]))))  
  out_2 <- as.data.frame(cbind(out_2[,1:3], q_val_2))
  
  colnames(out) <- c("Ref", "Com", "F", "P.value", "Adj. P.value")
  colnames(out_2) <- c("Ref", "Com", "F", "Adj. P.value")
  
  beta_list$download <- out
  beta_list$table <- out_2 
  
  return(beta_list)
}



permanova.paired.mult.united <- function(num_perm = 3000, beta_div, sam_dat, prim_id, block_id, method_adj = "BH", download = FALSE){
  total_list <- list()
  download_list <- list() 
  table_list <- list()
  
  for (i in 1:length(beta_div)){
    div_list <- beta.permanova.paired.multi(num_perm = 10000, i, beta_div, sam_dat, prim_id, block_id, method_adj)
    download_list[[i]] <- div_list$download 
    table_list[[i]] <- div_list$table
  }
  
  names(download_list) <- names(beta_div) 
  names(table_list) <- names(beta_div)
  
  total_list$download <- download_list
  total_list$table <- table_list
  return(total_list)
}

permanova.paired.mult.united_base <- function(num_perm = 3000, level, beta_div, sam_dat, prim_id, block_id, method_adj = "BH", download = FALSE){
  total_list <- list()
  download_list <- list() 
  table_list <- list()
  
  for (i in 1:length(beta_div)){
    div_list <- beta.permanova.paired.base (num_perm = 10000, level, i, beta_div, sam_dat, prim_id, block_id, method_adj)
    download_list[[i]] <- div_list$download 
    table_list[[i]] <- div_list$table
  }
  
  names(download_list) <- names(beta_div) 
  names(table_list) <- names(beta_div)
  
  total_list$download <- download_list
  total_list$table <- table_list
  return(total_list)
}


beta.permanovaFL.paired.multi <- function(num_perm = 3000, div_num, beta_div, sam_dat, prim_id, block_id, method_adj = "BH"){
  beta_list <- list() 
  
  ind <- list() 
  time.p.cat <- names(table(sam_dat[,prim_id]))
  for (i in 1:length(time.p.cat)){
    ind[[i]] <- which(sam_dat[,prim_id] == time.p.cat[i])
  }

  comb <- combn(length(time.p.cat), 2)
  out <- matrix(NA, ncol(comb), 4)
  out_2 <- matrix(NA, ncol(comb), 4)
  beta_div_ind <- beta_div[[div_num]]
  
  for (j in 1:ncol(comb)){
    combination <- comb[,j]
    combined_ind <- c(ind[[combination[[1]]]], ind[[combination[[2]]]])

    beta.dat.pair <<- beta_div_ind[combined_ind, combined_ind]
    sam.dat.pair <- sam_dat[combined_ind,]
    data_mani <- data.frame(prim.var = sam.dat.pair[[prim_id]], block.id = sam.dat.pair[[block_id]],  row.names = rownames(sam.dat.pair))
    
    fit_J <- LDM::permanovaFL(beta.dat.pair|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = num_perm)
    out[j,] <-c(time.p.cat[combination[1]], time.p.cat[combination[2]], fit_J$F.statistics, fit_J$p.permanova)
    out_2[j,] <- c(time.p.cat[combination[1]], time.p.cat[combination[2]], decimal_adjust(fit_J$F.statistics), fit_J$p.permanova)
  }
  
  q_val_2 <- p.value.0.1_char(decimal_adjust(p.adjust(out_2[,4], method = method_adj, n = length(out_2[,4]))))
  out_2 <- as.data.frame(cbind(out_2[,1:3], p.value.0.1_char(decimal_adjust(out_2[,4])), q_val_2))
  colnames(out_2) <- c("Ref", "Com", "F", "P.value", "Adj. P.value")
  
  q_val <- p.adjust(out[,4], method = method_adj)
  out <- as.data.frame(cbind(out[,1:3], out[,4], q_val))
  colnames(out) <- c("Ref", "Com", "F", "P.value", "Adj. P.value")
  
  beta_list$download <- out 
  beta_list$table <- out_2 
  
  return(beta_list)
}

permanovaFL.paired.mult.united <- function(num_perm = 3000, beta_div, sam_dat, prim_id, block_id, method_adj = "BH"){
  total_list <- list()
  download_list <- list() 
  table_list <- list()
  
  for (i in 1:length(beta_div)){
    div_list <- beta.permanovaFL.paired.multi(num_perm = 3000, i, beta_div, sam_dat, prim_id, block_id, method_adj)
    download_list[[i]] <- div_list$download 
    table_list[[i]] <- div_list$table
  }
  
  names(download_list) <- names(beta_div) 
  names(table_list) <- names(beta_div)
  
  total_list$download <- download_list
  total_list$table <- table_list
  
  return(total_list)
}

beta.bin.mult.cat.ref.func <- function(sel.bin.var, sam.dat, Ds.Ks){
  bin.var <- factor(unlist(sam.dat[,sel.bin.var]))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]]
  }
  return(list(bin.var = bin.var, Ds = Ds, Ks = Ks))
}

beta.permanova.mult.glob.p_val <- function(num_perm = 3000, beta_div, sam_dat, prim_id, block_id){
  p_val_result <- c() 
  for (i in 1:length(names(beta_div))){
    perm <- how(nperm = num_perm)
    setBlocks(perm) <- sam_dat[[block_id]]
    
    beta_div_ind <- beta_div[[i]]
    fit <- adonis2(beta_div_ind ~ sam_dat_pair, data = data.frame(sam_dat_pair = sam_dat[[prim_id]]), permutations = perm)
    
    each_p_result <- decimal_adjust(fit[1, "Pr(>F)"])
    p_val_result <- c(p_val_result, each_p_result)
  }

  return(p_val_result)
}

beta.permanova.paired.mult.glob <- function(num_perm = 3000, beta_div, sam_dat, prim_id, block_id, download = FALSE){
  
  global_result <- list() 
  p_val_result <- c() 
  
  out <- matrix(NA, length(beta_div), 3)
  
  for (i in 1:length(names(beta_div))){
    perm <- how(nperm = num_perm)
    setBlocks(perm) <- sam_dat[[block_id]]
    beta_div_ind <- beta_div[[i]]
    
    fit <- adonis2(beta_div_ind ~ sam_dat_pair, data = data.frame(sam_dat_pair = sam_dat[[prim_id]]), permutations = perm)
    
    if (download){
      out[i,] <- c(fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
      each_p_result <- decimal_adjust(fit[1, "Pr(>F)"])
      p_val_result <- c(p_val_result, each_p_result)
    }else{
      out[i,] <- c(decimal_adjust(fit[1, "F"]), decimal_adjust(fit[1,"R2"]), decimal_adjust(fit[1,"Pr(>F)"]))
      each_p_result <- decimal_adjust(fit[1, "Pr(>F)"])
      p_val_result <- c(p_val_result, each_p_result)
    }
  }
  
  colnames(out) <- c("F", "R2", "P.value")
  rownames(out) <- names(beta_div)
  global_result$table <- out
  global_result$p.val <- p_val_result 
  
  return(global_result)
}

beta.permanovaFL.paired.mult.glob <- function(num_perm = 3000, beta_div, sam_dat, prim_id, block_id, download = FALSE){
  
  beta_list <- list() 
  out <- matrix(NA, length(beta_div), 2)
  p_val_result <- c() 
  
  for (i in 1:length(names(beta_div))){
    beta.dat.ind <<- beta_div[[i]]
    data_mani <- data.frame(prim.var = sam_dat[[prim_id]], block.id = sam_dat[[block_id]], row.names = rownames(sam_dat))
    
    fit <- permanovaFL(beta.dat.ind|as.factor(block.id) ~ prim.var, dat = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none", cluster.id = block.id, n.perm.max = num_perm)
    
    if (download){
      out[i,] <- c(fit$F.statistics, fit$p.permanova)
      each_p_result <- fit$p.permanova
      p_val_result <- c(p_val_result, each_p_result)
    }else{
      out[i,] <- c(decimal_adjust(fit$F.statistics), decimal_adjust(fit$p.permanova))
      each_p_result <- fit$p.permanova
      p_val_result <- c(p_val_result, each_p_result)
    }
  }
  
  rownames(out) <- names(beta_div)
  colnames(out) <- c("F", "P.value")
  
  beta_list$table <- out 
  beta_list$p.val <- p_val_result 
  return(beta_list)
}

mirkat.bin.plot_2_mult <- function(p_val_result, beta.bin.out, sam_dat, sel.bin.var) {
  par(mfrow = c(3, 2), mar = c(5.5, 5, 2, 2))  
  n <- length(table(sam_dat[,sel.bin.var]))
  ind.p.sig <- which(p_val_result < 0.05)
  
  for (i in 1:length(beta.bin.out$Ds)) {
    if (is.element(i, ind.p.sig)){
      xlab.v <- paste("*p:", p.value.0.1(as.numeric(p_val_result[i])), sep = "")
    } else {
      xlab.v <- paste("p:", p.value.0.1(as.numeric(p_val_result[i])), sep = "")
    }
    
    if (is.element(i, ind.p.sig)){
      x.label = paste("PC 1 \n", "*p:", p.value.0.1(as.numeric(p_val_result[i])), sep = "")
    }else{
      x.label = paste("PC 1 \n", "p:", p.value.0.1(as.numeric(p_val_result[i])), sep = "")
    }
      
    mod <- betadisper(as.dist(beta.bin.out$Ds[[i]]), beta.bin.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i], xlab= x.label, ylab="PC 2",
         sub=NA, col = 1:n, label = FALSE, mgp=c(4,1,0), cex=1.5, label.cex = 0.5, cex.lab=1.2, cex.main=1.7)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, fill = c(1:n, cex=2.5, box.lty=0), legend = levels(beta.bin.out$bin.var), bty = "n", cex=1.6)
}


mirkat.bin.plot_2_base <- function(beta.bin.out, sam_dat, sel.bin.var) {
  par(mfrow = c(3, 2), mar = c(5.5, 5, 2, 2))  
  n <- length(table(sam_dat[,sel.bin.var]))
  
  for (i in 1:length(beta.bin.out$Ds)) {
    mod <- betadisper(as.dist(beta.bin.out$Ds[[i]]), beta.bin.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i],  xlab = "PC 1", ylab="PC 2",
         sub=NA, col = 1:n, label = FALSE, mgp=c(4,1,0), cex=1.5, label.cex = 0.5, cex.lab=1.2, cex.main=1.7)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, fill = c(1:n, cex=2.5, box.lty=0), legend = levels(beta.bin.out$bin.var), bty = "n", cex=1.6)
}


#################
######Shiny######
#################

prettyRadioButtons_beta_new <- function (
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
    inline = inline, type = "radio", choiceNames = args$choiceNames[1], 
    choiceValues = args$choiceValues[1], status = status, shape = shape, 
    outline = outline, fill = fill, thick = thick, animation = animation, 
    icon = icon, plain = plain, bigger = bigger
  )
  options2 <- shinyWidgets:::generatePretty(
    inputId = inputId, selected = selected, 
    inline = inline, type = "radio", choiceNames = args$choiceNames[2], 
    choiceValues = args$choiceValues[2], status = status, shape = shape, 
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

