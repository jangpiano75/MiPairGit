library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(quantreg)
library(fossil)
library(picante)
library(entropart)
library(lme4)
library(lmerTest)
#library(dirmult) 
#library(robustbase)
library(robCompositions) 
#library(BiasedUrn)
#library(CompQuadForm)
library(GUniFrac) 
library(ecodist) 
#library(MiRKAT)
#library(GLMMMiRKAT)
#library(proxy)

#####################
# Data manipulation #
#####################

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  return(bin.cat)
}

beta.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

beta.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, Ds.Ks) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  return(list(bin.var = bin.var, Ds = Ds, Ks = Ks))
}

beta.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  
  return(sam.dat)
}

##################
# Beta diversity #
##################

Ds.Ks.func <- function(rare.biom, biom.after.qc) {
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  no.rare.tree <- phy_tree(biom.after.qc)
  
  jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
  bc <- as.matrix(bcdist(t(rare.otu.tab)))
  unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
  u.unif <- unifs[, , "d_UW"]
  g.unif <- unifs[, , "d_0.5"]
  w.unif <- unifs[, , "d_1"]
  
  jac.k <- D2K(jac)
  bc.k <- D2K(bc)
  u.unif.k <- D2K(u.unif)
  g.unif.k <- D2K(g.unif)
  w.unif.k <- D2K(w.unif)
  
  rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
  rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
  rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
  rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
  rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
  
  return(
    list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
         Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k))
  )
}

beta.bin.cov.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var, sam.dat, Ds.Ks) {  
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  
  return(list(bin.var = bin.var, cov.var = cov.var, Ds = Ds, Ks = Ks))
}

beta.con.recode.func <- function(sam.dat, sel.con.var, rename.con.var, Ds.Ks) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][ind.nona, ind.nona]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][ind.nona, ind.nona]
  }
  
  return(list(con.var = con.var, Ds = Ds, Ks = Ks))
}

beta.con.cov.recode.func <- function(sam.dat, sel.con.var, sel.cov.var, rename.con.var, Ds.Ks) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  cov.var <- as.data.frame(sam.dat[ind.nona, sel.cov.var])
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][ind.nona, ind.nona]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][ind.nona, ind.nona]
  }
  
  return(list(con.var = con.var, cov.var = cov.var, Ds = Ds, Ks = Ks))
}

####################################################
# Comparative analysis for paired beta diversity  #
####################################################

beta.permanova.paired <- function(num_perm = 10000, beta_div, sam_dat, prim_id, block_id, download = FALSE){
  beta_list <- list() 
  
  out <- matrix(NA, length(names(beta_div)), 3)  # F, R2, P-value
  out_2 <- matrix(NA, length(names(beta_div)), 3)
  perm <- how(nperm = num_perm)
  setBlocks(perm) <- sam_dat[[block_id]]
  
  for (i in 1:length(names(beta_div))){
    fit <- adonis2(beta_div[[i]] ~ sam_dat, data = data.frame(sam_dat = sam_dat[[prim_id]]), permutations = perm)
    out[i,] <- c(fit[1, "F"], fit[1,"R2"], fit[1,"Pr(>F)"])
    out_2[i,] <- c(round(fit[1, "F"], 3), round(fit[1,"R2"], 3), round(fit[1,"Pr(>F)"], 3))
  }
  
  out <- as.data.frame(out)
  rownames(out) <- names(beta_div)
  colnames(out) <- c("F", "R2", "P.value")
  
  out_2 <- as.data.frame(out_2)
  rownames(out_2) <- names(beta_div)
  colnames(out_2) <- c("F", "R2", "P.value")
  
  beta_list$download <- out
  beta_list$table <- out_2 
  
  return(beta_list)
}

beta.permanovaFL.paired.ind <- function(sam_dat, prim, block, nperm = 3000){
  beta_list <- list() 
  
  out <- matrix(NA, 5, 2)
  out_2 <- matrix(NA, 5, 2)
  data_mani <- data.frame(prim.var = sam_dat[[prim]], block.id = sam_dat[[block]],  row.names = rownames(sam_dat))
  
  fit_J <- permanovaFL(div.Jac|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = nperm)
  fit_B <- permanovaFL(div.Bray|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = nperm)
  fit_U <- permanovaFL(div.U|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = nperm)
  fit_G <- permanovaFL(div.G|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = nperm)
  fit_W <- permanovaFL(div.W|as.factor(block.id) ~ prim.var , data = data_mani, dist.method = NULL, perm.within.type = "free", perm.between.type = "none",  cluster.id = block.id, n.perm.max = nperm)
  fit <- list(J = fit_J, B = fit_B, U = fit_U, G = fit_G, W = fit_W)
  
  for (i in 1:length(fit)){
    out[i,] <- c(round(fit[[i]]$F.statistics, 3), round(fit[[i]]$p.permanova,3)) 
    out_2[i,] <- c(fit[[i]]$F.statistics, fit[[i]]$p.permanova)
  }
  
  out <- as.data.frame(out)
  out_2 <- as.data.frame(out_2)
  colnames(out) <- c("F", "P.value")
  colnames(out_2) <- c("F", "P.value")
  
  beta_list$download <- out_2 
  beta_list$table <- out 
  
  return(beta_list)
}

###########################
# Visualization (Boxplot) #
###########################

mirkat.bin.plot_2 <- function(out, beta.bin.out) {
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.bin.out$Ds)) {
    if (out$P.value[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$P.value[i]), sep="")
    }
    if (out$P.value[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$P.value[i]), sep="")
    }
    mod <- betadisper(as.dist(beta.bin.out$Ds[[i]]), beta.bin.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub=NA, col = c("blue2", "red2"), mgp=c(2.5,1,0), cex=1.7, label.cex=1.3, cex.lab=1.2, cex.main=1.7)
    mtext(sub.tit, side=1, line=3.8, cex=1.0)
  }
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(beta.bin.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.6)
}

