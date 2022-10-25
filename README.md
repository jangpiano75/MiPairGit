# MiPair

Title: Integrative Web Cloud Computing and Analytics Using MiPair for Design-based Comparative Analysis with Paired Microbiome Data

Version: 1.0.0

Maintainer: Hyo Jung Jang <hyojung.jang@stonybrook.edu>

Description: MiPair is an integrative web cloud service for design-based comparative analysis with paired microbiome data. Pairing (or blocking) is a design technique that is widely used in comparative microbiome studies to efficiently control for the effects of potential confounders (e.g., genetic, environmental, or behavioral factors). Some typical paired (block) designs for human microbiome studies are repeated measures designs that profile each subject's microbiome twice (or more than twice) 1) for pre and post treatments to see the effects of a treatment on microbiome, or 2) for different organs of the body (e.g., gut, mouth, skin) to see the disparity in microbiome between (or across) organs. MiPair enables comprehensive comparative analysis in sequence for such paired microbiome studies on user-friendly web environments. Detailed features are as follows.

* A variety of data uploading, quality controlling, analytic and graphical procedures that produce publishable data, tables, and plots
* Comparative analysis between (or across) groups
* Comparative analysis between baseline (or reference) and other groups
* Parametric or non-parametric tests for incomplete or complete block designs
* Both ecological (alpha- and beta-diversity) and taxonomic (e.g., phylum, class, order, family, genus, species) analysis

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'gridExtra', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip'); GitHub ('LDM', 'volcano3D')

License: GPL 1, GPL 2 (GPL ≥ 2)

## URLs

* Web application (online implementation): http://mipair.micloud.kr
* GitHub repository (local implementation): https://github.com/yj7599/MiPairGit

## References

* Jang HJ, Koh H, Gu W, Kang B.  Integrative Web Cloud Computing and Analytics Using MiPair for Design-based Comparative Analysis with Paired Microbiome Data (*_in revision_*). 

## Prerequites

shiny
```
install.packages('shiny')
```
BiocManager
```
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
    
BiocManager::install(version = '3.14')
```
CRAN
```
cran.pkgs <- c('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'gridExtra', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip')

install.packages(cran.pkgs)
```
Bioconductor
```
library(BiocManager)

bio.pkgs <- c('BiocParallel', 'biomformat', 'phyloseq')

BiocManager::install(bio.pkgs)
```
GitHub
```
library(devtools)

remotes::install_github('KatrionaGoldmann/volcano3D')
```

# Launch App

```
library(shiny)

runGitHub('MiPairGit', 'yj7599', ref = 'main')
```

# Troubleshooting Tips

If you have any problems for using MiPair, please report in Issues (https://github.com/YJ7599/MiPairGit/issues) or email Hyo Jung Jang (hyojung.jang@stonybrook.edu).

* Tip 1. For the local implementation, depending on your pre-installed R libraries, your may need to some additional R packages. 
* Tip 2. For the local implementation, please make sure if you have the most recent package version for the local implementation
