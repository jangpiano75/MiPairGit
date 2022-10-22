# MiPair

Title: MiPair: An Integrative Web Cloud Service for Design-based Comparative Analysis with Paired Microbiome Data

Version: 1.0.0

Maintainer: Hyo Jung Jang <hyojung.jang@stonybrook.edu>

Description: MiPair is an integrative web cloud service for design-based comparative analysis with paired microbiome data. Pairing (or blocking) is a design technique that is widely used in comparative microbiome studies to efficiently control for the effects of potential confounders (e.g., genetic, environmental, or behavioral factors). Some typical paired (block) designs for human microbiome studies are repeated measures designs that profile each subjectâ€™s microbiome twice (or more than twice) 1) for pre and post treatments to see the effects of a treatment on microbiome, or 2) for different organs of the body (e.g., gut, mouse, skin) to see the disparity in microbiome between (or across) organs. MiPair enables comprehensive comparative analysis in sequence for such paired microbiome studies on user-friendly web environments. Detailed features are as follows.

* A variety of data uploading, quality controlling, analytic and graphical procedures that produce publishable data, tables, and plots
* Comparative analysis between (or across) groups
* Comparative analysis between baseline (or reference) and other groups
* Parametric or non-parametric tests for incomplete or complete block designs
* Both ecological (alpha- and beta-diversity) and taxonomic (e.g., phylum, class, order, family, genus, species) analysis

![workflow_git](https://user-images.githubusercontent.com/109124970/188030505-b6dcb1ad-a4bb-47ab-a9c5-75deb96e556a.png)

NeedsCompilation: No

Depends: R(>= 4.1.0), 'BiocManager', 'BiocParallel', 'biomformat', 'bios2mds', 'dashboardthemes', 'DESeq2', 'devtools', 'DiagrammeR', 'dplyr', 'DT', 'ecodist', 'entropart', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'google Vis', 'gridExtra', 'gridGraphics', 'GUniFrac', 'htmltools', 'ICSNP', 'LDM', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'phyloseq', 'picante', 'plotly', 'PMCMRplus', 'quantret', 'remotes', 'reticulate', 'rgl', 'robCompositions', 'seqinr', shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'volcano3D', 'xtable', 'zCompositions', 'zip'

License: GPL (≥ 2)

## URLs

* Web application (online implementation): http://223.194.200.100:3838
* GitHub repository (local implementation): https://github.com/yj7599/MiPairGit

## References

* Jang HJ, Koh H, Gu W, Kang B. MiPair: An integrative web cloud service for design-based comparative analysis with paired microbiome data (*_Submitted_*)

# Prerequites

shiny
```
install.packages("shiny")
```

# Launch App

```
library(shiny)

runGitHub("MiPairGit", "yj7599", ref = "main")
```

# Troubleshooting Tips

If you have any problems for using MiPair, please report in Issues (https://github.com/yj7599/MiPair/issues) or email Hyo Jung Jang (hyojung.jang@stonybrook.edu).
