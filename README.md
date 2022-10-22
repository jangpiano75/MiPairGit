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

Depends: R(>= 4.1.0), 'seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'bslib', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis',                       'xtable', 'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'dplyr', 'forestplot', 'zCompositions', 'fossil', 'picante',
                    'entropart', 'lme4', 'lmerTest', 'robCompositions', 'GUniFrac', 'ecodist', 'gridExtra', 'ggplot2', 'patchwork',                               'ggthemes', 'DiagrammeR', 'stringr','devtools', 'reticulate', 'nlme', 'remotes', 'gridGraphics', 'compositions', 'ICSNP',                     'xtable', 'rgl', 'BiocManager', 'PMCMRplus', 'vegan', 'lme4', 'GUniFrac', 'MiRKAT', 'BiocParallel', 'phyloseq',                               'biomformat', 'LDM', 'volcano3D'

NeedsCompilation: No

License: GPL-2

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
