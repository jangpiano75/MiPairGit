rm(list = ls())

install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')

install.packages(c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable'), repos='https://cloud.r-project.org/')
install.packages(c('DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante' ), repos='https://cloud.r-project.org/')
install.packages(c('entropart', 'lme4', 'lmerTest', 'dirmult', 'robustbase', 'BiasedUrn'), repos='https://cloud.r-project.org/')
install.packages(c('CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr'), repos='https://cloud.r-project.org/')
install.packages(c('devtools', 'betareg', 'remotes'), repos='https://cloud.r-project.org/')
install.packages(c('ICSNP', 'BiocManager', 'PMCMRplus', 'vegan'), repos='https://cloud.r-project.org/')
install.packages('robCompositions', repos='https://cloud.r-project.org/')

BiocManager::install('BiocParallel')
BiocManager::install('RcppGSL')

#remotes::install_github('KatrionaGoldmann/volcano3D', force = TRUE)
remotes::install_github('joey711/phyloseq')
remotes::install_github('joey711/biomformat')
remotes::install_github('hk1785/GLMM-MiRKAT')
remotes::install_github('nyiuab/NBZIMM')
remotes::install_github('https://github.com/karoliskoncevicius/matrixTests', repos = NULL)

#RUN R -e "install.packages('https://cran.rstudio.com/bin/macosx/contrib/4.1/Rfast_2.0.6.tgz', repos = NULL)
install.packages(c('gridGraphics', 'compositions'), repos = 'https://cloud.r-project.org/')
install.packages(c('rgl', 'vegan3d', 'pca3d', 'jpeg', 'splitTools', 'survival', 'survminer', 'coin'), repos = 'https://cloud.r-project.org/')
install.packages(c('randomForestSRC', 'kableExtra', 'caret', 'randomForest', 'glmnet'), repos = 'https://cloud.r-project.org/')
#RUN R -e "remotes::install_version('RcppGSL', version = '0.3.11')
remotes::install_version('RcppZiggurat', version = '0.1.6')
remotes::install_version('Rfast', version = '2.0.6')

library(seqinr)
library(shiny)
library(bslib)
library(shinydashboard)
library(dashboardthemes)
library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(MiRKAT)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(zip)
library(ICSNP)
library(xtable)
library(rgl)
library(BiocParallel)    
library(PMCMRplus)
library(vegan)
library(lme4)
library(GUniFrac)
library(ecodist) 
library(DiagrammeR)
library(phangorn)
library(quantreg)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(fossil)
library(picante)
library(entropart)
library(lmerTest)
library(ICSNP)
options(rgl.useNULL = TRUE)

install.packages('/root/Package/LDM_5.0.tar.gz', repos=NULL, type = 'source')
install.packages('/root/Package/volcano3D_2.0.0.tar.gz', repos=NULL, type = 'source', dependencies = TRUE)

library(LDM)
library(volcano3D)
source("Source/MiDataProc.Data.Paired.R")
source("Source/MiDataProc.Alpha.Paired.R")
source("Source/MiDataProc.Alpha.Paired.Mult.R")
source("Source/MiDataProc.Beta.Paired.R")
source("Source/MiDataProc.Beta.Paired.Mult.R")
source("Source/MiDataProc.Taxa.Paired.R")
source("Source/MiDataProc.Taxa.Paired.Mult.R")

# COMMENTS
{      
  TITLE = p("MiPair: An Integrative Web Cloud Service for Design-based Comparative Analysis with Paired Microbiome Data", style = "font-size:18pt")
  HOME_COMMENT = p(strong("MiPair", style = "font-size:15pt"), "is an integrative web cloud service for design-based comparative analysis with paired microbiome data. Pairing (or blocking) is a design technique that is widely used in comparative microbiome studies to efficiently control for the effects of potential confounders (e.g., genetic, environmental, or behavioral factors). Some typical paired (block) designs for human  microbiome  studies  are  repeated  measures  designs  that  profile  each  subject's microbiome twice (or more than twice) 1) for pre and post treatments to see the effects of a treatment on microbiome, or 2) for different organs of the body (e.g., gut, mouse, skin) to see the disparity in microbiome between (or across) organs. MiPair enables comprehensive comparative analysis in sequence for such paired microbiome studies on user-friendly web environments. Detailed features are as follows.", style = "font-size:13pt")
  
  HOME_COMMENT1 = ("A variety of data uploading, quality controlling, analytic and graphical procedures that produce publishable data, tables, and plots")
  HOME_COMMENT2 = ("Comparative analysis between(or across) groups")
  HOME_COMMENT3 = ("Comparative analysis between baseline (or reference) and other groups")
  HOME_COMMENT4 = ("Parametric or non-parametric tests for incomplete or complete block designs")
  HOME_COMMENT5 = ("Both ecological (alpha-and beta-diversity) and taxonomic (e.g., phylum, class, order, family, genus, species) analysis")
  HOME_COMMENT6 = p("Reference: Jang HJ, Koh H, Gu W, Kang B. MiPair: An integrative web cloud service for design-based comparative analysis with paired microbiome data (in review).", style = "font-size:13pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, 
                              metadata/sample information, and phylogenetic tree.", 
                              br(), br(), "Details:", br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                              (row names are feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, 
                              disease status or environmental/behavioral factors, where rows are subjects and columns are variables 
                              (row names are subject IDs, and column names are variable names).", br(), 
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiPair automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiPair will analyze only the matched features and subjects."
                              , style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download three example  data  sets(Zhang  et  al., 2018), ‘Pairs.Rdata', ‘CBD.3Groups.Rdata', and ‘IBD.3Groups.Rdata', for",  br(),
                              "1) a two-group comparison (a baseline group at the time of antibiotic  administration and 2 weeks afterwards)",  br(), 
                              "2) a  three-group comparison (a baseline group at the time of antibiotic administration and 2 weeks and 4 weeks afterwards) based on a complete block design, where every subject contains all possible three levels of baseline, 2 weeks and 4 weeks afterwards", br(), 
                              "3) a three-group comparison (a baseline group at the time of antibiotic administration and 2 weeks and 4 weeks afterwards) based on an  incomplete  block  design, where  not  every  subject  contains all possible three levels of baseline, 2 weeks and 4 weeks afterwards.", br(), br(), "The data are in the unified format, called phyloseq. For more details about 'phyloseq', see ", a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(),
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              
                              strong("1) Two-group comparison.", style = "color:black"), br(), br(), 
                              
                              "> load(file = ‘Pairs.Rdata')", br(), br(), 
                              "> otu.tab <- otu_table(Pairs)", br(),
                              "> tax.tab <- tax_table(Pairs)", br(), 
                              "> sam.dat <- sample_data(Pairs)", br(),
                              "> tree <- phy_tree(Pairs)", br(), br(), 
                              "> sam.dat$Antibiotic # Treatments", br(), 
                              "> sam.dat$MouseID # Block ID", br(), br(), 
                              
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              
                              "> identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              "> identical(rownames(otu.tab), tree$tip.label)", br(), 
                              "> identical(colnames(otu.tab), rownames(sam.dat))", br(), br(), 
                              
                              strong("2) Three-group comparison based on complete block design.", style = "color:black"), br(), br(), 
                              "> load(file = ‘CBD.3Groups.Rdata')", br(), br(),
                              "> otu.tab <- otu_table(CBD.3Groups)", br(), 
                              "> tax.tab <- tax_table(CBD.3Groups)", br(), 
                              "> sam.dat <- sample_data(CBD.3Groups)", br(), 
                              "> tree <- phy_tree(CBD.3Groups)", br(), br(),
                              
                              "> sam.dat$Antibiotic # Treatments", br(),
                              "> sam.dat$MouseID # Block ID", br(), br(), 
                              
                              "You can check if the features are matched and identical across feature table, taxonomic tableand phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              
                              "> identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              "> identical(rownames(otu.tab), tree$tip.label)", br(),
                              "> identical(colnames(otu.tab), rownames(sam.dat))", br(), br(), 
                              
                              strong("3) Three-group comparison based on incomplete block design.", style = "color:black"), br(), br(),
                              "> load(file = ‘IBD.3Groups.Rdata')", br(), br(), 
                              "> otu.tab <- otu_table(IBD.3Groups)", br(), 
                              "> tax.tab <- tax_table(IBD.3Groups)", br(), 
                              "> sam.dat <- sample_data(IBD.3Groups)", br(), 
                              "> tree <- phy_tree(IBD.3Groups)", br(), br(), 
                              
                              "> sam.dat$Antibiotic # Treatments", br(), 
                              "> sam.dat$MouseID # Block ID", br(), br(), 
                              
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              
                              "> identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              "> identical(rownames(otu.tab), tree$tip.label)", br(), 
                              "> identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                                    (row names are feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME.", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                                    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME.", br(), 
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, 
                                    disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and 
                                    column names are variable names).", br(), 
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiPair automatically roots the tree through midpoint 
                                    rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiPair will analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, feature table (otu.tab.txt), 
                                     taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     "> sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A kingdom of interest. Bacteria (default) for 16S data, Fungi for ITS data, or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove samples that have low library sizes (total read counts). ", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per sample", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative B-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Transform the data into five different formats 1) count, 2) count (rarefied), 3) proportion, 
                             4) CLR (centered log ratio) (Aitchison, 1982) for each taxonomic rank (phylum, class, order, familiy, genus, species), 5) arcsine (arcsine square root) () for each taxonomic rank (phylum, class, order, familiy, genus, species)." )
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
}


# UI
{
  ui = dashboardPage(
    skin = "yellow",
    title = "MiPair",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      tags$head(tags$style(HTML('.progress-bar {background-color: Tomato;}'))),
      setSliderColor(rep("Orange", 100), seq(1, 100)),
      
      chooseSliderSkin("Flat"),
      
      sidebarMenu(
        id = "side_menu",
        menuItem("Home", tabName = "home", icon = icon("home")),
        menuItem("Data Processing",  icon = icon("file-text-o"),
                 menuSubItem(span("Data Input",style = "font-size: 13px"), tabName = "step1", icon = icon("mouse")),
                 menuSubItem(span("Quality Control",style = "font-size: 13px"), tabName = "step2", icon = icon("chart-bar"))),
        menuItem("Ecological Analysis",  icon = icon("chart-pie"),
                 menuSubItem(span("Diversity Calculation", style = "font-size: 13px"), tabName = "divCalculation", icon = icon("calculator")),
                 menuSubItem(span("Alpha Diversity", style = "font-size: 13px"), tabName = "alphaDivanalysis", icon = icon("font")),
                 menuSubItem(span("Beta Diversity", style = "font-size: 13px"), tabName = "betaDivanalysis", icon = icon("bold"))),
        menuItem("Taxonomic Analysis", icon = icon("disease"),
                 menuSubItem(span("Data Transformation", style = "font-size: 12.5px"), tabName = "dataTransform", icon = icon("th-large")),
                 menuSubItem(span("Differential Abundance Analysis", style = "font-size: 12.5px"), tabName = "taxaAnalysis", icon = icon("align-left"))))),
    dashboardBody(
      
      
      tags$style(
        
        ".pretty {
      white-space: normal;
      margin-right: 1px; 
      margin-bottom: 5px;
    }
    .pretty .state label {
      line-height: 1.5em;
      margin-top: -4px;
    }
    .pretty .state label::after, .pretty .state label::before {
      top: -2px;
    }"
      ),
      
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      
      
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      uiOutput("themes"),
      tabItems(
        ##### HOME ####
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT,
                    tags$ol(
                      tags$li(HOME_COMMENT1), tags$li(HOME_COMMENT2), tags$li(HOME_COMMENT3), tags$li(HOME_COMMENT4), tags$li(HOME_COMMENT5),
                      style = "font-size:13pt"), 
                    div(tags$img(src="mipair_workflow.png", height = 900, weight = 500), style =  "text-align: center;"),
                    br(), 
                    HOME_COMMENT6, 
                    )),
        
        ##### DATA INPUT ####
        tabItem(tabName = "step1", br(),
                fluidRow(
                  column(width = 6, style='padding-left:+20px',
                         box(
                           width = NULL, status = "warning", solidHeader = TRUE,
                           title = strong("Data Input", style = "color:black"),
                           selectInput("inputOption", h4(strong("Data Type?")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                           div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                           uiOutput("moreOptions"))),
                  column(width = 6, style='padding-left:0px', uiOutput("addDownloadinfo"), uiOutput("data_input_ref"))
                  
                )),
        
        ##### QC ####
        tabItem(tabName = "step2", br(), 
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                               QC_KINGDOM_COMMENT,
                               tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                               tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'), 
                               
                               sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 2000, step = 1000),
                               QC_LIBRARY_SIZE_COMMENT1,
                               QC_LIBRARY_SIZE_COMMENT2,
                               
                               sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1, value = 0.002, step = 0.001,  post  = " %"),
                               QC_MEAN_PROP_COMMENT1,
                               QC_MEAN_PROP_COMMENT2,
                               
                               br(),
                               p(" ", style = "margin-bottom: -20px;"),
                               
                               h4(strong("Erroneous Taxonomic Names?")),
                               textInput("rem.str", label = "Complete Match", value = ""),
                               QC_TAXA_NAME_COMMENT1,
                               
                               textInput("part.rem.str", label = "Partial Match", value = ""),
                               QC_TAXA_NAME_COMMENT2,
                               
                               actionButton("run", (strong("Run!")), class = "btn-warning btn-info"), br(), br(),
                               uiOutput("moreControls")),
                  mainPanel(width = 9,
                            fluidRow(width = 12,
                                     status = "warning", solidHeader = TRUE, 
                                     valueBoxOutput("sample_Size", width = 3),
                                     valueBoxOutput("OTUs_Size", width = 3),
                                     valueBoxOutput("phyla", width = 3),
                                     valueBoxOutput("classes", width = 3)),
                            fluidRow(width = 12, 
                                     status = "warning", solidHeader = TRUE,
                                     valueBoxOutput("orders", width = 3),
                                     valueBoxOutput("families", width = 3),
                                     valueBoxOutput("genera", width = 3),
                                     valueBoxOutput("species", width = 3)),
                            fluidRow(style = "position:relative",
                                     tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist"),
                                                     sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot", 
                                                     plotlyOutput("boxplot"))),
                                     tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist2"),
                                                     sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot",
                                                     plotlyOutput("boxplot2"))))))),
        
        ##### DIVERSITY Calculation ####
        tabItem(tabName = "divCalculation", br(),
                column(width = 6, style = 'padding-left:0px',
                       box(title = strong("Diversity Calculation", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                           ALPHA_COMMENT, 
                           BETA_COMMENT, 
                           actionButton("divCalcRun", (strong("Run!")), class = "btn-warning btn-info"),
                           p(" ", style = "margin-bottom: +10px;"),
                           p("You have to click this Run button to perform following ecological (alpha-and beta-diversity) analyses.")
                       ),
                       uiOutput("divCalcDownload")),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                           p("Alpha Diversity", style = "font-size:12pt"),
                           ALPHA_REFERENCES,
                           p("Beta Diversity", style = "font-size:12pt"),
                           BETA_REFERENCES))),
        
        ##### ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysis", br(),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               uiOutput("primvars"), 
                               uiOutput("select_lev_alpha"),
                               uiOutput("prim_vars_types"), 
                               p(" ", style = "margin-bottom: +20px;"),
                               uiOutput("blockid"),
                               p(" ", style = "margin-bottom: +10px;"),
                               uiOutput("chooseTest"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("adjust_method"), br(),  
                               uiOutput("alpha_downloadTable"),
                               uiOutput("alpha_references")),
                  mainPanel(width = 9,
                            fluidRow(width = 9, 
                                     uiOutput("alpha_display_results")),
                            
                            uiOutput("barPanel"), br(), br(), br()
                  ))),
        
        ##### BETA DIVERSITY ####
        tabItem(tabName = "betaDivanalysis", br(),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               uiOutput("beta_primvar"),
                               uiOutput("select_lev_beta"),
                               uiOutput("beta_prim_vars_types"), 
                               p(" ", style = "margin-bottom: +20px;"),
                               uiOutput("beta_blockid"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("beta_chooseTest"), 
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("beta_adjust_method"), br(), 
                               uiOutput("beta_downloadTable"),
                               uiOutput("beta_references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("beta_display_results_cross")),
                            uiOutput("beta_barPanel"), br(), br(), br()))),
        
        ##### Data Transformation ####
        tabItem(tabName = "dataTransform", br(),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("Data Transformation", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                           DATA_TRANSFORM_COMMENT,
                           actionButton("datTransRun", (strong("Run!")), class = "btn-warning btn-info"), #class = "btn-info"
                           p(" ", style = "margin-bottom: +10px;"),
                           p("You have to click this Run button to perform following taxonomic differential abundance analyses.")
                       ),
                       uiOutput("datTransDownload")),
                column(width = 6, style='padding-left:0px', 
                       box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                           DATA_TRANSFORM_REFERENCE))),
        
        
        
        ##### Taxa Analysis ####
        tabItem(tabName = "taxaAnalysis", br(),
                sidebarLayout( 
                  sidebarPanel(width = 3,
                               uiOutput("primvars_taxa"),
                               uiOutput("select_levels"), 
                               p(" ", style = "margin-bottom: +20px;"),
                               uiOutput("morePrimvar_opt_taxa"), 
                               p(" ", style = "margin-bottom: +20px;"),
                               uiOutput("blockid_taxa"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("chooseTest_taxa"), 
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("taxa_adjust"),
                               uiOutput("lmm_level"), br(),   
                               uiOutput("downloadTable_taxa"),
                               uiOutput("taxa_references")),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     div(style='height:600px;overflow-y: scroll;', uiOutput("taxa_display")), br(),br(),
                                     uiOutput("taxa_display_dend")),
                            fluidRow(width = 12, 
                                     div(style = 'height:600px;overflow-y: scroll;', uiOutput("taxa_display_pairwise"))), br(), br(), 
                            uiOutput("taxa_barPanel"), 
                            br(), br(), br())))
      )
    )
  )}


# Server
server = function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  load(file="Data/IBD.3Groups.Rdata")
  load(file="Data/CBD.3Groups.Rdata")
  load(file="Data/Pairs.Rdata")
  
  biom_1 <- Pairs
  biom_2 <- CBD.3Groups
  biom_3 <- IBD.3Groups
  
  output$downloadData_1 <- downloadHandler(
    filename = function() {
      paste("Pairs",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom_1, file = file1)
    })
  
  output$downloadData_2 <- downloadHandler(
    filename = function() {
      paste("CBD.3Groups",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom_2, file = file1)
    })
  
  output$downloadData_3 <- downloadHandler(
    filename = function() {
      paste("IBD.3Groups",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom_3, file = file1)
    })
  
  
  ori.biom_1 <- biom_1
  otu.tab_1 <- otu_table(ori.biom_1)
  tax.tab_1 <- tax_table(ori.biom_1)
  tree_1 <- phy_tree(ori.biom_1)
  sam.dat_1 <- sample_data(ori.biom_1)
  
  
  output$downloadZip_1 <- downloadHandler(
    filename = function() {
      paste("Pairs",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab_1, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab_1, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat_1, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree_1, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ori.biom_2 <- biom_2
  otu.tab_2 <- otu_table(ori.biom_2)
  tax.tab_2 <- tax_table(ori.biom_2)
  tree_2 <- phy_tree(ori.biom_2)
  sam.dat_2 <- sample_data(ori.biom_2)
  
  output$downloadZip_2 <- downloadHandler(
    filename = function() {
      paste("CBD.3groups",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab_2, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab_2, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat_2, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree_2, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ori.biom_3 <- biom_3
  otu.tab_3 <- otu_table(ori.biom_3)
  tax.tab_3 <- tax_table(ori.biom_3)
  tree_3 <- phy_tree(ori.biom_3)
  sam.dat_3 <- sample_data(ori.biom_3)
  
  
  output$downloadZip_3 <- downloadHandler(
    filename = function() {
      paste("IBD.3Groups",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab_3, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab_3, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat_3, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree_3, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ## variable define ####
  infile = reactiveValues(biom = NULL, qc_biom = NULL, qc_biomNA = NULL, rare_biom = NULL, rare_biomNA = NULL, is.mon = NULL, ori.var = NULL)
  #infile = reactiveValues(biom = NULL, qc_biom = NULL, rare_biom = NULL)
  ds.Ks = reactiveValues(res = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, prim_vars = NULL, alpha.div = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, taxa.outNA = NULL, tax.tab = NULL, tax.tabNA = NULL)
 
  is.results = reactiveValues(result = NULL)
  is.results.long = reactiveValues(result = NULL)
  multi.test = reactiveValues(boolval = FALSE)
  multi.test.long = reactiveValues(boolval = FALSE)
  alpha.categos = reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.categos.long = reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.data.results = reactiveValues(table.out = NULL, data.q.out = NULL, table.p.out = NULL)
  alpha.results = reactiveValues(bin.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.results.cont = reactiveValues(alpha.con.out = NULL, alpha.table.out = NULL)
  alpha.reg.results = reactiveValues(bin.var = NULL, cov.var = NULL, alpha.div = NULL)
  alpha.resultslong = reactiveValues(bin.var = NULL, id.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.noncovs_res = reactiveValues(con.var = NULL, id.var = NULL, alpha.div = NULL)
  data.alphaBin_res = reactiveValues(table.output = NULL, data.output = NULL, alpha.bin.sum.out = NULL, table.p_outbin = NULL)
  data.results.cont_long = reactiveValues(table.out = NULL, data.q.out = NULL, table_p.out = NULL, alpha.table.out = NULL)
  
  beta.data.results = reactiveValues(data.q.out = NULL)
  beta.results = reactiveValues(result = NULL)
  beta.resultscont = reactiveValues(beta.cont.out = NULL)  
  beta.data.results_long = reactiveValues(beta.bin.out = NULL)
  beta.resultscon_long = reactiveValues(beta.con.out = NULL)
  beta.categos = reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.categos.long = reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.down.results = reactiveValues(CS = NULL, LONG = NULL)
  
  taxa.categos = reactiveValues(cat1 = NULL, cat2 = NULL)
  taxa.data.results = reactiveValues(data.q.out = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL, con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  taxa.types = reactiveValues(dataType = NULL, regression = NULL)
  taxa.outputs = reactiveValues(DAoutput = NULL, DAoutput_or = NULL, DAoutputlong = NULL)
  
  rcol = reactiveValues(selected = "lightblue") 
  

  ## options to input ####
  observeEvent(input$inputOption,{
    observe({
      if (input$inputOption == "Phyloseq") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.Rdata, .rds)", style = "color:black"), 
                      accept = c(".Rdata", ".rds"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-warning btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a Rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                downloadButton("downloadData_1", "Pairs, 2 groups", width = '30%', style = "background-color: red2"),
                downloadButton("downloadData_2", "CBD, 3 groups", width = '30%', style = "background-color: red2"), 
                downloadButton("downloadData_3", "IBD, 3 groups", width = '30%', style = "background-color: red2"), 
                
                br(),br(),
                INPUT_PHYLOSEQ_COMMENT2
            )
          )
        })
        
        output$data_input_ref <- renderUI({
          tagList(
            box(title = strong("Reference"), width = NULL, status = "warning", solidHeader = TRUE, 
                p("Zhang XS, Li J, Krautkramer KA, Badri M, Battaglia T, Borbet TC, et al. Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity. elife. 2018;7:e37816.")
            )
          )
        })
        
        
      } else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Individual_Data', 'Upload', class = "btn-warning btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                downloadButton("downloadZip_1", "Pairs, 2 groups", width = '30%', style = "color:black; background-color: red2"),
                downloadButton("downloadZip_2", "CBD, 3 groups", width = '30%', style = "color:black; background-color: red2"), 
                downloadButton("downloadZip_3", "IBD, 3 groups", width = '30%', style = "color:black; background-color: red2"),
                
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2
            )
          )
        })
        
        output$data_input_ref <- renderUI({
          tagList(
            box(title = strong("Reference"), width = NULL, status = "warning", solidHeader = TRUE, 
                p("Zhang XS, Li J, Krautkramer KA, Badri M, Battaglia T, Borbet TC, et al. Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity. elife. 2018;7:e37816.")
            )
          )
        })
      }
    })
    
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  observe({
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$taxTable) | is.null(input$samData) | is.null(input$tree)))
    toggleState("run", !is.null(infile$biom))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    
    toggleState("divCalcRun", !is.null(infile$rare_biom))
    toggleState("datTransRun", !is.null(infile$rare_biom))
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "lightblue"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        if (!is.null(input$otuTable) && !is.null(input$taxTable) && !is.null(input$samData) && !is.null(input$tree)) {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, tax.table, sam.data, tree.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") && (ext2 == "txt" | ext2 == "tsv") &&
                (ext3 == "txt" | ext3 == "csv") && (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu(otu.tab, tax.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } else {
          return(NULL)
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "lightblue"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ######################################
  # Quality control and transformation #
  ######################################
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    num_tax.rank = reactive({
      tax.tab = tax_table(infile$qc_biom)
      num.tax.rank(tax.tab)
    })
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "fuchsia")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    output$phyla <- renderValueBox({
      num.phyla = num_tax.rank()[1]
      valueBox(
        value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
        "Number of Phyla", icon = icon("sitemap"), color = "orange")
    })
    
    output$classes <- renderValueBox({
      num.classes = num_tax.rank()[2]
      valueBox(
        value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
        "Number of Classes", icon = icon("sitemap"), color = "purple")
    })
    
    output$orders <- renderValueBox({
      num.orders = num_tax.rank()[3]
      valueBox(
        value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
        "Number of Orders", icon = icon("sitemap"), color = "blue")
    })
    
    output$families <- renderValueBox({
      num.families = num_tax.rank()[4]
      valueBox(
        value = tags$p(paste0(num.families), style = "font-size: 75%;"),
        "Number of Families", icon = icon("sitemap"), color = "red")
    })
    
    output$genera <- renderValueBox({
      num.genera = num_tax.rank()[5]
      valueBox(
        value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
        "Number of Genera", icon = icon("sitemap"), color = "lime")
    })
    
    output$species <- renderValueBox({
      num.species = num_tax.rank()[6]
      valueBox(
        value = tags$p(paste0(num.species), style = "font-size: 75%;"),
        "Number of Species", icon = icon("sitemap"), color = "teal" )
    })
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    max.mean.prop = as.numeric(mean.prop.func(infile$qc_biom)$mean.prop.sum["3rd quartile"])
    maxi.slider2 = round(max.mean.prop, digits = 6)
    
    if (maxi.slider2 < 2e-05) {
      maxi.slider2 = 2e-05
    }
    
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
    updateSliderInput(session, "slider2", min = 0, max = maxi.slider2*100)
  })
  
  ######################################
  ######################################
  #########  Data Analysis   ###########
  ######################################
  ######################################
  observeEvent(chooseData$alpha.div,{
    
    ######################################
    ######################################
    ######### Alpha diversity ############
    ######################################
    ######################################
    
    output$primvars <- renderUI({
      tagList(
        h4(strong("Primary Variable?", style = "color:black")),
        p("A factor variable that contains multiple groups/levels of treatments or body sites.", style = "font-size:10pt"), 
        
        selectInput("primvar", label = NULL,
                    c("Choose one" = "", sort(chooseData$prim_vars)), selected = sort(chooseData$prim_vars)[1], width = '70%'))
    })
    
    output$prim_vars_types <- renderUI({
      tagList(
        uiOutput("morePrimvar_opt"))
    })
    
    
    observeEvent(input$primvar, {
      
      prim.ind <- which(names(chooseData$sam.dat) == input$primvar)
      
      output$blockid <- renderUI({
        tagList(
          h4(strong("Pair/Block ID?", style = "color:black")),
          p("Pair/Block IDs are, for example, subjects IDs for pre and post treatments or body sites.", style = "font-size:10pt"),
          selectInput("blockid", label= NULL, 
                      choices = c("Choose one" = "", names(chooseData$sam.dat)[-prim.ind]), selected = names(chooseData$sam.dat)[-prim.ind][1], width = '70%'))   
      })
      
      
      if (length(names(table(chooseData$sam.dat[, input$primvar]))) == 2){
        shinyjs::hide("adjust_method")
        shinyjs::hide("select_lev_alpha")
        
        alpha.categos$cat1 =  alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[1]
        alpha.categos$cat2 = alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[2]
        
        output$morePrimvar_opt <- renderUI({
          tagList(
            h4(strong("Rename Categories?", style = "color:black")),
            p("You can rename categories of the primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"), 
            textInput("alphaCat1", label = (paste0("Reference: ", alpha.categos$cat1)), value = alpha.categos$cat1, width = '80%'),
            textInput("alphaCat2", label = (paste0("Comparison: ", alpha.categos$cat2)), value = alpha.categos$cat2, width = '80%'))
        }) 
        
        output$chooseTest <- renderUI({
          tagList(
            prettyRadioButtons("chooseMethod", label = h4(strong("Method?", style = "color:black")), 
                               c("Paired t-test", "Wilcoxon signed-rank test", "Multivariate Hotelling's t-squared test"), shape = c("round"), selected = "Paired t-test", width = '80%'),
            actionButton("runbtn_bin", (strong("Run!")), class = "btn-warning btn-info"))
        })
      }else{
        shinyjs::show("adjust_method")
        shinyjs::show("select_lev_alpha")
        
        output$select_lev_alpha <- renderUI({
          tagList(
            h4(strong("Groups/Levels?", style = "color:black")), 
            p("Select at least two groups/levels to be compared.", style = "font-size:10pt"),
            p(" ", style = "margin-bottom: +20px;"),
            checkboxGroupInput("alpha_levels", label = NULL, choices = names(table(chooseData$sam.dat[, input$primvar])), selected = c(names(table(chooseData$sam.dat[, input$primvar]))[1], names(table(chooseData$sam.dat[, input$primvar]))[2])
            )
          )
        })
        
        observeEvent(input$alpha_levels, {
          
          if (length(input$alpha_levels) == 2){
            shinyjs::hide("runbtn_bin_mult")
            
            alpha.categos$cat1 = try(alpha.bin.cat.func(reduced_data(chooseData$sam.dat, input$primvar, input$alpha_levels), input$primvar)[1], silent = TRUE)
            alpha.categos$cat2 = try(alpha.bin.cat.func(reduced_data(chooseData$sam.dat, input$primvar, input$alpha_levels), input$primvar)[2], silent = TRUE)
            
            validate(
              need(alpha.categos$cat1 %in% names(table(chooseData$sam.dat[,input$primvar])), "") 
            )
            
            categos_alpha <- c(alpha.categos$cat1, alpha.categos$cat2)
            
            output$morePrimvar_opt <- renderUI({
              tagList(
                h4(strong("Rename Categories?", style = "color:black")),
                p("You can rename groups/levels of the chosen primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"), 
                try(textInput("alphaCat1", label = (paste0("Reference: ", categos_alpha[1])), value = categos_alpha[1], width = '80%'), silent = TRUE),
                try(textInput("alphaCat2", label = (paste0("Comparison: ", categos_alpha[2])), value = categos_alpha[2], width = '80%'), silent = TRUE))
            }) 
            
            output$chooseTest <- renderUI({
              tagList(
                prettyRadioButtons("chooseMethod", label = h4(strong("Method?", style = "color:black")), 
                                   c("Paired t-test", "Wilcoxon signed-rank test", "Multivariate Hotelling's t-squared test"), shape = "round", selected = "Paired t-test", width = '80%'),
                actionButton("runbtn_bin", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }
          
          else if (length(input$alpha_levels) > 2){
            shinyjs::hide("runbtn_bin")
            
            prim_length_alpha <- length(input$alpha_levels)
            alpha.categos <- input$alpha_levels
            re_sam_alpha_now <- reduced_data(chooseData$sam.dat, input$primvar, input$alpha_levels)
            
            output$morePrimvar_opt <- renderUI({
              tagList(
                h4(strong("Rename Categories?", style = "color:black")),
                p("You can rename groups/levels of the chosen primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"),
                lapply(1:prim_length_alpha, function(i){
                  textInput(paste0("alphaCat_", i), label = (paste0("Group/Level ", i, ": ", alpha.categos[i])), value = alpha.categos[i], width = '80%')
                }))
            }) 
            
            
            ntselected.prim_vars_alpha = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar)
            
            
            if (length(unique(as.numeric(table(re_sam_alpha_now[, input$primvar])))) == 1) {  
              output$chooseTest <- renderUI({
                
                tags$div(
                  style = "width: 250px",
                  h4(strong("Method?", style = "color:black")), 
                  prettyRadioButtons_new("chooseMethod", label = NULL, 
                                         c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LMM (LRT (global) with t-test (pairwise))"), selected = "ANOVA F-test (global) with Tukey's HSD (pairwise)", status = "warning", width = '80%')
                )
              })
            } else {
              output$chooseTest <- renderUI({
                
                tagList(
                  h4(strong("Method?", style = "color:black")), 
                  p(" ", style = "margin-bottom: +10px;"),
                  prettyRadioButtons_new("chooseMethod", label = NULL, 
                                         c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LMM (LRT (global) with t-test (pairwise))"), selected = "ANOVA F-test (global) with Tukey's HSD (pairwise)", width = '80%')
                )
              })
            }
            
            observeEvent(input$chooseMethod,{
              output$adjust_method <- renderUI({
                if (input$chooseMethod == "Friedman's test (global) with Conover's test (pairwise)"  | input$chooseMethod == "ANOVA F-test (global) with Tukey's HSD (pairwise)" | input$chooseMethod == "Durbin's test (global) with Conover's test (pairwise)"){
                  tagList(
                    p(" ", style = "margin-bottom: +20px;"),
                    actionButton("runbtn_bin_mult", (strong("Run!")), class = "btn-warning btn-info")
                  )
                }else if(input$chooseMethod == "LMM (LRT (global) with t-test (pairwise))"){
                  tagList(
                    h4(strong("Baseline/Reference?", style = "color:black")),
                    p("You need to set up a baseline/reference group to be compared with the other groups/levels in the chosen primary variable.", style = "font-size:10pt"),
                    p(" ", style = "margin-bottom: +20px;"),
                    radioButtons("lmm_level_alpha", label = NULL, alpha.categos),
                    actionButton("runbtn_bin_mult", (strong("Run!")), class = "btn-warning btn-info")
                  )
                }
              })
            })
          }
        })
      }
    })
  })
  
  observeEvent(ds.Ks$res, {
  ######################################
  ######################################
  ########## Beta diversity ############
  ######################################
  ######################################
  output$beta_primvar <- renderUI({
    tagList(
      
      h4(strong("Primary Variable?", style = "color:black")),
      p("A factor variable that contains multiple groups/levels of treatments or body sites.", style = "font-size:10pt"), 
      
      selectInput("beta.primvar", label = NULL,
                  c("Choose one" = "", sort(chooseData$prim_vars)), selected = sort(chooseData$prim_vars)[1], width = '70%'))
  })
  
  output$beta_prim_vars_types <- renderUI({
    tagList(
      uiOutput("beta_morePrimvar"))
  })
  
  observeEvent(input$beta.primvar,{
    
    prim.ind <- which(names(chooseData$sam.dat) == input$beta.primvar)
    
    output$beta_blockid <- renderUI({
      tagList(
        h4(strong("Pair/Block ID?", style = "color:black")),
        p("Pair/Block IDs are, for example, subjects IDs for pre and post treatments or body sites. ", style = "font-size:10pt"),
        selectInput("beta.blockid", label = NULL,
                    choices = c("Choose one" = "", names(chooseData$sam.dat)[-prim.ind]), selected = names(chooseData$sam.dat)[-prim.ind][1], width = '70%'))
    })
    
    if (length(names(table(chooseData$sam.dat[, input$beta.primvar]))) == 2){
      
      shinyjs::hide("beta_adjustment")
      shinyjs::hide("beta_runbtn_mult")
      shinyjs::hide("select_lev_beta")
      shinyjs::hide("")
      
      beta.categos$cat1 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar)[1]
      beta.categos$cat2 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar)[2]
      
      output$beta_morePrimvar <- renderUI({
        tagList(
          h4(strong("Rename Categories?", style = "color:black")),
          p("You  can  rename  groups/levels  of  the  chosen  primary  variable.  MiPair  keeps  up  to  8 characters on graphs.", style = "font-size:10pt"), 
          textInput("betaCat1", label = (paste0("Reference: ", beta.categos$cat1)), value = beta.categos$cat1, width = '80%'),
          textInput("betaCat2", label = (paste0("Comparison: ", beta.categos$cat2)), value = beta.categos$cat2, width = '80%'))
      }) 
      
      output$beta_chooseTest <- renderUI({
        tagList(
          prettyRadioButtons("beta_chooseMethod", label = h4(strong("Method?", style = "color:black")),
                             c("PERMANOVA"), selected = "PERMANOVA", width = '80%'),
          p(" ", style = "margin-bottom: -10px;"),
          actionButton("beta_runbtn", (strong("Run!")), class = "btn-warning btn-info"))
      })
    }else{
      shinyjs::show("beta_adjustment")
      shinyjs::show("beta_runbtn_mult")
      shinyjs::show("select_lev_beta")
      
      output$select_lev_beta <- renderUI({
        tagList(
          h4(strong("Groups/Levels?", style = "color:black")), 
          p("Select at least two groups/levels to be compared.", style = "font-size:10pt"),
          p(" ", style = "margin-bottom: +20px;"),
          checkboxGroupInput("beta_levels", label = NULL, choices = names(table(chooseData$sam.dat[, input$beta.primvar])), selected = c(names(table(chooseData$sam.dat[, input$beta.primvar]))[1], names(table(chooseData$sam.dat[, input$beta.primvar]))[2])
          )
        )
      })
      
      observeEvent(input$beta_levels, {
        if (length(input$beta_levels) == 2){
          shinyjs::show("beta_runbtn")
          shinyjs::hide("beta_runbtn_mult")
          shinyjs::hide("base_level_beta")
          shinyjs::hide("beta_adjust_method")
          
          beta.categos$cat1 = try(beta.bin.cat.func(reduced_data(chooseData$sam.dat, input$beta.primvar ,input$beta_levels), input$beta.primvar)[1], silent = TRUE) 
          beta.categos$cat2 = try(beta.bin.cat.func(reduced_data(chooseData$sam.dat, input$beta.primvar ,input$beta_levels), input$beta.primvar)[2], silent = TRUE) 
          
          validate(
            need(beta.categos$cat1 %in% names(table(chooseData$sam.dat[,input$beta.primvar])), "") 
          )
          
          categos_beta <- c(beta.categos$cat1, beta.categos$cat2)
          
          output$beta_morePrimvar <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You can rename groups/levels of the chosen primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"), 
              textInput("betaCat1", label = (paste0("Reference: ", categos_beta[1])), value = categos_beta[1], width = '80%'),
              textInput("betaCat2", label = (paste0("Comparison: ", categos_beta[2])), value = categos_beta[2], width = '80%'))
          }) 
          
          output$beta_chooseTest <- renderUI({
            tagList(
              prettyRadioButtons("beta_chooseMethod", label = h4(strong("Method?", style = "color:black")),
                                 c("PERMANOVA"), selected = "PERMANOVA", width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              actionButton("beta_runbtn", (strong("Run!")), class = "btn-warning btn-info"))
          })
          
        } else if (length(input$beta_levels) > 2){
          shinyjs::hide("beta_runbtn")
          shinyjs::show("base_level_beta")
          shinyjs::show("beta_adjust_method")
          
          prim_length_beta = length(input$beta_levels) 
          beta.categos <- input$beta_levels
          
          output$beta_morePrimvar <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You  can  rename  groups/levels  of  the  chosen  primary  variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:prim_length_beta, function(i){
                textInput(paste0("betaCat_", i), label = (paste0("Group/Level ", i, ": ", beta.categos[i])), value = beta.categos[i], width = '80%')
              }))
          }) 
          
          output$beta_chooseTest <- renderUI({
            tagList(
              prettyRadioButtons_beta_new("beta_chooseMethod", label = h4(strong("Method?", style = "color:black")),
                                          choiceNames = c("PERMANOVA", "PERMANOVA"), choiceValues = c("PERMANOVA_ACROSS", "PERMANOVA_BASE")
                                          , width = '80%'))
          })
          
          observeEvent(input$beta_chooseMethod, {
            if(input$beta_chooseMethod == "PERMANOVA_ACROSS"){
              output$beta_adjust_method <- renderUI({
                tagList(
                  actionButton("beta_runbtn_mult", (strong("Run!")), class = "btn-warning btn-info"))
              })
            } else if(input$beta_chooseMethod == "PERMANOVA_BASE"){
              output$beta_adjust_method <- renderUI({
                tagList(
                  h4(strong("Baseline/Reference?", style = "color:black")), 
                  p("You need to set up a baseline/reference group to be compared with the other groups/levels in the chosen primary variable.", style = "font-size:10pt"),
                  p(" ", style = "margin-bottom: +20px;"),
                  radioButtons("base_level_beta", label = NULL, beta.categos),
                  actionButton("beta_runbtn_mult", (strong("Run!")), class = "btn-warning btn-info"))
              })
            }
          })
          
        }
      })
      
    }
    
    ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar)
    
  })
  })
  observeEvent(chooseData$taxa.out,{
    
    ######################################
    ######################################
    ######## Taxonomic Analysis ##########
    ######################################
    ###################################### 
    output$primvars_taxa <- renderUI({
      tagList(
        prettyRadioButtons("dataType_taxa", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR", "Proportion", "Arcsine-root",  "Count (Rarefied)"), selected = "CLR",width = '70%'),
        p(" ", style = "margin-bottom: -7px;"),
        h4(strong("Primary Variable?", style = "color:black")),
        p("A factor variable that contains multiple groups/levels of treatments or body sites."),
        selectInput("primvar_taxa", label = NULL,
                    choices = chooseData$prim_vars, selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    observeEvent(input$dataType_taxa, {
      
      
      if (input$dataType_taxa == "Count (Rarefied)") {
        taxa.types$dataType = "rare.count"
      } else if (input$dataType_taxa == "CLR") {
        taxa.types$dataType = "clr"
      } else if (input$dataType_taxa == "Arcsine-root") {
        taxa.types$dataType = "arcsin"  
      } else if (input$dataType_taxa == "Proportion"){
        taxa.types$dataType = "imp.prop"
      }
    })
    
    observeEvent(input$primvar_taxa,{
      
      prim.ind <- which(names(chooseData$sam.dat) == input$primvar_taxa)
      
      output$blockid_taxa <- renderUI({
        tagList(
          h4(strong("Pair/Block ID?", style = "color:black")),
          p("Pair/Block IDs are, for example, subjects IDs for pre and post treatments or body sites. ", style = "font-size:10pt"),
          selectInput("blockid_taxa", label = NULL,
                      choices = c("Choose one" = "", names(chooseData$sam.dat)[-prim.ind]), selected = names(chooseData$sam.dat)[-prim.ind][1], width = '70%'))
      })
      
      
      if (length(names(table(chooseData$sam.dat[,input$primvar_taxa]))) == 2){
        shinyjs::hide("select_levels")
        shinyjs::hide("runbtn_bin_taxa")
        shinyjs::hide("taxa_display_pairwise")
        
        taxa.categos$cat1 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[1]
        taxa.categos$cat2 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[2]
        
        output$morePrimvar_opt_taxa <- renderUI({
          tagList(
            h4(strong("Rename Categories?", style = "color:black")),
            p("You can rename categories of the primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"),
            textInput("taxaCat1", label = (paste0("Reference: ",taxa.categos$cat1)), value = taxa.categos$cat1, width = '80%'),
            textInput("taxaCat2", label = (paste0("Comparison: ",taxa.categos$cat2)), value = taxa.categos$cat2, width = '80%'))
        })
        ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
        
        observeEvent(input$dataType_taxa,{
          if (input$dataType_taxa == "Proportion" | input$dataType_taxa == "Arcsine-root"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test", "LDM (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test", "LDM"), selected = "LDM", width = '80%'),
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }else if (input$dataType_taxa == "CLR"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test"), selected = "Wilcoxon signed-rank test", width = '80%'),
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }else if (input$dataType_taxa == "Count (Rarefied)"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test", "LDM (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test", "LDM"), selected = "LDM", width = '80%'),
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }
        })
      }
      
      else if (length(names(table(chooseData$sam.dat[, input$primvar_taxa]))) > 2){
        
        shinyjs::show("select_levels")
        shinyjs::show("runbtn_bin_taxa")
        output$select_levels <- renderUI({
          tagList(
            h4(strong("Groups/Levels?", style = "color:black")), 
            p("Select at least two groups/levels to be compared.", style = "font-size:10pt"),
            p(" ", style = "margin-bottom: +20px;"),
            checkboxGroupInput("taxa_levels", label = NULL, choices = names(table(chooseData$sam.dat[, input$primvar_taxa])), selected = c(names(table(chooseData$sam.dat[, input$primvar_taxa]))[1], names(table(chooseData$sam.dat[, input$primvar_taxa]))[2])
            )
          )
        })
      } 
    })
    observeEvent(input$taxa_levels, {
      if (length(input$taxa_levels) == 2){
        shinyjs::hide("taxa_adjustment")
        shinyjs::hide("lmm_level")
        shinyjs::hide("runbtn_bin_taxa")
        
        #shinyjs::hide("include_species.dend")
        shinyjs::hide("taxa_display_pairwise")
        
        taxa.categos$cat1 = taxa.bin.cat.func(reduced_data(chooseData$sam.dat, input$primvar_taxa ,input$taxa_levels), input$primvar_taxa)[1]
        taxa.categos$cat2 = taxa.bin.cat.func(reduced_data(chooseData$sam.dat, input$primvar_taxa ,input$taxa_levels), input$primvar_taxa)[2]
        
        output$morePrimvar_opt_taxa <- renderUI({
          tagList(
            h4(strong("Rename Categories?", style = "color:black")),
            p("You can rename groups/levels of the chosen primary  variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"),
            textInput("taxaCat1", label = (paste0("Reference: ",taxa.categos$cat1)), value = taxa.categos$cat1, width = '80%'),
            textInput("taxaCat2", label = (paste0("Comparison: ",taxa.categos$cat2)), value = taxa.categos$cat2, width = '80%'))
        })
        ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
        
        
        observeEvent(input$dataType_taxa,{
          if (input$dataType_taxa == "Proportion" | input$dataType_taxa == "Arcsine-root"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test", "LDM (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test", "LDM"), selected = "LDM", width = '80%'),
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }else if (input$dataType_taxa == "CLR"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test"), selected = "Wilcoxon signed-rank test", width = '80%'),
                
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }else if (input$dataType_taxa == "Count (Rarefied)"){
            output$chooseTest_taxa <- renderUI({
              tagList(
                prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                   c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                   shape = c("round"), width = '80%'),
                p(" ", style = "margin-bottom: -10px;"),
                prettyRadioButtons("chooseMethod_taxa_base", label = h4(strong("Method?", style = "color:black")), choiceNames = c("Paired t-test", "Wilcoxon signed-rank test", "LDM (default)"), choiceValues = c("Paired t-test", "Wilcoxon signed-rank test", "LDM"), selected = "LDM", width = '80%'),
                actionButton("runbtn_bin_taxa_base", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }
        })
        
      }else if (length(input$taxa_levels) > 2){
        
        shinyjs::hide("chooseMethod_taxa_base")
        shinyjs::show("taxa_display_pairwise")
        shinyjs::show("runbtn_bin_taxa")
        
        prim_length = length(input$taxa_levels)
        taxa.categos_mult <- input$taxa_levels
        
        re_sam_now <- reduced_data(chooseData$sam.dat, input$primvar_taxa, input$taxa_levels)  
        
        output$morePrimvar_opt_taxa<- renderUI({
          tagList(
            h4(strong("Rename Categories?", style = "color:black")),
            p("You can rename groups/levels of the chosen primary variable. MiPair keeps up to 8 characters on graphs.", style = "font-size:10pt"),
            lapply(1:prim_length, function(i){
              textInput(paste0("mult_taxaCat", i), label = (paste0("Group/Level ", i, ": ", taxa.categos_mult[i])), value = taxa.categos_mult[i], width = '80%')
            }))
        }) 
        
        ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
        
        observeEvent(input$dataType_taxa,{
          if (input$dataType_taxa == "Proportion" | input$dataType_taxa == "Arcsine-root"){
            if(length(unique(as.numeric(table(re_sam_now[, input$primvar_taxa])))) == 1){
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_4("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LDM (default)", "LMM (LRT (global) with t-test (pairwise))"), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LDM", "LMM (LRT (global) with t-test (pairwise))"), selected = "LDM", width = '80%'), 
                  p(" ", style = "margin-bottom: -8px;")
                )})
            } else{
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_4("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LDM (default)", "LMM (LRT (global) with t-test (pairwise))"), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LDM", "LMM (LRT (global) with t-test (pairwise))"), selected = "LDM", width = '80%'), 
                  p(" ", style = "margin-bottom: -8px;")
                )})
            }
            
          } else if(input$dataType_taxa == "Count (Rarefied)"){
            if (length(unique(as.numeric(table(re_sam_now[, input$primvar_taxa])))) == 1){
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_4("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LDM", "LMM (LRT (global) with t-test (pairwise))"), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LDM (default)", "LMM (LRT (global) with t-test (pairwise))"), selected = "LDM", width = '80%'), 
                  p(" ", style = "margin-bottom: -8px;")
                  
                  # p(" ", style = "margin-bottom: -10px;"),
                  # prettyRadioButtons_new("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LDM"), selected = "ANOVA F-test (global) with Tukey's HSD (pairwise)", width = '80%'), p(" ", style = "margin-bottom: -8px;")
                )})
            } else{
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_4("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LDM", "LMM (LRT (global) with t-test (pairwise))"), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LDM (default)", "LMM (LRT (global) with t-test (pairwise))"), selected = "LDM", width = '80%'), 
                  p(" ", style = "margin-bottom: -8px;")
                  
                )})
            }
          } else if (input$dataType_taxa == "CLR"){
            if (length(unique(as.numeric(table(re_sam_now[, input$primvar_taxa])))) == 1){
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_new("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)(default)", "LMM (LRT (global) with t-test (pairwise))"), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Friedman's test (global) with Conover's test (pairwise)", "LMM (LRT (global) with t-test (pairwise))"), selected = "Friedman's test (global) with Conover's test (pairwise)", width = '80%'), p(" ", style = "margin-bottom: -8px;")
                )})
            } else{
              output$chooseTest_taxa <- renderUI({
                tagList(
                  prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                     c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                     shape = c("round"), width = '80%'),
                  p(" ", style = "margin-bottom: -10px;"),
                  prettyRadioButtons_new("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), choiceNames = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)(default)", "LMM (LRT (global) with t-test (pairwise))"), choiceValues = c("ANOVA F-test (global) with Tukey's HSD (pairwise)", "Durbin's test (global) with Conover's test (pairwise)", "LMM (LRT (global) with t-test (pairwise))"), selected = "Durbin's test (global) with Conover's test (pairwise)", width = '80%'), p(" ", style = "margin-bottom: -8px;")
                )})
            }
          }
        })
        
        observeEvent(input$chooseMethod_taxa,{
          
          if (input$chooseMethod_taxa == "ANOVA F-test (global) with Tukey's HSD (pairwise)" | input$chooseMethod_taxa == "Friedman's test (global) with Conover's test (pairwise)" | input$chooseMethod_taxa == "Durbin's test (global) with Conover's test (pairwise)" | input$chooseMethod_taxa == "LDM"){
            output$taxa_adjust <- renderUI({tagList(
              p(" ", style = "margin-bottom: +25px;"),
              actionButton("runbtn_bin_taxa", (strong("Run!")), class = "btn-warning btn-info"))
            })}
          
          else if (input$chooseMethod_taxa == "LMM (LRT (global) with t-test (pairwise))"){
            output$taxa_adjust <- renderUI({
              tagList(
                h4(strong("Baseline/Reference?", style = "color:black")),
                p("You need to set up a baseline/reference group to be compared with the other groups/levels in the chosen primary variable.", style = "font-size:10pt"),
                p(" ", style = "margin-bottom: +20px;"),
                radioButtons("lmm_level", label = NULL, taxa.categos_mult), 
                actionButton("runbtn_bin_taxa", (strong("Run!")), class = "btn-warning btn-info"))
            })
          }
        })
      }
    })
  })
  
  ######################################
  ######################################
  ##########  RUN BUTTONS   ############
  ######################################
  ######################################
  
  ##########################################
  ##                  QC                  ##
  ##########################################
  
  observeEvent(input$run, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Trimming in progress")
        
        if (nchar(input$part.rem.str) == 0) {
          rem.tax.complete <- rem.tax.d
          rem.tax.partial <- rem.tax.str.d
        } else {
          rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
          rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
        }
        
        tax.tab <- tax_table(infile$biom)
        
        if (input$kingdom != "all") {
          ind <- is.element(tax.tab[,1], input$kingdom)
          validate(
            if (sum(ind) == 0) {
              showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                        paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))),
                               type = "error")
            } else {
              NULL
            }
          )
        }
        
        shinyjs::disable("run")
        shinyjs::disable("slider1")
        shinyjs::disable("slider2")
        shinyjs::disable("kingdom")
        shinyjs::disable("skip")
        shinyjs::disable("binwidth")
        shinyjs::disable("binwidth2")
        
        
        rcol$selected = "rgba(255, 0, 0, 0.6)"
        
        infile$qc_biom = biom.clean(infile$biom,  
                                    input$kingdom, 
                                    lib.size.cut.off = input$slider1, 
                                    mean.prop.cut.off = input$slider2/100,
                                    rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
        
        infile$qc_biomNA = biom.cleanSNA(infile$biom,
                                         input$kingdom,
                                         lib.size.cut.off = input$slider1,
                                         mean.prop.cut.off = input$slider2/100,
                                         rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                         surv.time = chooseData$surv.Time.select,
                                         censor = chooseData$censor.select,
                                         follow = chooseData$follow.time)
        
        incProgress(3/10, message = "Rarefying in progress")
        lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
        infile$rare_biom = rarefy.func(infile$qc_biom, 
                                       cut.off = lib_size.sum["Minimum"],
                                       multi.rarefy = 1)
        
        infile$rare_biomNA = rarefy.func(infile$qc_biomNA, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1)
        
        incProgress(2/10, message = "Saving File in progress")
        
        chooseData$sam.dat = sample_data(infile$qc_biom)
        chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.combined(chooseData$sam.dat)
        chooseData$prim_vars = prim_vars_pair(chooseData$sam.dat) #pri.func.combined(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        chooseData$tax.tab = tax_table(infile$rare_biom)
        chooseData$tax.tabNA = tax_table(infile$rare_biomNA)
        
        output$moreControls <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:13pt"),
                h5("Data after Quality Control"),
                downloadButton("downloadData2", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                h5("Data after Quality Control and Rarefaction"),
                downloadButton("downloadData3", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                      (rare.biom.after.qc) after QC and rarefaction.",
                  style = "font-size:11pt")
            )
          )
        })
        
        output$text <- renderText({"You are all set! You can proceed to data analysis!"})
        
        qc_biom = infile$qc_biom
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(qc_biom, file = file1)
          })
        
        rare_biom = infile$rare_biom
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("rare.biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(rare_biom, file = file1)
          })
        incProgress(1/10, message = "Done")
        shinyjs::enable("run")
        shinyjs::enable("slider1")
        shinyjs::enable("slider2")
        shinyjs::enable("kingdom")
        shinyjs::enable("skip")
        shinyjs::enable("binwidth")
        shinyjs::enable("binwidth2")
      })
  })
  
  ##########################################
  ## Data Calculation(Alpha, Beta, Taxa)  ##
  ##########################################
  observeEvent(input$divCalcRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        
        incProgress(3/10, message = "Calculating Diversity")
        
        chooseData$alpha.div.rare = alpha.v1.func(infile$rare_biom)
        chooseData$alpha.div.qc = alpha.v1.func(infile$qc_biom)
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        incProgress(3/10, message = "Calculating Distance")
        
        ds.Ks$res = Ds.Ks.func(infile$rare_biom, infile$qc_biom)
        
        output$divCalcDownload <- renderUI({
          tagList(
            box(title = strong("Download Diversity Data", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download alpha- and beta-diversity data.",
                  style = "font-size:11pt"), 
                h5("Alpha Diversity"),
                downloadButton("alphaDiv", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                h5("Beta Diversity"),
                downloadButton("betaDiv", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        alpha.div = chooseData$alpha.div
        
        output$alphaDiv <- downloadHandler(
          filename = function() {
            paste("Alpha.Diversity.txt")
          },
          content = function(alpha.file) {
            write.table(chooseData$alpha.div, file = alpha.file, row.names = TRUE, col.names = TRUE, sep = "\t")
          })
        
        output$betaDiv <- downloadHandler(
          filename = function() {
            paste("Beta.Diversity.zip")
          },
          content <- function(fname) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Jaccard.txt", "Bray.Curtis.txt", "U.UniFrac.txt" ,"G.UniFrac.txt", "W.UniFrac.txt")
            
            write.table(as.data.frame(ds.Ks$res$Ds$Jaccard), file = "Jaccard.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$Bray.Curtis), file = "Bray.Curtis.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$U.UniFrac), file = "U.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$G.UniFrac), file = "G.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$W.UniFrac), file = "W.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=fname, files=dataFiles)
          }
        )
        incProgress(1/10, message = "Done")
        shinyjs::enable("divCalcRun")
      })
  })
  observeEvent(input$datTransRun, {
    
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        incProgress(3/10, message = "Transformation in progress")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        rare.otu.tabNA <- otu_table(infile$rare_biomNA)
        rare.tax.tabNA <- tax_table(infile$rare_biomNA)
        no.rare.otu.tabNA <- otu_table(infile$qc_biomNA)
        no.rare.tax.tabNA <- tax_table(infile$qc_biomNA)
        
        
        chooseData$taxa.out <- tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        chooseData$taxa.outNA = tax.trans.na(no.rare.otu.tabNA, no.rare.tax.tabNA, rare.otu.tabNA, rare.tax.tabNA)
        chooseData$taxa.names.out <- taxa.names.rank(chooseData$taxa.out[[1]])
        chooseData$tax.tab <- rare.tax.tab
        chooseData$tax.tabNA = rare.tax.tabNA
        
        chooseData$NAadded <- add_NA(chooseData$taxa.outNA, chooseData$tax.tabNA)
        
        
        incProgress(3/10, message = "Transformation in progress")
        
        output$datTransDownload <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download taxonomic abundance data.",
                  style = "font-size:11pt"), 
                h5("Count"),
                downloadButton("taxadataCount", "Download", width = '50%', style = "color:black; background-color: red3"), br(), 
                h5("Count (Rarefied)"),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = "color:black; background-color: red3"), br(), 
                h5("Proportion"),
                downloadButton("taxadataProp", "Download", width = '50%', style = "color:black; background-color: red3"), br(), 
                h5("CLR"),
                downloadButton("taxadataCLR", "Download", width = '50%', style = "color:black; background-color: red3"), br(), 
                h5("Arcsine-root"),
                downloadButton("taxadataArc", "Download", width = '50%', style = "background-color: red3"), br(), br(),
            )
          )
        })
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        arc_biom = chooseData$taxa.out$arcsin
        
        output$taxadataCount <- downloadHandler(
          
          filename = function() {
            paste("Count.Data.zip")
          },
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        incProgress(3/10, message = "Saving")
        
        output$taxadataRareCount <- downloadHandler(
          
          filename = function() {
            paste("Rarefied.Count.Data.zip")
          },
          content = function(rare.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(rare_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=rare.file, files=dataFiles)
          }
        )
        output$taxadataProp <- downloadHandler(
          filename = function() {
            paste("Proportion.Data.zip")
          },
          content = function(prop.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(prop_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=prop.file, files=dataFiles)
          }
        )
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        output$taxadataArc <- downloadHandler(
          filename = function() {
            paste("Arcsine.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(arc_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        incProgress(1/10, message = "Done")
        
        shinyjs::enable("divCalcRun")
      })
  })
  
  #########################
  ## Alpha Data Analysis ##
  #########################
  
  observeEvent(input$runbtn_bin,{
    shinyjs::disable("primvar")
    shinyjs::disable("blockid")
    
    length_prim <- length(names(table(chooseData$sam.dat[, input$primvar])))
    
    if (length_prim == 2 | length(input$alpha_levels) == 2){
      
      shinyjs::disable("runbtn_bin")
      shinyjs::disable("alphaCat1")
      shinyjs::disable("alphaCat2")
      shinyjs::disable("chooseData")
      shinyjs::disable("chooseMethod")
      
      withProgress(
        message = 'Calculation in progress',
        detail = 'This may take a while...', value = 0, {
          
          incProgress(1/10, message = "Renaming")
          chooseData$alpha.div = chooseData$alpha.div.rare
          
          alpha.categors = c(alpha.categos$cat1, alpha.categos$cat2)
          alpha.bin_categos = c(input$alphaCat1, input$alphaCat2)
          
          rename.cats_ref = try(alpha.bin_categos[which(alpha.categors == alpha.categos$cat1)], silent= TRUE) 
          rename.cats_com = try(alpha.bin_categos[which(alpha.categors != alpha.categos$cat1)], silent = TRUE) 
          
          alpha.bin.cat.ref.ori.out <- try(alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar), silent = TRUE) 
          sam_dat_alpha <- try(alpha.bin.cat.recode.func(chooseData$sam.dat, input$primvar, input$alpha_levels, rename.cats_ref, rename.cats_com), silent = TRUE) 
          
          incProgress(3/10, message = "Calculating")
          
          alpha.div.bin.out <- try(alpha.bin.cat.ref.func(input$primvar, rename.cats_ref,
                                                          rename.cats_com, sam_dat_alpha,
                                                          chooseData$alpha.div), silent = TRUE) 
          
          alpha.results$bin.var <- alpha.div.bin.out$bin.var
          alpha.results$alpha_div <- alpha.div.bin.out$alpha.div
          alpha.results$alpha.bin.sum.out <- try(alpha.bin.sum.func(bin.var = alpha.div.bin.out$bin.var,
                                                                    alpha.div = alpha.div.bin.out$alpha.div), silent = TRUE) 
          
          
          if (length(input$alpha_levels) == 2){
            re_sam <- try(reduced_data(sam_dat_alpha, input$primvar, c(input$alphaCat1, input$alphaCat2)), silent = TRUE) 
            re_alpha <- try(reduced_alpha(chooseData$alpha.div, sam_dat_alpha, input$primvar, input$blockid, c(input$alphaCat1, input$alphaCat2)), silent = TRUE) 
          } else if (length(names(table(chooseData$sam.dat[,input$primvar]))) == 2){
            re_sam <- sam_dat_alpha
            re_alpha <- chooseData$alpha.div
          }
          
          data.out <- try(data_mani_paired(re_alpha, re_sam, input$primvar, input$blockid), silent = TRUE) 
          
          if (input$chooseMethod == "Paired t-test" | input$chooseMethod == "Wilcoxon signed-rank test" | input$chooseMethod == "Multivariate Hotelling's t-squared test") {
            if (input$chooseMethod == "Paired t-test") {
              incProgress(3/10, message = "T test")
              
              t.test.out <- tryCatch(alpha.t.paired (chooseData$alpha.div, data.out), error = function(e) {  
                message("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)})
              
              table.out <- t.test.out
            }
            else if (input$chooseMethod == "Wilcoxon signed-rank test") {
              incProgress(3/10, message = "Wilcoxon test")
              
              wilcox.test.out <- tryCatch(alpha.wilcox.paired (chooseData$alpha.div, data.out), error = function(e) {  
                message("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)}) 
              
              table.out <- wilcox.test.out
            }
            else if (input$chooseMethod == "Multivariate Hotelling's t-squared test"){
              incProgress(3/10, message = "Hotelling test")
              
              alpha_stand <- tryCatch(hotel.data.standard(data.out), error = function(e) {
                message("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)})
              
              Htest_result <- tryCatch(alpha.hotel.paired(alpha_stand), error = function(e) {
                message("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)})
              
              table.out <- Htest_result
            } 
            
            alpha.data.results$table.out = table.out
            
            incProgress(3/10, message = "Graph")
            
            if (input$chooseMethod == "Paired t-test" | input$chooseMethod == "Wilcoxon signed-rank test"){
              output$alpha_display_results = renderUI({
                tagList(
                  box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                      align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                      plotOutput("box_plots", height = 750, width = 500)
                  )
                )
              })
              output$box_plots = renderPlot({alpha.bin.paired(re_alpha, data.out, table.out)})
            }
            else if (input$chooseMethod == "Multivariate Hotelling's t-squared test"){
              output$alpha_display_results = renderUI({
                tagList(
                  box(title = strong("Multivariate Hotelling's t-squared test (Omnibus test)", style = "color:black"), 
                      align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                      plotOutput("box_plots", height = 350, width = 500)
                  )
                )
              })
              
              ref_dat_stand = try(hotel.data.mani.for.CI(alpha_stand, par = "ref"), silent = TRUE) 
              com_dat_stand = try(hotel.data.mani.for.CI(alpha_stand, par = "com"), silent = TRUE)
              
              
              Simul_conf_Standard <- tryCatch(simul_CI(ref_dat_stand, com_dat_stand, 0.95, type = "Simul"), error = function(e) {  
                message("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)}) 
              
              output$box_plots = try(renderPlot({Hotell_CI_plot_3(Htest_result, Simul_conf_Standard, type = "Simul", level = 0.95)}), silent = TRUE) 
            }
          }
          
          output$alpha_downloadTable = renderUI({
            tagList(
              box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  p("You can download the summary statistics and data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Summary Statistics"),
                  downloadButton("downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                  h5("Data Analysis Outputs"),
                  downloadButton("downloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$downloadTabl1 <- downloadHandler(
            filename = function() {
              paste("Alpha.Sum.Table.txt")
            },
            content = function(file) {
              sink(file); print(alpha.results$alpha.bin.sum.out) ;sink() 
            }
          )
          
          output$downloadTabl2 <- downloadHandler(
            filename = function() {
              paste("Alpha.DA.Output.txt")
            },
            content = function(file) {
              write.table(alpha.data.results$table.out, file)
            }
          )
          
          ref_string = REFERENCE_CHECK_P(method_name = isolate(input$chooseMethod))
          if (is.null(ref_string)) {
            shinyjs::hide("alpha_references")
          } else {
            shinyjs::show("alpha_references")
            output$alpha_references = renderUI({
              tagList(
                box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                    HTML(paste(ref_string, collapse="<br/>"))
                )
              )
            })
          }
          shinyjs::enable("runbtn_bin")
          shinyjs::enable("primvar")
          shinyjs::enable("blockid")
          shinyjs::hide("barPanel")
          shinyjs::enable("chooseMethod")
          shinyjs::enable("alphaCat1")
          shinyjs::enable("alphaCat2")
          shinyjs::enable("chooseData")
        })
    }},ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(input$runbtn_bin_mult,{
    
    prim_length = length(names(table(chooseData$sam.dat[,input$primvar])))
    
    
    for (i in 1:prim_length){
      shinyjs::disable(paste0("alphaCat_", i))
    }
    shinyjs::disable("chooseData")
    shinyjs::disable("runbtn_bin_mult")
    shinyjs::disable("adjustment")
    shinyjs::disable("barPanel")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(1/10, message = "Renaming")
        
        chooseData$alpha.div = chooseData$alpha.div.rare
        alpha.categors = try(alpha.bin.cat.func(chooseData$sam.dat, input$primvar), silent = TRUE)
        
        alpha.bin.categos = c() 
        for (i in 1:prim_length){
          name_ind <- eval(parse(text = paste0("input$alphaCat_", i)))
          alpha.bin.categos <- c(alpha.bin.categos, name_ind)
        }
        
        sam_dat_alpha <- try(bin.cat.recode.func.mult(chooseData$sam.dat, input$primvar, alpha.categors, alpha.bin.categos), silent = TRUE)
        incProgress(3/10, message = "Calculating")
        
        alpha.div.bin.out <- try(alpha.mult.pair.cat.ref.func(input$primvar, alpha.bin.categos, sam_dat_alpha, chooseData$alpha.div), silent = TRUE)
        alpha.results$bin.var <- alpha.div.bin.out$bin.var
        alpha.results$alpha_div <- alpha.div.bin.out$alpha.div
        dat_list <- list(bin.var = alpha.div.bin.out$bin.var,
                         alpha.div = alpha.div.bin.out$alpha.div)
        
        alpha.results$alpha.bin.sum.out <- dat_list
        
        re_sam <- try(reduced_data(sam_dat_alpha, input$primvar, alpha.bin.categos), silent = TRUE) 
        re_alpha <- try(reduced_alpha(chooseData$alpha.div, sam_dat_alpha, input$primvar, input$blockid, alpha.bin.categos), silent = TRUE) 
        
        data.out <- try(data_mani_paired(re_alpha, re_sam, input$primvar, input$blockid), silent = TRUE) 
        pass_t_test <- FALSE
        
        if (input$chooseMethod == "Paired t-test" | input$chooseMethod == "Wilcoxon signed-rank test") {
          if (input$chooseMethod == "Paired t-test") {
            
            incProgress(3/10, message = "T test")
            t.test.out <- tryCatch(alpha.t.paired.mult.united(chooseData$alpha.div, data.out, method_adj = "BH"), error = function(e) {  
              message("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)})
            
            table.out <- t.test.out$table 
            table.out_2 <- t.test.out$download
          }
          
          else if (input$chooseMethod == "Wilcoxon signed-rank test") {
            
            incProgress(3/10, message = "Wilcoxon test")
            wilcox.test.out <- tryCatch(alpha.wilcox.paired.mult.united(chooseData$alpha.div, data.out, "BH"), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })  
            
            table.out <- wilcox.test.out$table 
            table.out_2 <- wilcox.test.out$download 
          }
          
          alpha.data.results$table.out <- table.out_2
          
          incProgress(3/10, message = "Graph")
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots = tryCatch(renderPlot({alpha.bin.multi.paired(chooseData$alpha.div, data.out)}), error = function(e) {  
            message("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)})
          
          shinyjs::show("barPanel")
          shinyjs::show("adjustment")
          output$barPanel = renderUI({
            navbarPage(
              title = "Pairwise Comparison",
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")), 
              tabPanel("Simpson", dataTableOutput("alpha_3")), 
              tabPanel("InvSimpson", dataTableOutput("alpha_4")), 
              tabPanel("Fisher", dataTableOutput("alpha_5")), 
              tabPanel("Chao1", dataTableOutput("alpha_6")), 
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")), 
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(table.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(table.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(table.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(table.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(table.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(table.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(table.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(table.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(table.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          
        }else if (input$chooseMethod == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
          shinyjs::hide("barPanel")
          shinyjs::hide("adjustment")
          
          incProgress(3/10, message = "F test")
          
          F.test.out <- tryCatch(alpha.f.pair.mult.overall(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          F.test.out.2 <- tryCatch(alpha.f.pair.mult.overall_download(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          tukey.out <- tryCatch(alpha.f.pair.mult.tukey(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          tukey.out.2 <- tryCatch(alpha.f.pair.mult.tukey_download(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          
          table.out <- list(ANOVA_F = F.test.out, Tukey_HSD = tukey.out)
          table.out.2 <- list(ANOVA_F = F.test.out.2, Tukey_HSD = tukey.out.2)
          
          alpha.data.results$table.out <- table.out.2
          
          incProgress(3/10, message = "Graph")
          p.val.f.test<- try(F.test.out[,"P.value"], silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("ANOVA F-test (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots = try(renderPlot({alpha.bin.f.mult.paired(chooseData$alpha.div, data.out, p.val.f.test)}), silent = TRUE)
          
          shinyjs::show("barPanel")
          
          output$barPanel = renderUI({
            navbarPage(
              title = "Tukey's HSD (Pairwise Comparison)",
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(tukey.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(tukey.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(tukey.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(tukey.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(tukey.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(tukey.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(tukey.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(tukey.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(tukey.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          shinyjs::hide("adjustment")
        }
        else if (input$chooseMethod == "Friedman's test (global) with Conover's test (pairwise)"){
          
          shinyjs::hide("barPanel")
          shinyjs::show("adjustment")
          incProgress(3/10, message = "Friedman test")
          
          test.out <- tryCatch(alpha.f.pair.mult.friedman(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          test.out.2 <- tryCatch(alpha.f.pair.mult.friedman_download(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          conover.out <- tryCatch(alpha.f.pair.mult.conover(data.out, p_adjustment = "BH"), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          conover.out.2 <- tryCatch(alpha.f.pair.mult.conover_download(data.out, p_adjustment = "BH"), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          table.out <- list(Friedman_test = test.out, Conover_test = conover.out)
          table.out.2 <- list(Friedman_test = test.out.2, Conover_test = conover.out.2)
          alpha.data.results$table.out <- table.out.2 
          
          incProgress(3/10, message = "Graph")
          p.val.test <- try(test.out[,"P.value"])
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Friedman's test (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots <- tryCatch(renderPlot({alpha.bin.f.mult.paired(chooseData$alpha.div, data.out, p.val.test)}), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          shinyjs::show("barPanel")
          output$barPanel = renderUI({
            navbarPage(
              title = "Conover's test (Pairwise Comparison)",
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(conover.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(conover.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(conover.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(conover.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(conover.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(conover.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(conover.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(conover.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(conover.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          
        }else if (input$chooseMethod == "Durbin's test (global) with Conover's test (pairwise)"){
          
          shinyjs::hide("barPanel")
          shinyjs::show("adjustment")
          incProgress(3/10, message = "Durbin's test (global) with Conover's test (pairwise)")
          
          data.out.durbin <- try(fill_missing_na(data.out), silent = TRUE)
          
          durbin.test <- tryCatch(alpha.pair.mult.durbin(data.out.durbin), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          test.out <- durbin.test$table 
          test.out.2 <- durbin.test$download 
          
          pairwise.out <- tryCatch(alpha.pairwise.durbin(data.out.durbin, p_adjustment = "BH"), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          pairwise.out.2 <- tryCatch(alpha.pairwise.durbin_download(data.out.durbin, p_adjustment = "BH"), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          table.out <- list(Durbin_test = test.out, Durbin_pairwise_test = pairwise.out)
          table.out.2 <- list(Durbin_test = test.out.2, Durbin_pairwise_test = pairwise.out.2)
          alpha.data.results$table.out <- table.out.2 
          
          incProgress(3/10, message = "Graph")
          
          p.val.test<- try(test.out.2[,"P.value"], silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Durbin's test (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots <- tryCatch(renderPlot({alpha.bin.f.mult.paired(chooseData$alpha.div, data.out, p.val.test)}), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          shinyjs::show("barPanel")
          output$barPanel = renderUI({
            navbarPage(
              title = "Conover's test (Pairwise Comparison)",   
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(pairwise.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(pairwise.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(pairwise.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(pairwise.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(pairwise.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(pairwise.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(pairwise.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(pairwise.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(pairwise.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          
        } else if (input$chooseMethod == "LMM (LRT (global) with t-test (pairwise))"){
          
          shinyjs::hide("barPanel")
          shinyjs::show("adjustment")
          
          incProgress(3/10, message = "LMM (LRT (global) with t-test (pairwise))")
          
          new_ind <- which(alpha.categors == input$lmm_level_alpha)
          new_lmm_level <- alpha.bin.categos[new_ind]
          
          test.out <- tryCatch(global.alpha.lmm.dat(data.out, new_lmm_level), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          test.out.2 <- tryCatch(global.alpha.lmm.dat(data.out, new_lmm_level, TRUE), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          pairwise.out <- tryCatch(alpha.pairwise.lmm (data.out, new_lmm_level, "BH"), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          pairwise.out.2  <- tryCatch(alpha.pairwise.lmm (data.out, new_lmm_level, "BH", TRUE), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          table.out <- list(LMM_test = test.out, LMM_pairwise_test = pairwise.out)
          table.out.2 <- list(LMM_test = test.out.2, LMM_pairwise_test = pairwise.out.2)
          alpha.data.results$table.out <- table.out.2 
          
          incProgress(3/10, message = "Graph")
          p.val.test<- try(test.out.2[,"P.value"], silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong(" LMM (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots <- tryCatch(renderPlot({alpha.bin.f.mult.paired(chooseData$alpha.div, data.out, p.val.test)}), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          
          shinyjs::show("barPanel")
          output$barPanel = renderUI({
            navbarPage(
              title = "LMM (Pairwise Comparison)",   
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(pairwise.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(pairwise.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(pairwise.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(pairwise.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(pairwise.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(pairwise.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(pairwise.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(pairwise.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(pairwise.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
        }
        
        output$alpha_downloadTable = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        output$downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            sink(file); print(alpha.results$alpha.bin.sum.out); sink()  
          }
        )
        
        output$downloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            sink(file); print(alpha.data.results$table.out); sink()    
          }
        )
        
        ref_string = REFERENCE_CHECK_P(method_name = isolate(input$chooseMethod))
        if (is.null(ref_string)) {
          shinyjs::hide("alpha_references")
        } else {
          shinyjs::show("alpha_references")
          output$alpha_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("runbtn_bin_mult")
        shinyjs::enable("primvar")
        shinyjs::enable("blockid")
        shinyjs::enable("chooseMethod")
        shinyjs::enable("adjustment")
        for (i in 1:prim_length){
          shinyjs::enable(paste0("alphaCat_", i))
        }
        shinyjs::enable("chooseData")
        shinyjs::enable("barPanel")
      })
  }
  
  ,ignoreNULL = TRUE, ignoreInit = TRUE)
  
  #########################
  ## Beta Data Analysis ###
  #########################
  
  observeEvent(input$beta_runbtn,{
    shinyjs::disable("beta_runbtn")
    shinyjs::disable("beta_primvar")
    shinyjs::disable("beta_blockid")
    shinyjs::disable("beta_chooseMethod")
    
    length_prim <- length(names(table(chooseData$sam.dat[, input$beta.primvar])))
    
    if (length_prim== 2 | length(input$beta_levels) == 2){
      shinyjs::hide("beta_barPanel")
      shinyjs::disable("betaCat1")
      shinyjs::disable("betaCat2")
      
      withProgress(
        message = 'Calculation in progress',
        detail = 'This may take a while...', value = 0, {
          
          incProgress(3/10, message = "Rename")
          
          beta.categors = c(beta.categos$cat1, beta.categos$cat2)
          beta.bin_categos = c(input$betaCat1, input$betaCat2)
          
          rename.catsbin_ref = beta.bin_categos[which(beta.categors == beta.categos$cat1)]
          rename.catsbin_com = beta.bin_categos[which(beta.categors != beta.categos$cat1)]
          
          beta.bin.cat.ref.ori.out <- tryCatch(beta.bin.cat.ref.ori.func(chooseData$sam.dat, input$beta.primvar), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          beta.sam.dat.bin <- tryCatch(beta.bin.cat.recode.func(chooseData$sam.dat, input$beta.primvar,      #beta.chooseData$sam.dat.bin = sam.dat
                                                                beta.bin.cat.ref.ori.out,
                                                                rename.catsbin_ref, rename.catsbin_com), error = function(e) {
                                                                  message ("No outcome is available!")
                                                                  showModal(modalDialog(div("No outcome is available!")))
                                                                  return(NULL)
                                                                })
          
          incProgress(3/10, message = "Calculating")
          
          beta.div.bin.out <- tryCatch(beta.bin.cat.ref.func(input$beta.primvar, rename.catsbin_ref,
                                                             rename.catsbin_com, beta.sam.dat.bin, ds.Ks$res), error = function(e) {
                                                               message ("No outcome is available!")
                                                               showModal(modalDialog(div("No outcome is available!")))
                                                               return(NULL)
                                                             })
          
          beta.results$beta_div <- beta.div.bin.out$Ds   
          
          if (length(input$beta_levels) == 2) {
            re_sam_beta <- try(reduced_data(beta.sam.dat.bin, input$beta.primvar, c(input$betaCat1, input$betaCat2)), silent = TRUE) 
          } else if (length(names(table(chooseData$sam.dat[,input$beta.primvar]))) == 2) {
            re_sam_beta <- beta.sam.dat.bin
          }
          
          
          if (input$beta_chooseMethod == "PERMANOVA"){
            incProgress(3/10, message = "PERMANOVA")
            set.seed(521)
            beta.bin.out <- tryCatch(beta.bin.cat.ref.func(input$beta.primvar,
                                                           rename.catsbin_ref, rename.catsbin_com,
                                                           beta.sam.dat.bin, Ds.Ks = ds.Ks$res), error = function(e) {
                                                             message ("No outcome is available!")
                                                             showModal(modalDialog(div("No outcome is available!")))
                                                             return(NULL)
                                                           })
            
            beta.data.results$data.q.out <- beta.bin.out
          } 
          
          incProgress(3/10, message = "Plotting Graph")
          output$beta_display_results_cross = renderUI({
            tagList(
              box(title = strong("Graphs for Beta Diversity", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("beta_graph_plots.bin", height = 750, width = 500)
              )
            )
          })
          
          if (input$beta_chooseMethod == "PERMANOVA"){
            set.seed(521)
            beta.down.results <- tryCatch(beta.permanova.paired(num_perm = 3000, beta.results$beta_div, re_sam_beta, input$beta.primvar, input$beta.blockid), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            beta.down.results$CS <- beta.down.results$download 
            table.out <- beta.down.results$table
            
            output$beta_graph_plots.bin = renderPlot({
              
              tryCatch(mirkat.bin.plot_2(table.out, beta.data.results$data.q.out), error = function(e) {
                message ("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)
              })
            })} 
          
          else if (input$beta_chooseMethod == "PERMANOVA-FL"){
            beta.div <<-  ds.Ks$res$Ds
            div.Jac <<- beta.div$Jaccard
            div.Bray <<- beta.div$Bray.Curtis
            div.U <<- beta.div$U.UniFrac
            div.G <<- beta.div$G.UniFrac
            div.W <<- beta.div$W.UniFrac
            
            set.seed(521)
            beta.down.results <- tryCatch(beta.permanovaFL.paired.ind(beta.sam.dat.bin, input$beta.primvar, input$beta.blockid), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            beta.down.results$CS <- beta.down.results$download 
            table.out <- beta.down.results$table 
            
            output$beta_graph_plots.bin = renderPlot({
              tryCatch(mirkat.bin.plot_2(table.out, beta.data.results$data.q.out), error = function(e) {
                message ("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)
              })})
          }
          
          
          output$beta_downloadTable = renderUI({
            tagList(
              box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Data Analysis Outputs"),
                  downloadButton("beta_downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          incProgress(3/10, message = "SAVE")
          
          output$beta_downloadTabl1 <- downloadHandler(
            filename = function() {
              paste("Beta.DA.Output.txt")
            },
            content = function(file) {
              out_temp = as.data.frame(beta.down.results$CS)
              rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac")
              beta.down.results$CS = out_temp
              write.table(beta.down.results$CS, file, sep="\t")
            }
          )
          
          ref_string = REFERENCE_CHECK_P(method_name = isolate(input$beta_chooseMethod))
          if (is.null(ref_string)) {
            shinyjs::hide("beta_references")
          } else {
            shinyjs::show("beta_references")
            output$beta_references = renderUI({
              tagList(
                box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                    HTML(paste(ref_string, collapse="<br/>"))
                )
              )
            })
          }
          
          shinyjs::enable("beta_runbtn")
          shinyjs::enable("beta_primvar")
          shinyjs::enable("beta_blockid")
          shinyjs::enable("beta_chooseMethod")
          shinyjs::enable("betaCat1")
          shinyjs::enable("betaCat2")
        })
    }}
    
    ,ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(input$beta_runbtn_mult,{
    shinyjs::disable("beta_runbtn_mult")
    shinyjs::disable("beta_primvar")
    shinyjs::disable("beta_blockid")
    shinyjs::disable("beta_chooseMethod")
    
    prim_length = length(names(table(chooseData$sam.dat[,input$beta.primvar])))
    for (i in 1:prim_length){
      shinyjs::disable(paste0("betaCat_", i))
    }
    shinyjs::disable("beta_adjustment")  
    shinyjs::disable("beta_barPanel")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(1/10, message = "Renaming")
        
        beta.div.bin.out <- try(beta.bin.mult.cat.ref.func(input$beta.primvar, chooseData$sam.dat, ds.Ks$res), silent = TRUE) 
        beta.results$beta_div <- beta.div.bin.out$Ds
        beta.categors = try(alpha.bin.cat.func(chooseData$sam.dat, input$beta.primvar), silent = TRUE) 
        
        beta.bin.categos = c() 
        for (i in 1:prim_length){
          name_ind <- eval(parse(text = paste0("input$betaCat_", i)))
          beta.bin.categos <- c(beta.bin.categos, name_ind)
        }
        
        beta_sam_mult <- tryCatch(bin.cat.recode.func.mult(chooseData$sam.dat, input$beta.primvar, beta.categors, beta.bin.categos), error = function(e) {
          message ("No outcome is available!")
          showModal(modalDialog(div("No outcome is available!")))
          return(NULL)
        })
        
        incProgress(3/10, message = "Calculating")
        
        re_sam_beta <- try(reduced_data(beta_sam_mult, input$beta.primvar, beta.bin.categos), silent = TRUE)   #요기 
        #beta.bin.2 <- try(beta.bin.mult.cat.ref.func(input$beta.primvar, re_sam_beta, ds.Ks$res), silent = TRUE) 
        beta.bin.2 <- try(beta.bin.mult.cat.ref.func_sub(beta_sam_mult, input$beta.primvar, input$beta.blockid, ds.Ks$res, beta.bin.categos), silent = TRUE)
        
        beta.results$beta_div <- try(reduced_beta(beta.results$beta_div, beta_sam_mult, input$beta.primvar, input$beta.blockid, beta.bin.categos), silent = TRUE)   #요기 
        
        if (input$beta_chooseMethod == "PERMANOVA_ACROSS" | input$beta_chooseMethod == "PERMANOVA_BASE") {
          
          if (input$beta_chooseMethod == "PERMANOVA_ACROSS") {
            incProgress(3/10, message = "PERMANOVA") 
            set.seed(521)
            permanova.out <- tryCatch(permanova.paired.mult.united(num_perm = 3000, beta.results$beta_div, re_sam_beta, input$beta.primvar, input$beta.blockid, method_adj = "BH"), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            permanova.pair.out <- permanova.out$table
            permanova.pair.out.2 <- permanova.out$download 
            
            set.seed(521)
            permanova.global <- tryCatch(beta.permanova.paired.mult.glob(num_perm = 3000, beta.results$beta_div, re_sam_beta, input$beta.primvar, input$beta.blockid, download = TRUE), 
                                         error = function(e) {
                                           message ("No outcome is available!")
                                           showModal(modalDialog(div("No outcome is available!")))
                                           return(NULL)})
            
            permanova.global.out <- permanova.global$table 
            
            table.out <- try(list(Global_PERMANOVA = permanova.global.out, Pairwise_PERMANOVA = permanova.pair.out.2))
            
            shinyjs::show("beta_barPanel")
            shinyjs::show("beta_adjustment")
            output$beta_barPanel = renderUI({
              navbarPage(
                title = "PERMANOVA (Pairwise Comparison)",
                tabPanel("Jaccard", dataTableOutput("beta_1")),
                tabPanel("Bray.Curtis", dataTableOutput("beta_2")), 
                tabPanel("U.UniFrac", dataTableOutput("beta_3")), 
                tabPanel("G.UniFrac", dataTableOutput("beta_4")), 
                tabPanel("W.UniFrac", dataTableOutput("beta_5")),
                position = "static-top"
              )
            })
            
            output$beta_1  = try(renderDataTable({data.table(permanova.pair.out[["Jaccard"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_2  = try(renderDataTable({data.table(permanova.pair.out[["Bray.Curtis"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_3  = try(renderDataTable({data.table(permanova.pair.out[["U.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_4  = try(renderDataTable({data.table(permanova.pair.out[["G.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_5  = try(renderDataTable({data.table(permanova.pair.out[["W.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          }
          
          else if(input$beta_chooseMethod == "PERMANOVA_BASE"){
            incProgress(3/10, message = "PERMANOVA") 
            set.seed(521)
            permanova.out <- tryCatch(permanova.paired.mult.united_base(num_perm = 3000, input$base_level_beta, beta.results$beta_div, re_sam_beta, input$beta.primvar, input$beta.blockid, method_adj = "BH"), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            permanova.pair.out <- permanova.out$table
            permanova.pair.out.2 <- permanova.out$download 
            
            table.out <- permanova.pair.out.2
            
            set.seed(521)
            permanova.global <- tryCatch(beta.permanova.paired.mult.glob(num_perm = 3000, beta.results$beta_div, re_sam_beta, input$beta.primvar, input$beta.blockid, download = TRUE), 
                                         error = function(e) {
                                           message ("No outcome is available!")
                                           showModal(modalDialog(div("No outcome is available!")))
                                           return(NULL)})
            
            permanova.global.out <- permanova.global$table 
            
            shinyjs::show("beta_barPanel")
            shinyjs::show("beta_adjustment")
            output$beta_barPanel = renderUI({
              navbarPage(
                title = "PERMANOVA (Pairwise Comparison)",
                tabPanel("Jaccard", dataTableOutput("beta_1")),
                tabPanel("Bray.Curtis", dataTableOutput("beta_2")), 
                tabPanel("U.UniFrac", dataTableOutput("beta_3")), 
                tabPanel("G.UniFrac", dataTableOutput("beta_4")), 
                tabPanel("W.UniFrac", dataTableOutput("beta_5")),
                position = "static-top"
              )
            })
            
            output$beta_1  = try(renderDataTable({data.table(permanova.pair.out[["Jaccard"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_2  = try(renderDataTable({data.table(permanova.pair.out[["Bray.Curtis"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_3  = try(renderDataTable({data.table(permanova.pair.out[["U.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_4  = try(renderDataTable({data.table(permanova.pair.out[["G.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_5  = try(renderDataTable({data.table(permanova.pair.out[["W.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          }
          
          else if (input$beta_chooseMethod == "PERMANOVA-FL") {
            
            incProgress(3/10, message = "PERMANOVA-FL")
            set.seed(521)
            
            permanovaFL.out <- tryCatch(permanovaFL.paired.mult.united (num_perm = 3000, beta.results$beta_div, beta_sam_mult, input$beta.primvar, input$beta.blockid, method_adj = "BH"), 
                                        error = function(e) {
                                          message ("No outcome is available!")
                                          showModal(modalDialog(div("No outcome is available!")))
                                          return(NULL)})
            
            permanovaFL.pair.out.2 <- permanovaFL.out$download 
            permanovaFL.pair.out <- permanovaFL.out$table 
            
            set.seed(521)
            permanovaFL.global<- tryCatch(beta.permanovaFL.paired.mult.glob(num_perm = 3000, beta.results$beta_div, beta_sam_mult, input$beta.primvar, input$beta.blockid, TRUE), 
                                          error = function(e) {
                                            message ("No outcome is available!")
                                            showModal(modalDialog(div("No outcome is available!")))
                                            return(NULL)})
            
            permanovaFL.global.out <- permanovaFL.global$table 
            
            table.out <- list(Global_PERMANOVA_FL = permanovaFL.global.out, Pairwise_PERMANOVA_FL = permanovaFL.pair.out.2)
            
            output$beta_barPanel = renderUI({
              navbarPage(
                title = "Pairwise Comparison",
                tabPanel("Jaccard", dataTableOutput("beta_1")),
                tabPanel("Bray.Curtis", dataTableOutput("beta_2")), 
                tabPanel("U.UniFrac", dataTableOutput("beta_3")), 
                tabPanel("G.UniFrac", dataTableOutput("beta_4")), 
                tabPanel("W.UniFrac", dataTableOutput("beta_5")),
                position = "static-top"
              )
            })
            output$beta_1  = try(renderDataTable({data.table(permanovaFL.pair.out[["Jaccard"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_2  = try(renderDataTable({data.table(permanovaFL.pair.out[["Bray.Curtis"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_3  = try(renderDataTable({data.table(permanovaFL.pair.out[["U.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_4  = try(renderDataTable({data.table(permanovaFL.pair.out[["G.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
            output$beta_5  = try(renderDataTable({data.table(permanovaFL.pair.out[["W.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          }
        }
        
        beta.data.results$data.q.out <- table.out
        
        if(input$beta_chooseMethod == "PERMANOVA_ACROSS"){
          incProgress(3/10, message = "Plotting Graph")
          output$beta_display_results_cross = renderUI({
            tagList(
              box(title = strong("PERMANOVA (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("beta_graph_plots.bin", height = 800, width = 700)   #0711 수정 
              )
            )
          })
          
          set.seed(521)
          beta.p.val <- permanova.global$p.val 
          
          output$beta_graph_plots.bin = try(renderPlot({
            mirkat.bin.plot_2_mult(beta.p.val, beta.bin.2, re_sam_beta, input$beta.primvar)
          }), silent = TRUE)
        }
        else if(input$beta_chooseMethod == "PERMANOVA_BASE"){
          
          incProgress(3/10, message = "Plotting Graph")
          output$beta_display_results_cross = renderUI({
            tagList(
              box(title = strong("PERMANOVA (Global Test)", style = "color:black"), 
                  align = "center", width = NULL, status = "warning", solidHeader = TRUE,
                  plotOutput("beta_graph_plots.bin", height = 800, width = 700)   
              )
            )
          })
          
          beta.p.val <- permanova.global$p.val 
          
          output$beta_graph_plots.bin = try(renderPlot({
            mirkat.bin.plot_2_mult(beta.p.val, beta.bin.2, re_sam_beta, input$beta.primvar)
          }), silent = TRUE)
          
        }
        
        output$beta_downloadTable = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("beta_downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        incProgress(3/10, message = "SAVE")
        
        output$beta_downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            sink(file);print(beta.data.results$data.q.out);sink()           
          }
        )
        
        ref_string = REFERENCE_CHECK_P(method_name = isolate(input$beta_chooseMethod))
        if (is.null(ref_string)) {
          shinyjs::hide("beta_references")
        } else {
          shinyjs::show("beta_references")
          output$beta_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        shinyjs::enable("beta_runbtn_mult")
        shinyjs::enable("beta_primvar")
        shinyjs::enable("beta_blockid")
        shinyjs::enable("beta_chooseMethod")
        for (i in 1:prim_length){
          shinyjs::enable(paste0("betaCat_", i))
        }
        shinyjs::enable("beta_adjustment")
        shinyjs::enable("beta_barPanel")
      })
    
  },ignoreNULL = TRUE, ignoreInit = TRUE)
  
  
  ##################################
  # Taxa Comparative Data Analysis #
  ##################################
  
  observeEvent(input$runbtn_bin_taxa_base,{
    shinyjs::disable("runbtn_bin_taxa_base")
    shinyjs::disable("dataType_taxa")
    shinyjs::disable("primvar_taxa")
    shinyjs::disable("blockid_taxa")
    shinyjs::disable("chooseMethod_taxa_base")
    shinyjs::disable("chooseRef_taxa")
    shinyjs::hide("taxa_barPanel")
    shinyjs::disable("taxaCat1")
    shinyjs::disable("taxaCat2")
    
    
    if(length(input$taxa_levels) == 2 | length(names(table(chooseData$sam.dat[,input$primvar_taxa]))) == 2){
      
      withProgress(
        message = 'Calculation in progress',
        detail = 'This may take a while...', value = 0, {
          
          taxa.categors = c(taxa.categos$cat1, taxa.categos$cat2)
          taxa.bin_categos = c(input$taxaCat1, input$taxaCat2)
          rename.cats_ref = taxa.bin_categos[which(taxa.categors == taxa.categos$cat1)]
          rename.cats_com = taxa.bin_categos[which(taxa.categors != taxa.categos$cat1)]
          taxa.bin.cat.ref.ori.out <- try(taxa.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar_taxa), silent = TRUE) 
          
          if (input$include_species.dend == "Phylum - Species") {
            include = TRUE
          } else {
            include = FALSE
          }
          
          incProgress(1/10, message = "Examining Data in progress")
          
          re_sam_tax <- try(taxa.bin.cat.recode.func(chooseData$sam.dat, input$primvar_taxa, c(taxa.categos$cat1, taxa.categos$cat2),
                                                     rename.cats_ref, rename.cats_com), silent = TRUE) 
          
          taxa.bin.out <- try(taxa.bin.cat.ref.united.func.mult.pairwise(input$primvar_taxa, re_sam_tax, taxa = chooseData$taxa.out[[taxa.types$dataType]]), silent = TRUE) 
          
          if (length(input$taxa_levels) == 2){
            re_sam <- try(reduced_data(re_sam_tax, input$primvar_taxa, c(input$taxaCat1, input$taxaCat2)), silent = TRUE) 
            re_tax <- try(reduced_taxa(taxa.bin.out$taxa, re_sam_tax, input$primvar_taxa, input$blockid_taxa, c(input$taxaCat1, input$taxaCat2)), silent = TRUE) 
          } else if (length(names(table(re_sam_tax[,input$primvar_taxa]))) == 2){
            re_sam <- re_sam_tax
            re_tax <- taxa.bin.out$taxa
          }
          
          taxa.bin.out_2 <- try(taxa.bin.cat.ref.united.func.mult.pairwise(input$primvar_taxa, re_sam, taxa = re_tax), silent = TRUE) 
          taxa.results$bin.var <- taxa.bin.out_2$bin.var
          taxa.results$taxa <- taxa.bin.out_2$taxa
          taxa.results$taxa.bin.sum.out <- try(taxa.bin.sum.united.func(taxa.results$bin.var, taxa.results$taxa), silent = TRUE) 
          taxa_dataBinvar <- taxa.results$bin.var
          taxa_dataTaxa <- taxa.results$taxa
          
          if (input$chooseMethod_taxa_base == "Paired t-test"){
            incProgress(5/10, message = "Paired t-test")
            
            taxa.t.test.out  <- tryCatch(taxa.bin.t.paired.united(re_tax, re_sam, input$primvar_taxa, input$blockid_taxa), 
                                         error = function(e) {
                                           message ("No outcome is available!")
                                           showModal(modalDialog(div("No outcome is available!")))
                                           return(NULL)})
            
            taxa.t.test.q.out  <- tryCatch(bin.q.united.func(taxa.t.test.out, method = "BH"), 
                                           error = function(e) {
                                             message ("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            
            taxa.outputs$DAoutput <- taxa.t.test.q.out
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.t.test.q.out[[r]]$Q.value[complete.cases(taxa.t.test.q.out[[r]]$Q.value)] < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
            
            
          }else if (input$chooseMethod_taxa_base == "Wilcoxon signed-rank test"){
            incProgress(5/10, message = "Wilcoxon signed-rank test")
            
            taxa.wilcox.test.out  <- tryCatch(taxa.bin.wilcox.paired.united(re_tax, re_sam, input$primvar_taxa, input$blockid_taxa), 
                                              error = function(e) {
                                                message ("No outcome is available!")
                                                showModal(modalDialog(div("No outcome is available!")))
                                                return(NULL)})
            
            taxa.wilcox.test.q.out  <- tryCatch(bin.q.united.func(taxa.wilcox.test.out, method = "BH"), 
                                                error = function(e) {
                                                  message ("No outcome is available!")
                                                  showModal(modalDialog(div("No outcome is available!")))
                                                  return(NULL)})
            
            taxa.wilcox.test.est.added <- taxa.wilcox.test.est.func(taxa.results$bin.var, re_tax, rename.cats_ref, rename.cats_com, taxa.wilcox.test.q.out)
            
            
            taxa.outputs$DAoutput <- taxa.wilcox.test.est.added
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.wilcox.test.q.out[[r]]$Q.value[complete.cases(taxa.wilcox.test.q.out[[r]]$Q.value)] < 0.05)/4) 
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          }
          else if (input$chooseMethod_taxa_base == "LDM"){
            incProgress(5/10, message = "LDM")
            
            taxa_ex <<- re_tax
            
            taxa_1 <<- taxa_ex[[1]]
            taxa_2 <<- taxa_ex[[2]]
            taxa_3 <<- taxa_ex[[3]]
            taxa_4 <<- taxa_ex[[4]]
            taxa_5 <<- taxa_ex[[5]]
            taxa_6 <<- taxa_ex[[6]]
            
            if (input$dataType_taxa == "Arcsine-root"){
              
              set.seed(521)
              taxa.ldm.test.q.out  <- tryCatch(taxa.ldm.paired.united_trans(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa), 
                                               error = function(e) {
                                                 message ("No outcome is available!")
                                                 showModal(modalDialog(div("No outcome is available!")))
                                                 return(NULL)})
              
              taxa.ldm.test.est.added <- taxa.wilcox.test.est.func(taxa.results$bin.var, re_tax, rename.cats_ref, rename.cats_com, taxa.ldm.test.q.out)
              
              taxa.outputs$DAoutput <- taxa.ldm.test.est.edded
              
            }else{
              set.seed(521)
              taxa.ldm.test.q.out  <- tryCatch(taxa.ldm.paired.united_freq(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa), 
                                               error = function(e) {
                                                 message ("No outcome is available!")
                                                 showModal(modalDialog(div("No outcome is available!")))
                                                 return(NULL)})
              
              
              taxa.ldm.test.est.added <- taxa.wilcox.test.est.func(taxa.results$bin.var, re_tax, rename.cats_ref, rename.cats_com, taxa.ldm.test.q.out)
              taxa.outputs$DAoutput <- taxa.ldm.test.est.added 
            }
            
            nrow = numeric()
            
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.ldm.test.q.out[[r]]$Q.value[complete.cases(taxa.ldm.test.q.out[[r]]$Q.value)] < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          } 
          
          else if (input$chooseMethod_taxa_base == "LMM (LRT (global) with t-test (pairwise))"){
            incProgress(5/10, message = "LMM (LRT (global) with t-test (pairwise))")
            
            taxa.lmm.test.out  <- tryCatch(taxa.pair.lmm.united(re_tax, re_sam, input$primvar_taxa, input$blockid_taxa), 
                                           error = function(e) {
                                             message ("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
            
            taxa.outputs$DAoutput <- taxa.lmm.test.out
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.lmm.test.out[[r]]$Q.value[complete.cases(taxa.lmm.test.out[[r]]$Q.value)] < 0.05)/4)  
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          }
          
          shinyjs::show("taxa_display_dend")
          output$taxa_display_dend = renderUI({
            box(title = strong("Dendrogram", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                grVizOutput("dendrogram", height = 1000))
          })
          
          if(input$include_species.dend == "Phylum - Species"){
            incProgress(3/10, message = "Displaying Results in progress")
            output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong("Box Plot", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = nrow[1]*250, width = 750),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = nrow[2]*250, width = 750),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = nrow[3]*250, width = 750),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = nrow[4]*250, width = 750),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = nrow[5]*250, width = 750),
                       )
                       ,tabPanel("Species", align = "center",
                                 plotOutput("rank6", height = nrow[6]*250, width = 750),
                       )
                )
              )
            })
            
            output$rank1 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 1, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 2, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 3, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 4, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 5, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE) 
            
            output$rank6 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 6, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            
            if(input$chooseMethod_taxa_base == "Wilcoxon signed-rank test"){
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.wilcox.test.est.added, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            else if(input$chooseMethod_taxa_base == "LDM"){
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.ldm.test.est.added, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            else{
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            
            
            incProgress(1/10, message = "Displaying Results in progress")
            output$downloadTable_taxa = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +5px;"),
                box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                    p("You can download the summary statistics and data analysis outputs.", style = "font-size:11pt"),
                    h5("Summary Statistics"),
                    downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3"),
                    h5("Dendrogram"),
                    downloadButton("gdownload", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl1 <- downloadHandler(
              filename = function() {
                paste("Taxa.Sum.Table.zip")
              },
              content = function(sum.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                zip(zipfile=sum.file, files=dataFiles)
              }
            )
            
            output$tdownloadTabl2 <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              },
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                write.table(as.data.frame(taxa.outputs$DAoutput$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
            if(input$chooseMethod_taxa_base == "Wilcoxon signed-rank test"){
              dend_result <- taxa.wilcox.test.est.added
            }
            else if(input$chooseMethod_taxa_base == "LDM"){
              dend_result <- taxa.ldm.test.est.added
            }
            else{
              dend_result <- taxa.outputs$DAoutput
            }

            output$gdownload <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
              },
              content = function(file) {
                htmlwidgets::saveWidget(as_widget(taxa.sig.dend(dend_result, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file)
              }
            )
            
            
          }
          else {
            incProgress(3/10, message = "Displaying Results in progress")
            incProgress(3/10, message = "Displaying Results in progress")
            output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong("Box Plot", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = nrow[1]*250, width = 750),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = nrow[2]*250, width = 750),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = nrow[3]*250, width = 750),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = nrow[4]*250, width = 750),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = nrow[5]*250, width = 750),
                       )
                       
                )
              )
            })
            
            
            output$rank1 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 1, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 2, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 3, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 4, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.pair.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 5, mult.test.cor = TRUE, re_sam, input$primvar_taxa, input$blockid_taxa)
            }), silent = TRUE) 
            
            
            
            if(input$chooseMethod_taxa_base == "Wilcoxon signed-rank test"){
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.wilcox.test.est.added, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            else if(input$chooseMethod_taxa_base == "LDM"){
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.ldm.test.est.added, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            else{
              output$dendrogram = try(renderGrViz({
                taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)$flow.text
              }), silent = TRUE)
            }
            
            incProgress(1/10, message = "Displaying Results in progress")
            output$downloadTable_taxa = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +5px;"),
                box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                    p("You can download the summary statistics and data analysis outputs.", style = "font-size:11pt"),
                    h5("Summary Statistics"),
                    downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3"),
                    h5("Dendrogram"),
                    downloadButton("gdownload", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl1 <- downloadHandler(
              filename = function() {
                paste("Taxa.Sum.Table.zip")
              },
              content = function(sum.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                zip(zipfile=sum.file, files=dataFiles)
              }
            )
            
            output$tdownloadTabl2 <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              },
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                write.table(as.data.frame(taxa.outputs$DAoutput$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                write.table(as.data.frame(taxa.outputs$DAoutput$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
            if(input$chooseMethod_taxa_base == "Wilcoxon signed-rank test"){
              dend_result <- taxa.wilcox.test.est.added
            }
            else if(input$chooseMethod_taxa_base == "LDM"){
              dend_result <- taxa.ldm.test.est.added
            }
            else{
              dend_result <- taxa.outputs$DAoutput
            }
            
            output$gdownload <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
              },
              content = function(file) {
                htmlwidgets::saveWidget(as_widget(taxa.sig.dend(dend_result, chooseData$NAadded$tax.tab, "twopi", include)$flow.text), file)
              }
            )
          }
          
          ref_string = REFERENCE_CHECK_P(data_transform = input$dataType_taxa, method_name = isolate(input$chooseMethod_taxa_base), FDR = "Yes")
          
          if (is.null(ref_string)) {
            shinyjs::hide("taxa_references")
          } else {
            shinyjs::show("taxa_references")
            output$taxa_references = renderUI({
              tagList(
                box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                    HTML(paste(ref_string, collapse="<br/>"))
                )
              )
            })
          }
          shinyjs::enable("runbtn_bin_taxa_base")
          shinyjs::enable("dataType_taxa")
          shinyjs::enable("primvar_taxa")
          shinyjs::enable("blockid_taxa")
          shinyjs::enable("chooseMethod_taxa")
          shinyjs::enable("chooseMethod_taxa_base")
          shinyjs::enable("covariates_taxa")
          shinyjs::enable("chooseRef_taxa")
          shinyjs::enable("taxaCat1")
          shinyjs::enable("taxaCat2")
        })
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(input$runbtn_bin_taxa,{
    prim_length = length(input$taxa_levels)
    for (i in 1:prim_length){
      shinyjs::disable(paste0("mult_taxaCat", i))
    }
    shinyjs::disable("runbtn_bin_taxa")
    shinyjs::disable("dataType_taxa")
    shinyjs::disable("primvar_taxa")
    shinyjs::disable("blockid_taxa")
    shinyjs::disable("chooseMethod_taxa")
    shinyjs::disable("chooseRef_taxa")
    shinyjs::disable("taxa_barPanel")
    shinyjs::hide("taxa_display_dend")
    shinyjs::show("taxa_adjustment")
    shinyjs::show("lmm_level")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        if (input$include_species.dend == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        
        taxa.bin.categos = c() 
        
        for (i in 1:prim_length){
          name_ind <- eval(parse(text = paste0("input$mult_taxaCat", i)))
          taxa.bin.categos <- c(taxa.bin.categos, name_ind)
        }
        
        sam.dat_mult <- try(bin.cat.recode.func.mult(chooseData$sam.dat, input$primvar_taxa, input$taxa_levels, taxa.bin.categos), silent = TRUE) 
        
        taxa.bin.out <- try(taxa.bin.cat.ref.united.func.mult.pairwise(input$primvar_taxa, sam.dat_mult, taxa = chooseData$taxa.out[[taxa.types$dataType]]), silent = TRUE) 
        taxa.results$bin.var <- taxa.bin.out$bin.var
        taxa.results$taxa <- taxa.bin.out$taxa
        taxa.results$taxa.bin.sum.out <- try(taxa.bin.sum.united.func(taxa.results$bin.var, taxa.results$taxa), silent = TRUE) 
        taxa_dataBinvar <- taxa.results$bin.var
        taxa_dataTaxa <- taxa.results$taxa
        
        re_sam <- try(reduced_data(sam.dat_mult, input$primvar_taxa, taxa.bin.categos), silent = TRUE)
        re_tax <- try(reduced_taxa(taxa.results$taxa, sam.dat_mult, input$primvar_taxa, input$blockid_taxa, taxa.bin.categos), silent = TRUE)
        
        if (input$chooseMethod_taxa == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
          incProgress(5/10, message = "ANOVA F-test (global) with Tukey's HSD (pairwise)")
          
          taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa, include), silent = TRUE)
          f.overall.result <- tryCatch(taxa.f.pair.mult.overall(taxa.mani), 
                                       error = function(e) {
                                         message ("No outcome is available!")
                                         showModal(modalDialog(div("No outcome is available!")))
                                         return(NULL)})
          
          taxa.rmanova.out.ori <- try(f.overall.result$global_test, silent = TRUE)
          global.p.val.only <- try(f.overall.result$pval, silent = TRUE)
          
          taxa.rmanova.out <- tryCatch(q_val_combined_table.anova(taxa.rmanova.out.ori, global.p.val.only), 
                                       error = function(e) {
                                         message ("No outcome is available!")
                                         showModal(modalDialog(div("No outcome is available!")))
                                         return(NULL)})
          
          taxa.outputs$DAoutput <- taxa.rmanova.out
          
          pairwise_result.ori <- tryCatch(taxa.f.pair.mult.tukey(taxa.mani), 
                                          error = function(e) {
                                            message ("No outcome is available!")
                                            showModal(modalDialog(div("No outcome is available!")))
                                            return(NULL)})
          
          tukey.p.list <- try(make_p_val_list(pairwise_result.ori), silent = TRUE)
          tukey.q.list <- try(make_q_val_list(tukey.p.list), silent = TRUE)
          pairwise_result <- try(q_val_dat(tukey.q.list, pairwise_result.ori, FALSE), silent = TRUE)
          pairwise_result_fc <- try(q_val_dat_fc(tukey.q.list, pairwise_result.ori), silent = TRUE) 
          
          log2_result <- try(log2_fold_change(taxa.mani), silent = TRUE) 
          sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
          
          volcanodat <- try(volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
          volcano.dat <- try(volcano_list_to_dat(volcanodat, species = include), silent = TRUE) 
          
          if(length(taxa.bin.categos) == 3){
            pvalue_list <- try(make_tripe_pvalues(taxa.rmanova.out, pairwise_result_fc, taxa.mani, include), silent = TRUE) 
            pvalue_dat_3d <- try(triple_list_to_dat(pvalue_list), silent = TRUE) 
          }
          
          num <- 5 + include
          
          nrow = numeric()
          for (r in 1:num) {
            row.num <- ceiling(sum(global.p.val.only[[r]][complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)  
            if (row.num > 0) {
              nrow[r] <- row.num
            } else {
              nrow[r] <- 1
            } 
          }
          
        }else if (input$chooseMethod_taxa == "Friedman's test (global) with Conover's test (pairwise)"){
          shinyjs::show("taxa_adjustment")
          incProgress(5/10, message = "Friedman's test (global) with Conover's test (pairwise)")
          
          taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa, include), silent = TRUE) 
          
          friedman.result <- tryCatch(taxa.friedman.pair.mult.overall(taxa.mani), 
                                      error = function(e) {
                                        message ("No outcome is available!")
                                        showModal(modalDialog(div("No outcome is available!")))
                                        return(NULL)})
          
          taxa.friedman.out.ori <- try(friedman.result$global_test, silent = TRUE) 
          global.p.val.only <- try(friedman.result$pval, silent = TRUE)
          
          pairwise_result.ori <- tryCatch(taxa.pair.mult.conover(taxa.mani, p_adjustment = "BH"),  
                                          error = function(e) {
                                            message ("No outcome is available!")
                                            showModal(modalDialog(div("No outcome is available!")))
                                            return(NULL)})
          
          conover.p.list <- try(make_p_val_list(pairwise_result.ori), silent = TRUE)
          
          conover.q.list <- try(make_q_val_list(conover.p.list), silent = TRUE)
          pairwise_result <- try(q_val_dat(conover.q.list, pairwise_result.ori), silent = TRUE)
          pairwise_result_fc <- try(q_val_dat_fc(conover.q.list, pairwise_result.ori), silent = TRUE) 
          
          taxa.friedman.out <- try(q_val_combined_table.friedman(taxa.friedman.out.ori, global.p.val.only), silent = TRUE)
          
          taxa.outputs$DAoutput <- taxa.friedman.out
          
          log2_result <- log2_fold_change(taxa.mani)
          sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
          
          volcanodat <- try(volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
          volcano.dat <- try(volcano_list_to_dat(volcanodat, species = include), silent = TRUE) 
          
          if(length(taxa.bin.categos) == 3){
            pvalue_list <- try(make_tripe_pvalues(taxa.friedman.out, pairwise_result_fc, taxa.mani, include), silent = TRUE) 
            pvalue_dat_3d <- try(triple_list_to_dat(pvalue_list), silent = TRUE) 
          }
          
          num <- 5 + include
          nrow = numeric()
          for (r in 1:num) {
            row.num <- ceiling(sum(global.p.val.only[[r]][complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)  
            if (row.num > 0) {
              nrow[r] <- row.num
            } else {
              nrow[r] <- 1
            } 
          }
        }
        
        else if(input$chooseMethod_taxa == "Durbin's test (global) with Conover's test (pairwise)"){
          shinyjs::show("taxa_adjustment")
          incProgress(5/10, message = "Durbin's test (global) with Conover's test (pairwise)")
          
          taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa, include), silent = TRUE) 
          
          taxa.durbin.result <- tryCatch(taxa.durbin.pair.mult.overall(taxa.mani),  
                                         error = function(e) {
                                           message ("No outcome is available!")
                                           showModal(modalDialog(div("No outcome is available!")))
                                           return(NULL)})
          
          taxa.durbin.out.ori <- try(taxa.durbin.result$global_test, silent = TRUE)
          global.p.val.only <- try(taxa.durbin.result$pval, silent = TRUE)
          taxa.durbin.out <- try(q_val_combined_table.friedman(taxa.durbin.out.ori, global.p.val.only), silent = TRUE)
          
          pairwise_result.ori <-tryCatch(taxa.pairwise.mult.durbin (taxa.mani, p_adjustment = "BH"),  
                                         error = function(e) {
                                           message ("No outcome is available!")
                                           showModal(modalDialog(div("No outcome is available!")))
                                           return(NULL)}) #here 
          
          durbin.p.list <- try(make_p_val_list(pairwise_result.ori), silent = TRUE)
          durbin.q.list <- try(make_q_val_list(durbin.p.list), silent = TRUE)
          pairwise_result <- try(q_val_dat(durbin.q.list, pairwise_result.ori), silent = TRUE)
          taxa.outputs$DAoutput <- taxa.durbin.out
          pairwise_result_fc <- try(q_val_dat_fc(durbin.q.list, pairwise_result.ori), silent = TRUE) 
          
          log2_result <- try(log2_fold_change(taxa.mani), silent = TRUE) 
          sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
          
          volcanodat <- try(volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
          volcano.dat <- try(volcano_list_to_dat(volcanodat, species = include), silent = TRUE) 
          
          if(length(taxa.bin.categos) == 3){
            pvalue_list <- try(make_tripe_pvalues(taxa.durbin.out, pairwise_result_fc, taxa.mani, include), silent = TRUE) 
            pvalue_dat_3d <- try(triple_list_to_dat(pvalue_list), silent = TRUE) 
          }
          
          num <- 5 + include
          nrow = numeric()
          for (r in 1:num) {
            row.num <- ceiling(sum(global.p.val.only[[r]][complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)    
            if (row.num > 0) {
              nrow[r] <- row.num
            } else {
              nrow[r] <- 1
            } 
          }
        }
        
        else if (input$chooseMethod_taxa == "LDM"){
          incProgress(5/10, message = "LDM")
          shinyjs::show("taxa_adjustment")
          
          if (input$dataType_taxa == "Arcsine-root"){
            
            if(include){
              re_tax_2 <- re_tax 
            }else{
              re_tax_2 <- list(re_tax[[1]], re_tax[[2]], re_tax[[3]], re_tax[[4]], re_tax[[5]])
              names(re_tax_2) <- c("phylum", "class", "order", "family", "genus")
            }
            
            
            set.seed(521)
            taxa.ldm.result <- tryCatch(taxa.ldm.mult.glob(re_sam, re_tax_2, input$primvar_taxa, input$blockid_taxa, n.perm = 3000, trans = TRUE), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            taxa.ldm.test.q.out <- try(taxa.ldm.result$global_test, silent = TRUE) 
            taxa.outputs$DAoutput = taxa.ldm.test.q.out
            global.p.val.only <- try(taxa.ldm.result$pval, silent = TRUE) 
            
            set.seed(521)
            
            pairwise_result.ori <- tryCatch(taxa.ldm.paired.mult.united_trans(re_sam, re_tax_2, input$primvar_taxa, input$blockid_taxa, n.perm = 3000), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            pairwise_result <- try(q_val_convert(pairwise_result.ori), silent = TRUE)
            
            pairwise_result_fc <- try(q_val_convert_qc(pairwise_result.ori), silent = TRUE) 
            
            taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa, include), silent = TRUE) 
            log2_result <- try(log2_fold_change(taxa.mani), silent = TRUE) 
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            
            volcanodat <- try(volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
            volcano.dat <- try(volcano_list_to_dat(volcanodat, species = include), silent = TRUE) 
            
            if(length(taxa.bin.categos) == 3){
              pvalue_list <- try(make_tripe_pvalues(taxa.ldm.test.q.out, pairwise_result_fc, taxa.mani, include), silent = TRUE) 
              pvalue_dat_3d <- try(triple_list_to_dat(pvalue_list), silent = TRUE) 
            }
            
          }else{
            if(include){
              re_tax_2 <- re_tax 
            }else{
              re_tax_2 <- list(re_tax[[1]], re_tax[[2]], re_tax[[3]], re_tax[[4]], re_tax[[5]])
              names(re_tax_2) <- c("phylum", "class", "order", "family", "genus")
            }
         
            set.seed(521)
            
            taxa.ldm.result <- tryCatch(taxa.ldm.mult.glob(re_sam, re_tax_2, input$primvar_taxa, input$blockid_taxa, n.perm = 3000, trans = FALSE), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            taxa.ldm.test.q.out <- try(taxa.ldm.result$global_test, silent = TRUE) 
            taxa.outputs$DAoutput = taxa.ldm.test.q.out
            global.p.val.only <- try(taxa.ldm.result$pval, silent = TRUE) 
            
            set.seed(521)
            pairwise_result.ori <- try(taxa.ldm.paired.mult.united(re_sam, re_tax_2, input$primvar_taxa, input$blockid_taxa, n.perm = 3000), silent = TRUE)
            pairwise_result <- try(q_val_convert(pairwise_result.ori), silent = TRUE)
            
            pairwise_result_fc <- try(q_val_convert_qc(pairwise_result.ori), silent = TRUE) 
            
            taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa, include), silent = TRUE) 
            log2_result <- try(log2_fold_change(taxa.mani), silent = TRUE) 
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            
          
            volcanodat <- try(volcano_dat(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
            volcano.dat <- try(volcano_list_to_dat(volcanodat, FALSE, species = include), silent = TRUE) 
            
            if(length(taxa.bin.categos) == 3){
              pvalue_list <- try(make_tripe_pvalues(taxa.ldm.test.q.out, pairwise_result_fc, taxa.mani, include), silent = TRUE)
              pvalue_dat_3d <- try(triple_list_to_dat(pvalue_list), silent = TRUE) 
            }
          }
          
          
          num <- 5 + include
          nrow = numeric()
          for (r in 1:num) {
            print(r)
            row.num <- ceiling(sum(taxa.ldm.test.q.out[[r]]$P.value[as.numeric(taxa.ldm.test.q.out[[r]]$P.value)] < 0.05)/4) 
            if (row.num > 0) {
              nrow[r] <- row.num
            } else {
              nrow[r] <- 1
            } 
          }
          
        }else if (input$chooseMethod_taxa == "LMM (LRT (global) with t-test (pairwise))"){   
          shinyjs::show("taxa_adjustment")
          shinyjs::show("lmm_level")
          
          incProgress(5/10, message = "LMM (LRT (global) with t-test (pairwise))")
          
          taxa.categos <- input$taxa_levels
          
          new_ind <- which(taxa.categos == input$lmm_level)
          new_lmm_level <- taxa.bin.categos[new_ind]
          taxa.mani_2 <- try(taxa_mani(re_tax, re_sam, input$primvar_taxa), silent = TRUE)
          
          if (input$dataType_taxa == "Arcsine-root"){
            taxa.mani_2 <- try(no_constant_acr(taxa.mani_2, sam.dat_mult))
            taxa.lmm.out.ori <- tryCatch(taxa_total(taxa.mani_2, sam.dat_mult, input$primvar_taxa, input$blockid_taxa, new_lmm_level, "BH"), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            taxa.lmm.out <- try(q_val_convert(taxa.lmm.out.ori), silent = TRUE)
            
            taxa.outputs$DAoutput = try(global_total_list(re_tax, taxa.mani_2, input$primvar_taxa, input$blockid_taxa, new_lmm_level, TRUE), silent = TRUE)
            global.p.val.only <- try(global_total_p(taxa.mani_2, sam.dat_mult,  input$primvar_taxa, input$blockid_taxa, new_lmm_level), silent = TRUE)
            pairwise_result <- taxa.lmm.out 
            
            pairwise_result_fc <- try(q_val_convert_qc(pairwise_result), silent = TRUE) 
            
            taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa), silent = TRUE) 
            log2_result <- try(log2_fold_change_lmm(taxa.mani_2, new_lmm_level), silent = TRUE) 
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            
            volcanodat <- try(volcano_dat_lmm(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
            volcano.dat <- try(volcano_list_to_dat(volcanodat, TRUE, include), silent = TRUE) 
            
            num <- 5 + include
            nrow = numeric()
            for (r in 1:num) {
              row.num <- ceiling(sum(global.p.val.only[[r]][complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
            
          }else{
            taxa.lmm.out.ori <<- tryCatch(taxa_total(taxa.mani_2, sam.dat_mult, input$primvar_taxa, input$blockid_taxa, new_lmm_level, "BH"), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            taxa.lmm.out <- try(q_val_convert(taxa.lmm.out.ori), silent = TRUE)
            taxa.outputs$DAoutput <- try(global_total_list(re_tax, taxa.mani_2, input$primvar_taxa, input$blockid_taxa, new_lmm_level), silent = TRUE)
            global.p.val.only <- try(global_total_p(taxa.mani_2,  sam.dat_mult,  input$primvar_taxa, input$blockid_taxa, new_lmm_level), silent = TRUE) #여기 살펴보렴 
            pairwise_result <- taxa.lmm.out 
            pairwise_result_fc <- try(q_val_convert_qc(pairwise_result), silent = TRUE) 
            
            taxa.mani <- try(taxa.data.mani(re_sam, re_tax, input$primvar_taxa, input$blockid_taxa), silent = TRUE) 
            log2_result <- try(log2_fold_change_lmm(taxa.mani, new_lmm_level), silent = TRUE) 
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            
            volcanodat <- try(volcano_dat_lmm(pairwise_result_fc, log2_result, sig_global_taxa), silent = TRUE) 
            volcano.dat <- try(volcano_list_to_dat(volcanodat, TRUE, include), silent = TRUE)     #여기 수정필요 
            
            
            num <- 5 + include
            nrow = numeric()
            for (r in 1:num) {
              row.num <- ceiling(sum(global.p.val.only[[r]][complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          }
        }
        
        if (length(input$taxa_levels) < 4){
          width_boxplot <- 750 
        }else if (4 <= length(input$taxa_levels) | length(input$taxa_levels) <= 6){
          width_boxplot <- 800 
        }else {
          width_boxplot <- 850 
        }
        
        if(input$chooseMethod_taxa == "LMM (LRT (global) with t-test (pairwise))"){
          box_title <- "LMM (Global Test)"
        }else if(input$chooseMethod_taxa == "LDM"){
          box_title <- "LDM (Global Test)"
        }else if(input$chooseMethod_taxa == "Durbin's test (global) with Conover's test (pairwise)"){
          box_title <- "Durbin's test (Global Test)"
        }else if(input$chooseMethod_taxa == "Friedman's test (global) with Conover's test (pairwise)"){
          box_title <- "Friedman's test (Global Test)"
        }else if(input$chooseMethod_taxa == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
          box_title <- "ANOVA F-test (Global Test)"
        }
        
        if(input$include_species.dend == "Phylum - Species"){
          
          incProgress(3/10, message = "Displaying Results in progress")
          output$taxa_display = renderUI({
            tagList(
              tabBox(title = strong(box_title, style = "color:black"), width = NULL,
                     tabPanel("Phylum", align = "center",
                              plotOutput("rank1", height = nrow[1]*250, width = width_boxplot),
                     )
                     ,
                     tabPanel("Class", align = "center",
                              plotOutput("rank2", height = nrow[2]*250, width = width_boxplot),
                     )
                     ,tabPanel("Order", align = "center",
                               plotOutput("rank3", height = nrow[3]*250, width = width_boxplot),
                     )
                     ,tabPanel("Family", align = "center",
                               plotOutput("rank4", height = nrow[4]*250, width = width_boxplot),
                     )
                     ,tabPanel("Genus", align = "center",
                               plotOutput("rank5", height = nrow[5]*250, width = width_boxplot),
                     )
                     ,tabPanel("Species", align = "center",
                               plotOutput("rank6", height = nrow[6]*250, width = width_boxplot),
                     )
              )
            )
          })
          
          output$rank1 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 1, global.p.val.only)
          }), silent = TRUE)
          
          output$rank2 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 2, global.p.val.only)
          }), silent = TRUE)
          
          output$rank3 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 3, global.p.val.only)
          }), silent = TRUE)
          
          output$rank4 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 4, global.p.val.only)
          }), silent = TRUE)
          
          output$rank5 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 5, global.p.val.only)
          }), silent = TRUE)
          
          output$rank6 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 6, global.p.val.only)
          }), silent = TRUE)
          
          sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
          sig_pair <- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
          sig_pair_dat <- try(list_to_dat(sig_pair), silent = TRUE)
          
          shinyjs::show("taxa_barPanel")
          
          if (input$chooseMethod_taxa == "LDM"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "LDM (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
          }else if(input$chooseMethod_taxa == "Friedman's test (global) with Conover's test (pairwise)"){
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Conover's test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
          } else if (input$chooseMethod_taxa == "Durbin's test (global) with Conover's test (pairwise)"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Conover's test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
          }
          else if(input$chooseMethod_taxa == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Tukey's HSD (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
          }else if(input$chooseMethod_taxa == "LMM (LRT (global) with t-test (pairwise))"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "LMM (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
          }
          
          length_each_page <- nrow(sig_pair_dat[[1]])/length(sig_pair[[1]])
          
          length_page_1 <- nrow(sig_pair_dat[[1]])/length(sig_pair[[1]])
          length_page_2 <-  nrow(sig_pair_dat[[2]])/length(sig_pair[[2]])
          length_page_3 <- nrow(sig_pair_dat[[3]])/length(sig_pair[[3]])
          length_page_4 <- nrow(sig_pair_dat[[4]])/length(sig_pair[[4]])
          length_page_5 <- nrow(sig_pair_dat[[5]])/length(sig_pair[[5]])
          length_page_6 <- nrow(sig_pair_dat[[6]])/length(sig_pair[[6]])
         
          
          for (i in 1:length(sig_pair_dat)){
            colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
            if("P.value" %in% colnames(sig_pair_dat[[i]])){
              ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
              sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
            }
          }
          
          
          output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
          output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
          output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
          output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
          output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
          output$species_list <- try(renderDataTable({datatable(sig_pair_dat[[6]], rownames = FALSE, options = list(pageLength = length_page_6, dom = '<"top" p>'))}), silent = TRUE)
          
          nrow_volcano <- ceiling(length(taxa.bin.categos)/3)
          
          if(length(taxa.bin.categos) == 3 & input$chooseMethod_taxa != "LMM (LRT (global) with t-test (pairwise))"){
            output$taxa_dislay_pairwise <- renderUI({
              tagList(
                tabBox(title = strong("Volcano Plot", style = "color:black"), width = NULL, height = 500, 
                       tabPanel("2D", plotlyOutput("volcanoplot")), 
                       tabPanel("3D", plotlyOutput("volcanoplot_3d"))
                ))})
            
            dat_combined <- cbind(taxa.mani[[1]][,3:ncol(taxa.mani[[1]])], taxa.mani[[2]][,3:ncol(taxa.mani[[2]])], taxa.mani[[3]][,3:ncol(taxa.mani[[3]])], taxa.mani[[4]][,3:ncol(taxa.mani[[4]])], taxa.mani[[5]][,3:ncol(taxa.mani[[5]])], taxa.mani[[6]][,3:ncol(taxa.mani[[6]])])
            
            labels_3d <- try(rename_volcano_3d(factor(taxa.mani[[1]][,1], levels = names(table(taxa.mani[[1]][,1])))), silent = TRUE)
            syn_polar <- try(polar_coords(factor(taxa.mani[[1]][,1], levels = names(table(taxa.mani[[1]][,1]))), data = dat_combined, pvals = as.matrix(pvalue_dat_3d), labs = labels_3d), silent = TRUE) 
          
            print("here")
            print(volcano.dat)
            
            volcano_visualize_individual(volcano.dat) 
            
            output$volcanoplot <- renderPlotly(volcano <- volcano_visualize_individual(volcano.dat, include))
            output$volcanoplot_3d <- renderPlotly(volcano3d <- volcano3D_RE(syn_polar, type = 2))  
          }else {
            output$taxa_display_pairwise <- renderUI({
              tagList(
                tabBox(title = strong("Volcano Plot", style = "color:black"), width = NULL, height = 400*nrow_volcano,
                       tabPanel("2D", plotlyOutput("volcanoplot", height = paste0(400*nrow_volcano, "px")))
                ))})
            output$volcanoplot <- renderPlotly(volcano <- volcano_visualize_individual(volcano.dat, include))
          }
          
          incProgress(1/10, message = "Displaying Results in progress")
          output$downloadTable_taxa = renderUI({
            tagList(
              p(" ", style = "margin-bottom: +5px;"),
              box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  p("You can download the summary statistics and data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Summary Statistics"),
                  downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                  h5("Data Analysis Outputs"),
                  downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$tdownloadTabl1 <- downloadHandler(
            filename = function() {
              paste("Taxa.Sum.Table.zip")
            },
            content = function(sum.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              zip(zipfile=sum.file, files=dataFiles)
            }
          )
          
          output$tdownloadTabl2 <- downloadHandler(
            filename = function() {
              paste("Taxa.Analysis.Output.zip")
            }, 
            content = function(DA.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
              
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$phylum), Pairwise = pairwise_result_fc[[1]]), file = "Phylum.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$class), Pairwise = pairwise_result_fc[[2]]), file = "Class.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$order), Pairwise = pairwise_result_fc[[3]]), file = "Order.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$family), Pairwise = pairwise_result_fc[[4]]), file = "Family.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$genus), Pairwise = pairwise_result_fc[[5]]), file = "Genus.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$species), Pairwise = pairwise_result_fc[[6]]), file = "Species.txt")
              zip(zipfile=DA.file, files=dataFiles)
            }
          )
          
        }else{
          
          incProgress(3/10, message = "Displaying Results in progress")
          output$taxa_display = renderUI({
            tagList(
              tabBox(title = strong(box_title, style = "color:black"), width = NULL,
                     tabPanel("Phylum", align = "center",
                              plotOutput("rank1", height = nrow[1]*250, width = width_boxplot),
                     )
                     ,
                     tabPanel("Class", align = "center",
                              plotOutput("rank2", height = nrow[2]*250, width = width_boxplot),
                     )
                     ,tabPanel("Order", align = "center",
                               plotOutput("rank3", height = nrow[3]*250, width = width_boxplot),
                     )
                     ,tabPanel("Family", align = "center",
                               plotOutput("rank4", height = nrow[4]*250, width = width_boxplot),
                     )
                     ,tabPanel("Genus", align = "center",
                               plotOutput("rank5", height = nrow[5]*250, width = width_boxplot),
                     )
                     
              )
            )
          })
          
          output$rank1 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 1, global.p.val.only)
          }), silent = TRUE)
          
          output$rank2 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 2, global.p.val.only)
          }), silent = TRUE)
          
          output$rank3 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 3, global.p.val.only)
          }), silent = TRUE)
          
          output$rank4 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 4, global.p.val.only)
          }), silent = TRUE)
          
          output$rank5 = try(renderPlot({ 
            taxa.bin.mult.paired (re_tax, re_sam, input$primvar_taxa, input$blockid_taxa, 5, global.p.val.only)
          }), silent = TRUE)
          
          sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
          sig_pair <- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
          sig_pair_dat <- try(list_to_dat(sig_pair), silent = TRUE)
          
          shinyjs::show("taxa_barPanel")
          
          if (input$chooseMethod_taxa == "LDM"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "LDM (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
          }else if(input$chooseMethod_taxa == "Friedman's test (global) with Conover's test (pairwise)"){
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Conover's test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
          } else if (input$chooseMethod_taxa == "Durbin's test (global) with Conover's test (pairwise)"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Conover's test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
          }
          else if(input$chooseMethod_taxa == "ANOVA F-test (global) with Tukey's HSD (pairwise)"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Tukey's HSD (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
          }else if(input$chooseMethod_taxa == "LMM (LRT (global) with t-test (pairwise))"){
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "LMM (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
          }
          
          length_each_page <- nrow(sig_pair_dat[[1]])/length(sig_pair[[1]])
          
          length_page_1 <- nrow(sig_pair_dat[[1]])/length(sig_pair[[1]])
          length_page_2 <-  nrow(sig_pair_dat[[2]])/length(sig_pair[[2]])
          length_page_3 <- nrow(sig_pair_dat[[3]])/length(sig_pair[[3]])
          length_page_4 <- nrow(sig_pair_dat[[4]])/length(sig_pair[[4]])
          length_page_5 <- nrow(sig_pair_dat[[5]])/length(sig_pair[[5]])
          
          
          for (i in 1:length(sig_pair_dat)){
            colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
            if("P.value" %in% colnames(sig_pair_dat[[i]])){
              ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
              sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
            }
          }
          
          output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
          output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
          output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
          output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
          output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
          
          nrow_volcano <- ceiling(length(taxa.bin.categos)/3)
          
          if(length(taxa.bin.categos) == 3 & input$chooseMethod_taxa != "LMM (LRT (global) with t-test (pairwise))"){
            output$taxa_display_pairwise <- renderUI({
              tagList(
                tabBox(title = strong("Volcano Plot", style = "color:black"), width = NULL, height = 500, 
                       tabPanel("2D", plotlyOutput("volcanoplot")), 
                       tabPanel("3D", plotlyOutput("volcanoplot_3d"))
                ))})
            
            dat_combined <- cbind(taxa.mani[[1]][,3:ncol(taxa.mani[[1]])], taxa.mani[[2]][,3:ncol(taxa.mani[[2]])], taxa.mani[[3]][,3:ncol(taxa.mani[[3]])], taxa.mani[[4]][,3:ncol(taxa.mani[[4]])], taxa.mani[[5]][,3:ncol(taxa.mani[[5]])])
            
            labels_3d <- try(rename_volcano_3d(factor(taxa.mani[[1]][,1], levels = names(table(taxa.mani[[1]][,1])))), silent = TRUE) 
            syn_polar <- try(polar_coords(factor(taxa.mani[[1]][,1], levels = names(table(taxa.mani[[1]][,1]))), data = dat_combined, pvals = as.matrix(pvalue_dat_3d), labs = labels_3d), silent = TRUE)
            
            output$volcanoplot <- renderPlotly(volcano <- volcano_visualize_individual(volcano.dat, include))
            output$volcanoplot_3d <- renderPlotly(volcano3d <- volcano3D_RE(syn_polar, type = 2))  
            
          }else {
            output$taxa_display_pairwise <- renderUI({
              tagList(
                tabBox(title = strong("Volcano Plot", style = "color:black"), width = NULL, height = 400*nrow_volcano,
                       tabPanel("2D", plotlyOutput("volcanoplot", height = paste0(400*nrow_volcano, "px")))
                ))})
            output$volcanoplot <- renderPlotly(volcano <- volcano_visualize_individual(volcano.dat, include))
          }
          
          print(names(volcano.dat)) 
                    
          incProgress(1/10, message = "Displaying Results in progress")
          output$downloadTable_taxa = renderUI({
            tagList(
              p(" ", style = "margin-bottom: +5px;"),
              box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  p("You can download the summary statistics and data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Summary Statistics"),
                  downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                  h5("Data Analysis Outputs"),
                  downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$tdownloadTabl1 <- downloadHandler(
            filename = function() {
              paste("Taxa.Sum.Table.zip")
            },
            content = function(sum.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
              zip(zipfile=sum.file, files=dataFiles)
            }
          )
          
          output$tdownloadTabl2 <- downloadHandler(
            filename = function() {
              paste("Taxa.Analysis.Output.zip")
            }, 
            content = function(DA.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
              
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$phylum), Pairwise = pairwise_result_fc[[1]]), file = "Phylum.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$class), Pairwise = pairwise_result_fc[[2]]), file = "Class.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$order), Pairwise = pairwise_result_fc[[3]]), file = "Order.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$family), Pairwise = pairwise_result_fc[[4]]), file = "Family.txt")
              capture.output(list(Global = as.data.frame(taxa.outputs$DAoutput$genus), Pairwise = pairwise_result_fc[[5]]), file = "Genus.txt")
              zip(zipfile=DA.file, files=dataFiles)
            }
          )
          
        }
        ref_string = REFERENCE_CHECK_P(data_transform = input$dataType_taxa, method_name = isolate(input$chooseMethod_taxa), FDR = "Yes")
        
        if (is.null(ref_string)) {
          shinyjs::hide("taxa_references")
        } else {
          shinyjs::show("taxa_references")
          output$taxa_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "warning", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("dataType_taxa")
        shinyjs::enable("runbtn_bin_taxa")
        shinyjs::enable("primvar_taxa")
        shinyjs::enable("blockid_taxa")
        shinyjs::enable("chooseMethod_taxa")
        shinyjs::enable("covariates_taxa")
        shinyjs::enable("chooseRef_taxa")
        for (i in 1:prim_length){
          shinyjs::enable(paste0("mult_taxaCat", i))
        }
        shinyjs::enable("taxa_barPanel")
        shinyjs::enable("taxa_adjustment")
        shinyjs::enable("lmm_level")
        
      })
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
}

shinyApp(ui = ui, server = server)

