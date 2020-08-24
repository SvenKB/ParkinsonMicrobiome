#### Specify needed packages ####
packages <- c("ggplot2","tidyverse","drake",'rjson','knitr','kableExtra',
              'readr','hablar','foreach','phyloseq','vegan','lme4','microbiome',
              'merTools',"sjPlot",'lmerTest','reshape2','sjmisc','sjmisc',
              'texreg','ape','microbiome','rmdformats','seqinr','hillR','metafor')

#### Install packages ####
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos = "http://cran.us.r-project.org")  
}

#### Load packages ####
lapply(packages, library, character.only = TRUE)
