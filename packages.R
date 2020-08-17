#### Specify needed packages ####
packages <- c("ggplot2","tidyverse","drake",'rjson','knitr','kableExtra',
              'readr','hablar','foreach','phyloseq','vegan','lme4','microbiome',
              'merTools',"sjPlot",'lmerTest','reshape2','sjmisc','sjmisc',
              'texreg','ape','microbiome')

#### Install packages ####
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#### Load packages ####
lapply(packages, library, character.only = TRUE)
