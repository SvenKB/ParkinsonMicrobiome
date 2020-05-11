#### Specify needed packages ####
packages <- c("ggplot2","tidyverse","drake",'rjson','knitr','kableExtra')

#### Install packages ####
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#### Load packages ####
lapply(packages, library, character.only = TRUE)
