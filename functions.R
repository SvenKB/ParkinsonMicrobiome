#### Specify custom functions ####


#### Helper functions 
`%notin%` <- Negate(`%in%`)
`not.na` <- Negate(is.na)

#### Data preperation functions ####
tax <- c('Kingdom','Phylum','Class','Order','Family','Genus',"sequence")

prepare_taxonomy <- function(data) {
  
  d <- data %>% select(ID,all_of(tax)) %>% mutate(ID = paste0("ASV",ID))
  return(d)
}

## Transform the data into a tidy format
tidy_up <- function(data,trim=FALSE) {
  d <- data %>%
    select(-all_of(tax)) %>%
    rownames_to_column %>% 
    tidyr::gather(var, value, -rowname) %>% 
    tidyr::spread(rowname, value) %>%
    filter(var!="ID") %>%
    setNames(paste0('ASV', names(.))) %>%
    rename(ID=ASVvar)
  
  if(trim==TRUE) {d$ID <- unlist(strsplit(d$ID, split='-', fixed=TRUE))[seq(2, nrow(d)*2, 2)]}
  
  return(d)
}


## Load JSON from ENA 
loadENAMeta <- function(path) {
  d <- rjson::fromJSON(file=path,simplify = F)
  col_name <- names(d[[1]])
  
  meta <- as_tibble(data.frame(matrix(unlist(d), nrow=length(d), byrow=T),stringsAsFactors=FALSE))
  colnames(meta) <- col_name
  return(meta)
}
  

##  Data inspection

inspect_depth <- function(d) {
  names <- d %>% select(ID) %>% pull
  d <- d %>% select(-ID)
  depth <- apply(d,1,sum)
  names(depth) <- names
  return(depth)
}

inspect_reads <- function(d) {
  d <- d %>% select(-ID)
  count <- apply(d,2,sum)
  return(count)
}

 inspect_sparsity <- function(d,sp=c("asv")) {
  names <- d %>% select(ID) %>% pull
  d <- d %>% select(-ID)
  
  if(sp=="obs")   {spar <- apply(d,1,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='asv')   {spar <- apply(d,2,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='total') {spar <- sum(d==0) / (nrow(d)*ncol(d))*100}
 
  return(spar)
}


## Data cleaning functions


clean <- function(data,tax,min_tax="Family",seq_length=c(), min_abund=2,quant=T) {
  
  ## Filter unwanted kingdoms
  tot_asv <- nrow(tax) # Record full number of ASVs
  ft_king <- tax %>% filter(Kingdom=="Bacteria") %>% select(ID) %>% pull # Create Filter
  ft_king_nr <- tot_asv - length(ft_king) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Kingdom=="Bacteria") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  d <- data %>% select(all_of(ft_king)) # Apply filter to data
  
  ## Filter Chloroplasts
  ft_chloroplast <- tax %>% filter(Class!="Chloroplast") %>% select(ID) %>% pull # Create Filter
  ft_chloroplast_nr <- new_asv - length(ft_chloroplast) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Class!="Chloroplast") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  d <- d %>% select(all_of(ft_chloroplast)) # Apply filter to data
  
  ## Filter based on low abundance AND low taxonomic information
  
  abund <- apply(d,2,sum) # Calculate Abundance
  if(quant==T) {thresh <- summary(abund)[[2]]} else {thresh <- min_abund}
  ft_abund <- tax %>% filter(abund<thresh) %>% select(ID) %>% pull
  ft_mintax <- tax %>% filter(ID%in%ft_abund) %>% filter(is.na(!!as.symbol(min_tax))) %>% select(ID) %>% pull
  ft_mintax_nr <- length(ft_mintax)
  
  tax <- tax %>% filter(ID %notin% ft_mintax)
  new_asv <- nrow(tax)
  d <- d %>% select(-all_of(ft_mintax))
  
  ## Filter minimum taxonomic information
  #ft_mintax <- tax %>% filter(not.na(!!as.symbol(min_tax))) %>% select(ID) %>% pull # Create Filter
  #ft_mintax_nr <- new_asv - length(ft_mintax) # Calculate nr. of filtered ASVs
  
  #tax <- tax %>% filter(not.na(!!as.symbol(min_tax)))
  #new_asv <- nrow(tax)
  #d <- d %>% select(all_of(ft_mintax))
  
  tot_filt <- ft_king_nr + ft_chloroplast_nr + ft_mintax_nr
  
  
  cat('Total number of ASVs:',tot_asv,"\n",
      ft_king_nr,'ASVs filtered due to wrong Kingdom\n',
      ft_chloroplast_nr,'ASVs of Class Chloroplast filtered\n',
      ft_mintax_nr,'ASVs filtered due to low abundance (<',thresh,') and missing information at rank',min_tax,'\n',
      'Number of filtered ASVs:',tot_filt,'\n',
      'New number of ASVs:',new_asv)
  return(d)
}

 
#### Analysis functions ####



#### Summary / plot functions ####
