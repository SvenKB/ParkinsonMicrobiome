#### Specify custom functions ####


#### Helper functions 
`%notin%` <- Negate(`%in%`)
`not.na` <- Negate(is.na)

#### Data preperation functions ####
tax <- c('Kingdom','Phylum','Class','Order','Family','Genus',"sequence")

prepare_taxonomy <- function(data) {
  
  d <- data %>% dplyr::select(ID,all_of(tax)) %>% mutate(ID = paste0("ASV",ID))
  return(d)
}

prepare_taxonomy_sameregion <- function(x) {
  x %>% dplyr::select("ID"=X1,"Kingdom"=X3,"Phylum"=X6,"Class"=X9,"Order"=X12,"Family"=X15,"Genus"=X18) %>% mutate(ID = paste0("ASV",ID))}

## Transform the data into a tidy format
tidy_up <- function(data,trim=FALSE) {
  d <- data %>%
    dplyr::select(-all_of(tax)) %>%
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
  names <- d %>% dplyr::select(ID) %>% pull
  d <- d %>% dplyr::select(-ID)
  depth <- apply(d,1,sum)
  names(depth) <- names
  return(depth)
}

inspect_reads <- function(d) {
  d <- d %>% dplyr::select(-ID)
  count <- apply(d,2,sum)
  return(count)
}

 inspect_sparsity <- function(d,sp=c("asv")) {
  names <- d %>% dplyr::select(ID) %>% pull
  d <- d %>% dplyr::select(-ID)
  
  if(sp=="obs")   {spar <- apply(d,1,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='asv')   {spar <- apply(d,2,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='total') {spar <- sum(d==0) / (nrow(d)*ncol(d))*100}
 
  return(spar)
}


## Data cleaning functions

## Currently, the function only filters based on relative abundance, not taxonomic information, as taxonomic information is present in almost all ASVs
 
 
clean <- function(data,tax,min_tax="Family",seq_length=c(), filter=c("abs","rel","quant"),abs=100,rel=0.001) {
  
  id <- data %>% dplyr::select(ID)
  
  ## Filter kingdoms
  tot_asv <- nrow(tax) # Track full number of ASVs
  
  ft_king <- tax %>% filter(Kingdom=="Bacteria") %>% dplyr::select(ID) %>% pull # Create Filter
  ft_king_nr <- tot_asv - length(ft_king) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Kingdom=="Bacteria") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  d <- data %>% dplyr::select(all_of(ft_king)) # Apply filter to data
  
  ## Filter Chloroplasts
  ft_chloroplast <- tax %>% filter(Class!="Chloroplast") %>% dplyr::select(ID) %>% pull # Create Filter
  ft_chloroplast_nr <- new_asv - length(ft_chloroplast) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Class!="Chloroplast") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  d <- d %>% dplyr::select(all_of(ft_chloroplast)) # Apply filter to data
  
  ## Filter based on low abundance AND low taxonomic information
  
  abund <- apply(d,2,sum) # Calculate Abundance
  
  if(filter=="quant") {ft_abund <- tax %>% filter(abund<summary(abund)[[3]]) %>% dplyr::select(ID) %>% pull # Create filter based on quantile
                       thresh <- summary(abund)[[3]]}
  if(filter == "abs") {ft_abund <- tax %>% filter(abund<abs) %>% dplyr::select(ID) %>% pull # Create filter based on absolute reads
                       thresh <- abs}
  if(filter == "rel") {ind <- (abund/sum(d))*100<rel                                  # Create filter based on relative abundance
                       ft_abund <- tax %>% filter(ind) %>% dplyr::select(ID) %>% pull
                       thresh <- paste0(rel,"%")}
  
  # Apply combined filter -> only if abundance threshold AND missing annotation are are present, ASV will be filtered !: MAYBE BAD IDEA?
  # %>% filter(is.na(!!as.symbol(min_tax)))
  
  ft_mintax <- tax %>% filter(ID%in%ft_abund|is.na(!!as.symbol(min_tax))) %>% dplyr::select(ID) %>% pull
  ft_mintax_nr <- length(ft_mintax)
  
  ft_abund_nr <- length(ft_abund)
  
  tax <- tax %>% filter(ID %notin% ft_mintax)
  new_asv <- nrow(tax)
  d <- d %>% dplyr::select(-all_of(ft_mintax))
  
  
  ## Filter based on sequence length
  
  seq <- tax %>% dplyr::select(Sequence) %>% convert(chr(Sequence)) %>% pull
  seq_length <- sapply(str_split(seq,''),length)
  ft_seqlength <- tax %>% filter(seq_length > 150) %>% dplyr::select(ID) %>% pull
  ft_seqlength_nr <- nrow(tax)-length(ft_seqlength)
  tax <- tax %>% filter(ID %in% ft_seqlength)
  new_asv <- nrow(tax)
  
  d <- d %>% dplyr::select(all_of(ft_seqlength))
  
  tot_filt <- ft_king_nr + ft_chloroplast_nr + ft_mintax_nr + ft_seqlength_nr
  
  d <- data.frame(id,d)
  
  cat('Total number of ASVs:',tot_asv,"\n",
      ft_king_nr,'ASVs filtered due to wrong Kingdom\n',
      ft_chloroplast_nr,'ASVs of Class Chloroplast filtered\n',
      ft_mintax_nr,'ASVs filtered due to low abundance (<',thresh,')\n',
      "(Number of ASVs below threshold, but not filtered:",ft_abund_nr,"\n",
      ft_seqlength_nr,'ASVs filtered due to low sequence length < 150 \n',
      'Number of filtered ASVs:',tot_filt,'\n',
      'New number of ASVs:',new_asv,"\n")
  return(d)
}



# Prepare sameregion taxonomy


 
#### Analysis functions ####



#### Summary / plot functions ####
