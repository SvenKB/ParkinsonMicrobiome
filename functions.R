#### Specify custom functions ####

##########################
#### Helper functions ####
##########################

`%notin%` <- Negate(`%in%`)
`not.na` <- Negate(is.na)

####################################
#### Data preperation functions ####
####################################

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
 
prepare_DA_data <- function(data,tax="Genus",meta=c("status","author")) {
   
   d <- lapply(dat, function(x) as.data.frame(t(phyloseq::otu_table(aggregate_taxa(x,tax)))) %>% rownames_to_column("ID"))
   full_d <- reduce(d,full_join)
   full_d[is.na(full_d)] <- 0
   
   full_ID <- full_d[,1] # extraxct patient IDs
   rownames(full_d) <- full_ID # set rownames to patient IDs
   d <- full_d[,-1] # filter ID variable
   d <- d[,apply(d,2,sum)!=0] # Filter empty taxa
   
   
   filt <- apply(d,1,sum) != 0
   d <- d[filt,] # Filter observations without any OTU count
   d <- d %>% rownames_to_column("ID") %>% mutate(N=rowSums(.[,-1])) 
   
   ## Extract metadata
   met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(meta) %>% rownames_to_column("ID"))
   met1 <- reduce(met,full_join) %>% filter(filt)
   
   out.dat <- reduce(list(met1,d),full_join)
   names(out.dat) <- make.names(names(out.dat)) # make tidy variable names
   
   return(out.dat)
 }

 
prepare_beta_data <- function(dat,meta=c("status","author"),tax="Genus",method="bray") {
   
   ## Combine data to a single dataset
   d <- lapply(dat, function(x) as.data.frame(t(otu_table(aggregate_taxa(x,tax)))) %>% rownames_to_column("ID"))
   full_d <- reduce(d,full_join)
   full_d[is.na(full_d)] <- 0 # Set NA OTU counts to 0, these are taxa which where not detected in other datasets
   
   full_ID <- full_d[,1] # extraxct patient IDs
   rownames(full_d) <- full_ID # set rownames to patient IDs
   d <- full_d[,-1] # filter ID variable
   filt <- apply(d,1,sum) != 0
   d <- d[filt,] # Filter observations without any OTU count
   
   
   ## Extract metadata
   met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(meta) %>% rownames_to_column("ID"))
   met1 <- reduce(met,full_join) %>% filter(filt)
   
   ## Calculate distance matrix
   dist <- vegdist(d,method=method)
   
   
   out <- list("dist"=dist,
               "meta"=met1,
               "data"=d)
   return(out)
   
 }
 
estimate_alpha <- function(dat,tax="Sequence",meta=c("author")) {
  
  dat <- lapply(dat,function(x) aggregate_taxa(x,tax))
  res <- lapply(dat,function(x) {data.frame("richness"=hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=0),
                                            "shannon"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=1),0),
                                            "inv.simpson"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=2)),0)})
  met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(author) %>% rownames_to_column("ID"))
  met1 <- reduce(met,full_join)
  ls <- unlist(lapply(dat,function(x)  rowSums(t(otu_table(x))))) %>% data.frame("N"=.) %>% rownames_to_column("ID")
  
  res <- foreach::foreach(i=1:length(dat)) %do%  {res[[i]] %>% dplyr::mutate(ID=rownames(sample_data(dat[[i]])))}
  
  res <- foreach::foreach(i=1:length(dat)) %do% {dat[[i]] %>% phyloseq::sample_data() %>% plyr::mutate(ID=rownames(.)) %>% data.frame() %>% dplyr::select(ID,status) %>% left_join(res[[i]]) %>% plyr::mutate(study=i)}
  
  out <- res %>% reduce(full_join)
  out <- out %>% left_join(ls,by="ID")
  out <- out %>% left_join(met1,by="ID")
  out <- out %>% mutate(author = fct_reorder(author,study))
  return(out)
  
}

#################################
#### Data cleaning functions ####
#################################
 
 
## Currently, the function only filters based on relative abundance, not taxonomic information, as taxonomic information is present in almost all ASVs
 
 
clean <- function(data,tax,min_tax="Family",seq_length=c(), filter=c("abs","rel","quant"),abs=100,rel=0.01) {
  
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
  
  ## Filter based on low abundance AND/OR low taxonomic information
  
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
      #"Number of ASVs below threshold, but not filtered:",ft_abund_nr,"\n",
      ft_seqlength_nr,'ASVs filtered due to low sequence length < 150 \n',
      'Number of filtered ASVs:',tot_filt,'\n',
      'New number of ASVs:',new_asv,"\n")
  return(d)
}


###################################
#### Data exploration function ####
###################################


## calculate number of ASVs assigned to one Phylum
features_per_phylum <- function(phyl) {table(tax_table(phyl)[, "Phylum"], exclude = NULL)}


## Calculate average abundance of all ASVs in a Phylum and total abundance
prevalence_per_phylum <- function(phyl) {
  ps <- subset_taxa(phyl,!is.na(Phylum))
  # Compute prevalence of each feature, store as data.frame
  prevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
  prevdf <- prevdf %>% arrange(Phylum)
  prev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
  colnames(prev) <- c("Phylum","Avg. Prevalence", "Total Prevalence")
  return(prev)
}


## Plot the fraction of samples in which a taxon appears against its abundance by Phylum
prevalence_by_abundance <- function(phyl,thrs=0.05) {
  
  n <- nsamples(phyl)
  threshold <- thrs*n
  ps <- subset_taxa(phyl,!is.na(Phylum))
  prevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
  
  keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= threshold)]
  
  ps0 <- prune_taxa(keepTaxa,ps)
  prevdf = apply(X = otu_table(ps0),
                 MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))
  
  
  prevdf <- prevdf %>%
    arrange(Phylum) %>%
    mutate(prev.frac = Prevalence/n)
  
  ggplot(aes(x=TotalAbundance,y=prev.frac,color=Phylum),data=prevdf) +
    geom_point() +
    facet_wrap(~Phylum,scales = "free") +
    theme_minimal() +
    scale_x_log10() +
    xlab("Total abundance[log]") +
    ylab("Prevalence (fraction of samples)") +
    theme(legend.position="none")
}



############################ 
#### Analysis functions ####
############################

### Negative binomial random effect model

fitZINBMM <- function(outcome) {
  
  fit <- tryCatch(NBZIMM::glmm.zinb(fixed = as.formula(paste0(outcome,"~status+offset(log(N))")),
                                    random = ~ 1|author,data=full.dat,family = zi.negative.binomial(),niter=100,method="ML"),
                  error=function(e) NULL)
  return(fit)
  
}


sensitivity_ls <- function(data,index="richness",stepsize=1000,length=50) {
  
  n.seq <- seq(from=range(data$N)[1], to=range(data$N)[2],by=stepsize)[1:length]
  
  data.list <- lapply(n.seq, function(x) {data[data$N>x,]} )
  names(data.list) <- n.seq
  
  fitModel <- function(x) {fit <- lme4::glmer(as.formula(paste0(index,"~1+(status|study)")),data=x,family=poisson(link="log"),offset = log(N))
  exp(ranef(fit)$study)}
  
  res <- vector(mode = "list", length = length)
  res <- foreach(i=seq_along(data.list)) %do% fitModel(data.list[[i]])
  res <- lapply(res,function(x) x %>% rownames_to_column(var = "study") )
  foreach(i=seq_along(res)) %do% {res[[i]][,ncol(res[[i]])+1] <- n.seq[[i]]}
  
  
  plot.dat <- res %>% reduce(full_join) %>% mutate(N = V4) %>% dplyr::select(-V4) %>% mutate(study= as.numeric(study))
  
  samp_size <- function(data) {data %>% group_by(study) %>% summarize(s.size = n(),
                                                                      ratio = sum(status=="PD")/sum(status=="HC"))}
  
  ss <- lapply(data.list,samp_size)
  
  ss1 <- foreach(i=seq_along(ss)) %do% {ss[[i]] %>% mutate(N=as.numeric(names(ss)[[i]]))} 
  ss2 <- ss1 %>% reduce(full_join)
  
  plot.dat <- plot.dat %>% left_join(ss2)
  
  ggplot(data=plot.dat,aes(x=N,y=statusPD,col=as.factor(study))) +
    geom_smooth() +
    annotate("text",x=min(plot.dat$N),y=plot.dat$statusPD[plot.dat$N==min(plot.dat$N)]+0.005,label=plot.dat$s.size[plot.dat$N==min(plot.dat$N)]) +
    annotate("text",x=max(plot.dat$N),y=plot.dat$statusPD[plot.dat$N==max(plot.dat$N)]+0.005,label=plot.dat$s.size[plot.dat$N==max(plot.dat$N)]) 
  
}

##################################
#### Summary / plot functions ####
##################################

## Select significant results

is.sig <- function(fit) {
  ind <- summary(fit)$tTable[2,5] < 0.05
  return(ind)
}
 
#### Summary of unviariate analysis results

summarize_univariate <- function(fit) {
  
  fixed <- fit$coefficients$fixed[[2]]
  fixed.p <- summary(fit)$tTable[2,5]
  #fixed.p <- anova(fit)[[2,4]]
  theta <- fit$theta[[1]]
  
  out <- c(fixed,fixed.p,theta)
  names(out) <- c("estimate","p.value","theta")
  return(out)
  
}

##########################
#### Helper functions ####
##########################


addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
