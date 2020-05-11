#### load packages ####
source('packages.R')

#### load functions ####
source("functions.R")


#### Load data ####

# Load ASV data
A <- readr::read_delim("Data/A.txt",delim="\t")
B <- readr::read_delim("Data/B.txt",delim="\t")
C <- readr::read_delim("Data/C.txt",delim="\t")
E <- readr::read_delim("Data/E.txt",delim="\t")
H <- readr::read_delim("Data/H.txt",delim="\t")

# Load Meta data
A_Meta <- loadENAMeta(path="Data/A_meta.txt")
B_Meta <- loadENAMeta(path="Data/B_meta.txt")
#C_Meta <- loadENAMeta(path="Data/C_meta.txt")
#D_Meta <- loadENAMeta(path="Data/D_meta.txt")
E_Meta <- loadENAMeta(path="Data/E_meta.txt")
F_meta <- loadENAMeta(path="Data/F_meta.txt")
G_Meta <- loadENAMeta(path="Data/G_meta.txt")
H_Meta <- loadENAMeta(path="Data/H_meta.txt")
#I_Meta <- loadENAMeta(path="Data/I_meta.txt")


#### Tidy up the Data ####

# prepare taxonomies 
A_tax <- prepare_taxonomy(A)
B_tax <- prepare_taxonomy(B)
C_tax <- prepare_taxonomy(C)
E_tax <- prepare_taxonomy(E)
H_tax <- prepare_taxonomy(H)

# transpose into tidy format
A_dat <- tidy_up(A,trim=F)
B_dat <- tidy_up(B,trim=T)
C_dat <- tidy_up(C)
E_dat <- tidy_up(E,trim=T)
H_dat <- tidy_up(H)


#### Data cleaning ####

## Inspect data
A_dat %>% inspect_sparsity(sp='asv') %>% summary
B_dat %>% inspect_sparsity(sp='asv') %>% summary
C_dat %>% inspect_sparsity(sp='asv') %>% summary
D_dat %>% inspect_sparsity(sp='asv') %>% summary
E_dat %>% inspect_sparsity(sp='asv') %>% summary
F_dat %>% inspect_sparsity(sp='asv') %>% summary
G_dat %>% inspect_sparsity(sp='asv') %>% summary
H_dat %>% inspect_sparsity(sp='asv') %>% summary
I_dat %>% inspect_sparsity(sp='asv') %>% summary

A_dat %>% inspect_reads() %>% summary
B_dat %>% inspect_reads() %>% summary
C_dat %>% inspect_reads() %>% summary
D_dat %>% inspect_reads() %>% summary
E_dat %>% inspect_reads() %>% summary
F_dat %>% inspect_reads() %>% summary
G_dat %>% inspect_reads() %>% summary
H_dat %>% inspect_reads() %>% summary
I_dat %>% inspect_reads() %>% summary

A_dat %>% inspect_depth() %>% hist
B_dat %>% inspect_depth() %>% hist
C_dat %>% inspect_depth() %>% hist
D_dat %>% inspect_depth() %>% hist
E_dat %>% inspect_depth() %>% hist
F_dat %>% inspect_depth() %>% hist
G_dat %>% inspect_depth() %>% hist
H_dat %>% inspect_depth() %>% hist
I_dat %>% inspect_depth() %>% hist


## Quality filtering
A_clean <- clean(A_dat,A_tax,min_tax = "Class",min_abund=100)
B_clean <- clean(B_dat,B_tax,min_tax = "Class",min_abund=100)
C_clean <- clean(C_dat,C_tax,min_tax = "Class",min_abund=100)
D_clean <- clean(D_dat,D_tax,min_tax = "Class",min_abund=100)
E_clean <- clean(E_dat,E_tax,min_tax = "Class",min_abund=100)
F_clean <- clean(F_dat,F_tax,min_tax = "Class",min_abund=100)
G_clean <- clean(G_dat,G_tax,min_tax = "Class",min_abund=100)
H_clean <- clean(H_dat,H_tax,min_tax = "Class",min_abund=100)
I_clean <- clean(I_dat,I_tax,min_tax = "Class",min_abund=100)


#### Run analyses ####




#### Prepare output ####
