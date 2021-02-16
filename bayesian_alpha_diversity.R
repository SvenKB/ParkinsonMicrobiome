######################################################################
#### Microbiome Meta-Analysis - Bayesian alpha diversity analyses ####
######################################################################

#### load packages ####
source('packages.R')

#### load functions ####
source("functions.R")


##########################
#### Data preparation ####
##########################

dat <- readRDS("phyl_list.RDS")
dat_FS <-  readRDS("phyl_list_FS.RDS")



alpha_dat <- estimate_alpha(dat) %>% mutate(type="Same region")
alpha_dat_FS <- estimate_alpha(dat_FS) %>% mutate(type="Full Seq.")


author.names <- c("Pietrucci et al., 2019",
                  "Qian et al., 2018",
                  "Keshavarzian et al., 2015",
                  "Aho et al., 2019",
                  "Heintz-Buschart et al., 2018",
                  "Weis et al., 2019",
                  "Barichella et al., 2019")

names <- cbind(study=paste0("Study",1:7),author.names)


###########################
#### Inspect diversity ####
###########################

#### Same region ####

ggplot(aes(x=as.factor(author),y=(richness)),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Richness")

ggplot(aes(x=as.factor(author),y=(shannon)),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Shannon")

ggplot(aes(x=as.factor(author),y=inv.simpson),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Inv. Simpson")


#### Full sequence ####

ggplot(aes(x=as.factor(author),y=(richness)),data=alpha_dat_FS) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Richness")

ggplot(aes(x=as.factor(author),y=(shannon)),data=alpha_dat_FS) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Shannon")

ggplot(aes(x=as.factor(author),y=inv.simpson),data=alpha_dat_FS) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Inv. Simpson")



#### Combined plots ####

alpha_df <- bind_rows(alpha_dat,alpha_dat_FS)

plot_alpha <- alpha_df %>%
  pivot_longer(cols=c("richness","shannon","inv.simpson"),names_to="measure") %>%
  mutate(measure = factor(measure,levels=c("richness","shannon","inv.simpson")))


ggplot(aes(x=measure,y=value,fill=type),data=plot_alpha) +
  geom_boxplot() +
  geom_jitter(aes(group=type,color=type),position = position_jitterdodge(),alpha=.5) +
  #coord_flip() +
  facet_wrap(~author,nrow=1) +
  theme_pubclean() +
  ylab("Effective number of taxa")


#### Correlation of measures with seq depth ####
alpha_dat %>%
  group_by(author) %>%
  dplyr::summarize("Richness" = round(cor(richness,N),2),
            "Shannon" = round(cor(shannon,N),2),
            "Inverse simpson" = round(cor(inv.simpson,N),2)) %>%
  data.frame(row.names=NULL) %>%
  column_to_rownames("author") %>%
  t %>% kable(caption = "Correlation of alpha diversity indices with library size - Same region sequences") %>%
  kableExtra::kable_styling("striped")

alpha_dat_FS %>%
  group_by(author) %>%
  dplyr::summarize("Richness" = round(cor(richness,N),2),
                   "Shannon" = round(cor(shannon,N),2),
                   "Inverse simpson" = round(cor(inv.simpson,N),2)) %>%
  data.frame(row.names=NULL) %>%
  column_to_rownames("author") %>%
  t %>% kable(caption = "Correlation of alpha diversity indices with library size - Full sequences") %>%
  kableExtra::kable_styling("striped")



#### Inspect relationships between library size and diversity
plot_alpha

ggplot(aes(x=N,y=value),data=plot_alpha) +
  geom_point(alpha=.5) +
  geom_smooth(colour="red") +
  facet_grid(author~measure,scales="free") +
  theme_hc()


cor(alpha_dat$richness,log(alpha_dat$N))


######################################
#### Bayesian random effect model ####
######################################

### Legend of fit names 
# RF = Random effect model
# SM = Simple model
# SR = Same region sequence
# FS = Full sequences
# SH = Shannon index
# IS = Inverse Simpson
# R = Richness

#### Same region sequences ####


disper_fit <- glm(inv.simpson~status+logN,data=alpha_dat,family = poisson())

AER::dispersiontest(disper_fit)

alpha_dat %>%
  mutate(author = case_when(author==names[[1,2]] ~ names[[1,1]],
                            author==names[[2,2]] ~ names[[2,1]],
                            author==names[[3,2]] ~ names[[3,1]],
                            author==names[[4,2]] ~ names[[4,1]],
                            author==names[[5,2]] ~ names[[5,1]],
                            author==names[[6,2]] ~ names[[6,1]],
                            author==names[[7,2]] ~ names[[7,1]]))

fit_SR__RF_SH <- brms::brm(richness~1+status+(1+status|author),data=alpha_dat,family=poisson())




#### Full sequences ####











