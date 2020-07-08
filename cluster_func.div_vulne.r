#########################################################################################
## This script is used to:
## 1) Cluster tree species into functional groups based on their functional traits 
##    using an aglomerative clustering algorithm
## 2) Compute tree communities functional diverisity as the exponent of the Shannon 
##    diversity index applied to the relative abundance of functional groups
## 3) Compute vulnerability to natural disturbances at the tree community level
#########################################################################################

## Load libraries
options(warn=-1)
suppressPackageStartupMessages({
library(cluster)
library(FD)
library(scales)
library(tidyverse)
library(vegan)  
})



########################################## CLUSTERING ###############################################

## Read functional trait values for a list of species
func.trait.spp <- read.table("https://raw.githubusercontent.com/nuaquilue/Simply-to-use_Functional_Network/master/functional.trait.txt", header=T, sep= "\t" )
rownames(func.trait.spp) <- func.trait.spp$code

## Create a list of data frames, one per each type of traits
## Quantitative traits
quant.trait <- select(func.trait.spp, seed.mass:waterlog.tol)
# Dichotomous trait: 0/FALSE/no or 1/TRUE/yes
dicho.trait <- select(func.trait.spp, ass.mycorrhiza)
# Nominal trait:	Angiosperme / Gymnosperme
nominal.trait <- select(func.trait.spp, phylogen)
# Binary traits:  Mode of dispersion
aux <- 1*data.frame(func.trait.spp$dispersion=="A", func.trait.spp$dispersion=="W",
                    func.trait.spp$dispersion=="H", func.trait.spp$dispersion=="U")
binary.trait <- prep.binary(aux, col.blocks=4, label="dispersion")
# List of dataframes
ktab <- ktab.list.df(list(quant.trait, nominal.trait, dicho.trait, binary.trait))


## Compute the species dissimilarity matrix based on the generalization of the Gower's distance.
distrait <- dist.ktab(ktab, c("Q", "N", "D", "B"), c("scaledBYrange")) 


## Conduct agglomerative clustering on the Gower's dissimilarity matrix
clust <- agnes(distrait, diss=TRUE, method="ward")


## Plot a tree dendogram with 'nclust' number of clusters
nclust <- 5
plot(clust, cex=0.5, which.plots=2, main="", xlab="")
rect.hclust(clust, k=nclust, border=brewer_pal(palette = "Set1")(nclust))



##################################### FUNCTIONAL DIVERSITY ##########################################

## Read list of species and relative abundance of species within tree communites
species <- read.table("https://raw.githubusercontent.com/nuaquilue/Simply-to-use_Functional_Network/master/species.txt", header=T, sep= "\t" )
rel.abund <- read.table("https://raw.githubusercontent.com/nuaquilue/Simply-to-use_Functional_Network/master/spp.rel.abund.stand.txt", header=T, sep= "\t" )


## Assign a functional group to every species
func.trait.spp$fg <- cutree(clust, k=nclust)
species.fg <- left_join(species, select(func.trait.spp, code, fg), "code")


## Per tree community, sum the relative abundance of each functional group 
## and compute functional diversity as the exponent of the Shannon index
func.diver <- gather(rel.abund, code, abund, ABBA:TSCA, factor_key=TRUE) %>% 
              left_join(select(species.fg, code, fg), "code") %>% 
              group_by(id, fg) %>% summarise(x=sum(abund)) %>% 
              group_by(id) %>% summarise(fd.exp=exp(diversity(x, "shannon"))) %>%
              mutate(fd=(fd.exp-1)/(nclust-1))



######################################## VULNERABILITY ##############################################

## Read list of scores (-3 negatively affected to 3 positively affected) per species and disturbances,
## and uncertainty and future relevance per disturbance 
scores <- read.table("https://raw.githubusercontent.com/nuaquilue/Simply-to-use_Functional_Network/master/scores.txt", header=T, sep= "\t" )
disturb <- read.table("https://raw.githubusercontent.com/nuaquilue/Simply-to-use_Functional_Network/master/disturbances.txt", header=T, sep= "\t" )


## Compute species vulnerability as the mean weighted scores by disturbance uncertainty and future relevance
## Vulnerability of a species goes from 0-less to 10-most vulnerable to disturbances
spp.vuln <-  gather(scores, disturbance, x, Browsing:Wind, factor_key=TRUE) %>% 
             left_join(disturb, "disturbance") %>% mutate(wx=x*future.relev*uncertainty) %>%
             group_by(code) %>% summarise(v=mean(wx)) %>% mutate(vindx=rescale(v, to=c(0, -10))+10) %>% select(-v)


## Compute tree communities vulnerability
vulner <- gather(rel.abund, code, abund, ABBA:TSCA, factor_key=TRUE) %>%
          left_join(spp.vuln, by="code") %>% mutate(wvuln=abund*vindx) %>%
          group_by(id) %>% summarise(wvuln=sum(wvuln, na.rm = T), abund=sum(abund)) %>%
          mutate(vindx=wvuln/abund) %>%  select(-wvuln, -abund)
