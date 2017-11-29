### BAJA CALIFORNIA PASSIFLORA STUDY ###
## Landmarks ##

# Rinse the environment
rm(list=ls())

# Load dependencies
library(Momocs)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggfortify)

# Domestic functions
links_outline <- function(nb){
  c(1, rep(2:nb, each=2), 1) %>% matrix(ncol=2, byrow=TRUE)
}


## Import landmarks ##

ldk <- read.csv("landmarks_import3.csv")

# grab the covariates/cofactors
ldk_fac <- ldk %>% select(id, species, label)

# grab the coordinates
ldk_coo_df <- ldk %>% select(matches("^x|y"))

ldk_coo <- ldk_coo_df %>% nrow() %>% seq_len() %>%
  lapply(function(.) ldk_coo_df[., ] %>%
           as.numeric() %>%
           matrix(ncol=2, byrow = TRUE))

# build the Ldk
pass_ldk <- Ldk(coo=ldk_coo, fac=ldk_fac, links=links_outline(6))

# Procrustes
pass_fg <- fgProcrustes(pass_ldk)


## Ordinations ##

# PCA
pass_fg %>% PCA %>%
  plot(~species, pch = c(17, 16, 15), chull=TRUE, chull.filled=TRUE, ellipses=TRUE, ellipsesax=FALSE, conf.ellipses=.95, eigen=FALSE, labelsgroups=TRUE, box=FALSE, nb.grids=FALSE, morphospace=FALSE, center.origin=FALSE, zoom=1.0, rug=FALSE, cex=0.9, title=" ")
# add "axisvar=FALSE" to get rid of percentages, and "axisnames=FALSE" to get rid of PC names

# Box plot
pass_fg %>% PCA %>%
  boxplot(~species, nax=1:2)

# LDA
l1 <- pass_fg %>% PCA %>% LDA(~species)
l1

l1$CV.tab
l1$CV.correct

plot(l1, pch = c(17, 16, 15), chull=TRUE, chull.filled=TRUE, ellipses=TRUE, ellipsesax=FALSE, conf.ellipses=.95, eigen=FALSE, labelsgroups=TRUE, box=FALSE, nb.grids=FALSE, morphospace=FALSE, center.origin=FALSE, zoom=1.0, rug=FALSE, cex=0.75, title=" ")
plot_CV(l1)
plot_CV2(l1)