## JSE-2017-12-115
## Leaf Morphometric Analysis

# Set working directory
setwd("~/Desktop")

# Load dependencies
library(Momocs)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggfortify)

# Reformat the data
setwd("~/Desktop")

data <- read.csv("~/Desktop/newLandmarks.csv")

len <- length(data$x)

overall.table <- matrix(nrow=len/6, ncol=(6*2)+2)

overall.length <- len/6

for(j in c(1:overall.length)) {
  
  print(j)
  
  sub.data <- as.matrix(data[ (1+6*(j-1)):((1+6*(j-1))+5), c(1,2,3:4)])
  
  overall.table[j,1] <- sub.data[1,1]
  overall.table[j,2] <- sub.data[2,2]
  overall.table[j,3:4] <- sub.data[1,3:4]
  overall.table[j,5:6] <- sub.data[2,3:4]
  overall.table[j,7:8] <- sub.data[3,3:4]
  overall.table[j,9:10] <- sub.data[4,3:4]
  overall.table[j,11:12] <- sub.data[5,3:4]
  overall.table[j,13:14] <- sub.data[6,3:4]
}

head(overall.table)
tail(overall.table)

colnames(overall.table) <- c("order", "label", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5", "x6", "y6")

write.table(overall.table, "reformatted.txt")

# Domestic functions
links_outline <- function(nb){
  c(1, rep(2:nb, each=2), 1) %>% matrix(ncol=2, byrow=TRUE)
}

# Import landmarks
ldk <- read.csv("BajaPassiflora_landmarks_aridaComplex_3groups.csv")

# Grab the covariates/cofactors
ldk_fac <- ldk %>% select(species, accession)

# Grab the coordinates
ldk_coo_df <- ldk %>% select(matches("^x|y"))

ldk_coo <- ldk_coo_df %>% nrow() %>% seq_len() %>%
  lapply(function(.) ldk_coo_df[., ] %>%
           as.numeric() %>%
           matrix(ncol=2, byrow = TRUE))

# Build the Ldk
pass_ldk <- Ldk(coo=ldk_coo, fac=ldk_fac, links=links_outline(6))

# Procrustes superimposition
pass_fg <- fgProcrustes(pass_ldk)

# PCA ordination
pass_fg %>% PCA %>%
  plot(~species, chull=TRUE, chull.filled=TRUE, eigen=FALSE, labelsgroups=TRUE, box=FALSE, nb.grids=FALSE, morphospace=FALSE, center.origin=FALSE, zoom=1.0, rug=FALSE, cex=0.9, title=" ")

# PCA box plot
pass_fg %>% PCA %>%
  boxplot(~species, nax=1:2)

# Build the classifier
l1 <- pass_fg %>% PCA %>% LDA(~species)
l1$CV.tab
l1$CV.correct

# LDA ordination
plot(l1, chull=TRUE, chull.filled=TRUE, ellipsesax=FALSE, eigen=FALSE, labelsgroups=TRUE, box=FALSE, nb.grids=FALSE, morphospace=FALSE, center.origin=FALSE, zoom=1.0, rug=FALSE, cex=0.75, title=" ")
plot_CV(l1)
plot_CV2(l1)