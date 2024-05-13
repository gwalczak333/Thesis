library(dplyr)

HDL_X2 <- read.csv("~/Documents/BLAST/Thesis/Mice/X2.csv")
HDL_X2 <- na.omit(HDL_X2)
X2 <- HDL_X2[order(HDL_X2[, 5], decreasing = TRUE), ]
X2_10 <- X2[1:10,c(2,3,5)]
sprintf(X2_10$P, fmt = '%#.5f') 

HDL_X1 <- read.csv("~/Documents/BLAST/Thesis/Mice/X1.csv")
X1 <- X1[X1$SNP %in% X2$SNP,]
X1 <- X1[order(X1[, 5], decreasing = TRUE), ]
X1_10 <- X1[1:10,c(2,3,5)]
sprintf(X1_10$P, fmt = '%#.5f')

HDL_GOALS <- read.csv("~/Documents/BLAST/Thesis/Mice/HDL_GOALS.csv")
GOALS <- GOALS[GOALS$SNP %in% X2$SNP,]
GOALS <- GOALS[order(GOALS[, 5], decreasing = TRUE), ]
GOALS_10 <- GOALS[1:10,c(2,3,5)]

HDL_RATE <- read.csv("~/Documents/BLAST/Thesis/Mice/HDL_RATE.csv")
RATE <- RATE[RATE$SNP %in% X2$SNP,]
RATE <- RATE[order(RATE[, 5], decreasing = TRUE), ]
RATE_10 <- RATE[1:10,c(2,3,5)]

