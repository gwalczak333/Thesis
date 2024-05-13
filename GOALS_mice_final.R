### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(BGLR)
library(deepgp)
library(fields)
library(LearnPCA)
library(ggplot2)
library(dplyr)

### Read in the data set ###
data(mice)

### Load in the phenotypes ###
load("~/Documents/BLAST/Thesis/Mice/MousePheno.RData")

### Remove ID/Discrete Phenotypes ###
Pheno = Pheno[-c(5,7)]

### Extract the Phenotype Names ###
pnames = c()
for(i in 1:length(Pheno)){
  pnames = c(pnames,colnames(Pheno[[i]]))
}

### Create a new phenotype matrix for a subset of traits ###
Y = list()
Y[[1]] = Pheno$Glucose[,8]; names(Y[[1]]) = rownames(Pheno$Glucose)
Y[[2]] = Pheno$Obesity[,1]; names(Y[[2]]) = rownames(Pheno$Obesity)
Y[[3]] = Pheno$Immunology[,14]; names(Y[[3]]) = rownames(Pheno$Immunology)
Y[[4]] = Pheno$Haematology[,7]; names(Y[[4]]) = rownames(Pheno$Haematology)
Y[[5]] = Pheno$Biochemistry[,9]; names(Y[[5]]) = rownames(Pheno$Biochemistry)
Y[[6]] = Pheno$Biochemistry[,10]; names(Y[[6]]) = rownames(Pheno$Biochemistry)

names(Y) = pnames[c(42,60,92,102,121,122)]

### Create a y vector and X matrix for the j-th trait ### 
j = 5
X = mice.X

y = as.numeric(Y[[j]]); names(y) = names(Y[[j]])
y = y[which(names(y)%in%rownames(X))]

y = y[!is.na(y)]

X = X[which(rownames(X)%in%names(y)==TRUE),]
X = X[match(names(y),rownames(X)),]


p <- ncol(X)

y=(y-mean(y))/(sd(y))
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)


n <- length(y)


Chr.Positions1 <- read.csv("~/Documents/BLAST/Thesis/Mice/mice_map2.csv", stringsAsFactors = FALSE)
Chr.Positions1 <- Chr.Positions1[,c(1,3,4)]

Chr.Positions2 <- read.table("~/Documents/BLAST/Thesis/Mice/mice_mapfile.txt", header = TRUE)

Chr.Positions <- rbind(Chr.Positions2,
                       Chr.Positions1[which(Chr.Positions1[,1]%in%Chr.Positions2[,1]==FALSE),])

colnames(Chr.Positions) = c("SNP","CHR","BP")
Chr.Positions$SNP = as.character(Chr.Positions$SNP)
Chr.Positions$CHR = as.numeric(Chr.Positions$CHR)
Chr.Positions$CHR[is.na(Chr.Positions$CHR)] = 20

head(Chr.Positions)
head(colnames(X))

snps <- colnames(X)
snps <- unlist(strsplit(snps,"_A"))
snps <- unlist(strsplit(snps,"_C"))
snps <- unlist(strsplit(snps,"_G"))
snps <- unlist(strsplit(snps,"_T"))
colnames(X) <- snps

sum(!snps %in% Chr.Positions$SNP)

gene_1_20 <- Chr.Positions[Chr.Positions[2] == 1 | Chr.Positions[2] == 20]

gene_1_20  <- gene_1_20[gene_1_20 %in% colnames(X)]


pca_mice <- prcomp(X)

test <- PCAtoXhat(pca_mice)
all.equal(test, X, check.attributes = FALSE)

# scree plot
# this shows the amount of variance in the dataset explained by each of the 6 PCs
plot(pca_mice, main = "Scree Plot from PCA on Mice Dataset", xlab = "Principal Components")



# setting nmcmc
nmcmc <- 1000

# number of Principal Components (PCs)
numPC <- 3

v <- matrix(nrow = n, ncol = numPC)
mins <- ranges <- rep(NA, numPC)
for (i in 1:numPC) {
  v[,i] <- X%*%pca_mice$rotation[, i] # Same as pca_mice$x[,i] but we'll use this formula later
  mins[i] <- min(v[,i])
  ranges[i] <- max(v[,i] - min(v[,i]))
  v[,i] <- (v[,i] - mins[i])/ranges[i]
}


# one-layer
fit1 <- fit_one_layer(v, y, nmcmc = nmcmc, verb = TRUE, vecchia = TRUE)
muhat1 <- predict(fit1, v, lite = FALSE)

plot(v[, 1], y)
lines(v[, 1][order(v[, 1])], muhat1$mean[order(v[, 1])])

# two-layer
fit2 <- fit_two_layer(v, y, nmcmc = nmcmc, verb = TRUE, vecchia = TRUE)
muhat2 <- predict(fit2, v, lite = FALSE)

lines(v[, 1][order(v[, 1])], muhat2$mean[order(v[, 1])], col = "blue")


# Constructing MSE table
MSE <- matrix(data = NA, nrow = 10, ncol = 2)
for (PC in 1:10) {
  # setting up v (PCs)
  v <- matrix(nrow = n, ncol = PC)
  for (j in 1:PC) {
    v[,j] <- pca_mice$x[,j] 
    v[,j] <- (v[,j] - min(v[,j]))/max(v[,j] - min(v[,j]))
  }
  
  # one-layer
  fit1 <- fit_one_layer(v, y, nmcmc = nmcmc, verb = TRUE, vecchia = TRUE)
  muhat1 <- predict(fit1, v, lite = FALSE)
  
  MSE[PC, 1] <- round(mean((y - muhat1$mean)^2), 3)
  
  # two-layer
  fit2 <- fit_two_layer(v, y, nmcmc = nmcmc, verb = TRUE, vecchia = TRUE)
  muhat2 <- predict(fit2, v, lite = FALSE)
  
  MSE[PC, 2] <- round(mean((y - muhat2$mean)^2), 3)
}





# one-layer GOALS delta values
delta1 <- matrix(nrow = nrow(v), ncol = length(gene_1_20))

# If loop stops and needs to restart
# load("~/Documents/BLAST/Thesis/Mice/delta1.RData")
# delta_i <- 101
# for (i in gene_1_20[101:lengtH(gene_1_20)]) {
delta_i <- 1
for (i in gene_1_20) {
  cat("i=", i, "\n")
  newX <- X
  newX[, i] <- X[, i] + sort(unique(X[, i]))[3] - sort(unique(X[, i]))[2]
    
  vnew <- matrix(nrow = n, ncol = numPC)
  for (j in 1:ncol(vnew)) {
    vnew[,j] <-  newX%*%pca_mice$rotation[, j]
    vnew[,j] <- (vnew[,j] - mins[j])/ranges[j]
  }
  muhat1new <- predict(fit1, vnew, lite = FALSE)
    
  delta1[,delta_i] <- (muhat1$mean - muhat1new$mean)/(sort(unique(X[, i]))[3] - sort(unique(X[, i]))[2])
  plot(v[, 1], y, main = i)
  lines(v[, 1][order(v[, 1])], muhat1$mean[order(v[, 1])])
  lines(v[, 1][order(v[, 1])], muhat1new$mean[order(v[, 1])],
        col = "red")
  
  delta_i <- delta_i + 1
  
  save(delta1, file = "~/Documents/BLAST/Thesis/Mice/delta1.RData")

}

for (j in 1:length(gene_1_20)) {
  hist(delta1[, j])
}


plot(abs(colMeans(delta1)))

deltameans <- abs(colMeans(delta1))
delta1df <- data.frame(gene_1_20, deltameans)

ggplot(delta1df, aes(x = gene_1_20, y = deltameans)) +
  geom_point()



important_10 <- delta1df[order(delta1df[, 2], decreasing = TRUE), ]
important_10[1:10,]

# Manhattan plot function
manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray60", "gray80"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, novel=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...)
  
# ONE-LAYER MANHATTAN PLOT
load("~/Documents/BLAST/Thesis/Mice/delta1.RData")
names(delta1) = gene_1_20 # subset of only the gene variants I used
v = abs(colMeans(delta1))
names(v) <- gene_1_20
X1 = Chr.Positions[-which(Chr.Positions$SNP%in%names(v)==FALSE),]
v = v[which(names(v)%in%X1$SNP==TRUE)]
X1$"P"[match(names(v),X1$SNP)] = v

manhattan(X1, main = "HDL",cex = 0.5, cex.axis = 0.52, logp = FALSE,ylab = "GOALS One Layer", 
          xlab = "Chromosome",ylim = c(0,max(v, na.rm = TRUE)+0.001),chrlabs = c(1,"X"),
          col = c("blue", "dimgrey"), suggestiveline=FALSE, genomewideline=FALSE, annotateTop = TRUE)



# TWO-LAYER MANHATTAN PLOT
load("~/Documents/BLAST/Thesis/Mice/delta2.RData")
names(delta2) = gene_1_20 # subset of only the gene variants I used
v = abs(colMeans(delta2))
names(v) <- gene_1_20
X2 = Chr.Positions[-which(Chr.Positions$SNP%in%names(v)==FALSE),]
v = v[which(names(v)%in%X2$SNP==TRUE)]
X2$"P"[match(names(v),X2$SNP)] = v

manhattan(X2, main = "HDL",cex = 0.5, cex.axis = 0.52, logp = FALSE,ylab = "GOALS One Layer", 
          xlab = "Chromosome",ylim = c(0,max(v, na.rm = TRUE)+0.001),chrlabs = c(1,"X"),
          col = c("blue", "dimgrey"), suggestiveline=FALSE, genomewideline=FALSE)


one.two.scores.df <- data.frame(X1$P, X2$P)
ggplot(one.two.scores.df, aes(x = X1.P, y = X2.P)) +
  geom_point() +
  geom_abline() +
  xlab("One-Layer GOALS Scores") +
  ylab("Two-Layer GOALS Scores") +
  ggtitle("One-Layer vs. Two-Layer Importance Scores") +
  theme(plot.title = element_text(size = 13, face = "bold"))  +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13,face="bold"))

plot(X1$P, X2$P)
abline(a = 0, b = 1)

# two-layer GOALS
delta2 <- matrix(nrow = nrow(v), ncol = length(gene_1_20))

delta_i2 <- 1
for (i in gene_1_20) {
    cat("i=", i, "\n")
    newX <- X
    newX[, i] <- X[, i] + sort(unique(X[, i]))[3] - sort(unique(X[, i]))[2]
    
    vnew <- matrix(nrow = n, ncol = 1)
    for (j in 1:ncol(vnew)) {
      vnew[,j] <-  newX%*%pca_mice$rotation[, j]
      vnew[,j] <- (vnew[,j] - mins[j])/ranges[j]
    }
    muhat2new <- predict(fit2, vnew[, 1, drop = FALSE], lite = FALSE)
      
    delta2[,delta_i2] <- muhat2$mean - muhat2new$mean
    
    plot(v[, 1], y, main = i)
    lines(v[, 1][order(v[, 1])], muhat2$mean[order(v[, 1])])
    lines(v[, 1][order(v[, 1])], muhat2new$mean[order(v[, 1])],
          col = "red")
    
    delta_i2 <- delta_i2 + 1
    save(delta2, file = "~/Documents/BLAST/Thesis/Mice/delta2.RData")

}


for (j in 1:length(gene_index)) {
  hist(delta2[, j], breaks = seq(-0.003, 0.003, by = 0.0005))
}
plot(abs(colMeans(delta2)))



ggplot(delt2df, aes(x = gene_names, y = deltameans2, colour = deltameans2 >1.659434e-07)) +
  scale_colour_manual(name = 'Important', values = setNames(c('blue','orange'),c(T,F))) +
  geom_point(size = 5) +
  xlab("Gene Varient") +
  ylab("GOALS") +
  ggtitle("Gene Variant-Level Association Mapping for HDL content", subtitle = "*HDL = High-Density Lipoprotein") +
  theme(plot.title = element_text(size = 15, face = "bold"))  +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

delta2df$gene_names <- factor(delta2df$gene_names, levels = delta2df$gene_names)

ggplot(delt2df3) +
  geom_col(aes(x = gene_names, y = deltameans2)) +
  xlab("Gene Varient") +
  ylab("GOALS (importance)") +
  ggtitle("Gene Variant-Level Association Mapping for HDL content") +
  theme(plot.title = element_text(size = 15, face = "bold"))  +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=13,face="bold")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5))
  