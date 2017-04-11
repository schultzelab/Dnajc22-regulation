####################################################################################
#
#     Install and Load all relevant packages
#
####################################################################################
#install.packages("stringr")
library(stringr)

#install.packages("reshape2")
library(reshape2)

#install.packages("data.table")
library(data.table)

#install.packages("ggplot2")
library(ggplot2)

#source("https://bioconductor.org/biocLite.R")
#biocLite("kohonen")
library(kohonen)

#biocLite("DESeq2")
library('DESeq2')

# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# biocLite("pcaGoPromoter")
library(pcaGoPromoter)

# biocLite("pcaGoPromoter.Hs.hg19")
library(pcaGoPromoter.Hs.hg19)

#install.packages("dplyr")
library(dplyr)

#install.packages("ggdendro") 
library(ggdendro)

#install.packages("Rmisc")
library(Rmisc)

##########################################################################################
#
#             Loading data
#
##########################################################################################
# set directory containing expression tables of respective files
setwd("path/to/expression tables")

#Merging of expression tables
files <- dir(pattern = "counted.txt")
df<- NULL
df <- read.table(files[1], sep = "\t", col.names = c("SYMBOL", str_sub(files[1],1,as.numeric((nchar(files[1]) -4)))))
df <- data.table(df, key = "SYMBOL")
for (i in 2:length(files)){
  tmp <- read.table(files[i], sep = "\t", col.names = c("SYMBOL", str_sub(files[i],1,as.numeric((nchar(files[i]) -4)))))
  tmp <- data.table(tmp, key = "SYMBOL")
  df <- merge(df, tmp)
}

### Remove rows containing no gene information #########################################################
samples <- c("__alignment_not_unique", "__ambiguous", "__no_feature", "not_aligned", "__too_low_aQual")
df <- df[grep(paste(samples, collapse = "|"), df$SYMBOL, value = T, invert = T),]

### Optional: cleaning up the colnames ###########################################################################
colnames(df) <- gsub("_counted", "", colnames (df))

#########################################################################################################
#
#
#                         DESeq-normalization
#
#
#########################################################################################################
df <- as.data.frame(df)
#read in table as factors
countData = df
rownames(countData) <- as.character(countData$SYMBOL)
countData$SYMBOL <- NULL
summary(countData)

##visualization of the data in a boxplot
melted.df <- replace(countData, countData < 1, 1)
melted.df <- as.matrix(melted.df)
melted <- reshape2::melt(melted.df, id = colnames(melted.df))
ggplot(melted, aes (x=Var2, y=log(value,10))) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=90, size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, size=1))


#numeric values to nearest integers, because the input must be a matrix of non-negative integers
countData_rounded <- round(countData) 

#sample annotation as an ordered vector
cond <- str_sub(colnames(df)[2:length(df)], 1, as.numeric(nchar(colnames(df)[2:length(df)])-1))
colData = data.frame(condition = factor(cond))
sample_annotation <- as.vector(str_sub(colnames(df)[2:length(df)], 1, as.numeric(nchar(colnames(df)[2:length(df)])-1)))

#make DESeq dataset and run DESeq main function to estimate size factors and dispersions
dds <- DESeqDataSetFromMatrix(countData_rounded, colData, formula(~ condition))

if(length(unique(sample_annotation)) >= 3){
  dds_de<-DESeq(dds, betaPrior=FALSE)
} else {
  dds_de<-DESeq(dds)
}

#plot PCA
vsd <- varianceStabilizingTransformation(dds_de)
x11()
plotPCA(vsd)

#Dealing with count outliers (optional)
ddsClean <- replaceOutliersWithTrimmedMean(dds_de)

if(length(unique(sample_annotation)) >= 3){
  ddsClean<-DESeq(ddsClean, betaPrior=FALSE)
} else {
  ddsClean<-DESeq(ddsClean)
}

tab <- table(initial = results(dds_de)$padj < .1,
             cleaned = results(ddsClean)$padj < .1)
addmargins(tab)

#plot PCA
vsd <- varianceStabilizingTransformation(ddsClean)
x11()
plotPCA(vsd)
dev.off()

#create normalized data counts by dividing sizeFactors for each sample
norm<-sweep(counts(ddsClean),MARGIN=2,colData(ddsClean)$sizeFactor,`/`)
colnames(norm) <- colnames(countData)

#Save normalized table
write.table(norm,"DESeq-norm.txt",sep="\t",quote=F, row.names=T,col.names=NA)

##############################################################################################
#
#           Further preprocessing of data
#
##############################################################################################
##flooring
norm <- replace(norm, norm < 1, 1)

##visualization of the data in a boxplot after normalization
melted <- melt(norm, id = colname(norm) )
ggplot(melted, aes (x=Var2, y=log(value,10))) +
geom_boxplot() +
theme(axis.text.x=element_text(angle=90, size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, size=1))

#rename the column names
colnames(norm) <- sample_annotation

##filtering
test <- t(norm)
df.reduced <- aggregate(test, list(row.names(test)), mean, na.rm=T)
df.minimal <- as.data.frame(apply(df.reduced[-1], 2, max))
colnames(df.minimal)[1] <- "values"
df.minimal$IDs <- row.names(df.minimal)
df.minimal <- df.minimal[order(df.minimal$values, decreasing = F),]
#get summary statistics of data table
summary(df.minimal)
quantile(df.minimal$values, prob = 0.33)
#defining background as the values under the 33th percentile ( here ~ 200)
df.ids <- as.vector(row.names(df.minimal[df.minimal$values >= 200,]))
df <- norm[row.names(norm) %in% df.ids,]

#cleaning up workspace
rm(colData, countData, countData_rounded, norm, test, cond, dds, dds_de, ddsClean, files, i, sample_annotation, vsd, samples)

###################################################################################################
#
#
#                   SOM-clustering
#
#
###################################################################################################
#Prepare input-file for SOM-clustering
data <- t(df)
data.annot <- aggregate(data, list(row.names(data)), mean, na.rm=T)
data.annot <- as.data.frame(data.annot)
rownames(data.annot) <- data.annot$Group.1
data.annot$Group.1 <- NULL
som.input <- t(data.annot)
som.input <- log(som.input,10)

#Performing SOM-clustering
set.seed(101)
training <- sample(nrow(som.input), 1000)
Xtraining <- scale(som.input[training,])
Xtest <- scale(som.input,
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))

som.wines <- kohonen::som(Xtraining, grid = somgrid(10, 10))
som.prediction <- predict(som.wines, newdata = Xtest,
                          trainX = Xtraining,
                          trainY = factor(som.input[training]))

#Identification of SOM-cluster containing DNAJC22
tab = cbind(rownames(som.input),som.prediction$unit.classif)
genename = "DNAJC22"
cluster = tab[tab[,1] %in% genename,2]

################################################################################################################
#
#                             TF prediction                           
#
##############################################################################################################
# identify potential TFs within this cluster
genes_cluster <- tab[tab[,2] %in% cluster,1]
TFtable <- primo(genes_cluster, inputType = "geneSymbol", org = "Hs")
head(TFtable$overRepresented)

################################################################################################################
#
#                           hierarchical clustering
#
################################################################################################################
#typing between \\b*\\b the number of cluster
tab.subset = tab[grep("\\b50\\b", tab[,2]),]
df.subset <- df[row.names(df) %in% tab.subset[,1],]
df.subset <- log(df.subset,10)
clusters <- hclust(dist((1-df.subset)/2))
plot(clusters)

################################################################################################
#
#                           Expression profile
#
################################################################################################
#Subsetting the data frame based on the GOIs
df.subset <- as.data.frame(som.input[row.names(som.input) %in% c("DNAJC22", "HNF4A", "HNF1B"),])
# df.subset <- as.data.frame(som.input[grep ("\\bDNAJC3\\b|\\bDNAJC7\\b|\\bDNAJC13\\b|\\bDNAJC22\\b", row.names(som.input)),])
# df.subset <- as.data.frame(som.input[grep ("\\bNR5A2\\b|\\bFOXA3\\b|\\bGRHL2\\b|\\bHNF4G\\b|\\bHNF4A\\b|\\bHNF1B\\b", row.names(som.input)),])
# df.subset <- as.data.frame(som.input[grep ("\\bEP300\\b|\\bZEB1\\b|\\bPBX1\\b|\\bHNF1A\\b", row.names(som.input)),])
df.subset <- as.data.frame(t(df.subset))
df.subset <- cbind(ID=row.names(df.subset), df.subset[,1:length(df.subset)])

#ordering by DNAJC22 expression
df.subset <- df.subset[order(df.subset$DNAJC22, decreasing = T),]
df.subset$ID <- factor(df.subset$ID, levels = df.subset$ID)

#function to determine the number of ticks in the plot
number_ticks <- function(n) {function(limits) pretty(limits, n)}

#Plotting the expression profile
ggplot(df.subset, aes(x = ID, group = 1)) +
geom_point(aes(y= DNAJC22), col = "red", size = 3) + geom_line (aes(y= DNAJC22), col = "red") +
geom_point(aes(y= HNF4A), col = "blue", size = 3) + geom_line (aes(y=HNF4A), col = "blue") +
# geom_point(aes(y= HNF1B), col = "green", size = 3) + geom_line (aes(y=HNF1B), col = "green") +
#
# geom_point(aes(y= NR5A2), col = "limegreen") + geom_line (aes(y=NR5A2), col = "limegreen") +
# geom_point(aes(y= FOXA3), col = "lightcyan4") + geom_line (aes(y=FOXA3), col = "lightcyan4") +
# geom_point(aes(y= GRHL2), col = "orange") + geom_line (aes(y=GRHL2), col = "orange") +
# geom_point(aes(y= HNF4G), col = "purple") + geom_line (aes(y=HNF4G), col = "purple") +
# geom_point(aes(y= HNF4A), col = "red") + geom_line (aes(y=HNF4A), col = "red") +
# geom_point(aes(y= HNF1B), col = "blue") + geom_line (aes(y=HNF1B), col = "blue") +
#
# geom_point(aes(y= EP300), col = "purple") + geom_line (aes(y=EP300), col = "purple") +
# geom_point(aes(y= HNF1A), col = "limegreen") + geom_line (aes(y=HNF1A), col = "limegreen") +
# geom_point(aes(y= PBX1), col = "blue") + geom_line (aes(y=PBX1), col = "blue") +
# geom_point(aes(y= ZEB1), col = "orange") + geom_line (aes(y=ZEB1), col = "orange") +
# 
# geom_point(aes(y= DNAJC3), col = "purple") + geom_line (aes(y=DNAJC3), col = "purple") +
# geom_point(aes(y= DNAJC7), col = "limegreen") + geom_line (aes(y=DNAJC7), col = "limegreen") +
# geom_point(aes(y= DNAJC13), col = "blue") + geom_line (aes(y=DNAJC13), col = "blue") +
# geom_point(aes(y= DNAJC22), col = "red") + geom_line (aes(y=DNAJC22), col = "red") +
# 
  scale_y_continuous(expand = c(0,0), breaks=number_ticks(5)) +
  coord_cartesian(ylim = c(0, 5)) +
  geom_hline(aes(yintercept=log(200,10)), color = "grey", linetype="dashed") + 
  ylab ("log10(expression)") +
  ggtitle("expression profile") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, size=1))


#################################################################################################
#
#                     Pearson's Correlation Coefficient Matrix
#
#################################################################################################

#prepare data; samples in columns + genes in rows
df.subset <- as.data.frame(som.input[grep ("\\bZEB1\\b|\\bDNAJC13\\b|\\bEP300\\b|\\bDNAJC7\\b|\\bPBX1\\b|\\bGRHL2\\b|\\bHNF1B\\b|\\bDNAJC3\\b|\\bHNF4G\\b|\\bDNAJC22\\b|\\bHNF4A\\b|\\FOXA3\\b|\\bNR5A2\\b|\\bNFIC\\b|\\bFOXA1\\b|\\bHDAC2\\b|\\bJAZF1\\b|\\bMBD4\\b", row.names(som.input)),])
df.subset <- df.subset[grep("-AS1", row.names(df.subset), invert = T),]
df.new <- as.data.frame(t(df.subset))

#compute the correlation matrix
cormat <- round(cor(df.new),2)

#create the correlation heatmap with ggplot2
melted_cormat <- melt(cormat)
ggplot(melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Get the upper triangles of the cormat
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)

#Finished correlation matrix heatmap
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  coord_fixed()

#Reorder the correlation matrix
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1))+
  coord_fixed()

#ggdendrogram(hc, rotate = FALSE, size = 2)

#Add correlation coefficients on the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

