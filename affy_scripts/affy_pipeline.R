####################################################################################
#
#     Install and Load all relevant packages
#
####################################################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
library("affy")

#install.packages("path/to/kohonen_2.0.19.tar.gz", repos = NULL, type="source")
library(kohonen)

#install.packages("data.table")
library(data.table)

#biocLite("pcaGoPromoter")
library("pcaGoPromoter")

#biocLite("pcaGoPromoter.Hs.hg19")
library(pcaGoPromoter.Hs.hg19)
#biocLite("pcaGoPromoter.Mm.mm9")
library(pcaGoPromoter.Mm.mm9)

#biocLite("AnnotationDbi")
library(AnnotationDbi)

#install.packages("path/to/org.Mm.eg.db_3.3.0.tar.gz", repos = NULL, type="source")
library(org.Mm.eg.db)

#biocLite("mouse4302.db")
library(mouse4302.db)

#install.packages("dplyr")
library(dplyr)

#install.packages("ggplot2")
library(ggplot2)

#install.packages("reshape2") 
library(reshape2)

#install.packages("ggdendro") 
library(ggdendro)

#install.packages("Rmisc")
library(Rmisc)

####################################################################################
#
#                   Data preprocessing
#
####################################################################################
####################################################################################
#The following steps are only optional!
#We provide also a data frame which has already been RMA normalized
####################################################################################

# set directory containing all CEL files
dir.cel <- "path/to/GSE10246_RAW"

# load sample info file (should contain filenames + group info: headers "Filename" + "Tissue" +...
sample.info <- read.delim("path/to/sample_info.txt", header = T) 
subset <- read.table("path/to/samples_Fig1-dataset.txt", header = T)
rownames(sample.info) <- sample.info$Filename  

#Load CEL-files
print("reading CEL files")
setwd(dir.cel) 
data <- ReadAffy(phenoData = sample.info)

# RMA normalization
data <- rma(data)
data.pheno <- subset
row.names(data.pheno) <- data.pheno$Filename
data.pheno$Filename <- NULL

# get expression values 
data.ex <- exprs(data) # genes in rows, samples in columns
#######################################################################################
#Only execude the following steps when the steps before have been omitted
#######################################################################################
data.ex <- read.table("path/to/norm_ex_tbl.txt", sep="\t", row.names=1, header=T)
sample.info <- read.delim("path/to/sample_info.txt", header = T) 
subset <- read.table("path/to/samples_Fig1-dataset.txt", header = T)
rownames(sample.info) <- sample.info$Filename
#######################################################################################

data.pheno <- subset[order(subset$Filename),]
row.names(data.pheno) <- data.pheno$Filename 
data.pheno$Filename <- NULL
data.ex <- data.ex[,order(colnames(data.ex))]
data.ex<- data.ex[grep("^AFFX-",row.names(data.ex),invert=T),]

#subsetting the dataframe based on tissue of interest
interm <- t(data.ex)
data.ex <- t(interm[row.names(interm) %in% subset$Filename, ])

#annotation of array ids
keytypes(mouse4302.db)
ids <- row.names(data.ex)
ids.annot <- AnnotationDbi::select(mouse4302.db, ids, "SYMBOL", "PROBEID")
ids.annot <- ids.annot[complete.cases(ids.annot),]
df <- data.ex[row.names(data.ex) %in% unique(ids.annot$PROBEID),]
summary(df)

#defining a background value
df.trsp <- t(df)
df.trsp <- merge(data.pheno, df.trsp, all = TRUE, by = 0)
rownames(df.trsp) <- df.trsp$Row.names
df.trsp$Row.names <- NULL
df.reduced <- aggregate(df.trsp[,-1], list(df.trsp$Tissue), mean, na.rm=T)
df.reduced$Tissue <- NULL
colnames(df.reduced)[1] <- "Tissue"
df.minimal <- as.data.frame(apply(df.reduced[-1], 2, max))
colnames(df.minimal)[1] <- "values"
df.minimal$IDs <- row.names(df.minimal)
df.minimal <- df.minimal[order(df.minimal$values, decreasing = F),]
#defining background as an expression value of 6
df.ids <- as.vector(row.names(df.minimal[df.minimal$values > 6.0,]))
df <- df[row.names(df) %in% df.ids,]


#making gene symbols unique by selecting probes with the highest variance
df.reduced <- as.data.frame(apply(df,1,var))
df.reduced <- as.data.frame(cbind(row.names(df.reduced), df.reduced$`apply(df, 1, var)`))
row.names(df.reduced) <- df.reduced$V1
ids.annot <- AnnotationDbi::select(mouse4302.db, row.names(df.reduced), "SYMBOL", "PROBEID")
#prepare annotation file
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
#ids.annot <- ids.annot[ids.annot$PROBEID %in% row.names(df),]
probes <- unique(ids.annot[duplicated(ids.annot$PROBEID),][1])
probes <- as.character(probes$PROBEID)
probe.genes <- NULL
for (probe in 1:length(probes)) {
  tmp <- ids.annot[grep(probes[probe], ids.annot$PROBEID),][2]
  tmp <- as.character(tmp$SYMBOL)
  tmp <- sort(tmp)
  probe.genes <- c(probe.genes, paste(tmp,collapse = " /// "))
}
ids.annot1 <- ids.annot[ids.annot$PROBEID %in% outersect(ids.annot$PROBEID, probes),]
ids.annot2 <- as.data.frame(cbind(PROBEID = probes, SYMBOL = probe.genes))
ids.annot.new <- as.data.frame(rbind(ids.annot1, ids.annot2))
row.names(ids.annot.new) <- ids.annot.new$PROBEID
ids.annot.tmp <- ids.annot.new
ids.annot.tmp$PROBEID <- NULL
df.reduced <- merge(df.reduced, ids.annot.tmp, all = TRUE, by = 0)
summar_var <- df.reduced %>% group_by(SYMBOL) %>% top_n(1,V2)
df <- df[row.names(df) %in% summar_var$V1,]


# annotate df
ids.annot.new <- as.data.frame(ids.annot.new[row.names(ids.annot.new) %in% row.names(df),])
ids.annot.new$PROBEID <- NULL
df <- merge(ids.annot.new, df, all = TRUE, by = 0)

# save dataframe
# write.table(df, "destination path/dataframe_filtered_unique-genes_annotated.txt", row.names = F, col.names = T, quote = F, sep="\t")

#cleaning up workspace
rm(ids, df.trsp, df.reduced, df.minimal, df.ids, summar_var, probe, tmp, probes, probe.genes, ids.annot, ids.annot1, ids.annot2, ids.annot.new)

##########################################################################################################################################################
#
#                                 SOM clustering
#
##########################################################################################################################################################

workingDir = "destination path";
setwd (workingDir);

# df <- read.delim("path/to/dataframe_filtered_unique-genes_annotated.txt", check.names = FALSE, row.names = 1)
df <- df[,-1]
row.names(df) <- df$SYMBOL
df$SYMBOL <- NULL
data <- t(df)

#Prepare input-file for SOM-clustering
data <- merge(data.pheno, data, all = TRUE, by = 0)
data$Row.names <- NULL
data.annot <- aggregate(data[,-1], list(data$Tissue), mean, na.rm=T)
data.annot$Tissue <- NULL
colnames(data.annot)[1] <- "Tissue"
row.names(data.annot) <- data.annot$Tissue
data.annot$Tissue <- NULL
som.input <- t(data.annot)

#Performing SOM-clustering
set.seed(101)
training <- sample(nrow(som.input2), 1000)
Xtraining <- scale(som.input2[training,])
Xtest <- scale(som.input2,
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))

som.wines <- kohonen::som(Xtraining, grid = somgrid(10, 10))
som.prediction <- predict(som.wines, newdata = Xtest,
                          trainX = Xtraining,
                          trainY = factor(som.input2[training]))


#Identification of SOM-cluster containing Dnajc22
tab = cbind(rownames(som.input),som.prediction$unit.classif)
genename = "Creld1"
cluster = tab[tab[,1] %in% genename,2]

################################################################################################################
#
#                             TF prediction                           
#
##############################################################################################################
# identify overrepresented genes within this cluster
genes_cluster <- tab[tab[,2] %in% cluster,1]
TFtable <- primo(genes_cluster, inputType = "geneSymbol", org = "Mm")
head(TFtable$overRepresented)

################################################################################################################
#
#                           hierarchical clustering
#
################################################################################################################
#typing between \\b*\\b the number of cluster
tab.subset = tab[grep("\\b3\\b", tab[,2]),]
df.subset <- df[row.names(df) %in% tab.subset[,1],]
clusters <- hclust(dist(df.subset))

plot(clusters, labels = F)
abline(h=15, col="red")

#determine number of clusters under the abline
clusterCut <- cutree(clusters, h=15)
hclusters <- as.data.frame(clusterCut)
hclusters$SYMBOL <- row.names(hclusters)
clust <- as.character(hclusters[grep("Dnajc22", row.names(hclusters)),])

#subsetting based on branch of interest
clust.genes <- hclusters[grep("\\b7\\b", hclusters$clusterCut),]
tab.subset <- tab.subset[tab.subset[,1] %in% clust.genes$SYMBOL,]
df.subset <- df[row.names(df) %in% tab.subset[,1],]
clusters <- hclust(dist(df.subset))

plot(clusters)


################################################################################################
#
#                           Expression profile
#
################################################################################################
#Subsetting the data frame based on the GOIs
df.subset <- as.data.frame(som.input[row.names(som.input) %in% c("Dnajc22", "Hnf4a"),])
#df.subset <- as.data.frame(som.input[grep ("\\bDnajc3\\b|\\bDnajc7\\b|\\bDnajc13\\b|\\bDnajc22\\b", row.names(som.input)),])
#df.subset <- as.data.frame(som.input[grep ("\\bNr5a2\\b|\\bFoxa2\\b|\\bFoxa3\\b|\\bHnf4a\\b|\\bNr1h4\\b", row.names(som.input)),])
#df.subset <- as.data.frame(som.input[grep ("\\bDnajc22\\b|\\bCdx2\\b|\\bGata1\\b|\\bHnf1a\\b|\\bHnf1b\\b|\\bNr1h2\\b|\\bHnf4a\\b|\\bNr2f1\\b|\\bPparg\\b", row.names(som.input)),])
df.subset <- as.data.frame(t(df.subset))
df.subset <- cbind(ID=row.names(df.subset), df.subset[,1:length(df.subset)])

#ordering by Dnajc22 expression
df.subset <- df.subset[order(df.subset$Dnajc22, decreasing = T),]
df.subset$ID <- factor(df.subset$ID, levels = df.subset$ID)

#function to determine the number of ticks in the plot
number_ticks <- function(n) {function(limits) pretty(limits, n)}

#Plotting the expression profile
ggplot(df.subset, aes(x = ID, group = 1)) +
  # 
  geom_point(aes(y= Dnajc22), col = "red") + geom_line (aes(y= Dnajc22), col = "red") +
  #
  # geom_point(aes(y= Foxa3), col = "limegreen") + geom_line (aes(y=Foxa3), col = "limegreen") +
  # geom_point(aes(y= Nr5a2), col = "lightcyan4") + geom_line (aes(y=Nr5a2), col = "lightcyan4") +
  # geom_point(aes(y= Nr1h4), col = "orange") + geom_line (aes(y=Nr1h4), col = "orange") +
  # geom_point(aes(y= Foxa2), col = "purple") + geom_line (aes(y=Foxa2), col = "purple") +
  geom_point(aes(y= Hnf4a), col = "blue") + geom_line (aes(y=Hnf4a), col = "blue") +
  
  # geom_point(aes(y= Cdx2), col = "purple") + geom_line (aes(y=Cdx2), col = "purple") +
  # geom_point(aes(y= Gata1), col = "limegreen") + geom_line (aes(y=Gata1), col = "limegreen") +
  # geom_point(aes(y= Hnf1a), col = "orangered2") + geom_line (aes(y=Hnf1a), col = "orangered2") +
  # geom_point(aes(y= Hnf1b), col = "lightcyan4") + geom_line (aes(y=Hnf1b), col = "lightcyan4") +
  # geom_point(aes(y= Nr1h2), col = "darkred") + geom_line (aes(y= Nr1h2), col = "darkred") +
  # geom_point(aes(y= Nr2f1), col = "black") + geom_line (aes(y= Nr2f1), col = "black") +
  # geom_point(aes(y= Pparg), col = "gold") + geom_line (aes(y= Pparg), col = "gold") +
  
  # geom_point(aes(y= Dnajc3), col = "purple") + geom_line (aes(y=Dnajc3), col = "purple") +
  # geom_point(aes(y= Dnajc7), col = "limegreen") + geom_line (aes(y=Dnajc7), col = "limegreen") +
  # geom_point(aes(y= Dnajc13), col = "blue") + geom_line (aes(y=Dnajc13), col = "blue") +
  # geom_point(aes(y= Dnajc22), col = "red") + geom_line (aes(y=Dnajc22), col = "red") +
  # 
  scale_y_continuous(expand = c(0,0), breaks=number_ticks(5)) +
  coord_cartesian(ylim = c(3.5, 12)) +
  geom_hline(aes(yintercept=6.0), color = "grey", linetype="dashed") + 
  ylab ("log2(expression)") +
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
df.subset <- as.data.frame(som.input[grep ("\\bNr1h2\\b|\\bDnajc3\\b|\\bHnf1b\\b|\\bNr1h4\\b|\\bHnf4a\\b|\\bDnajc22\\b|\\bCdx2\\b|\\bNr5a2\\b|\\bFoxa2\\b|\\bFoxa3\\b|\\bHnf1a\\b|\\Dnajc13\\b|\\bDnajc7\\b|\\bGata1\\b|\\bPparg\\b|\\bNr2f1\\b", row.names(som.input)),])
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

