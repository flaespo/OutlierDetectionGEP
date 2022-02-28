# Ensemble method 
# Author: Laura Selicato Università degli studi di Bari Aldo Moro

# Citation: 
# Selicato, L., Esposito, F., Gargano, G., Vegliante, M. C., Opinto, G., Zaccaria, G. M., ... 
# Del Buono, N. (2021). A new ensemble method for detecting anomalies in gene expression matrices.
# Mathematics, 9(8), 882.

#Packages required
require(fpc);
require(NbClust);
require(cluster);
require(limma);
require(ggraph);
require(igraph);
require(RColorBrewer);
require(ape);
require(tidyverse);
require(factoextra);
require(dendextend);
require(stats);
require(gplots);

A <- as.matrix(read.table("C:/Users/Dataset.txt",header=TRUE,row.names=1,sep="\t"))
Af <- data.frame(A)

# Dataset comes from real tumor gene expression profiling data, with samples on the
# columns and genes on the rows.

################## Preliminary Analysis ########################################

# Plot density on samples:
dev.new()
plotDensities(A,legend=FALSE) 

# Dissimilarity matrix
# Euclidean distance:
d_a <- dist(t(Af),method = "euclidean")
# Hierarchical clustering using Average Hierarchical clustering
hc_a <- hclust(d_a, method = "average" ) 

# Pearson distance:
c_2 = cor(Af,method="pearson")
d_a_p <- as.dist(1-c_2)
# Hierarchical clustering using Average Hierarchical clustering
hc_a_p <- hclust(d_a_p, method = "average" ) 


# Distance evaluation 
# Cophenetic Coefficient:
res.coph <- cophenetic(hc_a)
res.coph_p <- cophenetic(hc_a_p)
# Correlation
cor(d_a,res.coph)
cor(d_a_p,res.coph_p)


# Plot the obtained dendrogram
dev.new()
plot(hc_a, hang=-1) 
plot(hc_a_p, hang=-1)


################# Hierarchical Clustering ######################################
hc.res_a <- eclust(t(Af), "hclust",k = 7, hc_metric = "pearson", 
                   hc_method = "average",  graph = FALSE)  
fviz_silhouette(hc.res_a)


# Optimal number of clusters using gap statistics
hc.res_a$nbclust

# Print result
hc.res_a

# Silhouette information
silinfo_a <- hc.res_a$silinfo
names(silinfo_a)
# Silhouette widths of each observation
head(silinfo_a$widths[, 1:3], 10)
# Average silhouette width of each cluster
silinfo_a$clus.avg.widths
# The total average (mean of all individual silhouette widths)
silinfo_a$avg.width
# The size of each clusters
hc.res_a$size



# Silhouette width of observation
sil_a <- hc.res_a$silinfo$widths[, 1:3]
# Objects with low silhouette
neg_sil_index_a <- which(sil_a[, 'sil_width'] < 0.2)  
sil_a[neg_sil_index_a, , drop = FALSE]

################# Robust PCA ###################################################

library(rospca)
A_t = t(A)
At= as.data.frame(A_t)

#Plot
Dataset_info = data.frame(names = colnames(Af))
outlier_name = c('name_1','name_2','name_3','name_4','name_5','name_6','name_7', 'name_8')
# name_i represents the outputs of Hierarchical clustering

index = which(Dataset_info$names%in%outlier_name)
Dataset_info$Regular = rep("Regular observarion",dim(Dataset_info)[1])
Dataset_info$Regular[index] ="Regular"

outlier1 = c('name_1','name_2')
outlier2 = c('name_3','name_4') 
outlier3 = c('name_5', 'name_6')
outlier4 = c('name_7', 'name_8')

index_1 = which(Dataset_info$names%in%outlier1)
index_2 = which(Dataset_info$names%in%outlier2)
index_3 = which(Dataset_info$names%in%outlier3)
index_4 = which(Dataset_info$names%in%outlier4)

Dataset_info$Regular = rep("Regular observation",dim(Dataset_info)[1])
Dataset_info$Regular[index_1] = "Bad Leverage Point"
Dataset_info$Regular[index_2] = "Orthogonal Point"
Dataset_info$Regular[index_3] = "Mislabeled sample"
Dataset_info$Regular[index_3] = "Alert sample"

group = Dataset_info$outlier_name

colors = rep("#70AD47",length(colnames(Af)))
pos_out_KS = which(Dataset_info$Regular=="Bad Leverage Point")
colors[pos_out_KS] = "#FF0000"

pos_out_SUM = which(Dataset_info$Regular=="Mislabeled sample")
colors[pos_out_SUM] = "#4472C4"

pos_out = which(Dataset_info$Regular=="Orthogonal Point")
colors[pos_out] = "#FA9106"

pos_out_UP = which(Dataset_info$Regular=="Alert sample")
colors[pos_out_UP] = "dark green"

dev.new()
diagPlot(resR0,title = "Robust PCA", col = colors, pch = 16, labelOut = TRUE,id = 8 )#id = 10
legend("topleft",c("Regular Observation","Bad Leverage Point", "Orthogonal Point","Mislabeled sample", "Alert Sample"),lty=1,col=c("#70AD47","#FF0000","#4472C4","#FA9106","dark green"),lwd=2)


