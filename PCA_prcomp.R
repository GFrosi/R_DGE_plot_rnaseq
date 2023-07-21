#PCA script 
#Acknowledgment: Fred White

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Installing tidyverse
install.packages("tidyverse")
install.packages("dplyr") 
install.packages('GGally')
install.packages('ggpubr')
install.packages('ggforce')



library(tidyverse)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggrepel)
library(ggforce)


#setting the workdir
setwd("/usr/rnaseq-2023/CT26_4BR/")

#Loading file
count_data = read.delim('dge_results.csv', sep='\t')
# count_data = read.csv('dge_results.csv') #HT29
# count_data_filter <- count_data[,c(1,9,10,11,12,13,14)]

##Count data filter normalized counts for 4BR
count_data_filter <- count_data[,c(1,9,10,11,12)] 

#replacing .raw in the header to get the same design names
colnames(count_data_filter) = gsub(".norm", "", colnames(count_data_filter))
count_data_filter <- count_data_filter %>% column_to_rownames(., var = 'id')
count_data_filter

#loading design
sample_info <- read.csv("R_PCA/design.csv")
sample_info <- sample_info %>% column_to_rownames(., var='Sample')
sample_info

#checking if the names in count and design are the same
all(colnames(count_data_filter) %in% rownames(sample_info))

#checking if the names are in the same orther
all(colnames(count_data_filter) == rownames(sample_info))


#new script
# count_data - rows are genes and columns are samples. Transpose by t() to create PCA object
transposed_df = t(count_data_filter)

#Create PCA object
PCA_obj = prcomp(transposed_df)
rownames(PCA_obj$x)
print(PCA_obj)


#calculate PCs values
PCs = data.frame(PCA_obj$x, sample = rownames(PCA_obj$x))
#Add condition to plot the legends (e.g treated and untreated)
PCs$treatment <- sample_info$Treatment 
PCs


## extract the proportion of variance explained by each PC. Transforming values into percent  
var = summary(PCA_obj)$importance[2,]
var

Var = data.frame("x"=1:length(var), "var" = as.vector(var)*100)
Var


# You can also extract PC loadings (gene contribution to each PC): rows are genes and columns are PCs 
loadings = PCA_obj$rotation
#saving the loadings
write.csv(as.data.frame(loadings), file="CT26_4br_gene_contribution.csv" )


#plot pca using ggplot2 + ggpubr
#scree plot
scree_plot = ggplot(Var, aes(x=x, y=var)) + geom_bar(stat = "identity", fill = "blue") + ggtitle("Scree plot CT26") + xlab("PCs") + ylab("Variance (%)")
scree_plot

#PCA
#PC1_PC2
pc1_pc2 = ggplot(PCs,aes(x=PC1, y=PC2, color=treatment, label = sample)) + geom_point(alpha=2, size=3) + ggtitle("PCA CT26 4BR PC1 vs PC2") + xlab(paste0("PC1 (",Var$var[1],")")) + ylab(paste0("PC2 (",Var$var[2],")")) + geom_label_repel(show.legend = FALSE) 
pc1_pc2
ggsave("PCA_CT26_4BR_pc1_pc2_ggplot.png",pc1_pc2,dpi=300)


#PC1_PC3
pc1_pc3 = ggplot(PCs,aes(x=PC1, y=PC3, color=treatment, label = sample)) + geom_point(alpha=2, size=3) + ggtitle("PCA CT26 4BR PC1 vs PC3") + xlab(paste0("PC1 (",Var$var[1],")")) + ylab(paste0("PC3 (",Var$var[3],")"))  + geom_label_repel(show.legend = FALSE)
pc1_pc3
ggsave("PCA_CT26_4BR_pc1_pc3_ggplot.png",pc1_pc3,dpi=300)


#Pairs plot ggforce
pairs = ggplot(PCs, aes(x = .panel_x, y = .panel_y, color=treatment)) + geom_point() + facet_matrix(vars(PC1, PC2, PC3, PC4)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle("PCA CT26 4BR all pairs")
pairs
ggsave("PCA_CT26_4BR_pairs_ggplot.png", pairs ,dpi=300)

