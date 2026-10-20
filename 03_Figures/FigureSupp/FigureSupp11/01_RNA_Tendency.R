#########################################################################
#Analyzes expression trajectories of CoM+, CoM−, and random gene sets across age 
#in postmortem brain RNA-seq datasets (GSE36192, GSE33000, GSE15222). For each gene set, expression 
#values are centered relative to a young reference age (35 years) and visualized per-sample using LOESS smoothing.
#########################################################################
library(dplyr)  
library(ggplot2)
library(ggpubr)
library(readr)
library(data.table)
library(sva)
library(factoextra)



data<-fread("/home/epi-lab/Escritorio/RNA_GSE/3_datasets_updated.txt", header=TRUE)
data<-as.data.frame(data)
rownames(data)<-data$V1
data<-data[,-1]

data<-as.data.frame(t(data))

data$Gender<-gsub("female", "F", data$Gender)
data$Gender<-gsub("Female", "F", data$Gender)
data$Gender<-gsub("male", "M", data$Gender)
data$Gender<-gsub("Male", "M", data$Gender)

table(data$Array, data$Gender)                           #To check the sex ratios of the studies

hist(as.numeric(data$Age))
max(as.numeric(data$Age))


#We split data and data-info
info<-as.data.frame(data$Age)
row.names(info)<-rownames(data)
colnames(info)<-"Age"
info$Gender<-data$Gender
info$Region<-data$Region
info$Array<-data$Array

info$Age<-as.numeric(data$Age)
hist(info$Age)

data<-subset(data, select = -c(Age,Gender,Region,Array))                               #We remove non-numeric columns (kept in info)
sapply(data, class)                                                                    #To check the class -> all them factors
data<- mutate_all(data, function(x) as.numeric(as.character(x)))                          #To convert to numeric objects
sapply(data, class)  






table(info$Region, info$Array)

aggregate(info$Age ~ info$Array,FUN = quantile,
          probs = c(0,0.25,0.50, 0.75,1))


hist(info$Age) 
table(info$Gender, info$Array)



#To load CoM+ and CoM- genes

genes_pos <- read_delim("Escritorio/Resultados/CHH/genes/genes_COM_pos_aggregated.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
genes_neg<- read_delim("Escritorio/Resultados/CHH/genes/genes_COM_neg_aggregated.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
min(info$Age)
max(info$Age)
x<-data
x$Age<-info$Age






#######################################################################################
#####  COM +  (direct on samples)
#######################################################################################

control_genes_pos<-subset(data, select= intersect(colnames(data),genes_pos$gene_symbol))
control_genes_pos$Mean<- apply(control_genes_pos, 1, median, na.rm=T)
control_genes_pos$SD<-apply(control_genes_pos[,1:(ncol(control_genes_pos)-1)], 1, sd, na.rm = TRUE)
control_genes_pos$Age<-as.numeric(info$Age)


set.zero.young2<-aggregate(. ~ control_genes_pos$Age==35, data = control_genes_pos, FUN = mean)
set.zero.young2<-set.zero.young2[2,2:ncol(set.zero.young2)]
set.zero.young2<- as.numeric(unlist(set.zero.young2))
control_genes_pos2 <- sweep(control_genes_pos, 2, set.zero.young2, "-")  #"/"


################################################################################# we select based on "X"

control_genes_pos3<-control_genes_pos2
control_genes_pos3$Age<-control_genes_pos$Age
control_genes_pos3$Mean<-control_genes_pos2$Mean
control_genes_pos3$SD<-control_genes_pos2$SD
control_genes_pos3$Mean.post<- apply(control_genes_pos3[,1:(ncol(control_genes_pos3)-3)], 1, median, na.rm=T)



################################################################################# large df to plot all indiv. loess
df_long <- control_genes_pos3 |>
  tidyr::pivot_longer(
    cols = -c(Age,Mean,SD,Mean.post),     
    names_to = "variable",
    values_to = "value")



ggplot(df_long, aes(x = Age, y = (2/pi)*(value), group = variable)) +  
  geom_smooth(method="loess", span=0.7, show.legend = FALSE, se=FALSE, linewidth=3, color = scales::alpha("#6b0d8c", 0.01))+  
  geom_hline(yintercept=0, color="white",linewidth=0.5,linetype=2)+
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1),  
    axis.ticks = element_line(linewidth = 1), 
    axis.text = element_text(size = 16, color = "black"),   
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold") 
  )+
  ylab("Z-score")+
  xlab("Age")+
  ggtitle("CoM+ genes per sample (GSE36192+GSE33000+GSE15222)")


#######################################################################################
#####  COM -  (direct on samples)
#######################################################################################

control_genes_neg<-subset(data, select= intersect(colnames(data),genes_neg$gene_symbol))
control_genes_neg$Mean<- apply(control_genes_neg, 1, median, na.rm=T)
control_genes_neg$SD<-apply(control_genes_neg[,1:(ncol(control_genes_neg)-1)], 1, sd, na.rm = TRUE)
control_genes_neg$Age<-info$Age


set.zero.young<-aggregate(. ~ control_genes_neg$Age==35, data = control_genes_neg, FUN = mean)
set.zero.young<-set.zero.young[2,2:ncol(set.zero.young)]
set.zero.young<- as.numeric(unlist(set.zero.young))
control_genes_neg2 <- sweep(control_genes_neg, 2, set.zero.young, "-")  


control_genes_neg3<-control_genes_neg2
control_genes_neg3$Age<-control_genes_neg$Age
control_genes_neg3$Mean<-control_genes_neg2$Mean
control_genes_neg3$SD<-control_genes_neg2$SD
control_genes_neg3$Mean.post<- apply(control_genes_neg3[,1:(ncol(control_genes_neg3)-3)], 1, median, na.rm=T)


################################################################################# large df to plot all indiv. loess
df_long2 <- control_genes_neg3 |>
  tidyr::pivot_longer(
    cols = -c(Age,Mean,SD,Mean.post),      
    names_to = "variable",
    values_to = "value")


ggplot(df_long2, aes(x = Age, y = (2/pi)*(value), group = variable)) +  # (2/pi)*atan
  geom_smooth( method="loess", span=0.7, show.legend = FALSE, se=FALSE, linewidth=3,color = scales::alpha("gold4", 0.01))+  #loess, lm, span 0-1 (defautl 0.75)
  geom_hline(yintercept=0, color="white",linewidth=0.5,linetype=2)+
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1),  
    axis.ticks = element_line(linewidth = 1), 
    axis.text = element_text(size = 16, color = "black"),   
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold") 
  )+
  ylab("Z-score")+
  xlab("Age")+
  ggtitle("CoM- genes per sample (GSE36192+GSE33000+GSE15222)")

#######################################################################################
#####  random genes  (direct on samples)
#######################################################################################


random.pos<-subset(data, select= setdiff(colnames(data),genes_pos$gene_symbol))
random.pos<-subset(data, select= setdiff(colnames(random.pos),genes_neg$gene_symbol))

numeric_cols <- function(df) {subset(df, select = sapply(df, is.numeric))}
random.pos<-numeric_cols(random.pos)

intersect(colnames(random.pos),colnames(genes_pos$gene_symbol))

set.seed(4)
POS.genes_random <- as.data.frame(sample(random.pos,nrow(control_genes_pos)))


POS.genes_random$Mean<-apply(POS.genes_random, 1, median, na.rm=T)
POS.genes_random$SD<-apply(POS.genes_random[,1:(ncol(POS.genes_random)-1)], 1, sd, na.rm = TRUE)
POS.genes_random$Age<-info$Age


table(info$Array)

set.zero.young.random<-aggregate(. ~ POS.genes_random$Age==35, data = POS.genes_random, FUN = mean)
set.zero.young.random<-set.zero.young.random[2,2:ncol(set.zero.young.random)]
set.zero.young.random<- as.numeric(unlist(set.zero.young.random))
POS.genes_random2 <- sweep(POS.genes_random, 2, set.zero.young.random, "-")  #"/"
POS.genes_random2$Age<-POS.genes_random$Age



df_long3 <- POS.genes_random2 |>
  tidyr::pivot_longer(
    cols = -c(Age,Mean,SD),     
    names_to = "variable",
    values_to = "value")


ggplot(df_long3, aes(x = Age, y = (2/pi)*atan(value), group = variable)) +  # (2/pi)*
  geom_smooth( method="loess", span=0.7, show.legend = FALSE, se=FALSE, linewidth=3,color = scales::alpha("gray35", 0.01))+  #loess, lm, span 0-1 (defautl 0.75)
  geom_hline(yintercept=0, color="white",linewidth=0.5,linetype=2)+
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 1),  
    axis.ticks = element_line(linewidth = 1), 
    axis.text = element_text(size = 16, color = "black"),   
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold") 
  )+
 ylab("Z-score")+
  xlab("Age")+
  ggtitle("Random genes per sample (GSE36192+GSE33000+GSE15222)")

