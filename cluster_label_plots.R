library(ggplot2)
library(tidyr)
library(imager)
#Loading the cluster labels (generated from different methods)
cluster_labels = read.csv('Threshold_Labels_Full_Data.csv')
cluster_labels = cluster_labels$X0

image1 = matrix(cluster_labels[1:3551],nrow = 67,ncol = 53,byrow = FALSE)
image2 = matrix(cluster_labels[3552:7519],nrow = 62,ncol = 64,byrow = FALSE)
image3 = matrix(cluster_labels[7520:11094],nrow = 65,ncol = 55,byrow = FALSE)
image4 = matrix(cluster_labels[11095:13168],nrow = 61,ncol = 34,byrow = FALSE)
image5 = matrix(cluster_labels[13169:15859],nrow = 39,ncol = 69,byrow = FALSE)
image6 = matrix(cluster_labels[15860:18883],nrow = 56,ncol = 54,byrow = FALSE)
image7 = matrix(cluster_labels[18884:21627],nrow = 56,ncol = 49,byrow = FALSE)
image8 = matrix(cluster_labels[21628:25839],nrow = 78,ncol = 54,byrow = FALSE)
image9 = matrix(cluster_labels[25840:28039],nrow = 50,ncol = 44,byrow = FALSE)

# Convert the matrix to a data frame
img1 <- as.data.frame(as.table(image1))
colnames(img1) <- c("Column", "Row", "Value")
img2 <- as.data.frame(as.table(image2))
colnames(img2) <- c("Column", "Row", "Value")
img3 <- as.data.frame(as.table(image3))
colnames(img3) <- c("Column", "Row", "Value")
img4 <- as.data.frame(as.table(image4))
colnames(img4) <- c("Column", "Row", "Value")
img5 <- as.data.frame(as.table(image5))
colnames(img5) <- c("Column", "Row", "Value")
img6 <- as.data.frame(as.table(image6))
colnames(img6) <- c("Column", "Row", "Value")
img7 <- as.data.frame(as.table(image7))
colnames(img7) <- c("Column", "Row", "Value")
img8 <- as.data.frame(as.table(image8))
colnames(img8) <- c("Column", "Row", "Value")
img9 <- as.data.frame(as.table(image9))
colnames(img9) <- c("Column", "Row", "Value")

# Define individual colors for each value in the matrix
value_colors <- c("#c6dbef", "#6baed6","#2171b5", "#08306b")

# Plot the matrix as a heatmap with individual colors for each value
#Image 1
ggplot(img1, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 2
ggplot(img2, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 3
ggplot(img3, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 4
ggplot(img4, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 5
ggplot(img5, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 6
ggplot(img6, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 7
ggplot(img7, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 8
ggplot(img8, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)

#Image 9
ggplot(img9, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(values = value_colors) +
  theme_minimal()+theme(axis.text.x = element_blank(), 
                        axis.ticks.x = element_blank(), 
                        axis.title.x = element_blank(),axis.text.y = element_blank(),                                                                                                                                                                                                axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank(),legend.position="none"
  )+coord_flip()+scale_x_discrete(limits=rev)
