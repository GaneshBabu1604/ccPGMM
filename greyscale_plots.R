library(ggplot2)
library(tidyr)
library(imager)
#Loading the data
Data = read.csv('Data.csv')
Data$X <-NULL
mean_Data<-apply(Data,1,mean)

image1 = matrix(mean_Data[1:3551],nrow = 67,ncol = 53,byrow = FALSE)
image2 = matrix(mean_Data[3552:7519],nrow = 62,ncol = 64,byrow = FALSE)
image3 = matrix(mean_Data[7520:11094],nrow = 65,ncol = 55,byrow = FALSE)
image4 = matrix(mean_Data[11095:13168],nrow = 61,ncol = 34,byrow = FALSE)
image5 = matrix(mean_Data[13169:15859],nrow = 39,ncol = 69,byrow = FALSE)
image6 = matrix(mean_Data[15860:18883],nrow = 56,ncol = 54,byrow = FALSE)
image7 = matrix(mean_Data[18884:21627],nrow = 56,ncol = 49,byrow = FALSE)
image8 = matrix(mean_Data[21628:25839],nrow = 78,ncol = 54,byrow = FALSE)
image9 = matrix(mean_Data[25840:28039],nrow = 50,ncol = 44,byrow = FALSE)

#Plotting greyscale
#Image 1
as.cimg(t(image1)) %>% plot(axes=FALSE)
#Image 2
as.cimg(t(image2)) %>% plot(axes=FALSE)
#Image 3
as.cimg(t(image3)) %>% plot(axes=FALSE)
#Image 4
as.cimg(t(image4)) %>% plot(axes=FALSE)
#Image 5
as.cimg(t(image5)) %>% plot(axes=FALSE)
#Image 6
as.cimg(t(image6)) %>% plot(axes=FALSE)
#Image 7
as.cimg(t(image7)) %>% plot(axes=FALSE)
#Image 8
as.cimg(t(image8)) %>% plot(axes=FALSE)
#Image 9
as.cimg(t(image9)) %>% plot(axes=FALSE)
