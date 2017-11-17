library("deldir") #Install those !!
library("R.matlab")

#########################################################################
#This is an exemple of the Kmeans algorithm (with some variations) with a fixed cluster number of 10.
#########################################################################


setwd("C:/Users/Lucien/Desktop/DataAnalysis/Classification (Kmeans)")
mat=readMat("dataset_1.mat")
matD=as.data.frame(mat)
matD = data.matrix(matD)
stockdistNcen = 0
stockdistNmat = 0
kWinner = 0
stockdistL = 100000
plot(matD)
vect = c(1:length(matD[,1]))
datF = as.data.frame(matD)

centroids = matrix(rnorm(n = 20, mean = 1.25, sd = 0.5),nrow=10, ncol=2, byrow = TRUE)


matD = data.matrix(matD)
centroids = data.matrix(centroids)
#points(centroids, col = "green", lwd = 2)
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

alpha = 0.9
beta = 1

for(i in 1 : 100){
  for(j in 1 : length(matD[,1])){
    for(k in 1 : length(centroids[,1])){
      if(euc.dist(matD[vect[j],],centroids[k,]) < stockdistL){
         stockdistL = euc.dist(matD[vect[j],],centroids[k,])
         kWinner = k
      }
    }
      centroidOld = centroids[kWinner,]
      centroids[kWinner,] = centroids[kWinner,] + alpha * (matD[vect[j],] - centroids[kWinner,])
      stockdistL = 1000
      if(i%%10 == 0 & j%%100 == 0){
      #segments(col = "blue",centroidOld[1], centroids[kWinner,1], centroidOld[2], centroids[kWinner,2])
      }
  }
  alpha = (alpha*beta)/(alpha+beta)
  #shuffle
  vect = sample(vect,length(matD[,1]))
}
 
points(centroids,pch = 17,  col="red")
centroidsDat = as.data.frame(centroids)
voronoi = deldir(centroidsDat$V1, centroidsDat$V2)
plot(voronoi, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1, col = "blue")

###################################################################
#Quantization

Codebook = matrix(c(1,1.5,1.5,1), byrow = TRUE, ncol = 2)
library(R.matlab)
Data = data.matrix(as.data.frame(readMat("dataset_1.mat")))
Concat = NULL

for(i in 1:length(Data[,1])){
  CentroidWinner = 1
  for(j in 1 : length(Codebook[,1])){
    if(euc.dist( Codebook[CentroidWinner,], Data[i,] ) > euc.dist( Codebook[j,], Data[i,]) ){
      CentroidWinner = j
    }
  }
  Concat = c(Concat, CentroidWinner)
}

QuantizedData = cbind(Data, Concat)

#for(i in 1 : length(Data[,1])){
#  segments(Data[i,1], Data[i,2], Codebook[QuantizedData[i,3],1], Codebook[QuantizedData[i,3],2], col = "green")
#}







