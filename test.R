library(FactoMineR)
library(Rcmdr)
library(RcmdrPlugin.FactoMineR)
library(RCurl)

#There is no order, use for test only.

#Data loading from github
file = getURL("https://raw.githubusercontent.com/dfghdf89/DataAnalysis/master/eurojob.txt")
Data = read.csv(text = file, sep = "\t")

#PCA ANALYSIS, ncp = 9 (all dimensions kept)
PCA(Data, quali.sup = c(1), ncp = 9)

#HCPC ON PCA, without consolidatio, (hierarchical clustering)
#Click to cut the tree where wanted on R plot !
HC = HCPC(res, consol = FALSE, nb.clust = 3)

#Cluster group 
HC$data.clust

#HCPC with consolidation 
HCPC(res, consol = TRUE)

#HCPC 4C¨, consol = true
PCA(Data, quali.sup = c(1), ncp = 4)
HCPC(res, consol = TRUE)


#new data creation 
matrix = matrix(rnorm(5*3), nrow = 5, ncol = 3)
#Ward algotihm : function hclust
dist = dist(matrix)
Delta = dist^2 / (2 * nrow(matrix))

###signif(Delta)

res.ward = hclust(Delta, method = "ward.D")
plot(res.ward)


#Preparing HCPC data for analysis 
res.hcpc = 