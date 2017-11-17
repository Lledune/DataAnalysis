library(RCurl)
library(plotrix)
#Getting data
file = getURL("https://raw.githubusercontent.com/dfghdf89/DataAnalysis/master/eurojob.txt")
Data = read.csv(text = file, sep = "\t")

R = Data

#Preparing Data
R = R[2:ncol(R)]
rownames(R) = Data[,1]

#R centered 
Rcs = scale(R)
Rcs = as.matrix(Rcs)

#X for factorial analysis
X = Rcs/sqrt(nrow(R))
X = as.matrix(X)

#C correlation
C = t(X) %*% X

#Eigenvectors + values 
eig = eigen(C)
U = eig$vectors 
U[,2] = U[,2] *-1    #Manual axis correction
eigenvalues = eig$values

psi = X %*% U
gamma = Rcs %*% U #On multiplie Rcs par les vecteurs propres ce qui nous donne les coordonn?es dans le nouveau plan.


Ul = matrix(ncol = ncol(R),nrow = nrow(R))


Ul = t(X) %*% psi

for(i in 1:ncol(R)){
  Ul[,i] = Ul[,i] / sqrt(eigenvalues[i])
}

plot(gamma[,1], gamma[,2], asp=1, xlim = c(-7, 7), main = "PCA implementation", xlab = round(eig$values[1]/ncol(R),2), ylab = round(eig$values[2]/ncol(R),2))
text(gamma[,1], gamma[,2], labels=rownames(R), cex= 0.7, pos = 4)
abline(h=0)
abline(v=0)


plot(Ul[,1], Ul[,2], xlim = c(-1, 1), ylim = c(-1, 1),  asp=1)
text(Ul[,1], Ul[,2], labels=colnames(R), cex= 0.7, pos = 3)

draw.circle(0,0,1)

abline(h=0)
abline(v=0)

for(i in 1:ncol(X)){
  segments(0, 0, Ul[i,1], Ul[i,2])
}

#To add new individuals you need to get them to the format of X, so R-mean/sqrt(n)*s
#Then X+ * Ualpha (eigenvector alpha)


colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

