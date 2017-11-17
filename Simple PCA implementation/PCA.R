#generating sample
one = rnorm(n = 10, mean = 5, sd = 2)
two = rnorm(n = 10, mean = 10, sd = 1)
three = rnorm(n = 10, mean = 1, sd = 2)

Data = cbind(one,two,three)

R = Data

#Mean vector 
mv = matrix(ncol = 3, nrow = 1)
for(i in 1:3){
  mv[i] = mean(R[,i])
}
mv = matrix(mv, ncol = 3, nrow = 10, byrow = T)

#Sd vector
sv = matrix(ncol = 3, nrow = 1)
for(i in 1:3){
  sv[i] = sd(R[,i])
}
sv = matrix(sv, ncol = 3, nrow = 10, byrow = T)

#R centered 
Rcs = (R-mv)/sv

#X for factorial analysis
X = Rcs/sqrt(10)

#C correlation
C = t(X) %*% X

#Eigenvectors + values 
eig = eigen(C)
U = eig$vectors 
eigenvalues = eig$values

psi = X %*% U
gamma = Rcs %*% U #On multiplie Rcs par les vecteurs propres ce qui nous donne les coordonnées dans le nouveau plan.


Ul = matrix(ncol = 3,nrow = 3)


Ul = t(X) %*% gamma

for(i in 1:3){
  Ul[,i] = Ul[,i] / sqrt(eigenvalues[i])
}


plot(psi[,1], psi[,2], main = "PCA implementation", xlab = round(eig$values[1]/3,2), ylab = round(eig$values[2]/3,2))
abline(h=0)
abline(v=0)

plot(Ul[,1], Ul[,2])
abline(h=0)
abline(v=0)

#To add new individuals you need to get them to the format of X, so R-mean/sqrt(n)*s
#Then X+ * Ualpha (eigenvector alpha)
