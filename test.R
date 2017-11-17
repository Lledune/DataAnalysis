library(FactoMineR)
#library(Rcmdr)
#library(RcmdrPlugin.FactoMineR)
library(RCurl)

#There is no order, use for test only.

#Data loading from github
file = getURL("https://raw.githubusercontent.com/dfghdf89/DataAnalysis/master/eurojob.txt")
Data = read.csv(text = file, sep = "\t")
Eurojob = Data
#PCA ANALYSIS, ncp = 9 (all dimensions kept)
res.pca = PCA(Data, quali.sup = c(1), ncp = 9)

#HCPC ON PCA, without consolidatio, (hierarchical clustering)
#Click to cut the tree where wanted on R plot !
HC = HCPC(res.pca, consol = FALSE, nb.clust = 3)

#Cluster group 
HC$data.clust

#HCPC with consolidation 
HCPC(res.pca, consol = TRUE)

#HCPC 4C¨, consol = true
res.pca2 = PCA(Data, quali.sup = c(1), ncp = 4)
HCPC(res.pca, consol = TRUE)

#################################################################################
#new data creation 
matrix = matrix(rnorm(5*3), nrow = 5, ncol = 3)

#Ward algotihm : function hclust
dist = dist(matrix)
Delta = dist^2 / (2 * nrow(matrix))

###signif(Delta)

res.ward = hclust(Delta, method = "ward.D")
plot(res.ward)


#Preparing HCPC data for analysis 
res.hcpc = HCPC(res.pca, consol = FALSE)
Dataclust = res.hcpc$data.clust

library(MASS)

#discriminant analysis
res.lda <- FDA(EuroClust,"clust")

############################################################
#More convenient fonction (used in TP) 



# Factorial Discriminant Analysis function
# Created by Cedric Taverne (on a Johan Segers' Code)
# from the Institut de Statistique, Biostatistique et Sciences Actuarielles (ISBA)
# of the Universit? Catholique de Louvain (UCL), Belgium
# Contact : Cedric Taverne : cedric.taverne@uclouvain.be
#
FDA <- function( Data , groups , X=NULL , groups.na.omit=TRUE , var.na.omit=FALSE ,
                 HuygensDetails=FALSE , IndPlot=TRUE , RowNamesPrint=TRUE , CorrPlot=TRUE ){
  #
  # Arguments
  # ---------
  # Data : A dataset (matrix or data.frame) with the variables used to discriminate the groups. The variable containing the groups can be also included. If the dataset contains other variables which are not used, the "X" argument have to appear. 
  # groups : If the variable containing the groups is not included in "Data", "groups" must be a vector indicating the group of each row in "Data". If the variable containing the groups is included in "Data", "groups" could be either the name or the position of the variable containing the information on the groups in "Data".
  # X : If "Data" contains non used variables, the variables used have to be declared in a vector form containing either the names of the variables or the position of the variables.
  # groups.na.omit : if TRUE, the individuals with NA on the "groups" variable are excluded. The default is TRUE.
  # var.na.omit : if TRUE, the individuals with NA on one of the "X" variables are excluded. If FALSE, missing values on quantitative variables are imputed by the mean of the variable. Missing values and empty levels on the qualitative variables are imputed by a "NA" level with a warning only for empty level. The default is FALSE.
  # HuygensDetails : if TRUE, the matrices T, B and W are returned. The default is FALSE.
  # IndPlot : if TRUE, a plot of the individuals on the facotrial axes is produced (with markers, colors and legend). The default is TRUE.
  # RowNamesPrint : if TRUE, the row names are printed on the individuals plot. The default is FALSE.
  # CorrPlot : if TRUE, a correlation circle plot is produced with the original variables and the factorial axes. The default is TRUE.
  #
  # Values
  # ------
  # Returns a list including: 
  # R2 : A pseudo coefficient of determination corresponding to 100 times the trace of B over the trace of T. It is an indicator of the proportion of the general dispersion which could be linked to the groups.
  # eig : The eigenvalues of T-1B. It is the ?discrimination power? of z = Xu and represents the fraction of the total dispersion of X which is due to the difference between the groups.
  # U : The eigenvectors of T-1B. Useful to calculate the position of a new observation.
  # Z : The centered discriminant variables.
  # R : The correlation matrix between the original and the discriminant variables.
  # T : The total variance matrix. Only returned if HuygensDetails=TRUE or if the Huygens control check exceed 1e-10.
  # W : The overall within group variance matrix. Only returned if HuygensDetails=TRUE or if the Huygens control check exceed 1e-10.
  # B : The variance matrix between the groups (between the centers of gravity of each group. Only returned if HuygensDetails=TRUE or if the Huygens control check exceed 1e-10.
  # HuygensError : The maximum computation error in the Huygens Theorem check : T-B-W. Only returned if the Huygens control check exceed 1e-10.
  # Call : A list containg : $G : a factor containg the group variable ; $DesignSummary : a table containing the information on the variable nature, position, number of levels, reference level and the presence/absence of NA ; $X : the complete design matrix. Be carefull, the reference level dummy of each factor variable is also included. Then, this design matrix is exactly singular.
  #
  ########################
  #
  # Begin of the function
  #
  # Separating the discriminant variables and the groups variable
  if( is.null(X)==FALSE ){
    if( length(groups)>1 ){
      G <- groups
      X <- Data[, X ]
      X <- X[, apply( X , 2 , function(x){ sum( x == G , na.rm=TRUE ) == length(x) } ) == FALSE ]
    } else if( is.character(groups) ){
      if( sum( colnames(Data)==groups ) == 1 ){
        G <- eval(parse(text=paste("Data$",groups,sep="")))
      } else stop("The groups variable is misspecified")
      X <- Data[, X ]
      X <- X[, (colnames(X)==groups) == FALSE ]
    } else if( is.numeric(groups) ){
      if( length(groups) == 1 & groups <= dim(Data)[2] ){
        G <- Data[,groups]
      } else stop("The groups variable is misspecified")
      X <- Data[, X ]
      X <- X[, 1:dim(X)[2] != groups ]
    } else stop("The groups variable is misspecified")
  } else if( length(groups)>1 ){
    G <- groups
    X <- Data[, apply( Data , 2 , function(x){ sum( x == groups , na.rm=TRUE ) == length(x) } ) == FALSE ]
  } else if( is.character(groups) ){
    if( sum(colnames(Data)==groups)==1 ){
      G <- eval(parse(text=paste("Data$",groups,sep="")))
    } else stop("The groups variable is misspecified")
    X <- Data[, (colnames(Data)==groups) == FALSE ]
  } else if( is.numeric(groups) ){
    if( length(groups) == 1 & groups <= dim(Data)[2] ){
      G <- Data[,groups]
    } else stop("The groups variable is misspecified") 
    X <- Data[, 1:dim(Data)[2] != groups ]
  } else stop("The groups variable is misspecified")
  #
  # Check if the dimensions matchs
  if( dim(X)[1] != length(G) ){
    stop( paste("The group variable has", length(G), "observations whereas the dataset has",
                dim(X)[1], "lines.") )
  }
  #
  # Remove (groups == NA)'s lines if necessary (based on groups.na.omit argument)
  if( groups.na.omit ){
    X <- X[which(is.na(G)==FALSE),]
    G <- na.omit(G)
  } else {
    if( ("NA" %in% levels(G))==FALSE ){
      G <- factor( G , levels=c(levels(G),"NA") )
    }
    G[is.na(G)] <- "NA"
  }
  #
  # NA imputation for the other variables
  if( var.na.omit ){
    for(i in 1:dim(X)[2] ){
      if( is.factor(X[,i]) ){ 
        if( "" %in% levels(X[,i]) ){
          is.na(X[,i]) <- (X[,i] == "") 
          X[,i] <- factor( X[,i] , levels=levels(X[,i])[levels(X[,i])!=""] )
          warning( paste("The empty level on the qualitative variable", colnames(X)[i],
                         "has been considered as NA") )
        } } }
    G <- G[ which( rowSums( is.na(X) ) == 0 ) ]
    X <- na.omit(X)
  } else { 
    for(i in 1:dim(X)[2] ){
      if( is.factor(X[,i]) ){ 
        if( "" %in% levels(X[,i]) ){
          is.na(X[,i]) <- (X[,i] == "") 
          X[,i] <- factor( X[,i] , levels=levels(X[,i])[levels(X[,i])!=""] )
          warning( paste("The empty level on the qualitative variable", colnames(X)[i],
                         "has been considered as NA") )
        }
        if( sum(is.na(X[,i]))>0 ){
          if( ("NA" %in% levels(X[,i]))==FALSE ){
            X[,i] <- factor( X[,i] , levels=c(levels(X[,i]),"NA") )
          }
          X[is.na(X[,i]),i] <- "NA"
        }
      } else {
        X[is.na(X[,i]),i] <- mean(X[,i],na.rm=TRUE)
      } } }
  #
  # Tracing zero variance variables and the empty levels on qualitative variable 
  Comment <- NA
  i <- 1
  while(i <= dim(X)[2]){
    if( is.factor(X[,i]) ){ 
      if( length(levels(X[,i])) != length(unique(X[,i])) ){
        EmptyLevel <- as.data.frame( table(X[,i]) != 0 )
        X[,i] <- factor( X[,i] , levels=rownames(EmptyLevel)[EmptyLevel[,1]==TRUE] )
        j <- 1
        while( j <= length(rownames(EmptyLevel)[EmptyLevel[,1]==FALSE]) ){
          if( is.list(Comment)==FALSE ){
            Comment <- list( paste("The level",rownames(EmptyLevel)[EmptyLevel[,1]==FALSE][j],
                                   "of the qualitative variable",colnames(X)[i],
                                   "has been deleted because it has no observation.") )
            list.indent <- 1
          } else {
            list.indent <- list.indent + 1
            Comment[[list.indent]] <- paste("The level",rownames(EmptyLevel)[EmptyLevel[,1]==FALSE][j],
                                            "of the qualitative variable",colnames(X)[i],
                                            "has been deleted because it has no observation.")
          }
          j <- j + 1				
        }
      }
    } else {
      if( is.na(sd(X[,i])) | sd(X[,i])==0 ){
        if( is.na(Comment) ){
          Comment <- paste("The quantitative variable",colnames(X)[i],
                           "has been deleted because its standard deviation is zero.")
        } else {
          Comment <- list( Comment , paste("The quantitative variable",colnames(X)[i],
                                           "has been deleted because its standard deviation is zero.") )
        }
        X <- X[,colnames(X)!=colnames(X)[i]]
      } 
    }	
    i <- i + 1 
  }
  #
  # Design matrix (1/2)
  DesignSummary <- cbind.data.frame( sapply( X , function(x){ is.factor(x) } ) , 0 ,
                                     sapply( X , function(x){ length( levels(x) ) } ) ,
                                     sapply( X , function(x){ if(is.factor(x)){ levels(x)[1] } else { NA } } ) ,
                                     sapply( X , function(x){ if(is.factor(x)){ length(levels(x))-length(unique(x)) } else { sum(sd(x)==0) } } ) )
  DesignSummary[,4] <-  paste( rownames(DesignSummary) , DesignSummary[,4] , sep="" )
  colnames(DesignSummary) <- c("Factor","Position","Level","Reference","EmptyLevel")
  X <- X[, order(DesignSummary[,1]) ]
  DesignSummary <- DesignSummary[ order(DesignSummary[,1]) , ]
  #
  # Remove the variables with no variance
  if( sum( DesignSummary$Level==1 ) > 0 ){
    warning( paste( "The variable" , rownames(DesignSummary)[DesignSummary$Level==1] , 
                    "has been removed since it presents no variance. \n " ) )
  }
  X <- X[, DesignSummary$Level!=1 ]
  DesignSummary <- DesignSummary[ DesignSummary$Level!=1 , ]
  #
  # Design matrix (2/2)
  X <- model.matrix( ~ . , X )
  X <- X[,2:dim(X)[2]]
  for(i in 1:dim(DesignSummary)[1]){
    if( DesignSummary[i,1]==FALSE ){
      if( i==1 ){ counter = i }
      DesignSummary[i,2] <- counter
      counter = counter + 1 
    } else {
      if( i==1 ){ counter = i }
      DesignSummary[i,2] <- counter
      counter = counter + DesignSummary[i,3] - 1
    } }
  #
  # Determining n and p
  n <- dim(X)[1]
  p <- dim(X)[2]
  #
  # Groups names and corresponding internal number
  g.names <- unique(G)
  g <- 1:length(g.names)
  #
  # Mean with check function of groups of size 1
  mymean <- function(X) {
    if (length(dim(X)) > 1) { 
      return( apply( X, MARGIN = 2, FUN = "mean" ) )
    } else { 
      return(X) 
    } }
  #
  # Creating X1 to Xg, ni, x1 to xg, S1 to Sg
  ni <- rep(NA,length(g))
  Si <- array(NA,c(p,p,length(g)))
  for(i in 1:length(g)){
    eval(parse(text=paste("X",g[i]," <- X[which(G==g.names[i]),]",sep=""))) 
    if( is.matrix( eval(parse(text=paste("X",g[i],sep=""))) ) ) { 
      eval(parse(text=paste("ni[",g[i],"] <- dim(X",g[i],")[1]",sep="")))
    } else { 
      eval(parse(text=paste("ni[",g[i],"] <- 1",sep="")))
    }
    eval(parse(text=paste("x",g[i]," <- mymean(X", g[i], ")", sep = "")))
    if( is.matrix( eval(parse(text=paste("X",g[i],sep=""))) ) ) { 
      eval(parse(text=paste("Si[,,",g[i],"] <- (1/ni[",g[i],"]) * t(X",g[i],") %*% X",g[i],
                            " - x",g[i]," %o% x",g[i],sep="")))
    } else { 
      eval(parse(text=paste("Si[,,",g[i],"] <- (1/ni[",g[i],"]) * X",g[i]," %o% X",g[i],
                            " - x",g[i]," %o% x",g[i],sep="")))
    } }
  #
  # Creating W
  W <- matrix(data = 0, nrow = p, ncol = p)
  for (i in 1:length(g)) {
    W <- W + ni[i] * Si[,,i]
  }
  W <- W / n
  rownames(W) <- colnames(W) <- colnames(X)
  #
  # Creating x
  x <- apply( X , 2, 'mean' )
  #
  # Creating T
  T <- (1/n) * t(X) %*% X - x %o% x
  #
  # Creating xi and Ci
  xi <- matrix(NA,p,length(g))
  for(i in 1:length(g)){ 
    xi[,i] <- eval(parse(text=paste("x",g[i],sep=""))) 
  }
  rownames(xi) <- colnames(X)
  colnames(xi) <- g.names
  Ci <- (xi-x) %*% diag(sqrt(ni/n))
  colnames(Ci) <- g.names
  #
  # Creating B
  B <- Ci %*% t(Ci)
  #
  # This is a recent function, so we pluged in some validation test. 
  # 	Here is a validation test for the Huygens' theorem
  # HuygensTolerance = Computation error tolerance accepted for the 
  #	Huygens Theorem check : T - B - W <= 1e-10 
  HuygensTolerance <- 1e-10
  HuygensError <- max(abs(T - W - B))
  if( HuygensError > HuygensTolerance ){
    warning( paste("The maximum Huygens' computation error is",round(HuygensError,7),
                   "It could remain an error in our code. Please contact cedric.taverne@uclouvain.be") )
  }
  #
  # Factorial Discriminant Analysis
  require(MASS)
  out <- eigen( t(Ci) %*% ginv(T) %*% Ci , symmetric=TRUE)
  #
  # Pseudo coefficient of determination
  R2 <- sum(diag(B))/sum(diag(T))*100
  #
  # Extracting the Lambdas and the Vectors
  lambda <- out$values
  Lambda <- cbind(lambda[1:(length(g)-1)],
                  lambda[1:(length(g)-1)]/sum(lambda[1:(length(g)-1)])*100,
                  lambda[1:(length(g)-1)]/sum(lambda[1:(length(g)-1)])*100)
  if(length(lambda)>2){
    for(i in 2:(length(g)-1)){
      Lambda[i,3] <- Lambda[(i-1),3] + Lambda[i,2]
    } }
  colnames(Lambda) <- c("eigenvalues","pct.of.between.inertia",
                        "cum.pct.of.between.inertia")
  rownames(Lambda) <- paste("Dim", g[1:(length(g)-1)] , sep=".")
  v <- out$vectors
  #
  # Creating ui and zi
  if( length(g)==2 ){
    U <- ginv(T) %*% Ci %*% v[,1:(length(g)-1)] %*% diag( 1 )
  } else {
    U <- ginv(T) %*% Ci %*% v[,1:(length(g)-1)] %*% diag( 1/ sqrt(lambda[1:(length(g)-1)]) )
  }
  Z <- X %*% U
  Z <- apply( Z , 2 , function(x) x-mean(x) )
  colnames(Z) <- colnames(U) <- paste("Dim", g[1:(length(g)-1)] , sep=".")
  #
  # Correlation with the original variables
  if(sum(DesignSummary[,1]==FALSE)>0){
    XFull <- X[,1:sum(DesignSummary[,1]==FALSE)]
  }
  if( sum(DesignSummary[,1]==FALSE)<p ){
    j <- sum(DesignSummary[,1]==FALSE) + 1
    jj <- 0
    while( j <= dim(DesignSummary)[1] ){
      if( j==1 ){
        XFull <- cbind( X[,DesignSummary[j,2]:(DesignSummary[j,2]+DesignSummary[j,3]-2)] , 1 )
      } else {
        XFull <- cbind( XFull , X[,DesignSummary[j,2]:(DesignSummary[j,2]+DesignSummary[j,3]-2)] , 1 )
      }
      for(i in (DesignSummary[j,2]+jj):(DesignSummary[j,2]+DesignSummary[j,3]-2+jj) ){
        XFull[,(DesignSummary[j,2]+DesignSummary[j,3]-1+jj)] <- XFull[,(DesignSummary[j,2]+DesignSummary[j,3]-1+jj)] -
          XFull[,i]
      }
      colnames(XFull)[(dim(XFull)[2]-DesignSummary[j,3]+1):dim(XFull)[2]] <- 
        c( colnames(X)[DesignSummary[j,2]:(DesignSummary[j,2]+DesignSummary[j,3]-2)], DesignSummary[j,4] )
      j <- j + 1
      jj <- jj + 1
    } }
  R <- cor(cbind(XFull,Z))
  #
  # Individual Plot
  if(IndPlot){
    if(length(lambda)>2){
      dev.new()
      for(i in 1:length(g)){
        eval(parse(text=paste("I",g[i]," <- which(G==g.names[i])",sep=""))) 
      }
      plot(Z[I1,1], Z[I1,2], pch=1,
           xlim = c(min(Z[,1]), max(Z[,1])),
           ylim = c(min(Z[,2]), max(Z[,2])),
           main = "Individual Plot",
           xlab = paste("Dim 1 (",round(lambda[1]/sum(lambda)*100,2),"%)",sep=""), 
           ylab = paste("Dim 2 (",round(lambda[2]/sum(lambda)*100,2),"%)",sep="") )
      abline(h = 0, col = "gray")
      abline(v = 0, col = "gray")
      for(i in 2:length(g)){	
        eval(parse(text=paste("points(Z[I",g[i],",1], Z[I",g[i],",2], pch=",g[i],
                              ", col=",g[i],")",sep=""))) 
      }
      for(i in 1:length(g)){	
        eval(parse(text=paste("points(mean(Z[I",g[i],",1]), mean(Z[I",g[i],
                              ",2]), pch=0, cex=2.2, col=",g[i],")",sep=""))) 
        eval(parse(text=paste("text(x = mean(Z[I",g[i],",1]), y = mean(Z[I",g[i],
                              ",2]), labels = g.names[",i,"], cex=.8, col=",g[i],")",sep=""))) 
      }
      if(RowNamesPrint){
        for(i in 1:length(g)){	
          eval(parse(text=paste("text(x = Z[I",g[i],",1], y = Z[I",g[i],
                                ",2], labels = rownames(Z)[I",g[i],"], pos = 2, col=",g[i],")",sep=""))) 
        }
      }
      legend( "topright", legend=g.names, col=g, pch=g)
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } else {
      dev.new()
      for(i in 1:length(g)){
        eval(parse(text=paste("I",g[i]," <- which(G==g.names[i])",sep=""))) 
      }
      stripchart(	Z[I1] , pch=1, xlim = c(min(Z), max(Z)), main = "Individual Plot", method="jitter",
                  xlab = "Factorial Discriminant Axis", ylab="Points are artificially exploded vertically" )
      stripchart(	Z[I2] , pch=2, col=2, method="jitter", add=TRUE)
      points( mean(Z[I1]) , 0.85 , pch=0, col=1, cex=2.2)
      points( mean(Z[I2]) , 0.85 , pch=0, col=2, cex=2.2)
      text( mean(Z[I1]),  0.8 , labels = g.names[1], cex=.8, col=1)
      text( mean(Z[I2]),  0.8 , labels = g.names[2], cex=.8, col=2)
      legend( "bottomright", legend=g.names, col=g, pch=g)
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } }
  #
  # Correlation Plot
  if(CorrPlot){
    if(length(lambda)>2){
      dev.new()
      thetas <- 2 * pi * (0:100)/100
      plot(cos(thetas), sin(thetas), col = "gray", type = "l", asp = 1,
           xlab = paste("Dim 1 (",round(lambda[1]/sum(lambda)*100,2),"%)",sep=""), 
           ylab = paste("Dim 2 (",round(lambda[2]/sum(lambda)*100,2),"%)",sep=""),
           main = "Correlation Circle")
      abline(h = 0, col = "gray")
      abline(v = 0, col = "gray")
      pFull <- dim(XFull)[2]
      if(sum(DesignSummary[,1]==FALSE)>0){
        arrows( 0, 0 , R[1:sum(DesignSummary[,1]==FALSE),(pFull+1)] , 
                R[1:sum(DesignSummary[,1]==FALSE),(pFull+2)] , length=0.1 )
      }
      if( sum(DesignSummary[,1]==FALSE)<pFull ){
        points( R[(sum(DesignSummary[,1]==FALSE)+1):pFull,(pFull+1)] , 
                R[(sum(DesignSummary[,1]==FALSE)+1):pFull,(pFull+2)] , pch=0 )
      }
      for(i in 1:pFull){
        if( abs(R[i,(pFull+1)]) >= abs(R[i,(pFull+2)]) ){
          if( R[i,(pFull+1)] >=0 ){
            text(x = R[i,(pFull+1)], y = R[i,(pFull+2)], labels = colnames(XFull)[i], pos = 4)
          } else { 
            text(x = R[i,(pFull+1)], y = R[i,(pFull+2)], labels = colnames(XFull)[i], pos = 2)
          }
        } else if( R[i,(pFull+2)] >=0 ){
          text(x = R[i,(pFull+1)], y = R[i,(pFull+2)], labels = colnames(XFull)[i], pos = 3)
        } else { 
          text(x = R[i,(pFull+1)], y = R[i,(pFull+2)], labels = colnames(XFull)[i], pos = 1)
        } }
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } else {
      dev.new()
      plot(NA, asp = 1, xlim=c(-1,1) , ylim=c(-1,1) , xlab="Correlation with the Factorial Discriminant Axis", 
           ylab="Variables are artificially exploded vertically", main="Correlation Plot", axes=FALSE)
      axis(1) 
      abline(v = c(-1,0,1), col = "gray")
      pFull <- dim(XFull)[2]
      RtoPlot <- sort(R[1:(dim(R)[1]-1),dim(R)[2]])
      y <- seq(-0.9,0.9,length.out=pFull)
      for(i in 1:pFull){	
        arrows( 0 , y[i] , RtoPlot[i] , y[i] , length=0.1 )
        text( RtoPlot[i] , y[i] , labels = names(RtoPlot)[i], pos=(if(RtoPlot[i]<0){2}else{4}) )
      }
      legend( -1, 1 , legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } }
  #
  # Return 
  if( HuygensError > HuygensTolerance ){
    resultats <- list( R2=R2 , eig=Lambda , U=U , Z=Z , R=R , T=T , W=W , B=B , 
                       HuygensError=HuygensError, Call=list( G=G , DesignSummary=DesignSummary , X=XFull ) )
    class(resultats) <- c("FDA", "list")
    return(resultats)
  } else if( HuygensDetails ){ 
    resultats <- list( R2=R2 , eig=Lambda , U=U , Z=Z , R=R , T=T , W=W , B=B, 
                       Call=list( G=G , DesignSummary=DesignSummary , X=XFull ) ) 
    class(resultats) <- c("FDA", "list")
    return(resultats)
  } else { 
    resultats <- list( R2=R2 , eig=Lambda , U=U , Z=Z , R=R, 
                       Call=list( G=G , DesignSummary=DesignSummary , X=XFull ) ) 
    class(resultats) <- c("FDA", "list")
    return(resultats)
  }
  #
  # End of the function
}
########



# Flexible Plot Function for the Output of the Factorial Discriminant Analysis function
# Created by Cedric Taverne 
# from Universit? Catholique de Louvain (UCL), Belgium
# Contact : Cedric Taverne : cedric.taverne@uclouvain.be
#
plot.FDA <- function( res , type='var' , axes=c(1,2) , RowNamesPrint=TRUE , CorrLim=0 , CorrProp=1 ){
  #
  # Arguments
  # ---------
  # res : A list produced by the FDA function. 
  # type : A character corresponding to the type of graph expected. Type could be 'var' for the 
  #	correlation plot or 'ind' for the plot of the individuals.
  # axes : a vector of size two or a scalar corresponding to the axes chosen.
  # RowNamesPrint : if TRUE, the row names are printed on the individuals plot. The 
  #	default is FALSE.
  # CorrLim : The absolute correlation limit to plot the variables. The sum of the absolute correlation with 
  #	both axes is considered. The default is 0 so that any variable is plotted.
  # CorrProp : The proportion of the variables to be ploted based on the quantile of the sum of the 
  #	absolutes correlations with both axes. The default is 1 so that any variable is plotted.
  #	If CorrLim is non-zero, CorrProp is ignored.
  #
  # Values
  # ------
  # Returns a graph of the factorial discriminant analysis.
  #
  #
  # Begin of the function
  #
  # Extracting the output from "res"
  lambda <- res$eig[,1]
  R2 <- res$R2 
  U <- res$U 
  Z <- res$Z
  R <- res$R
  DesignSummary <- res$Call$DesignSummary 
  G <- res$Call$G
  XFull <- res$Call$X
  g <- 1:length(levels(G))
  g.names <- levels(res$Call$G)
  #
  # Check the axes length
  if( length(axes)>2 ) stop("The axes argument is misspecified. It should have a maximum length of 2.")
  #
  # Reverse axes order if necessary
  if( axes[1]>axes[2] ) axes = c(axes[2],axes[1])
  #
  # Check the existence of the demanded axes
  if( missing(axes)==FALSE ){ if( max(axes)>(length(g)-1) ){
    axes = c((length(g)-2),(length(g)-1))
    warning(paste("At least one of the axes you expect to plot does not exist.",
                  "The last couple of axes is plotted instead."))
  }}
  #
  # Individual Plot
  if(type=='ind'){
    if(dim(Z)[2]>1){
      dev.new()
      for(i in 1:length(g)){
        eval(parse(text=paste("I",g[i]," <- which(G==g.names[i])",sep=""))) 
      }
      plot(Z[I1,axes[1]], Z[I1,axes[2]], pch=1,
           xlim = c(min(Z[,axes[1]]), max(Z[,axes[1]])),
           ylim = c(min(Z[,axes[2]]), max(Z[,axes[2]])),
           main = "Individual Plot",
           xlab = paste("Dim ",axes[1]," (",round(lambda[axes[1]]/sum(lambda)*100,2),"%)",sep=""), 
           ylab = paste("Dim ",axes[2]," (",round(lambda[axes[2]]/sum(lambda)*100,2),"%)",sep="") )
      abline(h = 0, col = "gray")
      abline(v = 0, col = "gray")
      for(i in 2:length(g)){	
        eval(parse(text=paste("points(Z[I",g[i],",",axes[1],"], Z[I",g[i],",",axes[2],"], pch=",g[i],
                              ", col=",g[i],")",sep=""))) 
      }
      for(i in 1:length(g)){	
        eval(parse(text=paste("points(mean(Z[I",g[i],",",axes[1],"]), mean(Z[I",g[i],
                              ",",axes[2],"]), pch=0, cex=2.2, col=",g[i],")",sep=""))) 
        eval(parse(text=paste("text(x = mean(Z[I",g[i],",",axes[1],"]), y = mean(Z[I",g[i],
                              ",",axes[2],"]), labels = g.names[",i,"], cex=.8, col=",g[i],")",sep=""))) 
      }
      if(RowNamesPrint){
        for(i in 1:length(g)){	
          eval(parse(text=paste("text(x = Z[I",g[i],",",axes[1],"], y = Z[I",g[i],
                                ",",axes[2],"], labels = rownames(Z)[I",g[i],"], pos = 2, col=",g[i],")",sep=""))) 
        }
      }
      legend( "bottomright", legend=g.names, col=g, pch=g)
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } else {
      dev.new()
      if(missing(axes)==FALSE){ if(length(axes)>1){ 
        warning("It exists only one factorial discriminant axis. So, the 'axes' argument has been ignored.")
      }}	
      for(i in 1:length(g)){
        eval(parse(text=paste("I",g[i]," <- which(G==g.names[i])",sep=""))) 
      }
      stripchart(	Z[I1] , pch=1, xlim = c(min(Z), max(Z)), main = "Individual Plot", method="jitter",
                  xlab = "Factorial Discriminant Axis", ylab="Points are artificially exploded vertically" )
      stripchart(	Z[I2] , pch=2, col=2, method="jitter", add=TRUE)
      points( mean(Z[I1]) , 0.85 , pch=0, col=1, cex=2.2)
      points( mean(Z[I2]) , 0.85 , pch=0, col=2, cex=2.2)
      text( mean(Z[I1]),  0.8 , labels = g.names[1], cex=.8, col=1)
      text( mean(Z[I2]),  0.8 , labels = g.names[2], cex=.8, col=2)
      legend( "bottomright", legend=g.names, col=g, pch=g)
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } }
  #
  # Correlation Plot
  if(type=='var'){
    if(dim(Z)[2]>1){
      dev.new()
      thetas <- 2 * pi * (0:100)/100
      plot(cos(thetas), sin(thetas), col = "gray", type = "l", asp = 1,
           xlab = paste("Dim ",axes[1]," (",round(lambda[axes[1]]/sum(lambda)*100,2),"%)",sep=""), 
           ylab = paste("Dim ",axes[2]," (",round(lambda[axes[2]]/sum(lambda)*100,2),"%)",sep=""),
           main = "Correlation Circle")
      abline(h = 0, col = "gray")
      abline(v = 0, col = "gray")
      pFull <- dim(XFull)[2]
      if(CorrLim<0|CorrLim>1) stop("CorrLim must be in the unit interval")
      if(CorrProp<0|CorrProp>1) stop("CorrLim must be in the unit interval")
      if(CorrLim>0 & CorrProp<1) warning("Since CorrLim is specified, CorrProp is ignored")
      if(CorrLim==0 & CorrProp<1){
        CorrLim <- quantile( ( abs( R[1:pFull,(pFull+axes[1])] ) + 
                                 abs( R[1:pFull,(pFull+axes[2])] ) ) , (1-CorrProp) ) 
      } 		
      DesignSummary <- cbind.data.frame( DesignSummary , 
                                         "CorrLim"=( abs( R[1:pFull,(pFull+axes[1])] ) + abs( R[1:pFull,(pFull+axes[2])] ) ) ,
                                         "CorrLimLog"=TRUE )
      for( i in 1:pFull ){ 
        if(DesignSummary$CorrLim[i]<CorrLim) DesignSummary$CorrLimLog[i] <- FALSE 
      }
      RPrintQt <- R[1:sum(DesignSummary[,1]==FALSE),c((pFull+axes[1]),(pFull+axes[2]))]
      RPrintQt <- RPrintQt[ DesignSummary$CorrLimLog[DesignSummary[,1]==FALSE] , ]
      arrows( 0, 0 , RPrintQt[,1] , RPrintQt[,2] , length=0.1 )
      for(i in 1:dim(RPrintQt)[1]){
        if( abs(RPrintQt[i,1]) >= abs(RPrintQt[i,2]) ){
          if( RPrintQt[i,1] >=0 ){
            text(x = RPrintQt[i,1], y = RPrintQt[i,2], 
                 labels = rownames(RPrintQt)[i], pos = 4)
          } else { 
            text(x = RPrintQt[i,1], y = RPrintQt[i,2], 
                 labels = rownames(RPrintQt)[i], pos = 2)
          }
        } else if( RPrintQt[i,2] >=0 ){
          text(x = RPrintQt[i,1], y = RPrintQt[i,2], 
               labels = rownames(RPrintQt)[i], pos = 3)
        } else { 
          text(x = RPrintQt[i,1], y = RPrintQt[i,2], 
               labels = rownames(RPrintQt)[i], pos = 1)
        } }
      if( sum(DesignSummary[,1]==FALSE)<pFull ){
        RPrintQl <- R[(sum(DesignSummary[,1]==TRUE)+1):pFull,c((pFull+axes[1]),(pFull+axes[2]))]
        RPrintQl <- RPrintQt[ DesignSummary$CorrLimLog[DesignSummary[,1]==TRUE] , ]
        points( RPrintQl[,1] , RPrintQl[,1] , pch=0 )
        for(i in 1:dim(RPrintQl)[1]){
          if( abs(RPrintQl[i,1]) >= abs(RPrintQl[i,2]) ){
            if( RPrintQl[i,1] >=0 ){
              text(x = RPrintQl[i,1], y = RPrintQl[i,2], 
                   labels = rownames(RPrintQl)[i], pos = 4)
            } else { 
              text(x = RPrintQl[i,1], y = RPrintQl[i,2], 
                   labels = rownames(RPrintQl)[i], pos = 2)
            }
          } else if( RPrintQl[i,2] >=0 ){
            text(x = RPrintQl[i,1], y = RPrintQl[i,2], 
                 labels = rownames(RPrintQl)[i], pos = 3)
          } else { 
            text(x = RPrintQl[i,1], y = RPrintQl[i,2], 
                 labels = rownames(RPrintQl)[i], pos = 1)
          } }
      }
      legend( "topleft", legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } else {
      dev.new()
      if(missing(axes)==FALSE){ if(length(axes)>1){ 
        warning("It exists only one factorial discriminant axis. So, the 'axes' argument has been ignored.")
      }}
      plot(NA, asp = 1, xlim=c(-1,1) , ylim=c(-1,1) , xlab="Correlation with the Factorial Discriminant Axis", 
           ylab="Variables are artificially exploded vertically", main="Correlation Plot", axes=FALSE)
      axis(1) 
      abline(v = c(-1,0,1), col = "gray")
      pFull <- dim(XFull)[2]
      RtoPlot <- sort(R[1:(dim(R)[1]-1),dim(R)[2]])
      y <- seq(-0.9,0.9,length.out=pFull)
      for(i in 1:pFull){	
        arrows( 0 , y[i] , RtoPlot[i] , y[i] , length=0.1 )
        text( RtoPlot[i] , y[i] , labels = names(RtoPlot)[i], pos=(if(RtoPlot[i]<0){2}else{4}) )
      }
      legend( -1, 1 , legend=paste("Pseudo-R2 : ",round(R2,2),"%",sep=""), pch=NULL)
    } }
  #
  # End of the function
}
########

################################################
#Assistant code for testing of previous functions
################################################

# rm(list=ls()) 

# Chargement des fonctions
source("FDA_function.R")
source("plot.FDA_function.R")
#source("dimdesc.FDA_function.R")

# Chargement des donn?es
Eurojob <- read.table("eurojob.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
row.names(Eurojob) <- as.character(Eurojob$Country)
Eurojob$Country <- NULL

# Chargement du package FactoMineR
require(FactoMineR)

# PCA + Clustering
res.PCA.4 <- PCA(Eurojob, scale.unit=TRUE, ncp=4, graph=FALSE)
res.HCPC <- HCPC(res.PCA.4, consol=TRUE, nb.clust=-1, graph=FALSE)
EuroClust <- res.HCPC$data.clust

# FDA
EuroFDA <- FDA(EuroClust,"clust")

# Valeurs propres
EuroFDA$eig

# Vecteurs propres
EuroFDA$U

# Variables discriminantes
EuroFDA$Z

# Graphe des individus
plot.FDA(EuroFDA,type='ind',RowNamesPrint=TRUE)

# Corr?lation des variables initiales avec les axes factoriels
round( EuroFDA$R[, dim(EuroClust)[2]:dim(EuroFDA$R)[2] ], 3)

# Graphe des variables
plot.FDA(EuroFDA,type='var')



