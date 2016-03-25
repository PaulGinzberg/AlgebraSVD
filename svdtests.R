



## The two Main functions:
AQRBC <- function(A,algebra=NULL,beta=getbeta(algebra),epsilon=1e-16,usenorm=function(x)max(abs(x)),QHinit=NULL) {
  # This function uses a "by columns" givens rotation algorithm to compute the QR decomposition.
  # ARGUMENTS
  # A is an mxnxd array, representing an mxn matrix with entries from a d-dimensional algebra.
  # algebra is the name of the algebra.
  # beta is a function which takes a d-dimensional vector as input (representing an element from the algebra)
  # and outputs a dxd matrix representation of a unitary element of the algebra, such that (beta(a)%*%a)[1] is large.
  # epsilon is the error tolerance
  # usenorm is a function which takes as input a d-dimensional vector (representing an element from the algebra) and outputs its norm
  # QHinit is an initialisation for QH. The default (NULL) is equivalent to QHinit being an identity matrix.
  # VALUE
  # The output is a mxmxd array QH (representing a unitary mxm matrix with entries from the algebra) and
  # a mxnxd array R (representing an mxn upper triangular matrix with entries from the algebra).
  # Note that R is only approximately upper triangular, in the sense that every entry below the main diagonal
  # has norm at most epsilon
  # in the algebra, QH times A is equal to R.
  # For non-NULL values of QHinit, QH is replaced by QH times QHinit.
  # sweeps is the number of sweeps performed
  # rotations is the number of Givens rotations performed
  m <- dim(A)[1]
  n <- dim(A)[2]
  d <- dim(A)[3]
  if(is.null(QHinit)) {
    QH <- array(0,dim=c(m,m,d))
    QH[cbind(1:m,1:m,1)] <- 1
  } else {
    QH <- QHinit
    rm(QHinit)
  }
  R <- A
  rm(A)
  sweeps <- 0
  rotations <- 0
  norms <- matrix(apply(R,c(1,2),usenorm),nrow=m,ncol=n)
  #g1 <- Inf
  g1 <- max(norms[lower.tri(norms,diag=FALSE)])
  if(m>1) {
    while(g1 > epsilon) {
      sweeps <- sweeps+1
      for(k in 1:min(m,n)) {
        #g2 <- Inf
        # Improve conditioning by making R[k,k,1] larger
        b <- beta(R[k,k,])
        R[k,,] <- t(t(b)%*%t(matrix(R[k,,],nrow=n,ncol=d)))
        QH[k,,] <- t(t(b)%*%t(matrix(QH[k,,],nrow=m,ncol=d)))
        norms[k,] <- apply(R[k,,,drop=FALSE],2,usenorm)
        if(k==m) {break}
        repeat {
          i <- which.max(norms[(k+1):m,k])+k
          g2 <- norms[i,k]
          #print(list(k=k,i=i,g2=g2,g1=g1))
          #print(R)
          #print(norms)
          if(g2 <= epsilon) {break}
          rotations <- rotations+1
          b <- beta(R[i,k,])
          R[i,,] <- t(t(b)%*%t(matrix(R[i,,],nrow=n,ncol=d)))
          # Note that here we have Q^H instead of Q for simplicity
          QH[i,,] <- t(t(b)%*%t(matrix(QH[i,,],nrow=m,ncol=d)))
          theta <- atan2(R[i,k,1],R[k,k,1])
          ctheta <- cos(theta)
          stheta <- sin(theta)
          R[c(k,i),,] <- rbind(c(ctheta,stheta),c(-stheta,ctheta)) %*% matrix(R[c(k,i),,],nrow=2,ncol=n*d)
          QH[c(k,i),,] <- rbind(c(ctheta,stheta),c(-stheta,ctheta)) %*% matrix(QH[c(k,i),,],nrow=2,ncol=m*d)
          #R[k,,] <- ctheta*R[k,,] + stheta*R[i,,]
          #R[i,,] <- -stheta*R[k,,] + ctheta*R[i,,]
          #QH[k,,] <- ctheta*QH[k,,] + stheta*QH[i,,]
          #QH[i,,] <- -stheta*QH[k,,] + ctheta*QH[i,,]
          R[i,,] <- t(b%*%t(matrix(R[i,,],nrow=n,ncol=d)))
          QH[i,,] <- t(b%*%t(matrix(QH[i,,],nrow=m,ncol=d)))
          norms[c(k,i),] <- apply(R[c(k,i),,,drop=FALSE],c(1,2),usenorm)
        }
      }
      g1 <- max(norms[lower.tri(norms,diag=FALSE)])
    }
  } else {
    b <- beta(R[1,1,])
    R[1,,] <- t(t(b)%*%t(matrix(R[1,,],nrow=n,ncol=d)))
    QH[1,,] <- t(t(b)%*%t(matrix(QH[1,,],nrow=m,ncol=d)))
  }
  signchange <- ifelse(diag(matrix(R[,,1],nrow=m,ncol=n))>=0,1,-1)
  if(n<m) {signchange <- c(signchange,rep(1,m-n))}
  R <- R*signchange
  QH <- QH*signchange
  return(list(QH=QH,R=R,sweeps=sweeps,rotations=rotations))
}
ASVD <- function(A,algebra=NULL,beta=getbeta(algebra),Aconj=getAconj(algebra),epsilon=1e-16,usenorm=function(x)max(abs(x))) {
  # This function computes the SVD of A by iterative application of the function AQRBC
  # ARGUMENTS
  # A is an mxnxd array, representing an mxn matrix with entries from a d-dimensional algebra.
  # algebra is the name of the algebra.
  # beta is a function which takes a d-dimensional vector as input (representing an element from the algebra)
  # and outputs a dxd matrix representation of a unitary element of the algebra, such that (beta(a)%*%a)[1] is large.
  # epsilon is the error tolerance
  # usenorm is a function which takes as input a d-dimensional vector (representing an element from the algebra) and outputs its norm
  # VALUE
  # The output is a mxmxd array U (representing a unitary mxm matrix with entries from the algebra).
  # a nxnxd array V (representing a unitary mxm matrix with entries from the algebra), and
  # a mxnxd array D (repersenting an mxn diagonal matrix with entries from the algebra).
  # Note that D is only approximately diagonal, in the sense that every off-diagonal entry
  # has norm at most epsilon
  # in the algebra, A time V is equal to U times D.
  # sweeps is the number of sweeps performed
  # rotations is the number of Givens rotations performed
  # qrds is the number of QR decompositions performed
  m <- dim(A)[1]
  n <- dim(A)[2]
  d <- dim(A)[3]
  UH <- array(0,dim=c(m,m,d))
  UH[cbind(1:m,1:m,1)] <- 1
  VH <- array(0,dim=c(n,n,d))
  VH[cbind(1:n,1:n,1)] <- 1
  D <- A
  rm(A)
  qrds <- 0
  sweeps <- 0
  rotations <- 0
  repeat {
    norms <- matrix(apply(D,c(1,2),usenorm),nrow=m,ncol=n)
    diag(norms) <- 0
    g <- max(norms)
    if(g<=epsilon) {break} # Exit condition
    # Note that here we use U^H and V^H instead of U and V for simplicity
    QRD <- AQRBC(A=D,beta=beta,epsilon=epsilon,usenorm=usenorm,QHinit=UH)
    qrds <- qrds+(QRD$rotations>0)
    sweeps <- sweeps + QRD$sweeps
    rotations <- rotations + QRD$rotations
    D <- aperm(array(apply(QRD$R,c(1,2),Aconj),dim=c(d,m,n)),c(3,2,1))
    UH <- QRD$QH
    QRD <- AQRBC(A=D,beta=beta,epsilon=epsilon,usenorm=usenorm,QHinit=VH)
    qrds <- qrds+(QRD$rotations>0)
    sweeps <- sweeps + QRD$sweeps
    rotations <- rotations + QRD$rotations
    D <- aperm(array(apply(QRD$R,c(1,2),Aconj),dim=c(d,n,m)),c(3,2,1))
    VH <- QRD$QH
  }
  # Now convert U^H and V^H to U and V for output
  UH <- aperm(array(apply(UH,c(1,2),Aconj),c(d,m,m)),c(3,2,1))
  VH <- aperm(array(apply(VH,c(1,2),Aconj),c(d,n,n)),c(3,2,1))
  return(list(U=UH,D=D,V=VH,qrds=qrds,sweeps=sweeps,rotations=rotations))
}

## Utility functions used by the two main functions

L2norm <- function(x)sqrt(sum(x^2))
LInfinitynorm <- function(x)max(abs(x))

# Cl_4_1_binary_basis is a binary representation of the basis of Cl(4,1).
# Multiplication in Cl(4,1) is then a twisting of the xor operation on pairs of binary vectors
Cl_4_1_binary_basis <- matrix(FALSE,nrow=32,ncol=5)
Cl_4_1_binary_basis[cbind((1:5)+1,1:5)] <- TRUE
Cl_4_1_binary_basis <- cbind(
  c(FALSE,
    TRUE,FALSE,FALSE,FALSE,FALSE,
    rep(TRUE,4),rep(FALSE,6),
    rep(TRUE,6),rep(FALSE,4),
    rep(TRUE,4),FALSE,
    TRUE),
  c(FALSE,
    FALSE,TRUE,FALSE,FALSE,FALSE,
    TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,
    TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,
    TRUE,TRUE,TRUE,FALSE,TRUE,
    TRUE),
  c(FALSE,
    FALSE,FALSE,TRUE,FALSE,FALSE,
    FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,
    TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,
    TRUE,TRUE,FALSE,TRUE,TRUE,
    TRUE),
  c(FALSE,
    FALSE,FALSE,FALSE,TRUE,FALSE,
    FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,
    FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,
    TRUE,FALSE,TRUE,TRUE,TRUE,
    TRUE),
  c(FALSE,
    rep(FALSE,4),TRUE,
    FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,
    FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,
    FALSE,TRUE,TRUE,TRUE,TRUE,
    TRUE)
  )
# Some preliminary work to write out an explicit basis for the real matrix representation of Cl(4,1)
Cl_4_1_basis <- list(diag(32))
for(k in 2:6) {
  ei <- matrix(0,32,32)
  for(j in 1:32) {
    indktimesj <- which(colSums(t(Cl_4_1_binary_basis)==xor(Cl_4_1_binary_basis[k,],Cl_4_1_binary_basis[j,]))==5)
    signktimesj <- (-2*(sum(Cl_4_1_binary_basis[j,(0:(k-2))[-1]]) %% 2)+1)*(2*(k!=6||Cl_4_1_binary_basis[j,k-1]==0)-1)
    ei[indktimesj,j] <- signktimesj
  }
  Cl_4_1_basis[[k]] <- ei
}
for(k in 7:32) {
  prodof <- Cl_4_1_basis[which(Cl_4_1_binary_basis[k,])+1]
  Cl_4_1_basis[[k]] <- Reduce(f="%*%",prodof)
}

Cl_4_1_matmult <- function(A,B) {
  # A is an m x n x 32 array representing an m x n matrix in Cl(4,1)
  # B is an n x p x 32 array representing an n x p matrix in Cl(4,1)
  # This function outputs the matrix product of A and B as an m x p x 32 array.
  m <- dim(A)[1]
  n <- dim(A)[2]
  d <- dim(A)[3]
  n2 <- dim(B)[1]
  p <- dim(B)[2]
  d2 <- dim(B)[3]
  if(d!=32||d2!=32) { #||n!=n2) {
    stop("A should be a m x n x 32 array, and B should be an n x p x 32 array.")
  }
  C <- array(0,dim=c(m,p,d))
  for(k1 in 1:d) {
    for(k2 in 1:d) {
      C <- C + ((matrix(A[,,k1],nrow=m,ncol=n) %*% matrix(B[,,k2],nrow=n,ncol=p)) %o% Cl_4_1_basis[[k1]][,k2])
    }
  }
  return(C)
}

Cl_4_1_representation <- function(x) {
  # If x is a vector of length 32 representing an element of Cl(4,1),
  # then this function computes the representation of x as a 4x4 complex matrix.
  # if x is a 4x4 complex matrix, then this function computes
  # the coresponding element of Cl(4,1) as a vector of length 32
  # This function also accepts as input an m by n by 32 array
  # or a 4m by 4n complex matrix

  # Equations taken from Tian, Yongge. "Universal Similarity Factorization Equalities over Real Clifford Algebras." Advances in Applied Clifford Algebras 8, no. 2 (1998): 365-402. doi:10.1007/BF03043105.
  # (2.2.2)
  P_2_0 <- array(0,dim=c(2,2,32))
  P_2_0[1,1,1:2] <- 0.5
  P_2_0[2,2,1] <- 0.5
  P_2_0[2,2,2] <- -0.5
  P_2_0[1,2,3] <- 0.5
  P_2_0[1,2,3] <- 0.5
  P_2_0[1,2,7] <- -0.5
  P_2_0[2,1,c(3,7)] <- 0.5
  # (2.4.4)
  P_3_1 <- array(0,dim=c(4,4,32))
  tmp <- array(0,dim=c(2,2,32))
  tmp[1,1,c(1,19)] <- 0.5
  tmp[2,2,c(1,19)] <- 0.5
  P_3_1[1:2,1:2,] <- Cl_4_1_matmult(tmp,P_2_0)
  #P_3_1[3:4,3:4,] <- P_3_1[1:2,1:2,]
  tmp[1,1,19] <- -0.5
  tmp[2,2,19] <- -0.5
  P_3_1[3:4,3:4,] <- Cl_4_1_matmult(tmp,P_2_0)
  tmp <- array(0,dim=c(2,2,32))
  tmp[1,1,17] <- 0.5
  tmp[2,2,17] <- 0.5
  tmp[1,1,15] <- -0.5
  tmp[2,2,15] <- -0.5
  P_3_1[1:2,3:4,] <- Cl_4_1_matmult(tmp,P_2_0)
  tmp <- array(0,dim=c(2,2,32))
  tmp[1,1,c(17,15)] <- -0.5
  tmp[2,2,c(17,15)] <- -0.5
  P_3_1[3:4,1:2,] <- Cl_4_1_matmult(tmp,P_2_0)
  # Cl_4_1_matmult(array(1,dim=c(1,1,32)),array(4,dim=c(2,2,3))) # This should be the identity
  
  
  realhalf <- function(x) {
    x[Cl_4_1_binary_basis[,4]] <- 0
    return(x)
  }
  e_all <- array(c(rep(0,31),1),dim=c(1,1,32))
  if(length(x)==32) {
    # Convert from Cl(4,1) to C^{4 x 4}
    x <- realhalf(x) - 1i*realhalf(Cl_4_1_matmult(array(x,dim=c(1,1,32)),e_all))
    y <- array(0,dim=c(4,4,32))
    for(k in 1:4) {y[k,k,] <- x}
    y <- Cl_4_1_matmult(Cl_4_1_matmult(P_3_1,y),P_3_1)
    return(y[,,1])
  } else if(is.complex(x)&&length(dim(x))==2&&all(dim(x)==4)) {
    # Convert from C^{4 x 4} to Cl(4,1) 
    y <- array(0,dim=c(4,4,32))
    y[,,1] <- Re(x)
    y[,,32] <- Im(x)
    y <- Cl_4_1_matmult(Cl_4_1_matmult(P_3_1,y),P_3_1)
    return(y[1,1,])
  } else if(length(dim(x))==3&&dim(x)[3]==32) {
    y <- apply(x,c(1,2),Cl_4_1_representation)
    y <- array(y,dim=c(4,4,dim(x)[1:2]))
    y <- aperm(y,c(1,3,2,4))
    y <- array(y,dim=dim(x)[1:2]*4)    
  } else if(is.complex(x)&&length(dim(x))==2) {
    m4 <- dim(x)[1]
    n4 <- dim(x)[2]
    m <- m4/4
    n <- n4/4
    x <- array(x,dim=c(4,m,4,n))
    x <- aperm(x,c(1,3,2,4))
    y <- apply(x,c(3,4),Cl_4_1_representation)
    y <- aperm(y,c(2,3,1))
    return(y)
  } else {
    stop("Wrong argument type")
  }
}




getbeta <- function(algebra=c("real","complex","splitcomplex","quaternion","cyclicpolynomial","cl(4,1)")) {
  algebra <- tolower(algebra)
  algebra <- match.arg(algebra)
  if(algebra=="real") {
    beta <- function(b) {
      return(1)
    }
  } else if(algebra=="cyclicpolynomial") {
    beta <- function(b) {
      d <- length(b)
      i <- which.max(abs(b))
      out <- matrix(0,nrow=d,ncol=d)
      out[cbind(1:d,((1:d - i) %% d)+1L)] <- 1
      return(out)
    }
  } else if(algebra=="complex") {
    beta <- function(b) {
      if(all(b==0)) {return(diag(2))}
      out <- b/sqrt(sum(b^2))
      return(cbind(out[1:2],c(-out[2],out[1])))
    }
  } else if(algebra=="quaternion") {
    beta <- function(b) {
      if(all(b==0)) {return(diag(4))}
      out <- b/sqrt(sum(b^2))
      return(cbind(out[1:4],c(-out[2],out[1],out[4],-out[3]),c(-out[3],-out[4],out[1:2]),c(-out[4],out[3],-out[2],out[1])))
    }
  } else if(algebra=="splitcomplex") {
    beta <- function(b) {
      if(abs(b[1])<abs(b[2])) {
        return(matrix(c(0,1,1,0),nrow=2,ncol=2))
      } else {
        return(diag(2))
      }
    }
  } else if(algebra=="cl(4,1)") {
    beta <- function(b) {
      i <- which.max(abs(b))
      return(Cl_4_1_basis[[i]])
    }
  }
  return(beta)
}
getAconj <- function(algebra=c("real","complex","splitcomplex","quaternion","cyclicpolynomial","cl(4,1)")) {
  algebra <- tolower(algebra)
  algebra <- match.arg(algebra)
  if(algebra=="real"||algebra=="splitcomplex") {
    Aconj <- identity
  } else if(algebra=="cyclicpolynomial") {
    Aconj <- function(a) {
      d <- length(a)
      return(a[(-(0:(d-1)) %% d)+1L])
    }
  } else if(algebra=="complex") {
    Aconj <- function(a) {
      return(a*c(1,-1))
    }
  } else if(algebra=="quaternion") {
    Aconj <- function(a) {
      return(a*c(1,-1,-1,-1))
    }
  } else if(algebra=="cl(4,1)") {
    Aconj <- function(a) {
      return(a*c(1,
                 1,1,1,1,-1,
                 ifelse(Cl_4_1_binary_basis[7:16,5],1,-1),
                 ifelse(Cl_4_1_binary_basis[17:26,5],1,-1),
                 ifelse(Cl_4_1_binary_basis[27:31,5],-1,1),
                 -1
                 ))
    }
  }
  return(Aconj)
}
Ahermconj <- function(A,algebra=NULL,Aconj=getAconj(algebra)) {
  return(aperm(array(apply(A,c(1,2),Aconj),dim=dim(A)[c(3,1,2)]),c(3,2,1)))
}





#### Some examples of use, including the plots created for the article:
## Paul Ginzberg and Christiana Mavroyiakoumou "The QRD and SVD of matrices over a real algebra", Linear Algebra and ist Applications, 2016.


oldpar <- par()
A <- array(rnorm(3*5*100),dim=c(3,5,100))
timing <- proc.time()
QR <- AQRBC(A,"CyclicPolynomial",epsilon=1e-10)
proc.time()-timing
QR$R

A <- array(rnorm(3*4*10),dim=c(3,4,10))
SVD <- ASVD(A,"CyclicPolynomial",epsilon=1e-10)

# A <- array(rnorm(2*2*2),dim=c(2,2,2))
# epsilons <- 2^(-(1:100))
# SVDs <- lapply(epsilons,function(x)ASVD(A,"splitcomplex",epsilon=x))
# counts <- sapply(SVDs,function(x)unlist(x[4:6]))
# matplot(epsilons,t(counts),log="xy",type="l")
# plot(epsilons,counts[3,]/counts[1,],type="l",log="x")
set.seed(1)
A <- array(rnorm(3*2*32),dim=c(3,2,32)) ###3x2
epsilons <- 2^(-(0:60))
SVDs <- lapply(epsilons,function(x)ASVD(A,"Cl(4,1)",epsilon=x))
QRs <- lapply(epsilons,function(x)AQRBC(A,"Cl(4,1)",epsilon=x))
counts <- sapply(SVDs,function(x)unlist(x[4:6]))
countsQR <- sapply(QRs,function(x)unlist(x[3:4]))
#matplot(epsilons,t(counts),log="xy",type="l",xlab=expression(epsilon),ylab="",col=)
#plot(epsilons,counts[3,]/counts[1,],type="l",log="x")

# Redo calculations using 4x4 complex representation of Cl(4,1)
af <- 16
rA <- Cl_4_1_representation(A)
rQR <- AQRBC(array(c(Re(rA),Im(rA)),dim=c(dim(rA),2)),algebra="complex",epsilon=1e-16)
rQRt <- rQR
rQRt$QH <- Cl_4_1_representation(rQR$QH[,,1]+1i*rQR$QH[,,2])
rQRt$R <- Cl_4_1_representation(rQR$R[,,1]+1i*rQR$R[,,2])
rSVD <- ASVD(array(c(Re(rA),Im(rA)),dim=c(dim(rA),2)),algebra="complex",epsilon=1e-16)
rSVDt <- rSVD
rSVDt$U <- Cl_4_1_representation(rSVD$U[,,1]+1i*rSVD$U[,,2])
rSVDt$D <- Cl_4_1_representation(rSVD$D[,,1]+1i*rSVD$D[,,2])
rSVDt$V <- Cl_4_1_representation(rSVD$V[,,1]+1i*rSVD$V[,,2])
rSVDs <- lapply(epsilons/af,function(x)ASVD(array(c(Re(rA),Im(rA)),dim=c(dim(rA),2)),algebra="complex",epsilon=x))
rcounts <- sapply(rSVDs,function(x)unlist(x[4:6]))
rQRs <- lapply(epsilons/af,function(x)AQRBC(array(c(Re(rA),Im(rA)),dim=c(dim(rA),2)),algebra="complex",epsilon=x))
rcountsQR <- sapply(rQRs,function(x)unlist(x[3:4]))

allvals <- c(counts["rotations",],rcounts["rotations",],rcountsQR["rotations",],countsQR["rotations",])
width <- 3
height <- 3
# pdf("plotrotations.pdf",width=3,height=3)
# plot(epsilons,counts["rotations",],pch=15,lty="solid",
#      type="l",log="xy",xlab=expression(paste(epsilon,' (log scale)')),ylab="Rotations (log scale)",lwd=3,
#      ylim=c(min(allvals),max(allvals)))
# lines(epsilons,countsQR["rotations",],type="l",col="black",lty="dotted",lwd=3,pch=4)
# lines(epsilons,rcounts["rotations",],type="l",col="red",lty="dashed",lwd=3,pch=15)
# lines(epsilons,rcountsQR["rotations",],type="l",col="red",lty="dotdash",lwd=3,pch=4)
# #legend("bottomleft",
# #       c("Using Algorithm 1",expression(paste("using the ",C^{4 %*% 4}," representation")),"QRD","SVD"),
# #       lty=c(1,2,NA,NA),lwd=2,pch=c(NA,NA,4,15),col=c("black","red","grey","grey"))
# dev.off()
pdf("plotrotationsSVD.pdf",width=3,height=3)
ylim <- range(c(counts["rotations",],rcounts["rotations",]/af))
xlim <- range(epsilons)
plot(epsilons,counts["rotations",],pch=15,lty="solid",
     type="l",log="xy",xlab=expression(paste(epsilon,' (log scale)')),ylab="Givens Rotations (log scale)",lwd=2,
     ylim=ylim,xlim=xlim,bty="l",xaxt="n",yaxt="n")
at <- c(1e-18,1e-14,1e-10,1e-6,1e-2)
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(1,at=at,labels=labels)
at <- at/10
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(3,at=at*af,labels=labels,col="red",col.axis="red")
at <- c(1e2,1e3,1e4,1e5)
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(2,at=at,labels=labels)
at <- at*10
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(4,at=at/af,labels=labels,col="red",col.axis="red")
lines(epsilons,rcounts["rotations",]/af,type="l",col="red",lty="solid",lwd=2,pch=15)
text(rev(epsilons)[3],rev(counts["rotations",])[3],"1",pos=1)
text(rev(epsilons)[3],rev(rcounts["rotations",])[3]/af,"2",col="red",pos=3)
dev.off()
pdf("plotrotationsQR.pdf",width=3,height=3)
ylim <- range(c(countsQR["rotations",],rcountsQR["rotations",]/af))
xlim <- range(epsilons)
plot(epsilons,countsQR["rotations",],pch=15,lty="solid",
     type="l",log="xy",xlab=expression(paste(epsilon,' (log scale)')),ylab="Givens Rotations (log scale)",lwd=2,
     ylim=ylim,xlim=xlim,bty="l",xaxt="n",yaxt="n")
at <- c(1e-18,1e-14,1e-10,1e-6,1e-2)
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(1,at=at,labels=labels)
at <- at/10
labels <- parse(text=paste0("10^",log10(at)))
labels[c(2,4)] <- ""
axis(3,at=at*af,labels=labels,col="red",col.axis="red")
at <- c(1e1,1e2,1e3)
labels <- parse(text=paste0("10^",log10(at)))
labels[2] <- ""
axis(2,at=at,labels=labels)
at <- at*10
labels <- parse(text=paste0("10^",log10(at)))
labels[2] <- ""
axis(4,at=at/af,labels=labels,col="red",col.axis="red")
lines(epsilons,rcountsQR["rotations",]/af,type="l",col="red",lty="solid",lwd=2,pch=15)
text(rev(epsilons)[3],rev(countsQR["rotations",])[3],"1",pos=1)
text(rev(epsilons)[3],rev(rcountsQR["rotations",])[3]/af,"2",col="red",pos=3)
dev.off()

#A <- array(rbinom(3*2*32,1,1/16),dim=c(3,2,32))
set.seed(1)
A <- array(rnorm(3*2*32),dim=c(3,2,32))
QR <- AQRBC(A,"cl(4,1)",epsilon=1e-16)
SVD <- ASVD(A,"cl(4,1)",epsilon=1e-16)
str(QR)
str(SVD)
Aplot(A)
Aplot(QR$R,ylab="R")
Aplot(QR$QH,ylab="QH")



algebra <- "cl(4,1)"
Q <- Ahermconj(QR$QH,algebra="cl(4,1)")
rQt <- Ahermconj(rQRt$QH,algebra="cl(4,1)")
Arec <- Cl_4_1_matmult(Q,QR$R)
#A2 <- A; A2[,,c(-1,-2)] <- 0; A2[,,2] <- A2[,,2]/10; QR2 <- AQRBC(A2,"cl(4,1)",epsilon=1e-16)
#Aplot(Cl_4_1_matmult(QR2$QH,A2))
#Aplot(Cl_4_1_matmult(QR2$QH,A2)-QR2$R)
#sum((Cl_4_1_matmult(QR2$QH,A2)-QR2$R)^2)
sum((Arec-A)^2)
sum((Cl_4_1_matmult(QR$QH,A)-QR$R)^2)
sum((Cl_4_1_matmult(rQRt$QH,A)-rQRt$R)^2)
sum(Cl_4_1_matmult(Q,QR$QH)^2)

Rthreshfun <- function(R) {
  for(j in 1:dim(R)[2]) {
    for(i in 1:dim(R)[1]) {
      if(i>j){
        R[i,j,] <- 0
      }
    }
  }
  return(R)
}
R <- Rthreshfun(QR$R)
sqrt(sum((Cl_4_1_matmult(Q,R)-A)^2))
rRt <- Rthreshfun(rQRt$R)
sqrt(sum((Cl_4_1_matmult(rQt,rRt)-A)^2))
max(abs(Cl_4_1_matmult(Q,R)-A))
max(abs(Cl_4_1_matmult(rQt,rRt)-A))
Dthreshfun <- function(D) {
  for(j in 1:dim(D)[2]) {
    for(i in 1:dim(D)[1]) {
      if(i!=j){
        D[i,j,] <- 0
      }
    }
  }
  return(D)
}
D <- Dthreshfun(SVD$D)
sqrt(sum((Cl_4_1_matmult(SVD$U,Cl_4_1_matmult(D,Ahermconj(SVD$V,algebra="cl(4,1)")))-A)^2))
rDt <- Dthreshfun(rSVDt$D)
sqrt(sum((Cl_4_1_matmult(rSVDt$U,Cl_4_1_matmult(rDt,Ahermconj(rSVDt$V,algebra="cl(4,1)")))-A)^2))

VH <- aperm(array(apply(SVD$V,c(1,2),getAconj(algebra)),dim=dim(SVD$V)[c(3,1,2)]),c(3,2,1))
Cl_4_1_matmult(SVD$V,VH)
Arec <- Cl_4_1_matmult(Cl_4_1_matmult(SVD$U,SVD$D),VH)
sum((Arec-A)^2)


pdf("plotA.pdf",width=width,height=height)
Aplot(A,ylab="A")
dev.off()
pdf("plotQ.pdf",width=width,height=height)
Aplot(Q,ylab="Q")
dev.off()
pdf("plotR.pdf",width=width,height=height)
Aplot(QR$R,ylab="R")
dev.off()
pdf("plotU.pdf",width=width,height=height)
Aplot(SVD$U,ylab="U")
dev.off()
pdf("plotD.pdf",width=width,height=height)
Aplot(SVD$D,ylab="D")
dev.off()
pdf("plotV.pdf",width=width,height=height)
Aplot(SVD$V,ylab="V")
dev.off()
pdf("plotQrt.pdf",width=width,height=height)
Aplot(rQt,ylab="Q")
dev.off()
pdf("plotRrt.pdf",width=width,height=height)
Aplot(rQRt$R,ylab="R")
dev.off()
pdf("plotUrt.pdf",width=width,height=height)
Aplot(rSVDt$U,ylab="U")
dev.off()
pdf("plotDrt.pdf",width=width,height=height)
Aplot(rSVDt$D,ylab="D")
dev.off()
pdf("plotVrt.pdf",width=width,height=height)
Aplot(rSVDt$V,ylab="V")
dev.off()

## Bits and pieces, Ignore.
# shell.exec(paste0(getwd(),"/plotA.pdf"))
# pdf("plotQ.pdf",paper="a4")
# Aplot(QR$Q)
# dev.off()
# pdf("plotD.pdf",paper="a4")
# Aplot(SVD$D)
# dev.off()
# shell.exec(paste0(getwd(),"/plotD.pdf"))
# pdf("plotU.pdf",paper="a4")
# Aplot(SVD$U)
# dev.off()
# shell.exec(paste0(getwd(),"/plotU.pdf"))
# pdf("plotR.pdf",paper="a4")
# Aplot(QR$R)
# dev.off()
# shell.exec(paste0(getwd(),"/plotR.pdf"))


# Z <- matrix(complex(4,c(10,2,3,4),c(5,6,6,6)),ncol=2)
# i <- complex(1,0,1)
# ct <- function(x) Conj(t(x))
# S <- Z%*%ct(Z)
# v <- function(x) rbind(Re(x),Im(x))
# vinv <- function(x) x[1:(nrow(x)/2),]+1i*x[(nrow(x)/2+1):nrow(x),]
# tilde <- function(x) cbind(v(x),v(x*1i))
# svdinv <- function(x) x$u%*%diag(x$d)%*%ct(x$v)
# vZ <- cbind()
# 
# 
# t1 <- cbind(c(1,0),c(0,1))
# t2 <- cbind(c(1,0),c(0,-1))
# t3 <- cbind(c(0,0),c(1,0))
# t3 %*% t3
