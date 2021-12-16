#This is the answer for Problem 3 in the final.


#Finds the angle of rotation
find_theta=function(A_matrix, i, j){
  A=A_matrix
  the_angle=(1/2)*atan((2*A[j,i])/(A[i,i]-A[j,j]))
  return(the_angle)
}


#Finds the rotation matrix
matrix_R2=function(A_matrix, i, j){
  A=A_matrix
  B=diag(ncol(A))
  t=find_theta(A, i, j)
  B[i,i]=cos(t)
  B[i,j]=sin(t)
  B[j,i]=-sin(t)
  B[j,j]=cos(t)
  
  return(B)
}


#Function to compute eigenvalues and eigenvectors.
eigen_fun_v3 = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        if(i!=j){
          # test whether abs(X[i,j]) is too large
          if(abs(X[i,j])>tol){
            # if so, get the Givens Rotation matrix R
            R=matrix_R2(X, i, j)
            # and update X and U
            X=R%*%X%*%t(R)
            U=U%*%t(R)
            # track that in this sweep through we have done at least one update
            did_we_update = TRUE
          }
        }
        
      }
    }
    if(!did_we_update) break
  }
  values=diag(X)
  sorted_values = sort(values,index.return=TRUE)
  values = sorted_values$x
  vectors = U[,sorted_values$ix]
  return(list(values=values, vectors=vectors))
}



#creating a 2 x 2 matrix
set.seed(1234567890)
z=matrix(rnorm(4), nrow = 2)
X=z+t(z)
XX=X

d1=eigen_fun_v3(X)
d2=eigen(XX)

sv=sort(d2$values, index.return=TRUE)


#test eigenvalues for 2 x 2 matrix
library(testthat)
test_that("test eigenvalues for 2 x 2 matrix", {
  for(i in 1:length(X)){
    expect_equal(d1$values[i], d2$values[sv$ix[i]])
  }
})

#test eigenvectors for 2 x 2 matrix
test_that("test eigenvectors for 2x2 matrix", {
  for(i in 1:ncol(X)){
    if(sign(d1$vectors[,i][1])==sign(d2$vectors[,sv$ix[i]][1])){
      expect_equal(d1$vectors[,i], d2$vectors[,sv$ix[i]])
    }
    else{
      expect_equal(d1$vectors[,i], -d2$vectors[,sv$ix[i]])
    }
  }
})



#creating a 3 x 3 matrix
set.seed(1234567890)
z=matrix(rnorm(9), nrow = 3)
X=z+t(z)
XX=X

d1=eigen_fun_v3(X)
d2=eigen(XX)

sv=sort(d2$values, index.return=TRUE)


#test eigenvalues for 3 x 3 matrix
library(testthat)
test_that("test eigenvalues for 3 x 3 matrix", {
  for(i in 1:length(X)){
    expect_equal(d1$values[i], d2$values[sv$ix[i]])
  }
})

#test eigenvectors for 3 x 3 matrix
test_that("test eigenvectors for 3 x 3 matrix", {
  for(i in 1:ncol(X)){
    if(sign(d1$vectors[,i][1])==sign(d2$vectors[,sv$ix[i]][1])){
      expect_equal(d1$vectors[,i], d2$vectors[,sv$ix[i]])
    }
    else{
      expect_equal(d1$vectors[,i], -d2$vectors[,sv$ix[i]])
    }
  }
})



#creating a 4 x 4 matrix
set.seed(1234567890)
z=matrix(rnorm(16), nrow = 4)
X=z+t(z)
XX=X

d1=eigen_fun_v3(X)
d2=eigen(XX)

sv=sort(d2$values, index.return=TRUE)


#test eigenvalues for 4 x 4 matrix
library(testthat)
test_that("test eigenvalues for 4 x 4 matrix", {
  for(i in 1:length(X)){
    expect_equal(d1$values[i], d2$values[sv$ix[i]])
  }
})

#test eigenvectors for 4 x 4 matrix
test_that("test eigenvectors for 4 x 4 matrix", {
  for(i in 1:ncol(X)){
    if(sign(d1$vectors[,i][1])==sign(d2$vectors[,sv$ix[i]][1])){
      expect_equal(d1$vectors[,i], d2$vectors[,sv$ix[i]])
    }
    else{
      expect_equal(d1$vectors[,i], -d2$vectors[,sv$ix[i]])
    }
  }
})



#creating a 5 x 5 matrix
set.seed(1234567890)
z=matrix(rnorm(25), nrow = 5)
X=z+t(z)
XX=X

d1=eigen_fun_v3(X)
d2=eigen(XX)

sv=sort(d2$values, index.return=TRUE)


#test eigenvalues for 5 x 5 matrix
library(testthat)
test_that("test eigenvalues for 5 x 5 matrix", {
  for(i in 1:length(X)){
    expect_equal(d1$values[i], d2$values[sv$ix[i]])
  }
})

#test eigenvectors for 5 x 5 matrix
test_that("test eigenvectors for 5 x 5 matrix", {
  for(i in 1:ncol(X)){
    if(sign(d1$vectors[,i][1])==sign(d2$vectors[,sv$ix[i]][1])){
      expect_equal(d1$vectors[,i], d2$vectors[,sv$ix[i]])
    }
    else{
      expect_equal(d1$vectors[,i], -d2$vectors[,sv$ix[i]])
    }
  }
})
