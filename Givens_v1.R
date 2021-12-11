diag(1,3,3,FALSE)

x=matrix(c(1:9), 3, 3)

x

x+t(x)

asin(1)
asin(0)

pi/2

acos(1)

z=x+t(x)
z
z[2,1]


find_theta=function(A_matrix, i, j){
  A=A_matrix
  the_angle=(1/2)*atan((2*A[j,i])/(A[i,i]-A[j,j]))
  return(the_angle)
}

find_theta(z, 2, 1)

(1/2)*atan(2*z[2,1]/(z[2,2]-z[1,1]))


#this does not work
matrix_R=function(A_matrix, i, j){
  A=A_matrix
  B=diag(ncol(A))
  t=find_theta(A, i, j)
  B[i,i]=sin(t)
  B[i,j]=cos(t)
  B[j,i]=-cos(t)
  B[j,j]=sin(t)
  
  return(B)
}



R=matrix_R(z, 2, 1)
R%*%z%*%t(R)

h=matrix(c(1:4), 2, 2)
h=h+t(h)

RR=matrix_R(h,1,2)
RR%*%h%*%t(RR)
RR%*%t(RR)


h

#this most likely works.
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


R2=matrix_R2(h, 1, 2)
R2%*%t(R2)
R2%*%h%*%t(R2)

R3=matrix_R2(z,2,1)
R3%*%z%*%t(R3)


e=matrix(c(3:11),3,3)
e=e+t(e)
R4=matrix_R2(e, 1, 3)
R4%*%e%*%t(R4)

u=matrix(c(3:18),4,4)
u=u+t(u)
R5=matrix_R2(u, 2, 4)
R5%*%u%*%t(R5)

?eigen

eigen(u, symmetric = TRUE)



eigen_fun_v1 = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in 2:d){
        if(i<j){
          # test whether abs(X[i,j]) is too large
          if(abs(X[i,j])>tol){
            # if so, get the Givens Rotation matrix R
            R=matrix_R2(X, i, j)
            # and update X and U
            X=R%*%X%*%t(R)
            U=U%*%R
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

s=matrix(c(1:9), 3, 3)
s=s+t(s)


d1=eigen_fun_v1(s)
d2=eigen(s, symmetric = TRUE)

sort(d2$values)

d1$values

d1$vectors

d2$vectors
d2$values

###
s=matrix(c(1:4),2,2)
s=s+t(s)
ss=s

R6=matrix_R2(s, 1, 2)
s=R6%*%s%*%t(R6)

U=diag(2)
U=U%*%R6

s
eigen(ss)$value

R6=matrix_R2(s, 1, 2)
s=R6%*%s%*%t(R6)
s
eigen(ss)$value

U=U%*%R6

eigen(ss)$vectors

U

str(eigen(ss)$vectors)


####
s=matrix(c(1:4),2,2)
s=s+t(s)
ss=s

R6=matrix_R2(s, 1, 2)
s=R6%*%s%*%t(R6)

U=diag(2)
U=R6%*%U

s
eigen(ss)$value

eigen(ss)$vectors
####


####
s=matrix(c(1:9), 3, 3)
s=s+t(s)
ss=s


eigen_fun_v1(s)

eigen(ss, symmetric = TRUE)
####


eigen_fun_v2 = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in 2:d){
        if(i!=j){
          # test whether abs(X[i,j]) is too large
          if(abs(X[i,j])>tol){
            # if so, get the Givens Rotation matrix R
            R=matrix_R2(X, i, j)
            # and update X and U
            X=R%*%X%*%t(R)
            U=U%*%R
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


s=matrix(c(1:9), 3, 3)
s=s+t(s)
ss=s

eigen_fun_v2(s)


s=matrix(c(1:9), 3, 3)
s=s+t(s)
eigen_fun_v1(s)

eigen(ss)$vectors

#this seems to work
eigen_fun_v3 = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in 2:d){
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

s=matrix(c(1:9), 3, 3)
s=s+t(s)
ss=s

eigen_fun_v3(s)
eigen(ss)$vectors



####this seems to work
####maybe not
eigen_fun_v4 = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in 2:d){
        if(i<j){
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

s=matrix(c(1:9), 3, 3)
s=s+t(s)
ss=s

eigen_fun_v4(s)
eigen(ss)$vectors


s=matrix(c(2:17), 4, 4)
s=s+t(s)
ss=s

eigen_fun_v3(s)
eigen(ss)$vectors

eigen_fun_v3(s)$values
eigen(ss)$values



s=matrix(c(5:13), 3, 3)
s=s+t(s)
ss=s

eigen_fun_v3(s)
eigen(ss)$vectors


s=matrix(rnorm(9), nrow = 3)
s
s=s+t(s)
ss=s

d1=eigen_fun_v3(s)
d2=eigen(ss)


###this seems to work
library(testthat)
test_that("test for 3 by 3", {
  expect_equal(d1$values, sort(d2$values))
})

str(d1$values)
str(d2$values)


vect1=c(1,2,3)
unlist(vect1)

vect2=c(2,1,3)

library(testthat)
test_that("test", {
  expect_equal(1, 1+10^-8)
})

####index return
sv=sort(d2$values, index.return=TRUE)
sv$ix

####this seems to work
####using index.return=TRUE in sort.
test_that("test for 3 by 3", {
  for(i in 1:3){
    expect_equal(d1$values[i], d2$values[sv$ix[i]])
  }
})

###does work
test_that("test for 3 by 3", {
  for(i in 1:3){
    expect_equal(abs(d1$vectors[,i]), abs(d2$vectors[,sv$ix[i]]))
  }
})

i=1
d1$vectors[,1]
d2$vectors[,1]
