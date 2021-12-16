#This is the answer for Problem 1 in the final.

setwd("~/Documents/Rstuff/final/final1")
set.seed(1234567890)

#Creating a list to store all the x's, where each x is 1000 X 500.
X=list()
for(i in 1:20){
  X[[i]]=matrix(rnorm(5*10^5), nrow = 1000)
}

#Creating a list to store all the y's, where each y is 1000 X 1.
Y=list()
for(i in 1:20){
  Y[[i]]=matrix(rnorm(1000), nrow = 1000)
}


#Defining the five functions
# direct inversion
r_squared_fun_1 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx)%*%xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# solve linear system
r_squared_fun_2 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx,xty)
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# invert using eigen system
r_squared_fun_3 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  eigen_xtx = eigen(xtx,symmetric=TRUE)
  u = eigen_xtx$vectors
  d_inv = diag(1/eigen_xtx$values)
  beta = u %*% d_inv %*% t(u) %*% xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# avoid inversion using svd
r_squared_fun_4 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# avoid inversion using svd
# ignore right singular vectors
r_squared_fun_5 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x,nv=0)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}


#Creating a vector of strings. Length=20. The i-th element is the name of the .out file
#from Rprof() for the i-th column of the 5 X 20 matrix.
output_filenames=c()
output_filenames[1]="f1_v4_1.out"
output_filenames[2]="f1_v4_2.out"
output_filenames[3]="f1_v4_3.out"
output_filenames[4]="f1_v4_4.out"
output_filenames[5]="f1_v4_5.out"
output_filenames[6]="f1_v4_6.out"
output_filenames[7]="f1_v4_7.out"
output_filenames[8]="f1_v4_8.out"
output_filenames[9]="f1_v4_9.out"
output_filenames[10]="f1_v4_10.out"
output_filenames[11]="f1_v4_11.out"
output_filenames[12]="f1_v4_12.out"
output_filenames[13]="f1_v4_13.out"
output_filenames[14]="f1_v4_14.out"
output_filenames[15]="f1_v4_15.out"
output_filenames[16]="f1_v4_16.out"
output_filenames[17]="f1_v4_17.out"
output_filenames[18]="f1_v4_18.out"
output_filenames[19]="f1_v4_19.out"
output_filenames[20]="f1_v4_20.out"


#giving the vector output_filenames a different name. 
#original name was too long to write
ofn=output_filenames


#Profiling for the i-th column for i in 1:20.
for(i in 1:20){
  Rprof(filename=ofn[i], interval=0.01, filter.callframe=TRUE)
  
  r_squared_fun_1(X[[i]], Y[[i]]) #Run the first function
  r_squared_fun_2(X[[i]], Y[[i]]) #Run the second function
  r_squared_fun_3(X[[i]], Y[[i]]) #Run the third function
  r_squared_fun_4(X[[i]], Y[[i]]) #Run the fourth function
  r_squared_fun_5(X[[i]], Y[[i]]) #Run the fifth function
  
  Rprof(NULL)
}
  
  

#Calling summaryRprof() for the 20 columns and storing the output in a list.
S=list()
T=list()
for(i in 1:20){
  S[[i]]=summaryRprof(ofn[i])
  T[[i]]=S[[i]]$by.total
}


#Creating a 5X20 matrix
matrix_5X20=matrix(rep(NA, 100), nrow = 5)


#creating a vector containing name of the 5 functions.
#to be used to find the time spent by each function for a given column in
#the 5X20 matrix
name_of_the_5_functions=c("\"r_squared_fun_1\"", "\"r_squared_fun_2\"", "\"r_squared_fun_3\"",
                          "\"r_squared_fun_4\"", "\"r_squared_fun_5\"")

#Renaming the vector name_of_the_5_functions. 
#Original name was too long to write. 
n5f=name_of_the_5_functions


#Filling up the 5X20 matrix
for(i in 1:5){
  for(j in 1:20){
    matrix_5X20[i,j]=T[[j]][n5f[i], ]$total.time
  }
}


#Setting up the rownames of the 5X20 matrix
row.names(matrix_5X20)=n5f



#Initialiazing a matrix to store the summary of each function
summary_of_functions=matrix(rep(NA, 30), nrow = 5)

#Setting up the rownames of the summary matrix
row.names(summary_of_functions)=n5f

#Setting up column names for the summary matrix
colnames(summary_of_functions)=c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")



#Filling up the summary matrix
for(i in 1:5){
  for(j in 1:6){
    summary_of_functions[i,j]=summary(matrix_5X20[n5f[i], ])[j]
  }
}

#This is the final answer
summary_of_functions 

#The second function is the fastest. 
#This is because solving a system of equations is faster than inverting a matrix 
#or finding singular value decomposition or finding eigenvalues.