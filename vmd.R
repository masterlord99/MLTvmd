library(vmd)
library(magrittr)
data = readRDS("../ProjectX/data/data/ETH/models_data.RDS")
ts = data$ETH$short_real$Close %>% tail(1000000) %>% na.omit()
ts = ts[c(T,rep(F,59))] %>% diff
y =ts[-(1)]
ts = ts[1:length(y)]
# 
holder = matrix(NA,length(ts),ncol = 5)
sig= vmd(ts[1:12999])
sig$calculate()
a = sig$as.data.frame()
holder[1:12999,] = as.matrix(a[,2:6])
i = 13000
while(i <= nrow(holder)){
  sig= vmd(ts[(i-12999):i])
  sig$calculate()
  a = sig$as.data.frame()
  holder[i,] = as.matrix(a[13000,2:6])
  i = i + 1
  print(i)
}


ind = complete.cases(holder)
holder = holder[ind,]
y = y[ind]

# b = data$ETH$short %>%  tail(100000)
# b = b[c(T,F,F,F,F),]
# b = b[2:19999,]

library(randomForest)

y = sign(y)
y[y==-1] = 0
table(y)

model = randomForest(x=holder[1:13000,],y = as.factor(y[1:13000]),ntree = 100)
model
a = predict(model)
sum(a == y[1:13000]) / length(a)

test = holder[13001:nrow(holder),]
test_y_real = y[13001:nrow(holder)]

test_y = predict(model,newdata = test,type = "prob")[,2]

plot(test_y_real[1:100],col="red")
lines(test_y[1:100],type="l")

message("train acc: ",round(sum(a == y[1:13000]) / length(a),3))
message("test acc: ",round(sum(round(test_y)==test_y_real ) / length(test_y_real),3))

sig$calculate



tau = 0
alpha = 2000
K = 3
DC = T
init = 0
tol = 1e-06
N = 500
orderModes =T

#Flip aswell as mirror?
flip       = FALSE #Mirror Only

# System Variables
eps        = .Machine$double.eps  #Smallest positive floating-point number

# Period and sampling frequency of input signal
lenOrg       = length(signal)
fs           = 1/lenOrg

# Extend the signal by mirroring
hw           = floor(lenOrg/2)                    #The Halfwidth
lhs          = rev(head(signal,0      + hw))      #First Half, Reversed
if(flip) lhs = tail(lhs,1) - c(lhs - tail(lhs,1)) #Flipped
rhs          = rev(tail(signal,lenOrg - hw))      #Last  Half, Reversed
if(flip) rhs = head(rhs,1) - c(rhs - head(rhs,1)) #Flipped
signalMir    = c(lhs,signal,rhs);                 #Mirrored Signal

# Time Domain 0 to T (of mirrored signal)
lenMir     = length(signalMir)  ##NB: Previously 'T' in original code, but T is reserved in R.
t          = seq_len(lenMir)/lenMir

# Spectral Domain discretization
freqs      = t -0.5 -(1.0/lenMir)

# For future generalizations: individual alpha for each mode
Alpha      = rep(alpha,K)

# Construct and center f_hat
f_hat      = private$fftshift(fft(signalMir))
f_hat_plus = f_hat
f_hat_plus[1:floor(lenMir/2)] = 0;

# Matrix keeping track of every iterant, could be discarded for mem
u_hat_plus = array(0,c(N,lenMir,K));

# Initialization of omega_k
omega_plus = array(0,c(N,K));
if(init == 1){
  omega_plus[1,] = (0.5/K)*((1:K) - 1)
}else if(init == 2){
  omega_plus[1,] = sort(exp(log(fs) + (log(0.5)-log(fs))*runif(K)));
}

# If DC mode imposed, set its omega to 0
if(DC)
  omega_plus[1,1] = 0

# Start with empty dual variables
lambda_hat = array(0,c(N,lenMir))

# Other inits
ix     = (floor(lenMir/2)+1):lenMir
uDiff  = Inf #update step
n      = 1   #loop counter
sum_uk = 0   #accumulator

# Main loop for iterative updates
while(uDiff > tol & n < N){
  
  # In the original matlab code, [A] The first mode is handled initially, and then [B] the subsequent
  # modes are looped (ie from 2:K), The following is a simplification, seeing as [A] and [B] largely
  # use the same code
  for(k in 1:K){
    
    # Accumulator
    sum_uk = u_hat_plus[`if`(k==1,n,n+1),,`if`(k==1,K,k-1)] + sum_uk - u_hat_plus[n,,k]
    
    # Mode spectrum
    u_hat_plus[n+1,,k] = (f_hat_plus - sum_uk - lambda_hat[n,]/2) / (1 + Alpha[k]*(freqs - omega_plus[n,k])^2)
    
    # Center frequencies
    if(!DC || k > 1)
      omega_plus[n+1,k] = (freqs[ix] %*% (abs(u_hat_plus[n+1,ix,k])^2)) / sum( abs(u_hat_plus[n+1,ix,k])^2 )
  }
  
  # Dual ascent
  lambda_hat[n+1,] = lambda_hat[n,] + tau*(rowSums(u_hat_plus[n+1,,]) - f_hat_plus)
  
  # Loop Counter
  n = n + 1
  
  # Converged Yet?
  uDiff = sapply(1:K,function(i){
    a = u_hat_plus[n,,i] - u_hat_plus[n-1,,i]
    b = Conj(a)
    (1/lenMir)*(a %*% b)
  })
  uDiff = abs(eps + sum(uDiff))
  
  #Reporting
  if(n > 0 && n %% 10 == 0)
    writeLines(sprintf("Iteration: %s, Diff: %.4g",n,uDiff))
  
  #Has it exploded?
  if(is.na(uDiff))
    stop("Problem converging, check parameters",call.=FALSE)
}

# Postprocessing and cleanup
N = min(N,n)
omega = omega_plus[1:N,]


#TUKAJ KONEC

# Signal reconstruction
u_hat            = array(0,c(lenMir,K));
u_hat[ix,]       =      u_hat_plus[N,ix,]
u_hat[ix[1]:2,]  = Conj(u_hat_plus[N,ix,])
u_hat[1,]        = Conj(u_hat[lenMir,]);

#NB: This Differs from original (it is transpose)
#    intentionally want consistency in having modes in columns
u = array(0,c(lenMir,K))
u[,1:K] = Reduce('cbind',lapply(1:K,function(k){
  Re( private$fftinv( private$fftshift(u_hat[,k],inverse = TRUE) ) )
}))

# Remove Mirror Part/s
ixRow = seq_len(length(signal)) + length(lhs)
u     = u[ixRow,]
u_hat = u_hat[ixRow,]
freqs = freqs[ixRow]
f_hat = f_hat[ixRow]

# Recompute spectrum
u_hat[,1:K] = Reduce('cbind',lapply(1:K,function(k){
  private$fftshift(fft(u[,k]))
}))

#Determine the ordering
ixCol = `if`(orderModes,order,seq_along)(tail(omega,1))

#Store the Result
private$varResult = list(signal = self$Signal,
                         u      = u[,    ixCol,drop=FALSE],
                         u_hat  = u_hat[,ixCol,drop=FALSE],
                         omega  = omega[,ixCol,drop=FALSE],
                         freqs  = freqs,
                         f_hat  = f_hat)

#Done