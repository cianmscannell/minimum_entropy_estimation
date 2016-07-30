## Initialisation
N = 20						
h = 0.1
X1 = numeric(N)
for(i in 1:N){
	X1[i] = 0.5 + (i-1)*0.95/(N-1)
}

## Implementation of the criterion for minimisation ##
##------------------------------------------------------------##

## function which returns residuals for given theta
e <- function(theta){			
	return(Y1 - model(theta,X1))
}


## function to calculate kernel, which we assume 
## to be gaussian with mean 0 and standard deviation 1
kern <- function(u){			
	return(1/sqrt(2*pi) * exp( - u^2 / 2))
}


# We created a function to return the values
# at which we evaluate the kernel, for each
# index i and j
# We do this for each side of our symmetrised estimator
a <- function(theta,ind1,ind2){	## component of kernel argument
	return((e(theta)[ind1] - e(theta)[ind2]) / h )
} 


b <- function(theta,ind1,ind2){	## component of kernel argument
	return((e(theta)[ind1] + e(theta)[ind2]) / h )
} 


# This function returns the f_n,h defined in our project
# This is 1/2Nh times the sum over j of the kernel evaluated at each
# value returned by a and b 
F <- function(theta,ind){		## function to return 1/2Nh * kernel
	s = 0
	for(i in 1:N){
		s = s + kern(a(theta,ind,i)) + kern(b(theta,ind,i)) 
	}
	s = 1/(2*N*h) * s
	return(s)
}


# This returns the actual criterion function to be
# minimised as defined in our report
J <- function(theta){		
	s1 = 0 
	for( j in 1:N){
		s1 = s1 + log(F(theta,j))
	}
	s1 = -1/N * s1
	return(s1)
}

## Example 1 ##

thbar = 1					## initial guess for theta
model <- function(theta,x){		## function to simulate model of interest
	return(exp(-theta*x))
}

Y1 = model(thbar,X) + rnorm(N,m=0,sd=0.1) ## first case where errors ~ N(0, sigma^2)

plot(X1, Y1, pch=20, main="data",col='indianred')	## visual assessment of data points
points(X1,model(X1,thbar),col='skyblue',pch=20)  ## model fit

## optimisation of minimum entropy criterion
th.start = 0.5		
opt.out <- optim(par=c(th.start),			
fn=J,method='Brent',lower =c(0.1),upper=c(2))
## tests whether values are reasonable
(opt.out$par)		

## Run a monte carlo simulation for entropy estimation
ans = numeric(1000)
for( i in 1:1000){
	Y1 = model(thbar,X) + rnorm(N,m=0,sd=0.1)
	opt.out <- optim(par=c(th.start),fn=J,method='Brent',lower =c(-0.5),upper=c(3))
	ans[i] = opt.out$par
}



## Least Squares Estimation
crit <- function(theta,x,y){
	# must have theta as 1st parameter and return a single value...
	return( sum( (y-model(theta,x))^2 ) )
}

## Run a monte carlo simulation for least squares
t = numeric(1000)
thinit = 0.005
for( i in 1:1000){
	Y1 = model(thbar,X1) + rnorm(N,m=0,sd=0.1)
	optim.out <- optim(par=c(thinit), fn=crit, x=X1, y=Y1,method='Brent',lower=c(0.1),upper=c(2))
	t[i] = optim.out$par
}


# Figure 1
plot(density(t),bty='n', ann='F',ylim=c(0,7))
lines(density(ans),lty=2)
legend("topright",'groups',c(expression(hat(theta)[ML]),expression(hat(theta)[e])),lty = c(1,2),col = c('black','black'),ncol=2)




#---------------------------------#
# Now try with uniform errors (Fig.2)
# Here we just repeat what we had already done.

Y2 = model(thbar,X1) + runif(N,-0.4,0.4)
e <- function(theta){			
	return(Y2 - model(theta,X1))
}

thu.start = 0.5
Y2 = model(thbar,X1) + runif(N,-0.4,0.4)
unif.out <- optim(par=c(thu.start),fn=J,method='Brent',lower =c(0.1),upper=c(2))
(unif.out$par)

# Monte Carlo
ansu = numeric(1000)
thu.start = 0.5
for( i in 1:1000){
	Y2= model(thbar,X1) + runif(N,-0.4,0.4)
	u.out <- optim(par=c(thu.start),fn=J,method='Brent',lower =c(-0.5),upper=c(3))
	ansu[i] = u.out$par
}

plot(density(ansu))

# LS compares to MLE under the wrong assumption that errors are normal

tr = numeric(1000)
thinit = 0.005
for( i in 1:1000){
	Y2 = model(thbar,X1) + runif(N,-0.4,0.4)
	o.out <- optim(par=c(thinit), fn=crit, x=X1, y=Y2,method='Brent',lower =c(-0.5),upper=c(3))
	tr[i] = o.out$par
}

# Figure 2
plot(density(ansu),bty='n', ann='F',ylim=c(0,4))
lines(density(tr),lty=2)
legend("topright",'groups',c(expression(hat(theta)[e]),expression(hat(theta)[LS])),lty = c(1,2),col = c('black','black'),ncol=2)

sd(ansu)
sd(tr)

#-----------------------------------------#

# Inverse transform method to generate a
# random sample from the laplace distribution

rlaplace <- function(N,sigma){
	u <- runif(N)
	uu<- runif(N)
	r <- numeric(N)
	for(i in 1:N){
		if(uu[i] < 0.5){
			r[i] = (sigma/(sqrt(2)) * log(2*u[i]))
		}
		else{
			r[i] = (- sigma/(sqrt(2)) * log(2*u[i]))
		}
	}
	return(r)
}

# we then added an outlier to our sample to 
# generate figure 3.

thl.start = 0.5
Y3 = model(thbar,X1) + rlaplace(N,0.15)
Y3[N] = 5

e <- function(theta){			
	return(Y3 - model(theta,X1))
}

lap.out <- optim(par=c(thl.start),fn=J,method='Brent',lower =c(0.1),upper=c(2))
(lap.out$par)

ansl = numeric(1000)
thl.start = 0.5
for( i in 1:1000){
	Y3 = model(thbar,X1) + rlaplace(N,0.15)
    Y3[N] = 5
	l.out <- 	optim(par=c(thl.start),fn=J,method='Brent',lower =c(-0.5),upper=c(3))
	ansl[i] = l.out$par
}
                                     
tl = numeric(1000)
thl.start = 0.5
for( i in 1:1000){
	Y3 = model(thbar,X1) + rlaplace(N,0.15)
    Y3[N] = 5
	o.out <- optim(par=c(thl.start), fn=crit, x=X1, y=Y3,method='Brent',lower =c(-0.5),upper=c(3))
	tl[i] = o.out$par
}

                                              
plot(X1,Y3)
lines(X1,model(l.out$par,X1))

plot(density(tl))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

var(ansl)
var(tl)

curve(exp(-x),0.5,1.45,pch=20,ylim=c(0,1),ann='F',lty=4,col='indianred',lwd=2)
points(X1,Y3,pch=4,col='skyblue')
points(sort(X1),model(sort(X1),o.out$par),type='l',lwd=2,col='darkorange')
points(sort(X1),model(sort(X1),l.out$par),type='l',lwd=2,col='forestgreen')

# Example 5

# Import data
dat <- read.csv("/Users/cianscannell/Desktop/Final_Year/Stats/group6-data.csv",header=FALSE)
head(dat)
X2 <- dat[,1]
Y4 <-dat[,2]
N<-length(X2)
h = 0.1

# linear regression model
lin <- function(th){
	return(th[1] + th[2]*X2)
}

# Implement our unsymmetrised criterion
eu <- function(theta){			
	return(Y4 - lin(theta))
}

au <- function(theta,ind1,ind2){	## component of kernel argument
	return((eu(theta)[ind1] - eu(theta)[ind2]) / h )
} 

F_unsym <- function(theta,ind){
	s = 0
	for(i in 1:N){
		s = s + kern(au(theta,ind,i))
	}
	s = 1/(N*h) * s
	return(s)
}


J_unsym <- function(theta){			## function to return criterion of interest
	s1 = 0 
	for( j in 1:N){
		s1 = s1 + log(F_unsym(theta,j))
	}
	s1 = -1/N * s1
	return(s1)
}

plot(X2,Y4)
abline(lm(Y4 ~ X2))

t.st <- c(0.0,0.0)
t.out <- optim(par=t.st,fn=J_unsym)
(t.out$par)

## Least Squares Estimation
crit_u <- function(theta){
	# must have theta as 1st parameter and return a single value...
	return( sum( (Y4-lin(theta))^2 ) )
}

e <- function(theta){			
  return(Y4 - lin(theta))
}


th.start = c(0,0)		
unsym.out <- optim(par=c(th.start),fn=J_unsym)
sym.out <- optim(par=c(th.start),fn=J)
ls.out <- optim(par=c(th.start),fn=crit_u)

(unsym.out$par)
(sym.out$par)
(ls.out$par)