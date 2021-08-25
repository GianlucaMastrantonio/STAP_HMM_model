library(moveHMM)

### ### ### ### ### ### ### ###
### Simulation 1
### ### ### ### ### ### ### ###

n 		= 20000

# simulation of the displacement coordinates



# bimodality
set.seed(100)
theta = rcauchy(n,pi-0.1,0.05)
r     = rgamma(n,1,1)

#skew
set.seed(100)
theta = rcauchy(n,pi-0.3,0.05)
r     = rgamma(n,7.5,1)

# the path
s         = matrix(NA, nrow=n+1, ncol=2)
s[1,]    	= c(0,0)

angle  		= cumsum(theta)

for(i in 2:(n+1))
{
	s[i,] = s[i-1,]+c(r[i-1]*cos(angle[i-1]), r[i-1]*sin(angle[i-1]))
}

#plot(s[1:10,], type="l")
colnames(s) = c("X","Y")
s = as.data.frame(s)

##

data = prepData(s,type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))




# displacement coordinates for time-difference 3, 5 and 20

par(mfrow=c(2,3))
for(i in 1:6)
{
	data = prepData(s[seq(1,n,by=i),],type="UTM",coordNames=c("X","Y"))

	plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
	abline(v=c(-pi,pi), col=2)
}




##################
############

data = prepData(s[seq(1,n,by=5),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)
par(mfrow=c(1,1))
plot(data$step,data$angle)


data = prepData(s[seq(1,n,by=10),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)
par(mfrow=c(1,1))
plot(data$step,data$angle)


##########
#########


library(moveHMM)

### ### ### ### ### ### ### ###
### Simulation 1
### ### ### ### ### ### ### ###

n 		= 20000

# simulation of the displacement coordinates
set.seed(100)
# theta = rcauchy(n,pi/2-0.1,0.1)
# r     = rgamma(n,5,1)

theta = rcauchy(n,pi/2,0.1)
r     = rgamma(n,1,1)

# the path
s         = matrix(NA, nrow=n+1, ncol=2)
s[1,]    	= c(0,0)

angle  		= cumsum(theta)

for(i in 2:(n+1))
{
	s[i,] = s[i-1,]+c(r[i-1]*cos(angle[i-1]), r[i-1]*sin(angle[i-1]))
}

colnames(s) = c("X","Y")
s = as.data.frame(s)
plot(s, type="l")

##

# data = prepData(s,type="UTM",coordNames=c("X","Y"))
#
#
# plot(density(data$step, na.rm=T))
# plot(density(data$angle, na.rm=T))

# displacement coordinates for time-difference 3, 5 and 20

data = prepData(s[seq(1,n,by=3),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)

par(mfrow=c(1,1))
plot(data$step,data$angle)


data = prepData(s[seq(1,n,by=5),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)
par(mfrow=c(1,1))
plot(data$step,data$angle)


data = prepData(s[seq(1,n,by=10),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)
par(mfrow=c(1,1))
plot(data$step,data$angle)






### ### ### ### ### ### ### ###
### ESEMPIO 2
### ### ### ### ### ### ### ###


n = 2000
theta = rcauchy(n,pi/2,1)
r     = rgamma(n,1,0.1)

# path

s         = matrix(NA, nrow=n+1, ncol=2)
s[1,]      = c(0,0)
angle  = cumsum(theta)

for(i in 2:(n+1))
{
	s[i,] = s[i-1,]+c(r[i-1]*cos(angle[i-1]), r[i-1]*sin(angle[i-1]))
}

colnames(s) = c("X","Y")
s = as.data.frame(s)
plot(s, type="l")

##
library(moveHMM)

data = prepData(s,type="UTM",coordNames=c("X","Y"))


plot(density(data$step, na.rm=T))
plot(density(data$angle, na.rm=T))

### sub samp

data = prepData(s[seq(1,n,by=3),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)


data = prepData(s[seq(1,n,by=5),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)

data = prepData(s[seq(1,n,by=20),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)



### ### ### ### ### ### ### ###
### ESEMPIO 3
### ### ### ### ### ### ### ###


n = 2000
theta = rcauchy(n,pi/2,1)
r     = rgamma(n,1,10)

# path

s         = matrix(NA, nrow=n+1, ncol=2)
s[1,]      = c(0,0)
angle  = cumsum(theta)

for(i in 2:(n+1))
{
	s[i,] = s[i-1,]+c(r[i-1]*cos(angle[i-1]), r[i-1]*sin(angle[i-1]))
}

colnames(s) = c("X","Y")
s = as.data.frame(s)
plot(s, type="l")

##
library(moveHMM)

data = prepData(s,type="UTM",coordNames=c("X","Y"))


plot(density(data$step, na.rm=T))
plot(density(data$angle, na.rm=T))

### sub samp

data = prepData(s[seq(1,n,by=3),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)


data = prepData(s[seq(1,n,by=5),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)

data = prepData(s[seq(1,n,by=20),],type="UTM",coordNames=c("X","Y"))

par(mfrow=c(1,2))
plot(density(data$step, na.rm=T))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)
