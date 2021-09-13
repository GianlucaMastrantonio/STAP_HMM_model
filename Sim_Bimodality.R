PLOT_DIRPLOT  = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/sheepdog/tex/"

library(moveHMM)

#### #### #### #### ####
#### TEST
#### #### #### ####


n 		= 40000

# simulation of the displacement coordinates
set.seed(100)
# bimodal 2
theta = rcauchy(n,pi,0.1)%%(2*pi)
r     = rweibull(n,1,1)


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
data = prepData(s,type="UTM",coordNames=c("X","Y"))

data = prepData(s[seq(1,n,by=2),],type="UTM",coordNames=c("X","Y"))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)


# par1 = c(pi-0.1, pi/12,pi/4-0.2, pi/4-0.015)
# par2 = c(0.01, 0.005,0.1, 0.01)
# par3 = c(1,0.3,5,10)
# par4 = c(1,1,10,2.3)
#asymmetric 7
n = 20000
set.seed(100)
# theta = rcauchy(n,pi/12-0.02,0.01)%%(2*pi)
# r     = rweibull(n,10,2.2) 24
theta = rcauchy(n,pi/4-0.1,0.05)%%(2*pi)
r     = rweibull(n,1.7,15)
summary(r)
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
par(mfrow=c(4,3))
for(i in 1:30)
{
	data = prepData(s[seq(1,n,by=i),],type="UTM",coordNames=c("X","Y"))

	plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi), main=i)
	abline(v=c(-pi,pi), col=2)
}



data = prepData(s[seq(1,n,by=3),],type="UTM",coordNames=c("X","Y"))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)

# asymmetrix 7
set.seed(100)
theta = rcauchy(n,pi/2,0.001)%%(2*pi)
r     = rweibull(n,1.4,1)
plot(density(r))
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
data = prepData(s,type="UTM",coordNames=c("X","Y"))

data = prepData(s[seq(1,n,by=7),],type="UTM",coordNames=c("X","Y"))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)



# bimodal
set.seed(100)
theta = rcauchy(n,pi/2,0.0001)%%(2*pi)
r     = rweibull(n,1.4,1)

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
data = prepData(s,type="UTM",coordNames=c("X","Y"))

data = prepData(s[seq(1,n,by=9),],type="UTM",coordNames=c("X","Y"))
plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)



# asymmetric 18
set.seed(100)
theta = rcauchy(n,pi/12,0.005)%%(2*pi)
r     = rweibull(n,0.3,1)

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


data = prepData(s[seq(1,n,by=16),],type="UTM",coordNames=c("X","Y"))

plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)


par(mfrow=c(4,3))
for(i in 1:36)
{
	data = prepData(s[seq(1,n,by=i),],type="UTM",coordNames=c("X","Y"))

	plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
	abline(v=c(-pi,pi), col=2)
}




n = 100000
#  bimodality 18
set.seed(100)
theta = rcauchy(n,(pi-0.5)/6,0.001)%%(2*pi)
r     = rweibull(n,0.1,1)


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
#data = prepData(s[seq(1,n,by=18),],type="UTM",coordNames=c("X","Y"))



par(mfrow=c(4,3))
for(i in 1:12)
{
	data = prepData(s[seq(1,n,by=i),],type="UTM",coordNames=c("X","Y"))

	plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
	abline(v=c(-pi,pi), col=2)
}




data = prepData(s[seq(1,n,by=10),],type="UTM",coordNames=c("X","Y"))

plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
abline(v=c(-pi,pi), col=2)





#  bimodality 18
set.seed(100)
theta = rcauchy(n,pi/5-0.1,0.001)%%(2*pi)
r     = rweibull(n,0.5,1)


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


par(mfrow=c(4,3))
for(i in 1:36)
{
	data = prepData(s[seq(1,n,by=i),],type="UTM",coordNames=c("X","Y"))

	plot(density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3), xlim=c(-pi,pi))
	abline(v=c(-pi,pi), col=2)
}



###### ###### ###### ###### ###### ######
###### ###### Bimodality
###### ###### ###### ###### ###### ######
library(ggplot2)
par1 = c(pi-0.1, pi/12,pi/4-0.2, pi-0.35)
par2 = c(0.1, 0.005,0.1, 0.1)
par3 = c(1,0.3,5,1.7)
par4 = c(1,1,10,15)

# theta = rcauchy(n,pi/4-0.2,0.1)%%(2*pi)
# r     = rweibull(n,5,10)
iset = c(2,16,9,3)
# theta = rcauchy(n,pi-0.35,0.1)%%(2*pi)
# r     = rweibull(n,1.7,15)
# theta = rcauchy(n,pi/3-0.5,0.1)%%(2*pi)
# r     = rweibull(n,10,10)
# theta = rcauchy(n,pi/4-0.015,0.01)%%(2*pi)
# r     = rweibull(n,10,2.3)
# theta = rcauchy(n,pi/12-0.01,0.01)%%(2*pi)
# r     = rweibull(n,10,2.1)
##### i = 1
i = 1
nn = 500
dens_circ = matrix(NA, ncol=4, nrow= nn)
dens_lin = matrix(NA, ncol=4, nrow= nn)
theta_seq = seq(-pi, pi, length.out=nn)
r_seq     = seq(0, 10, length.out=nn)

dens_circ[,i] = dcauchy(theta_seq,par1[i],par2[i])+dcauchy(theta_seq+2*pi,par1[i],par2[i])+dcauchy(theta_seq-2*pi,par1[i],par2[i])
dens_lin[,i]  = dweibull(r_seq,par3[i],par4[i])


Data = data.frame("X" = theta_seq, "Y"=dens_circ[,i] )
P1=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))

P1




 Data = data.frame("X" = r_seq , "Y"=dens_lin[,i] )
P2=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Step-length")

pdf(paste(PLOT_DIRPLOT ,"SimAng1.pdf",sep=""))
print(P1)
dev.off()

pdf(paste(PLOT_DIRPLOT ,"SimLength1.pdf",sep=""))
print(P2)
dev.off()


i = 2
nn = 500
theta_seq = seq(-pi, pi, length.out=nn)
r_seq     = seq(0, 5, length.out=nn)

dens_circ[,i] = dcauchy(theta_seq,par1[i],par2[i])+dcauchy(theta_seq+2*pi,par1[i],par2[i])+dcauchy(theta_seq-2*pi,par1[i],par2[i])
dens_lin[,i]  = dweibull(r_seq,par3[i],par4[i])


Data = data.frame("X" = theta_seq, "Y"=dens_circ[,i] )
P1=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-Angle")



Data = data.frame("X" = r_seq, "Y"=dens_lin[,i] )
P2=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
 axis.text.x = element_text(face="bold",size=25),
 axis.text.y = element_text(face="bold",size=25),
 axis.title.x = element_text(face="bold",size=25),
 axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Step-length")


pdf(paste(PLOT_DIRPLOT ,"SimAng2.pdf",sep=""))
print(P1)
dev.off()

pdf(paste(PLOT_DIRPLOT ,"SimLength2.pdf",sep=""))
print(P2)
dev.off()




i = 3
nn = 500
theta_seq = seq(-pi, pi, length.out=nn)
r_seq     = seq(0, 5, length.out=nn)

dens_circ[,i] = dcauchy(theta_seq,par1[i],par2[i])+dcauchy(theta_seq+2*pi,par1[i],par2[i])+dcauchy(theta_seq-2*pi,par1[i],par2[i])
dens_lin[,i]  = dweibull(r_seq,par3[i],par4[i])


Data = data.frame("X" = theta_seq, "Y"=dens_circ[,i] )
P1=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-Angle")
P1


Data = data.frame("X" = r_seq, "Y"=dens_lin[,i] )
P2=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
 axis.text.x = element_text(face="bold",size=25),
 axis.text.y = element_text(face="bold",size=25),
 axis.title.x = element_text(face="bold",size=25),
 axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Step-length")


pdf(paste(PLOT_DIRPLOT ,"SimAng3.pdf",sep=""))
print(P1)
dev.off()

pdf(paste(PLOT_DIRPLOT ,"SimLength3.pdf",sep=""))
print(P2)
dev.off()




i = 4
nn = 500
theta_seq = seq(-pi, pi, length.out=nn)
r_seq     = seq(0, 0.01, length.out=nn)

dens_circ[,i] = dcauchy(theta_seq,par1[i],par2[i])+dcauchy(theta_seq+2*pi,par1[i],par2[i])+dcauchy(theta_seq-2*pi,par1[i],par2[i])
dens_lin[,i]  = dweibull(r_seq,par3[i],par4[i])


Data = data.frame("X" = theta_seq, "Y"=dens_circ[,i] )
P1=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-Angle")



Data = data.frame("X" = r_seq, "Y"=dens_lin[,i] )
P2=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
 axis.text.x = element_text(face="bold",size=25),
 axis.text.y = element_text(face="bold",size=25),
 axis.title.x = element_text(face="bold",size=25),
 axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Step-length")


pdf(paste(PLOT_DIRPLOT ,"SimAng4.pdf",sep=""))
print(P1)
dev.off()

pdf(paste(PLOT_DIRPLOT ,"SimLength4.pdf",sep=""))
print(P2)
dev.off()

##### Sim
n = 100000


for(ii in 1:length(par1))
{
	#ii = 1
	set.seed(100)
	theta = rcauchy(n,par1[ii],par2[ii])%%(2*pi)
	r     = rweibull(n,par3[ii],par4[ii])

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
	data = prepData(s,type="UTM",coordNames=c("X","Y"))

	data = prepData(s[seq(1,n,by=iset[ii]),],type="UTM",coordNames=c("X","Y"))
	dcirc = density(c(data$angle,data$angle-2*pi,data$angle+2*pi), na.rm=T, adjust=1/3)
	dlin  = density(data$step[data$step<100], na.rm=T, from=0,adjust=1.3)

	W = dcirc$x>=-pi & dcirc$x<pi

	Data = data.frame("X" = dcirc$x[W], "Y"=dcirc$y[W] )
	P1=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
	  axis.text.x = element_text(face="bold",size=25),
	  axis.text.y = element_text(face="bold",size=25),
	  axis.title.x = element_text(face="bold",size=25),
	  axis.title.y = element_text(face="bold",size=25)
	)+ylab("Density")+xlab("Turning-Angle")+ylim(c(0, max(dcirc$y[W])))
	P1

	Data = data.frame("X" = dlin$x, "Y"=dlin$y )
	P2=ggplot(Data, aes(x=X, y=Y))+geom_line()+theme(
	 axis.text.x = element_text(face="bold",size=25),
	 axis.text.y = element_text(face="bold",size=25),
	 axis.title.x = element_text(face="bold",size=25),
	 axis.title.y = element_text(face="bold",size=25)
	)+ylab("Density")+xlab("Step-length")
	P2

	pdf(paste(PLOT_DIRPLOT ,"HSimAng",ii, ".pdf",sep=""))
	print(P1)
	dev.off()

	pdf(paste(PLOT_DIRPLOT ,"HSimLength",ii, ".pdf",sep=""))
	print(P2)
	dev.off()

}
