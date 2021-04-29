#### #### #### #### #### #### #### ####
####  Values that must be changed
#### #### #### #### #### #### #### ####
rm(list=ls())
PLOT_DIRDATA  = ""
PLOT_DIRPLOT  = ""
MOD_STAP_NAME = ""
MOD_BRW_NAME   = ""
MOD_CRW_NAME   = ""

#MOD_SIM_NAME = "Sim"


#### #### #### #### #### #### #### ####
####  Libraries
#### #### #### #### #### #### #### ####

library(ggplot2)
library(ggmosaic)
library(moveHMM)

### ### ### ### ### ###
### Functions
### ### ### ### ### ###

findmode = function(x){
        TT = table(as.vector(x))
        return(as.numeric(names(TT)[TT==max(TT)][1]))
}
dmnorm=function (x, mean = rep(0, d), varcov, log = FALSE)
{
    d <- if (is.matrix(varcov))
        ncol(varcov)
    else 1
    if (d > 1 & is.vector(x))
        x <- matrix(x, 1, d)
    n <- if (d == 1)
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(varcov) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    if (log)
        logPDF
    else exp(logPDF)
}
rmnorm=function(n = 1, mean = rep(0, d), varcov)
{
   d <- if (is.matrix(varcov))
       ncol(varcov)
   else 1
   z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
   y <- t(mean + t(z))
   return(y)
}

### ### ### ### ### ###
### Colors palette
### ### ### ### ### ###
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- rev(c(
  "#ffffbf",
"#d7191c",
"#fdae61",
"#abdda4",
"#2b83ba")
)



#### #### #### #### #### #### #### ####
#### Coordinate Systems
#### #### #### #### #### #### #### ####

angle = -1/8*pi
RotMat = matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),nrow=2)
#RotMat = diag(1,2)
CoordAx1 =  c(-1,0)%*%RotMat
CoordAx2 =  c(1,0)%*%RotMat
CoordAy1 =  c(0,1.)%*%RotMat
CoordAy2 =  c(0,-0.5)%*%RotMat

PT = qplot(c(-10),10,geom="line", xlab="",ylab="", ylim=c(-1,1), xlim=c(-1,1))+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")+ coord_fixed(ratio=1)

PT1 = PT + geom_segment(aes(x = CoordAx1[1], y = CoordAx1[2], xend = CoordAx2[1], yend = CoordAx2[2]), linetype =  2)+geom_segment(aes(x = CoordAy1[1], y = CoordAy1[2], xend = CoordAy2[1], yend = CoordAy2[2]), linetype =  2)+
geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), linetype =  2)+
geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype =  2)

OBS = matrix(c(-0.8,0,0.8,0,0,0.6),nrow=3,ncol=2)%*%RotMat

PT2 = PT1+geom_line(aes(x=OBS[,1],y=OBS[,2]))+geom_point(aes(x=OBS[,1],y=OBS[,2],size=1))
PT2
Desx1 = rbind(c(0.8,0.6),cbind(0.8,0))%*%RotMat
Desx2 = rbind(c(0.8,0.6),cbind(0,0.6))%*%RotMat
Desx1_v2 = Desx1
Desx1_v2[2,2] = 0
Desx1_v2[2,1] = Desx1_v2[1,1]
Desx2_v2 = Desx2
Desx2_v2[2,1] = 0
Desx2_v2[2,2] = Desx2_v2[1,2]

PT3 = PT2+geom_line(aes(x=Desx1[,1],y=Desx1[,2]), linetype = 3)+geom_line(aes(x=Desx2[,1],y=Desx2[,2]), linetype = 3)+geom_line(aes(x=Desx1_v2[,1],y=Desx1_v2[,2]), linetype = 3) + geom_line(aes(x=Desx2_v2[,1],y=Desx2_v2[,2]), linetype = 3)


Xc 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Yc 		= rbind(c(-0.02,0.6),c(0.02,0.6))
Xc 		= Xc%*%RotMat
Yc 		= Yc%*%RotMat
Xc_v2 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Xc_v2[1,]    =  OBS[3,]
Xc_v2[1,1]    =  -0.02
Xc_v2[2,1]    =  0.02
Xc_v2[2,2]  = Xc_v2[2,2]=  OBS[3,2]

Yc_v2 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Yc_v2[1,]    =  OBS[3,]
Yc_v2[1,2]    =  -0.02
Yc_v2[2,2]    =  0.02
Yc_v2[2,1]  = Yc_v2[2,1]=  OBS[3,1]

PT4 = PT3+geom_line(aes(x=Xc[,1],y=Xc[,2]), linetype = 1,size=1)+geom_line(aes(x=Yc[,1],y=Yc[,2]), linetype = 1,size=1)+geom_line(aes(x=Xc_v2[,1],y=Xc_v2[,2]), linetype = 1,size=1)+geom_line(aes(x=Yc_v2[,1],y=Yc_v2[,2]), linetype = 1,size=1)


PT4

PT5 =PT4+
geom_text(aes(x=Xc[2,1]+0.1,y=Xc[2,2]-0.06), label=deparse(bquote(paste('y'['t'['i']*','][1]))),parse=TRUE,size=7)+
geom_text(aes(x=Yc[2,1]-0.11,y=Yc[2,2]-0.04), label=deparse(bquote(paste('y'['t'['i']*','][2]))),parse=TRUE,size=7)
PT5




PT6 =PT5+
geom_text(aes(x=OBS[1,1]-0.1,y=OBS[1,2]+0.04), label=deparse(bquote(paste('s'['t'['i-1']]))),parse=TRUE,size=7)+
geom_text(aes(x=OBS[2,1]-0.12,y=OBS[2,2]+0.02), label=deparse(bquote(paste('s'['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=OBS[3,1]+0.12,y=OBS[3,2]-0.02), label=deparse(bquote(paste('s'['t'['i+1']]))),parse=TRUE,size=7)
PT6


cangle = seq(0.38,pi/2-0.5,length.out=100)
S1 = sin(cangle)*0.2
C1 = cos(cangle)*0.2

PT7 = PT6 +geom_line(aes(x=C1,y=S1), linetype = 1,size=0.4)+geom_text(aes(x=0.2,y=0.18), label=deparse(bquote(paste(theta['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=0.37,y=0.45), label=deparse(bquote(paste('r'['t'['i']]))),parse=TRUE,size=7)
PT7


PT8 = PT7+geom_segment(aes(x = OBS[1,1], y = OBS[1,2], xend = 1, yend = OBS[1,2]), linetype =  3)+
geom_segment(aes(x = OBS[2,1], y = OBS[2,2], xend = 1, yend = OBS[2,2]), linetype =  3)


cangle = seq(0.38,0,length.out=100)
S11 = sin(cangle)*0.2
C11 = cos(cangle)*0.2

cangle = seq(1.03,0,length.out=100)
S12 = sin(cangle)*0.35
C12 = cos(cangle)*0.35

PT9 = PT8+geom_line(aes(x=C11+OBS[1,1],y=S11+OBS[1,2]), linetype = 1,size=0.4)+
geom_line(aes(x=C12+OBS[2,1],y=S12+OBS[2,2]), linetype = 1,size=0.4)+
geom_text(aes(x=0.40,y=0.05), label=deparse(bquote(paste(phi['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=-0.45,y=-0.25), label=deparse(bquote(paste(phi['t'['i-1']]))),parse=TRUE,size=7)+
geom_text(aes(x=-0.1,y=0.9), label=deparse(bquote(paste('v'['t'['i']*','][2]))),parse=TRUE,size=7)+
geom_text(aes(x=0.5,y=-0.1), label=deparse(bquote(paste('v'['t'['i']*','][1]))),parse=TRUE,size=7)
PT9


library(car)
library(mvtnorm)

app = ellipse(c(0,0), matrix(c(1.9,-0.7,-0.7,0.5), ncol=2)*0.02, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+0.35
app[,2] = app[,2]+1

PT10 = PT9+geom_path(aes(x=app[,1],y=app[,2]),  size = 1)+
geom_segment(aes(x = 0, y = 0, xend = 0.35, yend = 1, color= "Prev. Dir."),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1, linetype=1, color= "black")+scale_x_continuous(breaks= c(-1,-0.5,0, 0.5,1), labels = c(5,6,7,8,9), limits=c(-1,1))+scale_y_continuous(breaks= c(-1,-0.5,0, 0.5,1), labels = c(5,6,7,8,9)+5, limits=c(-0.7,1.3))
PT10+ coord_fixed()
pdf(paste(PLOT_DIRPLOT,"TransX.pdf",sep=""),width=7, height=7)
PT10+ coord_fixed()
dev.off()

####################################################
#### Projected normal
####################################################
dProj= function( x,mu1, mu2, sigma2_1, sigma2_2=1,rho)
{
	a = 1/(sqrt(sigma2_1)*sqrt(sigma2_2)*sqrt(1-rho^2))
	C = a^2*(sigma2_2*cos(x)^2-rho*sqrt(sigma2_1)*sqrt(sigma2_2)*sin(2*x)+ sigma2_1*sin(x)^2)
	D = a^2*C^(-0.5)*(
mu1*sqrt(sigma2_2)*(sqrt(sigma2_2)*cos(x)-rho*sqrt(sigma2_1)*sin(x))+
mu2*sqrt(sigma2_1)*(sqrt(sigma2_1)*sin(x)-rho*sqrt(sigma2_2)*cos(x)   )
)


	f = dmnorm(c(mu1,mu2), c(0,0), matrix(c(sigma2_1,sqrt(sigma2_1)*sqrt(sigma2_2)*rho,sqrt(sigma2_1)*sqrt(sigma2_2)*rho,sigma2_2), ncol=2))+ a*D*pnorm(D,0,1)*dnorm(a*C^(-0.5)*(mu1*sin(x)-mu2*cos(x)))
	return(f/C)
}

theta = seq(-pi,pi, length.out=100)
dataD = data.frame(Angle=rep(theta,3),
	Density = c(
			dProj( x=theta,mu1=1.2, mu2=0, sigma2_1=0.5, sigma2_2=1,rho=0),
			dProj( x=theta,mu1=0.5, mu2=-0.8, sigma2_1=0.5, sigma2_2=1,rho=0),
			dProj( x=theta,mu1=0.2, mu2=0.2, sigma2_1=0.5, sigma2_2=1,rho=0)
		),
	Index = factor(rep(1:3, each=length(theta)))
)

D1 = ggplot(dataD, aes(x =Angle, y = Density, group = Index , color=Index))+geom_line(size=2)+  theme(axis.text.y = element_text(face="bold",size=20),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=20),
legend.title = element_text(face="bold",size=25) ,legend.text.align = 0)+scale_color_manual(values = c(cbPalette[c(4,2,1)]),name="",
labels=c(expression(mu ~"=(1.2,0)'"),expression(mu ~"=(0.5,-0.8)'"),expression(mu ~"=(0.2,0.2)'"))
)+theme(legend.position="bottom")

pdf(paste(PLOT_DIRPLOT,"ProjN.pdf",sep=""),width=7, height=7)
D1
dev.off()


####################################################
#### Plot CRW - DRW - STAP
####################################################

Mu = matrix(c(0,0), ncol=1)
Nu = matrix(c(0,6), ncol=1)
tau = 0.25
sigma21 = 0.2
sigma22 = 1
corr   = 0
rho = c(0.33,0.66)
Sigma = matrix(c(sigma21,(sigma21*sigma22)^0.5*corr,(sigma21*sigma22)^0.5*corr,sigma22), ncol=2)

angle = rev(rep(seq(0,-pi,by=-pi/3), each=4))
xseq = yseq = seq(-20,20, length.out=4)
xseq = yseq = seq(-20,20, length.out=4)
Grid = expand.grid(xseq, yseq)
colnames(Grid) = c("Longitude", "Latitude")

# BRW




P = qplot(c(210),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

meanPlot_2 = meanPlot

i = 1
SigmaPlot = Sigma
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+meanPlot[i,1]
app[,2] = app[,2]+meanPlot[i,2]
ell = app
for(i in 2:nrow(Grid))
{

    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+meanPlot[i,1]
    app[,2] = app[,2]+meanPlot[i,2]
    ell = rbind(ell,NA, app)

}
ell_BRW = ell
P_BWR = P1

# CRW


#P = qplot(c(0),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        #axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")+ coord_fixed(ratio=1)

P1 = P_BWR


i = 1
arr = Grid

ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot = R%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = R%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
    ell = rbind(ell,NA, app)

}

cc = factor(c(1,2,3))



col_ell_1 = rep(cbPalette[1], nrow(ell_BRW))
col_ell_2 = rep(cbPalette[4], nrow(ell))
ell = rbind(ell,NA, ell_BRW)

P4 = P1# + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="CRW"),
                  arrow = arrow(length = unit(0.3, "cm")), size = 1.2)
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot_2[,1], yend = meanPlot_2[,2], color="BRW"),
                arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2], color= "Prev. Dir."),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")

PTOT = P8+scale_color_manual(values = c(cbPalette[c(1,4)],"black"),name="")+  theme(axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25) )+
geom_path(aes(x=ell[1:length(col_ell_2),1],y=ell[1:length(col_ell_2),2]),  size = 1, color=c(col_ell_2))+
geom_path(aes(x=ell[-c(1:(length(col_ell_2)+1)),1],y=ell[-c(1:(length(col_ell_2)+1)),2]),  size = 1., color=c(col_ell_1))+
theme(legend.position="bottom")
PTOT+xlim(c(-30,30))+ylim(c(-30,30))

pdf(paste(PLOT_DIRPLOT,"EsMov.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-30,30))+ylim(c(-30,30))
dev.off()




### ### ### ### ###
### STAP
### ### ### ### ###



P = qplot(c(120),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

j=1
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}
ell_2 = ell
arr_2 = arr
muPlot_2 = arr
j=2
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}

col_ell_1 = rep(cbPalette[1], nrow(ell_2))
col_ell_2 = rep(cbPalette[4], nrow(ell))
ell = rbind(ell,NA, ell_2)

P4 = P1# + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                  arrow = arrow(length = unit(0.3, "cm")), size = 1.2)
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")

PTOT =P8+scale_color_manual(values = c(cbPalette[c(1,4)],"black"),name="", labels=c(expression(rho ~"=1/3"),expression(rho ~"=2/3")))+  theme(axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25),legend.text.align = 0 )+
geom_path(aes(x=ell[1:length(col_ell_2),1],y=ell[1:length(col_ell_2),2]),  size = 1, color=c(col_ell_2))+
geom_path(aes(x=ell[-c(1:(length(col_ell_2)+1)),1],y=ell[-c(1:(length(col_ell_2)+1)),2]),  size = 1, color=c(col_ell_1))+
theme(legend.position="bottom")
PTOT





pdf(paste(PLOT_DIRPLOT,"EsMov2.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-30,30))+ylim(c(-30,30))
dev.off()






####################################################
#### Data
####################################################

load(paste(PLOT_DIRDATA , MOD_STAP_NAME, ".Rdata",sep=""))

DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2] )



p = ggplot(DataZ, aes(Longitude, Latitude))
p = p+
#geom_path(size = 0.1)+
geom_point(size = 0.7)
p = p+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)
p
pdf(paste(PLOT_DIRPLOT ,"Data.pdf",sep=""))
print(p+xlim(c(-5,5))+ylim(c(-5,5)))
dev.off()

DataCoords_2 = as.data.frame(DataCoords)
colnames(DataCoords_2) = c("x","y")
dataDog  <- prepData(DataCoords_2,type="UTM")
dataDog$step = log(dataDog$step)
P = ggplot(dataDog, aes(x=step))+geom_density(aes(x=step))+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Log Step-length")
P
pdf(paste(PLOT_DIRPLOT ,"ObsStepLength.pdf",sep=""))
print(P)
dev.off()

dataDog2 = rbind(dataDog[,-1] -2*pi,dataDog[,-1],dataDog[,-1]+2*pi)
l <- density(dataDog2$angle,na.rm=T, bw=0.3)
ld =data.frame(x=l$x,y=l$y)
P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))
pdf(paste(PLOT_DIRPLOT ,"ObsTurning.pdf",sep=""))
print(P)
dev.off()


P = ggplot(dataDog, aes(x=x))+geom_density(aes(x=x))+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Coordinate X")
P
pdf(paste(PLOT_DIRPLOT ,"ObsCoordX.pdf",sep=""))
print(P)
dev.off()




P = ggplot(dataDog, aes(x=y))+geom_density(aes(x=y))+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Coordinate Y")
P
pdf(paste(PLOT_DIRPLOT ,"ObsCoordY.pdf",sep=""))
print(P)
dev.off()



nt = nrow(DataAn3)
DataPec = rbind(DataAn3[,1:2],DataAn3[,1:2+2],DataAn3[,1:2+4],DataAn3[,1:2+6])
DataPec = as.data.frame(DataPec)
DataPec$Sheep = factor(rep(1:4,nt))
colnames(DataPec) = c("Longitude", "Latitude","Sheep")


DataPec$Longitude = (DataPec$Longitude-Meantot1)/SDtot1
DataPec$Latitude = (DataPec$Latitude-Meantot2)/SDtot2



p = ggplot(DataPec, aes(Longitude, Latitude))
p = p +geom_point(size = 2.5)
p = p+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
)+xlim(min(DataZ[,1], na.rm=T), max(DataZ[,1], na.rm=T))+ylim(min(DataZ[,2], na.rm=T), max(DataZ[,2], na.rm=T))

p

pdf(paste(PLOT_DIRPLOT ,"DataSheeps.pdf",sep=""))
print(p)
dev.off()


##### Movement prediction

Mu = matrix(c(0,0), ncol=1)/5
Nu = matrix(c(4,0), ncol=1)/5
tau = 0.25
sigma21 = 0.05
sigma22 = 0.5
corr   = -0.8
rho = c(0.33,0.66)
Sigma = matrix(c(sigma21,(sigma21*sigma22)^0.5*corr,(sigma21*sigma22)^0.5*corr,sigma22), ncol=2)/5^2

angle = rev(rep(seq(0,-pi,by=-pi/4), each=5))
xseq = yseq = seq(-20,20, by=10)/5
xseq = yseq = seq(-20,20, by=10)/5
Grid = expand.grid(xseq, yseq)
colnames(Grid) = c("Longitude", "Latitude")

# BRW




P = qplot(c(0),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

meanPlot_2 = meanPlot

i = 1
SigmaPlot = Sigma
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+meanPlot[i,1]
app[,2] = app[,2]+meanPlot[i,2]
ell = app
for(i in 2:nrow(Grid))
{

    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+meanPlot[i,1]
    app[,2] = app[,2]+meanPlot[i,2]
    ell = rbind(ell,NA, app)

}
ell_BRW = ell
P_BWR = P1

# CRW


#P = qplot(c(0),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        #axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")+ coord_fixed(ratio=1)

P1 = P_BWR


i = 1
arr = Grid

ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot = R%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = R%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
    ell = rbind(ell,NA, app)

}

cc = factor(c(1,2,3))




col_ell_1 = rep(cbPalette[1], nrow(ell_BRW))
col_ell_2 = rep(cbPalette[3], nrow(ell))
ell = rbind(ell,NA, ell_BRW)

P4 = P1 + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="CRW"),
                  arrow = arrow(length = unit(0.2, "cm")))
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot_2[,1], yend = meanPlot_2[,2], color="BRW"),
                arrow = arrow(length = unit(0.2, "cm")))

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2], color= "Prev. Dir."),arrow = arrow(length = unit(0.2, "cm")))

PTOT =P8+scale_color_manual(values = cbPalette[c(1,3,4)],name="Velocity")+ theme(legend.text=element_text(size=13),legend.title=element_text(size=18))



pdf(paste(PLOT_DIRPLOT,"EsMovdd.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-26/5,26/5))+ylim(c(-26/5,26/5))
dev.off()


### ### ### ### ###
### STAP
### ### ### ### ###



P = qplot(c(0),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

j=1
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}
ell_2 = ell
arr_2 = arr
muPlot_2 = arr
j=2
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}

col_ell_1 = rep(cbPalette[1], nrow(ell_2))
col_ell_2 = rep(cbPalette[3], nrow(ell))
ell = rbind(ell,NA, ell_2)

P4 = P1 + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                  arrow = arrow(length = unit(0.2, "cm")))
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                arrow = arrow(length = unit(0.2, "cm")))

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2], color= "3"),arrow = arrow(length = unit(0.2, "cm")))

PTOT =P8+scale_color_manual(values = cbPalette[c(1,3,4)],name="Velocity", labels=c(expression(rho ~"=1/3"),expression(rho ~"=2/3"),"Prev. Dir"))+ theme(legend.text=element_text(size=13),legend.title=element_text(size=18))



pdf(paste(PLOT_DIRPLOT,"EsMov2ddd.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-26/5,26/5))+ylim(c(-26/5,26/5))
dev.off()
####################################################
#### STAP-MODEL
####################################################

load(paste(PLOT_DIRDATA , MOD_STAP_NAME, ".Rdata",sep=""))
#WM = c(14,164,71,23,26)
WM = c(14,189,120,99,114)
WMgen = WM

Z_MAP = apply(ModelOUT$zeta,2,findmode)
table(Z_MAP)

### Posterior estimates
print("")
print("")
print("")
print("Posterior Estimates - STAP MODEL")
print("")
for(i in 1:2)
{
  dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$mu0[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu0[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu0[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}
for(i in 1:2)
{
  dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$muC[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$muC[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$muC[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}
for(i in 1:1)
{
  dd = paste("$backslash tau_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$psi[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$psi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$psi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}

for(i in 1:1)
{
  dd = paste("$backslash rho_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$rho[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$rho[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$rho[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}


for(i in 1:2)
{
  for(j in i:2)
  {
    dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }


}

ik=1
ddpi =c()
for(k in WM)
{

  dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
  ff = "(CI) "
  ii = 1
  for(i in WM)
  {
      ddpi [(ik-1)*length(WM)+ii] = round(mean(ModelOUT$pi[,i,k]),3)
    dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    ii = ii+1
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
  ik = ik+1
}

ik = 1
for(i in 1:1)
{
  dd = paste("$backslash beta_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}
gammaC = ModelOUT$gammaDP
kappaC = ModelOUT$rhoDP*ModelOUT$akDP
alphaC = ModelOUT$akDP-kappaC
print("backslash hline backslash hline ")
dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
dd = paste("$backslash hat{}$ ",sep="")
ff = "(CI) "

dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd1 = paste(dd1, " \\")
dd = paste(dd, " \\")
ff = paste(ff, " \\")
print(dd1)
print("backslash hline  ")
print(dd)
print(ff)


print("Distribution of K")
table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)

print("Confusion Matrix")
sum(diag(table(LatentClassification[-length(LatentClassification)],Z_MAP)[order(WM),])/length(Z_MAP))

#
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbPalette <- rev(c(
#   "#ffffbf",
# "#d7191c",
# "#fdae61",
# "#abdda4",
# "#2b83ba")
# )

#plot(1:8,col=cbPalette,cex=2,pch=20)

AppWM = rep(0,max(WM))
AppWM[WM] = 1:5

Z_MAP2_gen = AppWM[Z_MAP]

cbPalette <- c(
"#e41a1c",
"#377eb8",
"#4daf4a",
"#984ea3",
"#ff7f00"
)

DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2],Cluster = as.factor(Z_MAP2_gen) )
for(i in 1:5)
{
  p = ggplot(subset(DataZ, Cluster %in% c(paste(i))), aes(Longitude, Latitude, shape = Cluster))
  p = p +geom_point(size = 2.5, shape=16,color=cbPalette[i])
  p = p+theme(
    axis.text.x = element_text(face="bold",size=25),
    axis.text.y = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )+xlim(min(DataZ[,1], na.rm=T), max(DataZ[,1], na.rm=T))+ylim(min(DataZ[,2], na.rm=T), max(DataZ[,2], na.rm=T))
  p
  if((i==4) | (i==5))
  {
    Xp = mean(ModelOUT$mu0[,1,WM[i]])
    Yp = mean(ModelOUT$mu0[,2,WM[i]])
    p = p +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
    #+
    #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
  }

  pdf(paste(PLOT_DIRPLOT ,"DataPost",i,".pdf",sep=""))
  print(p+xlim(c(-5,5))+ylim(c(-5,5)))
  dev.off()
}



###
nmcmc = dim(ModelOUT$psi)[1]
nseq = 70
thetaseq = seq(-pi,pi,length.out=nseq+1)[-(nseq+1)]
rseq = seq(-8,1.2,length.out=nseq)+0.0000000000001
diff1 = thetaseq[2]-thetaseq[1]
diff2 = rseq[2]-rseq[1]
Densmat = list()
Densmat[[1]] = matrix(0, ncol=nseq,nseq)
Densmat[[2]] = matrix(0, ncol=nseq,nseq)
Densmat[[3]] = matrix(0, ncol=nseq,nseq)
Densmat[[4]] = matrix(0, ncol=nseq,nseq)
Densmat[[5]] = matrix(0, ncol=nseq,nseq)
Vecy  = matrix(NA, ncol=2, nrow=nseq*nseq)
Rvec     = rep(exp(rseq), each=nseq)
Vecy[,1] = rep(exp(rseq), each=nseq)*cos(rep(thetaseq,times=nseq))
Vecy[,2] = rep(exp(rseq), each=nseq)*sin(rep(thetaseq,times=nseq))
for(imcmc in 1:nmcmc)
{
  k1 = 0
  for(k in WM)
  {
    k1 = k1+1
    mu = ModelOUT$muC[imcmc,,k]
    sigma = matrix(ModelOUT$sigma[imcmc,,k],2)

    Densmat[[k1]][,] = Densmat[[k1]][,]+matrix(dmnorm(Vecy,mu, sigma)*Rvec^2, ncol=nseq)/nmcmc

  }

}


DataTR = data.frame(Theta = c(rowSums(Densmat[[1]]*diff2),rowSums(Densmat[[2]]*diff2)),R = c(colSums(Densmat[[1]]*diff1),colSums(Densmat[[2]]*diff1)) , Gruppo = as.factor(rep(c(1,2),each =nseq )),thetaseq = rep(thetaseq, times=2),rseq = rep(rseq, times=2))


p = ggplot(DataTR,aes(x=thetaseq, y =Theta,group = Gruppo   ))+
geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
) +ylim(0,0.300)+xlab("Turning-Angle")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette[c(1,2)])+guides(shape = guide_legend(override.aes = list(size = 5)))
p

pdf(paste(PLOT_DIRPLOT ,"Turn1.pdf",sep=""))
print(p)
dev.off()

p = ggplot(DataTR,aes(x=rseq, y =R,group = Gruppo   ))+
geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
) +ylim(0,0.8)+xlab("Log Step-Length")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette[c(1,2)])+guides(shape = guide_legend(override.aes = list(size = 5)))
p

pdf(paste(PLOT_DIRPLOT ,"Step1.pdf",sep=""))
print(p)
dev.off()

#### #### #### #### #### #### #### ####
#### Prediction
#### #### #### #### #### #### #### ####


mu0List = list()
muCList = list()
psiList = list()
rhoList = list()
sigmaList = list()
piList = list()
kk = 1
for(k in WM)
{
    mu0List[[kk]] = matrix(colMeans(ModelOUT$mu0[,1:2,k]),ncol=1)
    muCList[[kk]] = matrix(colMeans(ModelOUT$muC[,1:2,k]),ncol=1)
    psiList[[kk]] = matrix(mean(ModelOUT$psi[,1,k]),ncol=1)
    rhoList[[kk]] = matrix(mean(ModelOUT$rho[,1,k]),ncol=1)
    sigmaList[[kk]] = matrix(colMeans(ModelOUT$sigma[,1:4,k]),nrow=2)
    piList[[kk]] = matrix(colMeans(ModelOUT$pi[,1:length(WM),k]),ncol=length(WM))
    kk = kk+1

}


DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2],Cluster = as.factor(Z_MAP2_gen) )

kk = 4
for(kk in 1:length(WM))
{




    Mu = mu0List[[kk]]
    Nu = muCList[[kk]]
    tau = psiList[[kk]]


    rho = rhoList[[kk]]
    Sigma = sigmaList[[kk]]

    xseq = yseq = seq(-4,4,length.out=4)
    xseq = yseq = seq(-4,4, length.out=4)
    nnnn = length(xseq)
    angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

    Grid = expand.grid(xseq, yseq)
    colnames(Grid) = c("Longitude", "Latitude")


    P =qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
            axis.text.y = element_text(face="bold",size=25),
            axis.title.x = element_text(face="bold",size=25),
            axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



    meanPlot = Grid[,]
    meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
    meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

    P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
    #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                     # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

    j=1
    i = 1
    arr = Grid
    ang = angle[i]
    R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
    arr[i,] = arr[i,]+muPlot
    SigmaPlot = R%*%Sigma%*%t(R)
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = app
    for(i in 2:nrow(Grid))
    {
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        SigmaPlot = R%*%Sigma%*%t(R)
        arr[i,] = arr[i,]+muPlot
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = rbind(ell,NA, app)

    }
    ell_2 = ell
    arr_2 = arr
    muPlot_2 = arr


    # col_ell_1 = rep(cbPalette[1], nrow(ell_2))
    # col_ell_2 = rep(cbPalette[3], nrow(ell))
    # ell = rbind(ell,NA, ell_2)

    P4 = P1


    i = 1
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    arr2 = Grid
    muPlot = Grid[i,]+(R)%*%matrix(c(-5/5,0), ncol=1)
    arr2[i,] = muPlot
    for(i in 2:nrow(Grid))
    {
        ang = angle[i]
        R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot = Grid[i,]+(R)%*%matrix(c(-5/5,0), ncol=1)
        arr2[i,] = muPlot
    }
    #P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                      #arrow = arrow(length = unit(0.2, "cm")))
    #P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                    #arrow = arrow(length = unit(0.2, "cm")))


    P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")
    P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(kk)],
                    arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

    PTOT =P7+
    #scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
     theme(
         axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 2,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

    PP = PTOT+xlim(c(-5,5))+ylim(c(-5,5))+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)

PP

	PP2 = PP

	# if((kk==4) | (kk==5))
	# {
	#   Xp = mean(ModelOUT$mu0[,1,WM[i]])
	#   Yp = mean(ModelOUT$mu0[,2,WM[i]])
	#   PP2 = PP2 +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
	#   #+
	#   #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
	# }
	#PP2


    pdf(paste(PLOT_DIRPLOT,kk,"Mov.pdf",sep=""),width=7, height=7)
    print(PP2+xlim(c(-5,5))+ylim(c(-5,5)))
    dev.off()

}






kk = 1
Mu = mu0List[[kk]]
Nu = muCList[[kk]]
tau = psiList[[kk]]


rho = rhoList[[kk]]
Sigma = sigmaList[[kk]]

xseq = yseq = seq(-0.1,0.1 ,length.out=2)
xseq = yseq = seq(-0.1,0.1, length.out=2)
nnnn = length(xseq)
angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))

Grid = expand.grid(xseq, yseq)
colnames(Grid) = c("Longitude", "Latitude")


P =qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
				axis.text.y = element_text(face="bold",size=25),
				axis.title.x = element_text(face="bold",size=25),
				axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
								 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

j=1
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
		ang = angle[i]
		R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
		R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
		muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
		SigmaPlot = R%*%Sigma%*%t(R)
		arr[i,] = arr[i,]+muPlot
		app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
		app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
		app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
		ell = rbind(ell,NA, app)

}
ell_2 = ell
arr_2 = arr
muPlot_2 = arr


# col_ell_1 = rep(cbPalette[1], nrow(ell_2))
# col_ell_2 = rep(cbPalette[3], nrow(ell))
# ell = rbind(ell,NA, ell_2)

P4 = P1


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-0.25/5,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
		ang = angle[i]
		R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
		muPlot = Grid[i,]+(R)%*%matrix(c(-0.25/5,0), ncol=1)
		arr2[i,] = muPlot
}
#P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
									#arrow = arrow(length = unit(0.2, "cm")))
#P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
								#arrow = arrow(length = unit(0.2, "cm")))


P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")
P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(kk)],
								arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

PTOT =P7+
#scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
 theme(
		 axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 2,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

PP = PTOT+xlim(c(-0.2,0.2))+ylim(c(-0.2,0.2))
#+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)

PP

PP2 = PP

# if((kk==4) | (kk==5))
# {
#   Xp = mean(ModelOUT$mu0[,1,WM[i]])
#   Yp = mean(ModelOUT$mu0[,2,WM[i]])
#   PP2 = PP2 +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
#   #+
#   #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
# }
#PP2


pdf(paste(PLOT_DIRPLOT,kk,"DetMov.pdf",sep=""),width=7, height=7)
print(PP2)
dev.off()




kk = 2
Mu = mu0List[[kk]]
Nu = muCList[[kk]]
tau = psiList[[kk]]


rho = rhoList[[kk]]
Sigma = sigmaList[[kk]]

xseq = yseq = seq(-0.6,0.6 ,length.out=2)
xseq = yseq = seq(-0.6,0.6, length.out=2)
nnnn = length(xseq)
angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))

Grid = expand.grid(xseq, yseq)
colnames(Grid) = c("Longitude", "Latitude")


P =qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
				axis.text.y = element_text(face="bold",size=25),
				axis.title.x = element_text(face="bold",size=25),
				axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
								 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

j=1
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
		ang = angle[i]
		R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
		R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
		muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
		SigmaPlot = R%*%Sigma%*%t(R)
		arr[i,] = arr[i,]+muPlot
		app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
		app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
		app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
		ell = rbind(ell,NA, app)

}
ell_2 = ell
arr_2 = arr
muPlot_2 = arr


# col_ell_1 = rep(cbPalette[1], nrow(ell_2))
# col_ell_2 = rep(cbPalette[3], nrow(ell))
# ell = rbind(ell,NA, ell_2)

P4 = P1


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-1/5,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
		ang = angle[i]
		R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
		muPlot = Grid[i,]+(R)%*%matrix(c(-1/5,0), ncol=1)
		arr2[i,] = muPlot
}
#P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
									#arrow = arrow(length = unit(0.2, "cm")))
#P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
								#arrow = arrow(length = unit(0.2, "cm")))


P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")
P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(kk)],
								arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

PTOT =P7+
#scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
 theme(
		 axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 2,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

PP = PTOT+xlim(c(-1.1,1.1))+ylim(c(-1.1,1.1))
#+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)

PP

PP2 = PP

# if((kk==4) | (kk==5))
# {
#   Xp = mean(ModelOUT$mu0[,1,WM[i]])
#   Yp = mean(ModelOUT$mu0[,2,WM[i]])
#   PP2 = PP2 +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
#   #+
#   #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
# }
#PP2


pdf(paste(PLOT_DIRPLOT,kk,"DetMov.pdf",sep=""),width=7, height=7)
print(PP2)
dev.off()




#### #### #### #### #### ####
#### MCMC
#### #### #### #### #### ####



kk = 2
for(kk in 1:length(WM))
{
    nmcmc = dim(ModelOUT$mu0)[1]

    imcmc = 1
    kkk =1
    for(k in WM)
    {
        mu0List[[kkk]] = matrix((ModelOUT$mu0[imcmc,1:2,k]),ncol=1)
        muCList[[kkk]] = matrix((ModelOUT$muC[imcmc,1:2,k]),ncol=1)
        psiList[[kkk]] = matrix((ModelOUT$psi[imcmc,1,k]),ncol=1)
        rhoList[[kkk]] = matrix((ModelOUT$rho[imcmc,1,k]),ncol=1)
        sigmaList[[kkk]] = matrix((ModelOUT$sigma[imcmc,1:4,k]),nrow=2)
        piList[[kkk]] = matrix((ModelOUT$pi[imcmc,1:length(WM),k]),ncol=length(WM))
        kkk = kkk+1

    }

    Mu = mu0List[[kk]]
    Nu = muCList[[kk]]
    tau = psiList[[kk]]


    rho = rhoList[[kk]]
    Sigma = sigmaList[[kk]]

    xseq = yseq = seq(-4,4,length.out=4)
    xseq = yseq = seq(-4,4, length.out=4)
    nnnn = length(xseq)
    angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

    Grid = expand.grid(xseq, yseq)
    colnames(Grid) = c("Longitude", "Latitude")


    P = qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
            axis.text.y = element_text(face="bold",size=25),
            axis.title.x = element_text(face="bold",size=25),
            axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



    meanPlot = Grid[,]
    meanPlot[,1] = Grid[,1]+tau[1]*(Mu[1,1]-Grid[,1])
    meanPlot[,2] = Grid[,2]+tau[1]*(Mu[2,1]-Grid[,2])

    P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
    #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                     # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

    j=1
    i = 1
    arr = Grid
    ang = angle[i]
    R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
    arr[i,] = arr[i,]+muPlot
    SigmaPlot = R%*%Sigma%*%t(R)
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = app
    for(i in 2:nrow(Grid))
    {
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        SigmaPlot = R%*%Sigma%*%t(R)
        arr[i,] = arr[i,]+muPlot
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = rbind(ell,NA, app)

    }
    ell_2 = ell
    arr_2 = arr
    muPlot_2 = arr
    for(imcmc in 2:nmcmc)
    {
        kkk =1
        for(k in WM)
        {
            mu0List[[kkk]] = matrix((ModelOUT$mu0[imcmc,1:2,k]),ncol=1)
            muCList[[kkk]] = matrix((ModelOUT$muC[imcmc,1:2,k]),ncol=1)
            psiList[[kkk]] = matrix((ModelOUT$psi[imcmc,1,k]),ncol=1)
            rhoList[[kkk]] = matrix((ModelOUT$rho[imcmc,1,k]),ncol=1)
            sigmaList[[kkk]] = matrix((ModelOUT$sigma[imcmc,1:4,k]),nrow=2)
            piList[[kkk]] = matrix((ModelOUT$pi[imcmc,1:length(WM),k]),ncol=length(WM))
            kkk = kkk+1

        }

        Mu = mu0List[[kk]]
        Nu = muCList[[kk]]
        tau = psiList[[kk]]


        rho = rhoList[[kk]]
        Sigma = sigmaList[[kk]]

        xseq = yseq = seq(-4,4,length.out=4)
        xseq = yseq = seq(-4,4, length.out=4)
        nnnn = length(xseq)
        angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

        Grid = expand.grid(xseq, yseq)
        colnames(Grid) = c("Longitude", "Latitude")


        # P = qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
        #         axis.text.y = element_text(face="bold",size=25),
        #         axis.title.x = element_text(face="bold",size=25),
        #         axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



        meanPlot = Grid[,]
        meanPlot[,1] = Grid[,1]+tau[1]*(Mu[1,1]-Grid[,1])
        meanPlot[,2] = Grid[,2]+tau[1]*(Mu[2,1]-Grid[,2])

        #P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
        #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                         # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

        j=1
        i = 1
        arr = Grid
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        arr[i,] = arr[i,]+muPlot
        SigmaPlot = R%*%Sigma%*%t(R)
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = app
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
            R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
            SigmaPlot = R%*%Sigma%*%t(R)
            arr[i,] = arr[i,]+muPlot
            app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
            app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
            app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
            ell = rbind(ell,NA, app)

        }
        ell_2 = ell_2+ell
        arr_2 = arr_2+arr
        muPlot_2 = muPlot_2+arr
    }
    ell_2 = ell_2/nmcmc
    arr_2 = arr_2/nmcmc
    muPlot_2 = muPlot_2/nmcmc




    # col_ell_1 = rep(cbPalette[1], nrow(ell_2))
    # col_ell_2 = rep(cbPalette[3], nrow(ell))
    # ell = rbind(ell,NA, ell_2)

    P4 = P1


    i = 1
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    arr2 = Grid
    muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
    arr2[i,] = muPlot
    for(i in 2:nrow(Grid))
    {
        ang = angle[i]
        R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot = Grid[i,]+(R)%*%matrix(c(-3/5,0), ncol=1)
        arr2[i,] = muPlot
    }
    #P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                      #arrow = arrow(length = unit(0.2, "cm")))
    #P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                    #arrow = arrow(length = unit(0.2, "cm")))


        P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2], color= "3"),arrow = arrow(length = unit(0.3, "cm")), size = 1.)
        P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                        arrow = arrow(length = unit(0.3, "cm")), size = 1.)

        PTOT =P7+scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+ theme(axis.text.y = element_text(face="bold",size=25),
        axis.text.x = element_text(face="bold",size=25),
        axis.title.x = element_text(face="bold",size=25),
        axis.title.y = element_text(face="bold",size=25),
        legend.text = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 1.5,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")
        PTOT


    pdf(paste(PLOT_DIRPLOT,kk,"MovMCMC.pdf",sep=""),width=7, height=7)
    print(PTOT+xlim(c(-5,5))+ylim(c(-5,5)))
    dev.off()




		#######
		P4 = P1


		i = 1
		ang = angle[i]
		R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
		arr2 = Grid
		muPlot = Grid[i,]+(R)%*%matrix(c(-5/5,0), ncol=1)
		arr2[i,] = muPlot
		for(i in 2:nrow(Grid))
		{
				ang = angle[i]
				R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
				muPlot = Grid[i,]+(R)%*%matrix(c(-5/5,0), ncol=1)
				arr2[i,] = muPlot
		}
		#P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
											#arrow = arrow(length = unit(0.2, "cm")))
		#P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
										#arrow = arrow(length = unit(0.2, "cm")))


		P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1, linetype=6, color= "black")
		P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(kk)],
										arrow = arrow(length = unit(0.3, "cm")), size = 1)

		PTOT =P7+
		#scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
		 theme(
				 axis.text.y = element_text(face="bold",size=25),
		axis.text.x = element_text(face="bold",size=25),
		axis.title.x = element_text(face="bold",size=25),
		axis.title.y = element_text(face="bold",size=25),
		legend.text = element_text(face="bold",size=25),
		legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 1.2,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

		PP = PTOT+xlim(c(-5,5))+ylim(c(-5,5))+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)
		#+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)

		PP

		PP2 = PP
		pdf(paste(PLOT_DIRPLOT,kk,"MovMCMC.pdf",sep=""),width=7, height=7)
    print(PP2)
    dev.off()

}
for(kk in 1:2)
{
    nmcmc = dim(ModelOUT$mu0)[1]

    imcmc = 1
    kkk =1
    for(k in WM)
    {
        mu0List[[kkk]] = matrix((ModelOUT$mu0[imcmc,1:2,k]),ncol=1)
        muCList[[kkk]] = matrix((ModelOUT$muC[imcmc,1:2,k]),ncol=1)
        psiList[[kkk]] = matrix((ModelOUT$psi[imcmc,1,k]),ncol=1)
        rhoList[[kkk]] = matrix((ModelOUT$rho[imcmc,1,k]),ncol=1)
        sigmaList[[kkk]] = matrix((ModelOUT$sigma[imcmc,1:4,k]),nrow=2)
        piList[[kkk]] = matrix((ModelOUT$pi[imcmc,1:length(WM),k]),ncol=length(WM))
        kkk = kkk+1

    }

    Mu = mu0List[[kk]]
    Nu = muCList[[kk]]
    tau = psiList[[kk]]


    rho = rhoList[[kk]]
    Sigma = sigmaList[[kk]]

    xseq = yseq = seq(-4,4,length.out=4)
    xseq = yseq = seq(-4,4, length.out=4)
    nnnn = length(xseq)
    angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

		if(kk==1)
		{
			xseq = yseq = seq(-0.1,0.1 ,length.out=2)
			xseq = yseq = seq(-0.1,0.1, length.out=2)
			nnnn = length(xseq)
			angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))
		}
		if(kk==2)
		{
			xseq = yseq = seq(-0.6,0.6 ,length.out=2)
			xseq = yseq = seq(-0.6,0.6, length.out=2)
			nnnn = length(xseq)
			angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))
		}
    Grid = expand.grid(xseq, yseq)
    colnames(Grid) = c("Longitude", "Latitude")


    P = qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
            axis.text.y = element_text(face="bold",size=25),
            axis.title.x = element_text(face="bold",size=25),
            axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



    meanPlot = Grid[,]
    meanPlot[,1] = Grid[,1]+tau[1]*(Mu[1,1]-Grid[,1])
    meanPlot[,2] = Grid[,2]+tau[1]*(Mu[2,1]-Grid[,2])

    P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
    #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                     # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

    j=1
    i = 1
    arr = Grid
    ang = angle[i]
    R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
    arr[i,] = arr[i,]+muPlot
    SigmaPlot = R%*%Sigma%*%t(R)
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = app
    for(i in 2:nrow(Grid))
    {
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        SigmaPlot = R%*%Sigma%*%t(R)
        arr[i,] = arr[i,]+muPlot
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = rbind(ell,NA, app)

    }
    ell_2 = ell
    arr_2 = arr
    muPlot_2 = arr
    for(imcmc in 2:nmcmc)
    {
        kkk =1
        for(k in WM)
        {
            mu0List[[kkk]] = matrix((ModelOUT$mu0[imcmc,1:2,k]),ncol=1)
            muCList[[kkk]] = matrix((ModelOUT$muC[imcmc,1:2,k]),ncol=1)
            psiList[[kkk]] = matrix((ModelOUT$psi[imcmc,1,k]),ncol=1)
            rhoList[[kkk]] = matrix((ModelOUT$rho[imcmc,1,k]),ncol=1)
            sigmaList[[kkk]] = matrix((ModelOUT$sigma[imcmc,1:4,k]),nrow=2)
            piList[[kkk]] = matrix((ModelOUT$pi[imcmc,1:length(WM),k]),ncol=length(WM))
            kkk = kkk+1

        }

        Mu = mu0List[[kk]]
        Nu = muCList[[kk]]
        tau = psiList[[kk]]


        rho = rhoList[[kk]]
        Sigma = sigmaList[[kk]]

        xseq = yseq = seq(-4,4,length.out=4)
        xseq = yseq = seq(-4,4, length.out=4)
        nnnn = length(xseq)
        angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

				if(kk==1)
				{
					xseq = yseq = seq(-0.1,0.1 ,length.out=2)
					xseq = yseq = seq(-0.1,0.1, length.out=2)
					nnnn = length(xseq)
					angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))
				}
				if(kk==2)
				{
					xseq = yseq = seq(-0.6,0.6 ,length.out=2)
					xseq = yseq = seq(-0.6,0.6, length.out=2)
					nnnn = length(xseq)
					angle = rev(rep(seq(0,-pi,by=-pi/(3)), each=1))
				}

        Grid = expand.grid(xseq, yseq)
        colnames(Grid) = c("Longitude", "Latitude")


        # P = qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
        #         axis.text.y = element_text(face="bold",size=25),
        #         axis.title.x = element_text(face="bold",size=25),
        #         axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



        meanPlot = Grid[,]
        meanPlot[,1] = Grid[,1]+tau[1]*(Mu[1,1]-Grid[,1])
        meanPlot[,2] = Grid[,2]+tau[1]*(Mu[2,1]-Grid[,2])

        #P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
        #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                         # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

        j=1
        i = 1
        arr = Grid
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        arr[i,] = arr[i,]+muPlot
        SigmaPlot = R%*%Sigma%*%t(R)
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = app
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
            R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot =  (1-rho[j])*tau[1]*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
            SigmaPlot = R%*%Sigma%*%t(R)
            arr[i,] = arr[i,]+muPlot
            app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
            app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
            app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
            ell = rbind(ell,NA, app)

        }
        ell_2 = ell_2+ell
        arr_2 = arr_2+arr
        muPlot_2 = muPlot_2+arr
    }
    ell_2 = ell_2/nmcmc
    arr_2 = arr_2/nmcmc
    muPlot_2 = muPlot_2/nmcmc




    # col_ell_1 = rep(cbPalette[1], nrow(ell_2))
    # col_ell_2 = rep(cbPalette[3], nrow(ell))
    # ell = rbind(ell,NA, ell_2)

    P4 = P1



    #P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                      #arrow = arrow(length = unit(0.2, "cm")))
    #P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                    #arrow = arrow(length = unit(0.2, "cm")))





		#######
		P4 = P1




		if(kk==1)
		{

			i = 1
			ang = angle[i]
			R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
			arr2 = Grid
			muPlot = Grid[i,]+(R)%*%matrix(c(-0.25/5,0), ncol=1)
			arr2[i,] = muPlot
			for(i in 2:nrow(Grid))
			{
					ang = angle[i]
					R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
					muPlot = Grid[i,]+(R)%*%matrix(c(-0.25/5,0), ncol=1)
					arr2[i,] = muPlot
			}
		}
		if(kk==2)
		{
			i = 1
	    ang = angle[i]
	    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
	    arr2 = Grid
	    muPlot = Grid[i,]+(R)%*%matrix(c(-1/5,0), ncol=1)
	    arr2[i,] = muPlot
	    for(i in 2:nrow(Grid))
	    {
	        ang = angle[i]
	        R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
	        muPlot = Grid[i,]+(R)%*%matrix(c(-1/5,0), ncol=1)
	        arr2[i,] = muPlot
	    }
		}
		#P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
											#arrow = arrow(length = unit(0.2, "cm")))
		#P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
										#arrow = arrow(length = unit(0.2, "cm")))


		P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1, linetype=6, color= "black")
		P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(kk)],
										arrow = arrow(length = unit(0.3, "cm")), size = 1)

		PTOT =P7+
		#scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
		 theme(
				 axis.text.y = element_text(face="bold",size=25),
		axis.text.x = element_text(face="bold",size=25),
		axis.title.x = element_text(face="bold",size=25),
		axis.title.y = element_text(face="bold",size=25),
		legend.text = element_text(face="bold",size=25),
		legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 1.2,color=cbPalette[c(kk)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")


		if(kk==1)
		{
			ccc = 0.2
		}
		if(kk==2)
		{
			ccc = 1.1
		}
		PP = PTOT+xlim(c(-ccc,ccc))+ylim(c(-ccc,ccc))
		#+geom_point(aes(DataZ$Longitude[as.numeric(DataZ$Cluster)==kk],DataZ$Latitude[as.numeric(DataZ$Cluster)==kk]),size=0.99)

		PP

		PP2 = PP
		pdf(paste(PLOT_DIRPLOT,kk,"DetMovMCMC.pdf",sep=""),width=7, height=7)
    print(PP2)
    dev.off()

}


#### #### #### #### #### #### ####
#### SIMULATIONS
#### #### #### #### #### #### ####


mu0List = list()
muCList = list()
psiList = list()
rhoList = list()
sigmaList = list()
piList = list()
kk = 1
for(k in WM)
{
    mu0List[[kk]] = matrix(colMeans(ModelOUT$mu0[,1:2,k]),ncol=1)
    muCList[[kk]] = matrix(colMeans(ModelOUT$muC[,1:2,k]),ncol=1)
    psiList[[kk]] = matrix(mean(ModelOUT$psi[,1,k]),ncol=1)
    rhoList[[kk]] = matrix(mean(ModelOUT$rho[,1,k]),ncol=1)
    sigmaList[[kk]] = matrix(colMeans(ModelOUT$sigma[,1:4,k]),nrow=2)
    piList[[kk]] = matrix(colMeans(ModelOUT$pi[,1:length(WM),k]),ncol=length(WM))
    kk = kk+1

}


set.seed(100)
for(iisim in 1:10)
{
    kprev = 1
    nsim_z = 1000
    z_sim = matrix(NA, ncol=2, nrow=nsim_z)
    z_sim[1,] = c(runif(1,-5,5),runif(1,-5,5))
    K = length(WM)
    theta = 0
    i = 1
    for(i in 1:(nsim_z-1))
    {
        k = sample(1:K, 1, prob=piList[[kprev]])
        R = matrix(c( cos(rhoList[[k]]*theta),sin(rhoList[[k]]*theta), -sin(rhoList[[k]]*theta), cos(rhoList[[k]]*theta) ), ncol=2)
        R2 = matrix(c( cos(theta),sin(theta), -sin(theta), cos(theta) ), ncol=2)

        mean = z_sim[i,]+psiList[[k]][1]*(1-rhoList[[k]][1])*(mu0List[[k]][]-z_sim[i,])
        mean = mean+rhoList[[k]][1]*R2%*%muCList[[k]][]

        z_sim[i+1,] = rmnorm(1,mean, R%*%sigmaList[[k]]%*%t(R))

        theta = atan2(z_sim[i+1,2]-z_sim[i,2],z_sim[i+1,1]-z_sim[i,1])
        kprev = k
    }
    z_sim = as.data.frame(z_sim)
    colnames(z_sim) = c("Longitude","Latitude")
    p = ggplot(z_sim,aes(x=Longitude,y=Latitude))+geom_point(color=cbPalette[iisim%%5+1])+ theme(axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_text(face="bold",size=25) )+ylab("Latitude")+xlab("Longitude")+xlim(c(-5,5))+ylim(c(-5,5))

    pdf(paste(PLOT_DIRPLOT,iisim,"sim.pdf",sep=""),width=7, height=7)
    print(p)
    dev.off()
}


#### #### #### #### #### #### #### #### ####
#### OU model
#### #### #### #### #### #### #### #### ####


load(paste(PLOT_DIRDATA , MOD_BRW_NAME, ".Rdata",sep=""))
WM = c(39,136,14,127,23)
WM0 = WM


Z_MAP = apply(ModelOUT$zeta,2,findmode)
DIFF = matrix(NA, ncol=2, nrow=nrow(DataCoords)-1)
DIFF[,1] = diff(DataCoords[,1])
DIFF[,2] = diff(DataCoords[,2])


print("")
print("")
print("")
print("Posterior Estimates - OU MODEL")
print("")

for(i in 1:2)
{
  dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$mu0[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu0[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu0[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}

for(i in 1:1)
{
  dd = paste("$backslash tau_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$psi[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$psi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$psi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}


for(i in 1:2)
{
  for(j in i:2)
  {
    dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }


}
ddc = c()
ik=1
for(k in WM)
{
  dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
  ff = "(CI) "
  ii = 1
  for(i in WM)
  {
      ddc[(ik-1)*length(WM)+ii] = round(mean(ModelOUT$pi[,i,k]),3)
    dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    ii = ii+1
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
  ik = ik+1
}

ik = 1
for(i in 1:1)
{
  dd = paste("$backslash beta_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}
gammaC = ModelOUT$gammaDP
kappaC = ModelOUT$rhoDP*ModelOUT$akDP
alphaC = ModelOUT$akDP-kappaC
print("backslash hline backslash hline ")
dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
dd = paste("$backslash hat{}$ ",sep="")
ff = "(CI) "

dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd1 = paste(dd1, " \\")
dd = paste(dd, " \\")
ff = paste(ff, " \\")
print(dd1)
print("backslash hline  ")
print(dd)
print(ff)

print("Distribution of K")
table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)


AppWM = rep(0,max(WM))
AppWM[WM] = 1:5

Z_MAP2_0 = AppWM[Z_MAP]




#### #### #### #### #### #### ####
#### SIMULATIONS
#### #### #### #### #### #### ####


mu0List = list()
muCList = list()
psiList = list()
rhoList = list()
sigmaList = list()
piList = list()
kk = 1
for(k in WM)
{
    mu0List[[kk]] = matrix(colMeans(ModelOUT$mu0[,1:2,k]),ncol=1)
    muCList[[kk]] = matrix(colMeans(ModelOUT$muC[,1:2,k]),ncol=1)
    psiList[[kk]] = matrix(mean(ModelOUT$psi[,1,k]),ncol=1)
    rhoList[[kk]] = matrix(mean(ModelOUT$rho[,1,k]),ncol=1)
    sigmaList[[kk]] = matrix(colMeans(ModelOUT$sigma[,1:4,k]),nrow=2)
    piList[[kk]] = matrix(colMeans(ModelOUT$pi[,1:length(WM),k]),ncol=length(WM))
    kk = kk+1

}


set.seed(100)
for(iisim in 1:10)
{
    kprev = 1
    nsim_z = 1000
    z_sim = matrix(NA, ncol=2, nrow=nsim_z)
    z_sim[1,] = c(runif(1,-5,5),runif(1,-5,5))
    K = length(WM)
    theta = 0
    i = 1
    for(i in 1:(nsim_z-1))
    {
        k = sample(1:K, 1, prob=piList[[kprev]])
        R = matrix(c( cos(rhoList[[k]]*theta),sin(rhoList[[k]]*theta), -sin(rhoList[[k]]*theta), cos(rhoList[[k]]*theta) ), ncol=2)
        R2 = matrix(c( cos(theta),sin(theta), -sin(theta), cos(theta) ), ncol=2)

        mean = z_sim[i,]+psiList[[k]][1]*(1-rhoList[[k]][1])*(mu0List[[k]][]-z_sim[i,])
        mean = mean+rhoList[[k]][1]*R2%*%muCList[[k]][]

        z_sim[i+1,] = rmnorm(1,mean, R%*%sigmaList[[k]]%*%t(R))

        theta = atan2(z_sim[i+1,2]-z_sim[i,2],z_sim[i+1,1]-z_sim[i,1])
        kprev = k
    }
    z_sim = as.data.frame(z_sim)
    colnames(z_sim) = c("Longitude","Latitude")
    p = ggplot(z_sim,aes(x=Longitude,y=Latitude))+geom_point(color=cbPalette[iisim%%5+1])+ theme(axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_text(face="bold",size=25) )+ylab("Latitude")+xlab("Longitude")+xlim(c(-5,5))+ylim(c(-5,5))

    pdf(paste(PLOT_DIRPLOT,iisim,"simBRW.pdf",sep=""),width=7, height=7)
    print(p)
    dev.off()
}




#### #### #### #### #### #### #### #### ####
#### ST
#### #### #### #### #### #### #### #### ####

load(paste(PLOT_DIRDATA , MOD_CRW_NAME, ".Rdata",sep=""))
WM = c(14,160,28)


WM1 = WM

Z_MAP = apply(ModelOUT$zeta,2,findmode)
DIFF = matrix(NA, ncol=2, nrow=nrow(DataCoords)-1)


print("")
print("")
print("")
print("Posterior Estimates - ST MODEL")
print("")

for(i in 1:2)
{
  dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$muC[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$muC[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$muC[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}



for(i in 1:2)
{
  for(j in i:2)
  {
    dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }


}

ik=1
for(k in WM)
{
  dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
  ff = "(CI) "
  ii = 1
  for(i in WM)
  {
      ddc[(ik-1)*length(WM)+ii] = round(mean(ModelOUT$pi[,i,k]),3)
    dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    ii = ii+1
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
  ik = ik+1
}


ik = 1
for(i in 1:1)
{
  dd = paste("$backslash beta_{j}$ ",sep="")
  ff = "(CI) "
  for(k in WM)
  {
    dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
    ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
  }
  dd = paste(dd, " \\")
  ff = paste(ff, " \\")
  print(dd)
  print(ff)
}
gammaC = ModelOUT$gammaDP
kappaC = ModelOUT$rhoDP*ModelOUT$akDP
alphaC = ModelOUT$akDP-kappaC
print("backslash hline backslash hline ")
dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
dd = paste("$backslash hat{}$ ",sep="")
ff = "(CI) "

dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
dd1 = paste(dd1, " \\")
dd = paste(dd, " \\")
ff = paste(ff, " \\")
print(dd1)
print("backslash hline  ")
print(dd)
print(ff)



print("Distribution of K")
table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)



AppWM = rep(0,max(WM))
AppWM[WM] = 1:3

Z_MAP2_1 = AppWM[Z_MAP]



DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2],Cluster = as.factor(Z_MAP2_1) )



for(i in 1:3)
{
  p = ggplot(subset(DataZ, Cluster %in% c(paste(i))), aes(Longitude, Latitude, shape = Cluster))
  p = p +geom_point(size = 2.5, shape=16,color=cbPalette[i])
  p = p+theme(
    axis.text.x = element_text(face="bold",size=25),
    axis.text.y = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )+xlim(min(DataZ[,1], na.rm=T), max(DataZ[,1], na.rm=T))+ylim(min(DataZ[,2], na.rm=T), max(DataZ[,2], na.rm=T))
  p
  pdf(paste(PLOT_DIRPLOT ,"DataPost_ST",i,".pdf",sep=""))
  print(p+xlim(c(-5,5))+ylim(c(-5,5)))
  dev.off()
}



Densmat_2 = list()
Densmat_2[[1]] = matrix(0, ncol=nseq,nseq)
Densmat_2[[2]] = matrix(0, ncol=nseq,nseq)
Densmat_2[[3]] = matrix(0, ncol=nseq,nseq)
Densmat_2[[4]] = matrix(0, ncol=nseq,nseq)
Densmat_2[[5]] = matrix(0, ncol=nseq,nseq)
for(imcmc in 1:nmcmc)
{
  k1 = 0
  for(k in WM)
  {
    k1 = k1+1
    mu = ModelOUT$muC[imcmc,,k]
    sigma = matrix(ModelOUT$sigma[imcmc,,k],2)

    Densmat_2[[k1]][,] = Densmat_2[[k1]][,]+matrix(dmnorm(Vecy,mu, sigma)*Rvec^2, ncol=nseq)/nmcmc

  }

}





DataTR = data.frame(Theta = c(rowSums(Densmat_2[[1]]*diff2),rowSums(Densmat_2[[2]]*diff2),rowSums(Densmat_2[[3]]*diff2)),R = c(colSums(Densmat_2[[1]]*diff1),colSums(Densmat_2[[2]]*diff1),colSums(Densmat_2[[3]]*diff1)) , Gruppo = as.factor(rep(c(1,2,3),each =nseq )),thetaseq = rep(thetaseq, times=3),rseq = rep(rseq, times=3))


p = ggplot(DataTR,aes(x=thetaseq, y =Theta,group = Gruppo   ))+
geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
) +ylim(0,0.300)+xlab("Turning-Angle")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
 p


pdf(paste(PLOT_DIRPLOT ,"Turn2.pdf",sep=""))
print(p)
dev.off()

p = ggplot(DataTR,aes(x=rseq, y =R,group = Gruppo   ))+
geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
) +ylim(0,0.8)+xlab("Log Step-Length")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p

pdf(paste(PLOT_DIRPLOT ,"Step2.pdf",sep=""))
print(p)
dev.off()

#### #### #### #### #### #### #### #### ####
#### TIme of the day and mosaic plots
#### #### #### #### #### #### #### #### ####

table(Z_MAP2_0,Z_MAP2_gen)
table(Z_MAP2_1,Z_MAP2_gen)

DataZ = data.frame("STAP"=as.factor(Z_MAP2_gen),"OU"=as.factor(Z_MAP2_0),"ST"=as.factor(Z_MAP2_1))

p = ggplot(data = DataZ) +
   geom_mosaic(aes(x = product(STAP,OU), fill=STAP))  +guides(fill=guide_legend(title="STAP-HMM"))+theme(
     axis.text.y = element_blank(),
     axis.text.x = element_text(face="bold",size=25),
     axis.title.x = element_text(face="bold",size=25),
     axis.title.y = element_text(face="bold",size=25),
     legend.text = element_text(face="bold",size=25),
     legend.title = element_text(face="bold",size=25)
 )+scale_x_productlist("BRW-HMM")+scale_y_productlist("STAP-HMM")+scale_fill_manual(values=cbPalette)
p

pdf(paste(PLOT_DIRPLOT ,"Behav_05_0.pdf",sep=""))
print(p)
dev.off()

p = ggplot(data = DataZ) +
   geom_mosaic(aes(x = product(STAP,ST), fill=STAP))  +guides(fill=guide_legend(title="STAP-HMM"))+theme(
     axis.text.y = element_blank(),
     axis.text.x = element_text(face="bold",size=25),
     axis.title.x = element_text(face="bold",size=25),
     axis.title.y = element_text(face="bold",size=25),
     legend.text = element_text(face="bold",size=25),
     legend.title = element_text(face="bold",size=25)
 )+scale_x_productlist("CRW-HMM")+scale_y_productlist("STAP-HMM")+scale_fill_manual(values=cbPalette)
p
pdf(paste(PLOT_DIRPLOT ,"Behav_05_1.pdf",sep=""))
print(p)
dev.off()




DataObs  = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2],"Prop"= DataZ[,1],"rho0"=DataZ[,2],"rho1"=DataZ[,3], Time = 1:(nrow(DataCoords)-1)+18*2 )
DataObs$Time2 = DataObs$Time%%48

TT1 = table(DataObs$Prop,DataObs$Time2)/matrix(colSums(table(DataObs$Prop,DataObs$Time2)),byrow=T,ncol=48, nrow=5)
TT2 = table(DataObs[,"rho0"],DataObs$Time2)/matrix(colSums(table(DataObs[,"rho0"],DataObs$Time2)),byrow=T,ncol=48, nrow=5)
TT3 = table(DataObs[,"rho1"],DataObs$Time2)/matrix(colSums(table(DataObs[,"rho1"],DataObs$Time2)),byrow=T,ncol=48, nrow=3)

DataT = as.data.frame(TT1)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_05.pdf",sep=""))
print(p)
dev.off()



DataT = as.data.frame(TT2)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_0.pdf",sep=""))
print(p)
dev.off()



DataT = as.data.frame(TT3)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_1.pdf",sep=""))
print(p)
dev.off()


DataObs  = data.frame(Longitude = DataCoords[-nrow(DataCoords),1], Latitude = DataCoords[-nrow(DataCoords),2],"Prop"= DataZ[,1],"rho0"=DataZ[,2],"rho1"=DataZ[,3], Time = 1:(nrow(DataCoords)-1)+18*2 )
DataObs$Time2 = DataObs$Time%%48

TT1 = table(DataObs$Prop,DataObs$Time2)/matrix(rowSums(table(DataObs$Prop,DataObs$Time2)),byrow=F,ncol=48, nrow=5)
TT2 = table(DataObs[,"rho0"],DataObs$Time2)/matrix(rowSums(table(DataObs[,"rho0"],DataObs$Time2)),byrow=F,ncol=48, nrow=5)
TT3 = table(DataObs[,"rho1"],DataObs$Time2)/matrix(rowSums(table(DataObs[,"rho1"],DataObs$Time2)),byrow=F,ncol=48, nrow=3)

DataT = as.data.frame(TT1)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_05_V2.pdf",sep=""))
print(p)
dev.off()



DataT = as.data.frame(TT2)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_0_V2.pdf",sep=""))
print(p)
dev.off()



DataT = as.data.frame(TT3)
p = ggplot(DataT, aes(x=Var2,y = Freq, group=Var1))+
geom_line(aes(color=Var1),size=1.5)+
theme(aspect.ratio = 1/2,
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10),
  axis.title.x = element_text(face="bold",size=10),
  axis.title.y = element_text(face="bold",size=10),
  legend.text = element_text(face="bold",size=10),
  legend.title = element_text(face="bold",size=10)
)+xlab("Time")+ylab("")+ scale_x_discrete(breaks=(0:4)*12, labels=c(
 "00:00",
 "06:00",
 "12:00",
 "18:00",
 "24:00"
))+ labs(color='Behavior')+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
p
pdf(paste(PLOT_DIRPLOT ,"TimeBehav_1_V2.pdf",sep=""))
print(p)
dev.off()




#### #### #### #### #### #### ####
#### SIMULATIONS
#### #### #### #### #### #### ####


mu0List = list()
muCList = list()
psiList = list()
rhoList = list()
sigmaList = list()
piList = list()
kk = 1
for(k in WM)
{
    mu0List[[kk]] = matrix(colMeans(ModelOUT$mu0[,1:2,k]),ncol=1)
    muCList[[kk]] = matrix(colMeans(ModelOUT$muC[,1:2,k]),ncol=1)
    psiList[[kk]] = matrix(mean(ModelOUT$psi[,1,k]),ncol=1)
    rhoList[[kk]] = matrix(mean(ModelOUT$rho[,1,k]),ncol=1)
    sigmaList[[kk]] = matrix(colMeans(ModelOUT$sigma[,1:4,k]),nrow=2)
    piList[[kk]] = matrix(colMeans(ModelOUT$pi[,1:length(WM),k]),ncol=length(WM))
    kk = kk+1

}


set.seed(100)
for(iisim in 1:10)
{
    kprev = 1
    nsim_z = 1000
    z_sim = matrix(NA, ncol=2, nrow=nsim_z)
    z_sim[1,] = c(runif(1,-5,5),runif(1,-5,5))
    K = length(WM)
    theta = 0
    i = 1
    for(i in 1:(nsim_z-1))
    {
        k = sample(1:K, 1, prob=piList[[kprev]])
        R = matrix(c( cos(rhoList[[k]]*theta),sin(rhoList[[k]]*theta), -sin(rhoList[[k]]*theta), cos(rhoList[[k]]*theta) ), ncol=2)
        R2 = matrix(c( cos(theta),sin(theta), -sin(theta), cos(theta) ), ncol=2)

        mean = z_sim[i,]+psiList[[k]][1]*(1-rhoList[[k]][1])*(mu0List[[k]][]-z_sim[i,])
        mean = mean+rhoList[[k]][1]*R2%*%muCList[[k]][]

        z_sim[i+1,] = rmnorm(1,mean, R%*%sigmaList[[k]]%*%t(R))

        theta = atan2(z_sim[i+1,2]-z_sim[i,2],z_sim[i+1,1]-z_sim[i,1])
        kprev = k
    }
    z_sim = as.data.frame(z_sim)
    colnames(z_sim) = c("Longitude","Latitude")
    p = ggplot(z_sim,aes(x=Longitude,y=Latitude))+geom_point(color=cbPalette[iisim%%5+1])+ theme(axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_text(face="bold",size=25) )+ylab("Latitude")+xlab("Longitude")

    pdf(paste(PLOT_DIRPLOT,iisim,"CRWsim.pdf",sep=""),width=7, height=7)
    print(p)
    dev.off()
}


#
# #### #### #### #### #### #### #### #### ####
# #### Observed turning-angle for the first two
# #### OU-HMM behaviors
# #### #### #### #### #### #### #### #### ####
#
# load(paste(PLOT_DIRDATA , MOD_STAP_NAME,sep=""))
# WM = c(14,164,71,23,26)
# WMgen = WM
#
# MAP = apply(ModelOUT$zeta,2,findmode)
# INDEX_FIRST = which(MAP==14)
# INDEX_SECOND = which(MAP==164)
# length(INDEX_FIRST)
#
#
#
# DataCoords_2 = as.data.frame(DataCoords)
# colnames(DataCoords_2) = c("x","y")
# dataDog  <- prepData(DataCoords_2,type="UTM")
#
# ###  usare l$y*3 anche nelle figure dei dati osservati
# dataDog2 = rbind(dataDog[INDEX_FIRST,-1] -2*pi,dataDog[INDEX_FIRST,-1],dataDog[INDEX_FIRST,-1]+2*pi)
# l <- density(dataDog2$angle,na.rm=T, bw=0.3)
# ld =data.frame(x=l$x,y=l$y*3)
# P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
#   axis.text.x = element_text(face="bold",size=25),
#   axis.text.y = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25)
# )+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
# P
# pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k1.pdf",sep=""))
# print(P)
# dev.off()
#
#
# dataDog2 = rbind(dataDog[INDEX_SECOND,-1] -2*pi,dataDog[INDEX_SECOND,-1],dataDog[INDEX_SECOND,-1]+2*pi)
# l <- density(dataDog2$angle,na.rm=T, bw=0.3)
# ld =data.frame(x=l$x,y=l$y*3)
# P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
#   axis.text.x = element_text(face="bold",size=25),
#   axis.text.y = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25)
# )+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
# P
# pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k2.pdf",sep=""))
# print(P)
# dev.off()
#
#
#
# load(paste(PLOT_DIRDATA , MOD_BRW_NAME,sep=""))
# WM = c(39,136,14,127,23)
#
# s1    = mean(ModelOUT$sigma[,1,WM[1]])
# s12   = mean(ModelOUT$sigma[,2,WM[1]])
# s2    = mean(ModelOUT$sigma[,4,WM[1]])
#
# X = rmnorm(50000, c(0,0), matrix(c(s1,s12,s12,s2), ncol=2))
# X = apply(X,2,cumsum)
# colnames(X) = c("x","y")
# X = data.frame(X)
# dataDog  = prepData(X,type="UTM", LLangle=F)
#
#
# # DataCoords_2 = as.data.frame(DataCoords)
# # colnames(DataCoords_2) = c("x","y")
# # dataDog  <- prepData(DataCoords_2,type="UTM",LLangle=F)
#
# ###  usare l$y*3 anche nelle figure dei dati osservati
# dataDog2 = rbind(dataDog[,-1] -2*pi,dataDog[,-1],dataDog[,-1]+2*pi)
# l <- density(dataDog2$angle,na.rm=T, bw=0.3)
# ld =data.frame(x=l$x,y=l$y*3)
# P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
#   axis.text.x = element_text(face="bold",size=25),
#   axis.text.y = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25)
# )+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
# P
# pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k1_OU.pdf",sep=""))
# print(P)
# dev.off()
#
#
#
#
# s1    = mean(ModelOUT$sigma[,1,WM[2]])
# s12   = mean(ModelOUT$sigma[,2,WM[2]])
# s2    = mean(ModelOUT$sigma[,4,WM[2]])
#
# X = rmnorm(50000, c(0,0), matrix(c(s1,s12,s12,s2), ncol=2))
# X = apply(X,2,cumsum)
# colnames(X) = c("x","y")
# X = data.frame(X)
# dataDog  = prepData(X,type="UTM", LLangle=F)
#
#
# # DataCoords_2 = as.data.frame(DataCoords)
# # colnames(DataCoords_2) = c("x","y")
# # dataDog  <- prepData(DataCoords_2,type="UTM",LLangle=F)
#
# ###  usare l$y*3 anche nelle figure dei dati osservati
# dataDog2 = rbind(dataDog[,-1] -2*pi,dataDog[,-1],dataDog[,-1]+2*pi)
# l <- density(dataDog2$angle,na.rm=T, bw=0.3)
# ld =data.frame(x=l$x,y=l$y*3)
# P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
#   axis.text.x = element_text(face="bold",size=25),
#   axis.text.y = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25)
# )+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
# P
# pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k2_OU.pdf",sep=""))
# print(P)
# dev.off()



# #### #### #### #### #### #### #### #### ####
# #### SIMULATIONS
# #### #### #### #### #### #### #### #### ####
#
# cbPalette <- c(
# "#e41a1c",
# "#377eb8",
# "#4daf4a",
# "#984ea3",
# "#ff7f00"
# )
#
#
# isim = 1
# load(paste(PLOT_DIRDATA ,MOD_SIM_NAME ,isim,"ModelOUT.Rdata",sep=""))
#
#
# DataZ = data.frame(Longitude = DataCoords[,1], Latitude = DataCoords[,2] , Behavior = as.factor(LatentClassification))
#
# p = ggplot(DataZ, aes(x=Longitude, y=Latitude,group=Behavior))
# p = p+    geom_point(aes(shape=Behavior, color=Behavior),size = 1.5)
#
# p = p+theme(
#   axis.text.y = element_text(face="bold",size=25),
#   axis.text.x = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25),
#   legend.text = element_text(face="bold",size=25),
#   legend.title = element_text(face="bold",size=25)
# )+scale_color_manual(values=cbPalette[c(1,2,3)])  +guides(shape = guide_legend(override.aes = list(size = 15)), colour = guide_legend(override.aes = list(size=5)))
# p
#
# pdf(paste(PLOT_DIRPLOT ,"Sim",isim,".pdf",sep=""))
# print(p)
# dev.off()
#
# WM = c(76,186,121)
# WMgen = WM
#
# Z_MAP = apply(ModelOUT$zeta,2,findmode)
#
# ### Posterior estimates
# print("")
# print("")
# print("")
# print("Posterior Estimates - SIM1")
# print("")
# for(i in 1:2)
# {
#   dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$mu0[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu0[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu0[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:2)
# {
#   dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$muC[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$muC[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$muC[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:1)
# {
#   dd = paste("$backslash tau_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$psi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$psi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$psi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
# for(i in 1:1)
# {
#   dd = paste("$backslash rho_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$rho[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$rho[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$rho[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
#
# for(i in 1:2)
# {
#   for(j in i:2)
#   {
#     dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
#     ff = "(CI) "
#     for(k in WM)
#     {
#       dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
#       ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     }
#     dd = paste(dd, " \\")
#     ff = paste(ff, " \\")
#     print(dd)
#     print(ff)
#   }
#
#
# }
#
# ik=1
# for(k in WM)
# {
#   dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
#   ff = "(CI) "
#   ii = 1
#   for(i in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     ii = ii+1
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
#   ik = ik+1
# }
# ik = 1
# for(i in 1:1)
# {
#   dd = paste("$backslash beta_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# gammaC = ModelOUT$gammaDP
# kappaC = ModelOUT$rhoDP*ModelOUT$akDP
# alphaC = ModelOUT$akDP-kappaC
# print("backslash hline backslash hline ")
# dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
# dd = paste("$backslash hat{}$ ",sep="")
# ff = "(CI) "
#
# dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd1 = paste(dd1, " \\")
# dd = paste(dd, " \\")
# ff = paste(ff, " \\")
# print(dd1)
# print("backslash hline  ")
# print(dd)
# print(ff)
#
#
#
#
# print("Distribution of K")
# table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)
#
# print("Confusion Matrix")
# sum(diag(table(LatentClassification[-length(LatentClassification)],Z_MAP)[order(WM),])/length(Z_MAP))
#
#
#
# isim = 2
# load(paste(PLOT_DIRDATA ,MOD_SIM_NAME ,isim,"ModelOUT.Rdata",sep=""))
#
#
# DataZ = data.frame(Longitude = DataCoords[,1], Latitude = DataCoords[,2] , Behavior = as.factor(LatentClassification))
#
# p = ggplot(DataZ, aes(x=Longitude, y=Latitude,group=Behavior))
# p = p+    geom_point(aes(shape=Behavior, color=Behavior),size = 1.5)
#
# p = p+theme(
#   axis.text.y = element_text(face="bold",size=25),
#   axis.text.x = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25),
#   legend.text = element_text(face="bold",size=25),
#   legend.title = element_text(face="bold",size=25)
# )+scale_color_manual(values=cbPalette[c(1,2,3)])  +guides(shape = guide_legend(override.aes = list(size = 15)), colour = guide_legend(override.aes = list(size=5)))
# p
#
# pdf(paste(PLOT_DIRPLOT ,"Sim",isim,".pdf",sep=""))
# print(p)
# dev.off()
#
# WM = c(27,173,190)
# WMgen = WM
#
# Z_MAP = apply(ModelOUT$zeta,2,findmode)
#
#
# ### Posterior estimates
# print("")
# print("")
# print("")
# print("Posterior Estimates - SIM2")
# print("")
# for(i in 1:2)
# {
#   dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$mu0[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu0[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu0[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:2)
# {
#   dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$muC[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$muC[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$muC[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:1)
# {
#   dd = paste("$backslash tau_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$psi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$psi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$psi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
# for(i in 1:1)
# {
#   dd = paste("$backslash rho_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$rho[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$rho[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$rho[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
#
# for(i in 1:2)
# {
#   for(j in i:2)
#   {
#     dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
#     ff = "(CI) "
#     for(k in WM)
#     {
#       dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
#       ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     }
#     dd = paste(dd, " \\")
#     ff = paste(ff, " \\")
#     print(dd)
#     print(ff)
#   }
#
#
# }
#
# ik=1
# for(k in WM)
# {
#   dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
#   ff = "(CI) "
#   ii = 1
#   for(i in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     ii = ii+1
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
#   ik = ik+1
# }
# ik = 1
# for(i in 1:1)
# {
#   dd = paste("$backslash beta_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# gammaC = ModelOUT$gammaDP
# kappaC = ModelOUT$rhoDP*ModelOUT$akDP
# alphaC = ModelOUT$akDP-kappaC
# print("backslash hline backslash hline ")
# dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
# dd = paste("$backslash hat{}$ ",sep="")
# ff = "(CI) "
#
# dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd1 = paste(dd1, " \\")
# dd = paste(dd, " \\")
# ff = paste(ff, " \\")
# print(dd1)
# print("backslash hline  ")
# print(dd)
# print(ff)
#
#
#
# dd = paste("$backslash beta_{j}$ ",sep="")
# ff = "(CI) "
#
#
# print("Distribution of K")
# table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)
#
# print("Confusion Matrix")
# sum(diag(table(LatentClassification[-length(LatentClassification)],Z_MAP)[order(WM),])/length(Z_MAP))
#
#
#
# str(ModelOUT)
#
#
#
# isim = 3
# load(paste(PLOT_DIRDATA ,MOD_SIM_NAME ,isim,"ModelOUT.Rdata",sep=""))
#
#
# DataZ = data.frame(Longitude = DataCoords[,1], Latitude = DataCoords[,2] , Behavior = as.factor(LatentClassification))
#
# p = ggplot(DataZ, aes(x=Longitude, y=Latitude,group=Behavior))
# p = p+    geom_point(aes(shape=Behavior, color=Behavior),size = 1.5)
#
# p = p+theme(
#   axis.text.y = element_text(face="bold",size=25),
#   axis.text.x = element_text(face="bold",size=25),
#   axis.title.x = element_text(face="bold",size=25),
#   axis.title.y = element_text(face="bold",size=25),
#   legend.text = element_text(face="bold",size=25),
#   legend.title = element_text(face="bold",size=25)
# )+scale_color_manual(values=cbPalette[c(1,2,3)])  +guides(shape = guide_legend(override.aes = list(size = 15)), colour = guide_legend(override.aes = list(size=5)))
# p
#
# pdf(paste(PLOT_DIRPLOT ,"Sim",isim,".pdf",sep=""))
# print(p)
# dev.off()
#
# WM = c(65,192,150)
# WMgen = WM
#
# Z_MAP = apply(ModelOUT$zeta,2,findmode)
#
#
# ### Posterior estimates
# print("")
# print("")
# print("")
# print("Posterior Estimates - SIM3")
# print("")
# for(i in 1:2)
# {
#   dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$mu0[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu0[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu0[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:2)
# {
#   dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$muC[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$muC[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$muC[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# for(i in 1:1)
# {
#   dd = paste("$backslash tau_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$psi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$psi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$psi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
# for(i in 1:1)
# {
#   dd = paste("$backslash rho_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$rho[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$rho[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$rho[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
#
#
# for(i in 1:2)
# {
#   for(j in i:2)
#   {
#     dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
#     ff = "(CI) "
#     for(k in WM)
#     {
#       dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
#       ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     }
#     dd = paste(dd, " \\")
#     ff = paste(ff, " \\")
#     print(dd)
#     print(ff)
#   }
#
#
# }
#
# ik=1
# for(k in WM)
# {
#   dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
#   ff = "(CI) "
#   ii = 1
#   for(i in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#     ii = ii+1
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
#   ik = ik+1
# }
# ik = 1
# for(i in 1:1)
# {
#   dd = paste("$backslash beta_{j}$ ",sep="")
#   ff = "(CI) "
#   for(k in WM)
#   {
#     dd = paste(dd, " & ", round(mean(ModelOUT$betaDP[,i,k]),3),sep = )
#     ff = paste(ff, " & ","(",round(quantile(ModelOUT$betaDP[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$betaDP[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
#   }
#   dd = paste(dd, " \\")
#   ff = paste(ff, " \\")
#   print(dd)
#   print(ff)
# }
# gammaC = ModelOUT$gammaDP
# kappaC = ModelOUT$rhoDP*ModelOUT$akDP
# alphaC = ModelOUT$akDP-kappaC
# print("backslash hline backslash hline ")
# dd1 = paste("& $backslash alpha$ & $backslash kappa$  & $backslash gamma$ ",sep="")
# dd = paste("$backslash hat{}$ ",sep="")
# ff = "(CI) "
#
# dd = paste(dd, " & ", round(mean(alphaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(alphaC[,1], probs=c(0.025)),3), " ",round(quantile(alphaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(kappaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(kappaC[,1], probs=c(0.025)),3), " ",round(quantile(kappaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd = paste(dd, " & ", round(mean(gammaC[,1]),3),sep = )
# ff = paste(ff, " & ","(",round(quantile(gammaC[,1], probs=c(0.025)),3), " ",round(quantile(gammaC[,1], probs=c(1-0.025)),3)  ,   ")",sep="")
# dd1 = paste(dd1, " \\")
# dd = paste(dd, " \\")
# ff = paste(ff, " \\")
# print(dd1)
# print("backslash hline  ")
# print(dd)
# print(ff)
#
#
# print("Distribution of K")
# table(apply(ModelOUT$zeta,1,function(x) length(unique(x))))/nrow(ModelOUT$zeta)
#
# print("Confusion Matrix")
# sum(diag(table(LatentClassification[-length(LatentClassification)],Z_MAP)[order(WM),])/length(Z_MAP))
#
#
# ### beta
#
# str(mod)
# #ret = function(mod = ModelOUT)
# #{
#     mod = ModelOUT
# #}
