#### #### #### #### #### #### #### ####
####  Values that must be changed
#### #### #### #### #### #### #### ####

PLOT_DIRDATA  = ""
PLOT_DIRPLOT  = ""
MOD_STAP_NAME = ""
MOD_OU_NAME   = ""
MOD_ST_NAME   = ""

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
CoordAy1 =  c(0,1)%*%RotMat
CoordAy2 =  c(0,-1)%*%RotMat

PT = qplot(c(-10),10,geom="line", xlab="",ylab="", ylim=c(-1,1), xlim=c(-1,1))+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")

PT1 = PT + geom_segment(aes(x = CoordAx1[1], y = CoordAx1[2], xend = CoordAx2[1], yend = CoordAx2[2]), linetype =  2)+geom_segment(aes(x = CoordAy1[1], y = CoordAy1[2], xend = CoordAy2[1], yend = CoordAy2[2]), linetype =  2)

OBS = matrix(c(-0.8,0,0.6,0,0,0.3),nrow=3,ncol=2)%*%RotMat

PT2 = PT1+geom_line(aes(x=OBS[,1],y=OBS[,2]))+geom_point(aes(x=OBS[,1],y=OBS[,2],size=1))

Desx1 = rbind(c(0.6,0.3),cbind(0.6,0))%*%RotMat
Desx2 = rbind(c(0.6,0.3),cbind(0,0.3))%*%RotMat
PT3 = PT2+geom_line(aes(x=Desx1[,1],y=Desx1[,2]), linetype = 3)+geom_line(aes(x=Desx2[,1],y=Desx2[,2]), linetype = 3)


Xc 		= rbind(c(0.6,0.02),c(0.6,-0.02))
Yc 		= rbind(c(-0.02,0.3),c(0.02,0.3))
Xc 		= Xc%*%RotMat
Yc 		= Yc%*%RotMat

PT4 = PT3+geom_line(aes(x=Xc[,1],y=Xc[,2]), linetype = 1,size=1)+geom_line(aes(x=Yc[,1],y=Yc[,2]), linetype = 1,size=1)


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


cangle = seq(0.38,pi/2-0.7,length.out=100)
S1 = sin(cangle)*0.2
C1 = cos(cangle)*0.2

PT7 = PT6 +geom_line(aes(x=C1,y=S1), linetype = 1,size=0.4)+geom_text(aes(x=0.23,y=0.15), label=deparse(bquote(paste(theta['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=0.2,y=0.32), label=deparse(bquote(paste('r'['t'['i']]))),parse=TRUE,size=7)
PT7
c(-0.8,0,0.6,0,0,0.3)

PT8 = PT7+geom_segment(aes(x = OBS[1,1], y = OBS[1,2], xend = 1, yend = OBS[1,2]), linetype =  3)+
geom_segment(aes(x = OBS[2,1], y = OBS[2,2], xend = 1, yend = OBS[2,2]), linetype =  3)


cangle = seq(0.38,0,length.out=100)
S11 = sin(cangle)*0.2
C11 = cos(cangle)*0.2

cangle = seq(0.85,0,length.out=100)
S12 = sin(cangle)*0.4
C12 = cos(cangle)*0.4

PT9 = PT8+geom_line(aes(x=C11+OBS[1,1],y=S11+OBS[1,2]), linetype = 1,size=0.4)+
geom_line(aes(x=C12+OBS[2,1],y=S12+OBS[2,2]), linetype = 1,size=0.4)+
geom_text(aes(x=0.5,y=0.05), label=deparse(bquote(paste(phi['t'['i+1']]))),parse=TRUE,size=7)+
geom_text(aes(x=-0.45,y=-0.25), label=deparse(bquote(paste(phi['t'['i']]))),parse=TRUE,size=7)

pdf(paste(PLOT_DIRPLOT,"TransX.pdf",sep=""),width=7, height=7)
PT9
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
print(p)
dev.off()

DataCoords_2 = as.data.frame(DataCoords)
colnames(DataCoords_2) = c("x","y")
dataDog  <- prepData(DataCoords_2,type="UTM")

P = ggplot(dataDog, aes(x=step))+geom_density(aes(x=step))+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Step-length")
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



####################################################
#### STAP-MODEL
####################################################

load(paste(PLOT_DIRDATA , MOD_STAP_NAME, ".Rdata",sep=""))
WM = c(14,164,71,23,26)
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
  dd = paste("$backslash nu_{j}$ ",sep="")
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
for(k in WM)
{
  dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
  ff = "(CI) "
  ii = 1
  for(i in WM)
  {
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


AppWM = rep(0,max(WM))
AppWM[WM] = 1:5

Z_MAP2_gen = AppWM[Z_MAP]


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
  print(p)
  dev.off()
}



###
nmcmc = dim(ModelOUT$psi)[1]
nseq = 200
thetaseq = seq(-pi,pi,length.out=nseq+1)[-(nseq+1)]
rseq = seq(0,3,length.out=nseq)+0.0000000000001
diff1 = thetaseq[2]-thetaseq[1]
diff2 = rseq[2]-rseq[1]
Densmat = list()
Densmat[[1]] = matrix(0, ncol=nseq,nseq)
Densmat[[2]] = matrix(0, ncol=nseq,nseq)
Densmat[[3]] = matrix(0, ncol=nseq,nseq)
Densmat[[4]] = matrix(0, ncol=nseq,nseq)
Densmat[[5]] = matrix(0, ncol=nseq,nseq)
Vecy  = matrix(NA, ncol=2, nrow=nseq*nseq)
Rvec     = rep(rseq, each=nseq)
Vecy[,1] = rep(rseq, each=nseq)*cos(rep(thetaseq,times=nseq))
Vecy[,2] = rep(rseq, each=nseq)*sin(rep(thetaseq,times=nseq))
for(imcmc in 1:nmcmc)
{
  k1 = 0
  for(k in WM)
  {
    k1 = k1+1
    mu = ModelOUT$muC[imcmc,,k]
    sigma = matrix(ModelOUT$sigma[imcmc,,k],2)

    Densmat[[k1]][,] = Densmat[[k1]][,]+matrix(dmnorm(Vecy,mu, sigma)*Rvec, ncol=nseq)/nmcmc

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
) +ylim(0,25.)+xlab("Step-Length")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette[c(1,2)])+guides(shape = guide_legend(override.aes = list(size = 5)))
p

pdf(paste(PLOT_DIRPLOT ,"Step1.pdf",sep=""))
print(p)
dev.off()



#### #### #### #### #### #### #### #### ####
#### OU model
#### #### #### #### #### #### #### #### ####


load(paste(PLOT_DIRDATA , MOD_OU_NAME, ".Rdata",sep=""))
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
  dd = paste("$backslash nu_{j}$ ",sep="")
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

ik=1
for(k in WM)
{
  dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
  ff = "(CI) "
  ii = 1
  for(i in WM)
  {
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

AppWM = rep(0,max(WM))
AppWM[WM] = 1:5

Z_MAP2_0 = AppWM[Z_MAP]



#### #### #### #### #### #### #### #### ####
#### ST
#### #### #### #### #### #### #### #### ####

load(paste(PLOT_DIRDATA , MOD_ST_NAME, ".Rdata",sep=""))
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
  print(p)
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

    Densmat_2[[k1]][,] = Densmat_2[[k1]][,]+matrix(dmnorm(Vecy,mu, sigma)*Rvec, ncol=nseq)/nmcmc
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
) +ylim(0,25)+xlab("Step-Length")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
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
   )+scale_x_productlist("OU-HMM")+scale_y_productlist("STAP-HMM")+scale_fill_manual(values=cbPalette)
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
   )+scale_x_productlist("ST-HMM")+scale_y_productlist("STAP-HMM")+scale_fill_manual(values=cbPalette)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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
theme(
  axis.text.y = element_text(face="bold",size=25),
  axis.text.x = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25),
  legend.text = element_text(face="bold",size=25),
  legend.title = element_text(face="bold",size=25)
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




#### #### #### #### #### #### #### #### ####
#### Observed turning-angle for the first two
#### OU-HMM behaviors
#### #### #### #### #### #### #### #### ####

load(paste(PLOT_DIRDATA , MOD_STAP_NAME,sep=""))
WM = c(14,164,71,23,26)
WMgen = WM

MAP = apply(ModelOUT$zeta,2,findmode)
INDEX_FIRST = which(MAP==14)
INDEX_SECOND = which(MAP==164)
length(INDEX_FIRST)



DataCoords_2 = as.data.frame(DataCoords)
colnames(DataCoords_2) = c("x","y")
dataDog  <- prepData(DataCoords_2,type="UTM")

###  usare l$y*3 anche nelle figure dei dati osservati
dataDog2 = rbind(dataDog[INDEX_FIRST,-1] -2*pi,dataDog[INDEX_FIRST,-1],dataDog[INDEX_FIRST,-1]+2*pi)
l <- density(dataDog2$angle,na.rm=T, bw=0.3)
ld =data.frame(x=l$x,y=l$y*3)
P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
P
pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k1.pdf",sep=""))
print(P)
dev.off()


dataDog2 = rbind(dataDog[INDEX_SECOND,-1] -2*pi,dataDog[INDEX_SECOND,-1],dataDog[INDEX_SECOND,-1]+2*pi)
l <- density(dataDog2$angle,na.rm=T, bw=0.3)
ld =data.frame(x=l$x,y=l$y*3)
P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
P
pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k2.pdf",sep=""))
print(P)
dev.off()



load(paste(PLOT_DIRDATA , MOD_OU_NAME,sep=""))
WM = c(39,136,14,127,23)

s1    = mean(ModelOUT$sigma[,1,WM[1]])
s12   = mean(ModelOUT$sigma[,2,WM[1]])
s2    = mean(ModelOUT$sigma[,4,WM[1]])

X = rmnorm(50000, c(0,0), matrix(c(s1,s12,s12,s2), ncol=2))
X = apply(X,2,cumsum)
colnames(X) = c("x","y")
X = data.frame(X)
dataDog  = prepData(X,type="UTM", LLangle=F)


# DataCoords_2 = as.data.frame(DataCoords)
# colnames(DataCoords_2) = c("x","y")
# dataDog  <- prepData(DataCoords_2,type="UTM",LLangle=F)

###  usare l$y*3 anche nelle figure dei dati osservati
dataDog2 = rbind(dataDog[,-1] -2*pi,dataDog[,-1],dataDog[,-1]+2*pi)
l <- density(dataDog2$angle,na.rm=T, bw=0.3)
ld =data.frame(x=l$x,y=l$y*3)
P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
P
pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k1_OU.pdf",sep=""))
print(P)
dev.off()




s1    = mean(ModelOUT$sigma[,1,WM[2]])
s12   = mean(ModelOUT$sigma[,2,WM[2]])
s2    = mean(ModelOUT$sigma[,4,WM[2]])

X = rmnorm(50000, c(0,0), matrix(c(s1,s12,s12,s2), ncol=2))
X = apply(X,2,cumsum)
colnames(X) = c("x","y")
X = data.frame(X)
dataDog  = prepData(X,type="UTM", LLangle=F)


# DataCoords_2 = as.data.frame(DataCoords)
# colnames(DataCoords_2) = c("x","y")
# dataDog  <- prepData(DataCoords_2,type="UTM",LLangle=F)

###  usare l$y*3 anche nelle figure dei dati osservati
dataDog2 = rbind(dataDog[,-1] -2*pi,dataDog[,-1],dataDog[,-1]+2*pi)
l <- density(dataDog2$angle,na.rm=T, bw=0.3)
ld =data.frame(x=l$x,y=l$y*3)
P = ggplot(ld, aes(x=x,y=y))+geom_line()+theme(
  axis.text.x = element_text(face="bold",size=25),
  axis.text.y = element_text(face="bold",size=25),
  axis.title.x = element_text(face="bold",size=25),
  axis.title.y = element_text(face="bold",size=25)
)+ylab("Density")+xlab("Turning-angle")+xlim(c(-pi,pi))+ylim(c(0,0.32))
P
pdf(paste(PLOT_DIRPLOT ,"ObsTurning_k2_OU.pdf",sep=""))
print(P)
dev.off()
