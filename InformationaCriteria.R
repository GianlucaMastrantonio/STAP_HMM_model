#### #### #### #### #### #### #### ####
####  Values that must be changed
#### #### #### #### #### #### #### ####

PLOT_DIRDATA  = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/sheepdog/analisi/plot/"
PLOT_DIRPLOT  = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/sheepdog/analisi/plot/"
MOD_STAP_NAME = "RealData_Dog13_14_MISSING_STAP"
MOD_OU_NAME   = "RealData_Dog13_14_MISSING_OU"
MOD_ST_NAME   = "RealData_Dog13_14_MISSING_ST"



ICL   = c()
#DIC_part1 = c()
ICL_DATA = c()
ICL_ZETA = c()

ICL_DATAlist = list()
ICL_ZETAlist = list()


# DIC2_part1 = c()
# DIC3_part1 = c()

for(idata in 1:3)
{
	if(idata==1)
	{
		load(paste(PLOT_DIRDATA, MOD_STAP_NAME,".Rdata",sep=""))
		WWW = c(14,164,71,23,26)
	}
	if(idata==2)
	{
		load(paste(PLOT_DIRDATA, MOD_OU_NAME,".Rdata",sep=""))
		WWW = 	c(39,136,14,127,23)
	}
	if(idata==3)
	{
		load(paste(PLOT_DIRDATA, MOD_ST_NAME,".Rdata",sep=""))
		WWW = c(14,160,28)
	}
	Z_MAP = apply(ModelOUT$zeta,2,findmode)


	data = DataCoords
	nsim = dim(ModelOUT$psi)[1]
	nobs = nrow(data)
	### ICL
	for(isim in 1:nsim)
	{
		psi 	= ModelOUT$psi[isim,,]
		sigma = ModelOUT$sigma[isim,,]
		muC   = ModelOUT$muC[isim,,]
		mu0   = ModelOUT$mu0[isim,,]
 		zeta  = ModelOUT$zeta[isim,]
 		Pi    = ModelOUT$pi[isim,,]
 		rho   = ModelOUT$rho[isim,,]

		missing  = ModelOUT$missing[i,,]
		IndexRow = ModelOUT$missingIndexRow

		data[IndexRow,1:2] =  missing[,1:2]

		DIC_part1[isim] = 0
		angle = 0
		R = matrix(NA, ncol=2,nrow=2)
		logdens = matrix(NA, length(WWW))
		zpre = 1
		DIC_part1[isim] = 0
		DIC_part1_app = rep(0,nobs)
		LogMatdens = matrix(0, ncol=length(WWW),nrow = nobs)
		sigmaInv = list()
		logdet   = list()
		for(k in 1:length(WWW))
		{
			kk = WWW[k]
			sigmaInv[[kk]] = solve(matrix(sigma[,kk], ncol=2))
			logdet[[kk]]   = as.numeric(determinant(matrix(sigma[,kk], ncol=2) ,log=T)[[1]])
		}
		app = c()
		logdens = 0
		logzeta = 0
		zprec = 1
		for(iobs in 2:nobs)
		{

			kks = Z_MAP[iobs-1]

			R[,] = c(cos(rho[kks]*angle) ,sin(rho[kks]*angle), -sin(rho[kks]*angle) ,cos(rho[kks]*angle))
			mu   = (1.0-rho[kks])*psi[kks]*(mu0[,kks, drop=F]-t(data[iobs-1,, drop=F]))+rho[kks]*R%*%(muC[,kks,drop=F])


			NormDesn = log(2*pi)-0.5*logdet[[kks]]-0.5*t(t(data[iobs,,drop=F])-mu)%*%(R%*%sigmaInv[[kks]]%*%t(R))%*%(t(data[iobs,,drop=F])-mu)

			logdens = logdens+NormDesn

			angle = atan2(data[iobs,][2]-data[iobs-1,][2],data[iobs,][1]-data[ iobs-1,][1])

			logzeta = logzeta+log(Pi[kks,zprec])

			zprec = kks

		}
		ICL_DATA[isim] = logdens
		ICL_ZETA[isim] = logzeta

	}

ICL[idata] = max(ICL_DATA)+( max(ICL_ZETA)+log( sum(exp(ICL_ZETA-max(ICL_ZETA))    ) ) -log(nsim) )


ICL_DATAlist[[idata]] = ICL_DATAlist
ICL_ZETAlist[[idata]] = ICL_ZETAlist

names(ICL) = c("STAP","OU","ST")
print("")
print("")
print("")
print("ICL")
print(round(ICL,3))
