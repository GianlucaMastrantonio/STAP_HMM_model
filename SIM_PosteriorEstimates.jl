#### #### #### #### #### #### #### ####
####  Values that must be changed
#### #### #### #### #### #### #### ####

IndModel = [1]# [1] for STAP, [2] for BRW, [3] for CRW
NAME 	= "Sim"
DIR  	= "/work/gmastrantonio/Sheep/"

##### Packages
using Pkg
Pkg.activate("/work/gmastrantonio/Sheep/")
using BayesianAnimalMovementModels
using Distributions, Random
using LinearAlgebra, PDMats
using RCall


## Directories
DIRDATA = string(DIR,"Data/")
DIROUT  = string(DIR,"out/")


for ipar1 in [1; 2; 3; 4]
	for ipar2 in [1; 2; 3; 4]

		## Seed
		Seed 	= 100
		rng     = MersenneTwister(Seed);
		Random.seed!(Seed)

		## Options
		IndAn   	= [9;10]
		rhoTT 		= [0.5;0.0;1.0]
		sampleTT 	= [true;false;false]
		nameTT   	= ["_STAP";"_BRW"; "_CRW"]

		## Load data from R
		@rput DIRDATA
		@rput DIROUT
		@rput IndAn
		R" load(paste(DIRDATA,'Data3.Rdata'  ,sep='')) "
		@rget DataAn3

		R"""

			DataCoords 	= DataAn3[,c(IndAn[1],IndAn[2])]

			SDtot1   	= sd(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
			SDtot2   	= sd(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

			Meantot1 	= mean(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
			Meantot2 	= mean(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

			SDtot    	= mean(c(SDtot1,SDtot2))
			Meantot  	= mean(c(Meantot1,Meantot2))

			DataCoords[,seq(1,ncol(DataCoords),by=2 )] = (DataCoords[,seq(1,ncol(DataCoords),by=2 )]-Meantot1)/SDtot
			DataCoords[,seq(2,ncol(DataCoords),by=2 )] = (DataCoords[,seq(2,ncol(DataCoords),by=2 )]-Meantot2)/SDtot


		"""
		@rget DataCoords
		# NANindex 				= isnan.(DataCoords)
		# DataCoords  			= Matrix{Union{Float64,Missing}}(DataCoords)
		#DataCoords[NANindex]   .= missing



		Dataset = Matrix{Union{Float64,Missing}}(DataCoords[:,:])

		## The model parameters
		kmax        = 200
		NAMETOT     = NAME
		nt          = size(Dataset,1)

		nc          = size(Dataset,2)
		nanim       = Int32(nc/2)


		InitSigma   = Vector{Matrix{Float64}}()
		for k in 1:kmax
		    push!(InitSigma ,Matrix{Float64}(I,nc,nc))
		end
		rand(InverseWishart(nc*3.0,Matrix{Float64}(I,nc,nc)*25))

		NAMETOTapp = deepcopy(NAMETOT)


		for irho in IndModel

			#global NAMETOTapp
			NAMETOT     = string(NAMETOTapp,nameTT[irho],"par1_",ipar1,"par2_",ipar2)
			MCMCLikelihood = OptionsLikelihood(
			    # numer of regimes
			    kmax                = Int16(kmax),
			    # data
			    data                = Dataset,
			    # model type
			    likelihood_type     = "OU_CircLin",

			    # Missing
			    update_missing      = true,
			    # mu0
			    update_mu0          = true,
			    inits_mu0           = [rand(nc,kmax)[:,i] for i in 1:kmax],
			    prior_mu0           = Dict("name"=>"MvNormal","μ"=>zeros(nc), "Σ"=>Matrix{Float64}(I,nc,nc)*1000.0),
			    # psi
			    update_psi          = true,
			    inits_psi           = [0.5*ones(Integer(nanim),kmax)[:,i] for i in 1:kmax],
			    prior_psi           = Dict("name"=>"MvUniform","a"=>zeros(nanim), "b"=>ones(nanim)),
				# muC
			    update_muC          = true,
			    inits_muC           = [zeros(nc,kmax)[:,i] for i in 1:size(zeros(nc,kmax),2)],
			    prior_muC           = Dict("name"=>"MvNormal","μ"=>zeros(nc), "Σ"=>Matrix{Float64}(I,nc,nc)*1000.0),
			    # psi
			    update_rho          = sampleTT[irho],
			    inits_rho           = [rhoTT[irho]*ones(Integer(nanim),kmax)[:,i] for i in 1:size(rhoTT[irho]*ones(Integer(nanim),kmax),2)],
			    prior_rho           = Dict("name"=>"MvUniform","a"=>zeros(nanim), "b"=>ones(nanim)),
			    # Sigma
			    update_sigma        = true,
			    inits_sigma         = InitSigma,
			    prior_sigma         = Dict("name"=>"InverseWishart","df"=>Float64(nc+1), "Ψ"=>Matrix{Float64}(I,nc,nc)),
			    # zeta
			    update_zeta         = true,
			    inits_zeta          = Int16.(sample(1:10, nt, replace = true))
			)

			InitPi   = Vector{Vector{Float64}}()
			for k in 1:kmax
			    push!(InitPi ,ones(kmax)/kmax)
			end

			MCMCClusterization = OptionsClusterization(
			    Likelihood =           MCMCLikelihood,
			    clustering_type     = "HDP-HMM"::String, # one of "HMM"
			    # pi
			    update_pi           = true::Bool,
			    inits_pi            = InitPi,
			    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>ones(kmax)),
				prior_ak 			= [[0.01; 0.1; 1; 10.0][ipar1]; 1.0],
			    prior_gamma 		= [[0.01; 0.1; 1; 10.0][ipar1]; 1.0],
			    prior_rho 			= [[10.0; 10.0; 1.0; 1.0][ipar2]; [1.0; 10.0; 1.0; 10.0][ipar2]]
			    )
			molt = 5
			MCMCout = OptionsMCMC(MCMCLikelihood,MCMCClusterization;
					iterations = Int64(molt*25000),
					burnin = Int64(molt*15000),
					thin = Int64(molt*2),

			        )

			## MCMC algorithm

			ModelOUT_TOT = MCMCalgorithm(
			    MCMCout,
			    MCMCLikelihood,
			    MCMCClusterization
			)


			ModelOUT = ModelOUT_TOT["PosteriorSamples"]


			@rput ModelOUT
			@rput kmax
			@rput NAME
			@rput NAMETOT
			@rput DIROUT
			#@rput
			#@rput

			## Save the results as an R object

			R"""
				save.image(paste(DIROUT,NAMETOT,".Rdata",sep=""))
			"""
			R"""
				DIR_PLOT=DIROUT

				save.image(paste(DIR_PLOT,NAMETOT,".Rdata",sep=""))

				#library(GenFunc)
			    library(lattice)
				findmode = function(x){
                    TT = table(as.vector(x))
                    return(as.numeric(names(TT)[TT==max(TT)][1]))
            }

			    Z_MAP = apply(ModelOUT$zeta,2,findmode)


			    pdf(paste(DIR_PLOT,NAMETOT,".pdf",sep=""))
				barplot(table(Z_MAP))
				app  = apply(ModelOUT$zeta,2,function(x) length(table(x)))
				barplot(app)
				barplot(table(app))
			    #mosaicplot(table(Z_MAP,LatentClassification[-length(LatentClassification)]))
				plot(Z_MAP, type="l")
			    for(i in 1:nanim)
			    {
			        plot(range(c(  DataCoords[,seq(1,nCov, by=2)] ),na.rm=T),range(c(  DataCoords[,1+seq(1,nCov, by=2)] ),na.rm=T), type="n")
			        points(DataCoords[,1:2+(i-1)*2], type="l")
			        points(DataCoords[,1:2+(i-1)*2], pch=20, cex=0.4, col=c(1,Z_MAP+1))

			    }
				par(mfrow=c(3,3))
			    for(i in 1:(nCov-1))
			    {
			        for(j in (i+1):nCov)
			        {
			            plot(DataCoords[,c(i,j)], pch=20, cex=0.4, col=c(1,Z_MAP+1))
			        }
			    }
				par(mfrow=c(3,3))
				Unz = unique(Z_MAP)
				Unz = Unz[order(Unz)]
				## mu0

				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in 1:dim(ModelOUT$mu0)[2])
				  {
				      plot(ModelOUT$mu0[,i,k], type="l", main=paste("mu0",i,k))
				  }
				}
				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in 1:dim(ModelOUT$muC)[2])
				  {
				      plot(ModelOUT$muC[,i,k], type="l", main=paste("muC",i,k))
				  }
				}
				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in 1:dim(ModelOUT$psi)[2])
				  {
				      plot(ModelOUT$psi[,i,k], type="l", main=paste("psi",i,k))
				  }
				}
				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in 1:dim(ModelOUT$rho)[2])
				  {
				      plot(ModelOUT$rho[,i,k], type="l", main=paste("rho",i,k))
				  }
				}
				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in 1:dim(ModelOUT$sigma)[2])
				  {
				      plot(ModelOUT$sigma[,i,k], type="l", main=paste("sigma",i,k))
				  }
				}
				for(k in Unz)
				{
				  par(mfrow=c(3,3))
				  for(i in Unz)
				  {
				      plot(ModelOUT$pi[,i,k], type="l", main=paste("pi",i,k))
				  }
				}
				for(k in Unz)
				{
				  print(levelplot(matrix(colMeans(ModelOUT$sigma[,,k]),ncol=nanim*2)[(nanim*2):1,], main=paste(k)))
				}
				dev.off()
			"""
		end
	end
end
