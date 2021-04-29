#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_muC!(Likelihood::AbstractLikelihood, muCPar::AbstractVecPar)
    error(string("sample_muC not defined for Likelihood ", typeof(Likelihood), " and muC ", typeof(mu0Par)) )
end



###
function sample_muC!(Likelihood::Likelihood_OU_CircLinmodel, muCPar::VecParMvNormal)

    # Likelihood  = MCMCLikelihood
    # muCPar      = MCMCLikelihood.muC
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.psi.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.muC.parameteracc
    mu0         = Likelihood.mu0.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    InvMat_P      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(InvMat_P, deepcopy(muCPar.prior_invmat.mat))
    end

    Mean_P      = Vector{Vector{Float64}}()
    for k in 1:kmax
        push!(Mean_P, deepcopy(muCPar.prior.μ))
    end

    Dpsirho =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        app  = repeat([1:nanim;], inner=2, outer=1)
        push!(Dpsirho, diagm(  psi[k][ app ].*(1.0 .-rho[k][ app ])   )   )
    end
    Drho =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        app  = repeat([1:nanim;], inner=2, outer=1)
        push!(Drho, deepcopy(diagm(rho[k][ app ]  )))
    end

    IdMatrix = Matrix(1.0*I, nc, nc)
    appVar =  Vector{Matrix{Float64}}()
    appMean =  Vector{Matrix{Float64}}()



    Cangle  = deepcopy(StartAngle[1])
    MatR = zeros(Float64,nc,nc)
    MatR2 = zeros(Float64,nc,nc)
    for i in 2:nt
        k       = zeta[i-1]
        yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
        yt      = Obs[i-1]

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
        end

        ##
        app         =  Drho[k] *transpose(MatR2)*MatR* sigmainv[k].mat*transpose(MatR)
        InvMat_P[k] += app*MatR2 * Drho[k]
        Mean_P[k]   += app*(yt1P-yt-Dpsirho[k]*(mu0[k]-yt))
        ##

        for j in 1:nanim
            WW = [1 2] .+ (j-1)*2
            Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
        end


    end

    for k in 1:kmax
        Covmat = inv(PDMat(Symmetric(InvMat_P[k])))
        #PDMat(Symmetric(InvMat_P[k])).mat-InvMat_P[k]
        #B = PDMat(Symmetric(InvMat_P[k])).mat
        Mean   = Covmat*Mean_P[k]

        muCPar.parameteracc[k] = rand(MvNormal(Mean,Covmat))
        muCPar.parameterprop[k] = deepcopy(muCPar.parameteracc[k])
    end
    #isposdef(Likelihood.sigma.parameteraccinv[k])

    return nothing

end
