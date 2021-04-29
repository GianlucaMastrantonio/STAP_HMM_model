function sample_sigma!(Likelihood::AbstractLikelihood, sigmaPar::AbstractPosDefMatPar)
    error(string("sample_sigma not defined for Likelihood ", typeof(Likelihood), " and sigma ", typeof(sigmaPar)) )
end


function sample_sigma!(Likelihood::Likelihood_OUmodel, sigmaPar::PosDefMatInverseWishart)

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    Dpsi =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        app_psi = repeat([1:nanim;], inner=2, outer=1)
        push!(Dpsi, deepcopy(diagm(Likelihood.psi.parameteracc[k][ app_psi ])))
    end

    Mat_P =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(Mat_P, deepcopy(sigmaPar.prior.Ψ.mat))
    end

    for i in 2:nt
        k       = Likelihood.clusterization.zeta[i-1]
        yt1P    = Likelihood.data.data[i]
        yt      = Likelihood.data.data[i-1]

        Obs     = yt1P-(Likelihood.mu0.parameteracc[k]+Dpsi[k]*(yt- Likelihood.mu0.parameteracc[k]))

        Mat_P[k] += Obs*transpose(Obs)

    end


    for k in 1:kmax
        par1 = sigmaPar.prior.df+typeof(sigmaPar.prior.df)(Likelihood.clusterization.n_obsC[k])
        sigmaPar.parameteracc[k]    = PDMat(rand(InverseWishart(par1, Mat_P[k] )))
        sigmaPar.parameteraccinv[k] = inv(sigmaPar.parameteracc[k])

        sigmaPar.parameterprop[k] = deepcopy(sigmaPar.parameteracc[k])
        sigmaPar.parameterpropinv[k] = deepcopy(sigmaPar.parameteraccinv[k])
    end

    Mat_P   = nothing
    Dpsi    = nothing
    yt1P    = nothing
    yt      = nothing

end


####


function sample_sigma!(Likelihood::Likelihood_OU_CircLinmodel, sigmaPar::PosDefMatInverseWishart)

    # Likelihood  = MCMCLikelihood
    # sigmaPar      = MCMCLikelihood.sigma

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.psi.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.muC.parameteracc
    mu0         = Likelihood.mu0.parameteracc
    #sigmainv    = Likelihood.sigma.parameteraccinv

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc


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


    Mat_P =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(Mat_P, deepcopy(sigmaPar.prior.Ψ.mat))
    end

    Cangle  = deepcopy(StartAngle[1])
    MatR    = zeros(Float64,nc,nc)
    MatR2    = zeros(Float64,nc,nc)
    for i in 2:nt
        k       = zeta[i-1]
        yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
        yt      = Obs[i-1]

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim])  cos(Cangle[ianim] )  ]
        end

        SumsObs   = yt1P-(yt+Dpsirho[k]*(mu0[k]-yt) +Drho[k]*MatR2*muC[k])

        Mat_P[k] += transpose(MatR)*SumsObs*transpose(SumsObs)*MatR

        for j in 1:nanim
            WW = [1 2] .+ (j-1)*2
            Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
        end

    end


    for k in 1:kmax
        par1 = sigmaPar.prior.df+typeof(sigmaPar.prior.df)(Likelihood.clusterization.n_obsC[k])
        sigmaPar.parameteracc[k]    = PDMat(rand(InverseWishart(par1,PDMat(Symmetric(Mat_P[k]))  )))
        sigmaPar.parameteraccinv[k] = inv(sigmaPar.parameteracc[k])

        sigmaPar.parameterprop[k] = deepcopy(sigmaPar.parameteracc[k])
        sigmaPar.parameterpropinv[k] = deepcopy(sigmaPar.parameteraccinv[k])
    end

    return nothing

end
