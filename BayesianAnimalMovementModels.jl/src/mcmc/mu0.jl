#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_mu0!(Likelihood::AbstractLikelihood, mu0Par::AbstractVecPar)
    error(string("sample_mu0 not defined for Likelihood ", typeof(Likelihood), " and mu0 ", typeof(mu0Par)) )
end


function sample_mu0!(Likelihood::Likelihood_OUmodel, mu0Par::VecParMvNormal)

    # Likelihood  = MCMCLikelihood
    # mu0Par      = MCMCLikelihood.mu0
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt


    InvMat_P      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(InvMat_P, deepcopy(mu0Par.prior_invmat.mat))
    end

    Mean_P      = Vector{Vector{Float64}}()
    for k in 1:kmax
        push!(Mean_P, deepcopy(mu0Par.prior.μ))
    end

    Dpsi =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        app_psi = repeat([1:nanim;], inner=2, outer=1)
        push!(Dpsi, deepcopy(diagm(Likelihood.psi.parameteracc[k][ app_psi ])))
    end

    IdMatrix = Matrix(1.0*I, nc, nc)
    appVar =  Vector{Matrix{Float64}}()
    appMean =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(appVar, (IdMatrix-Dpsi[k])*Likelihood.sigma.parameteraccinv[k].mat*(IdMatrix-Dpsi[k]) )
        push!(appMean, (IdMatrix-Dpsi[k])*Likelihood.sigma.parameteraccinv[k].mat )
    end


    for i in 2:nt
        k       = Likelihood.clusterization.zeta[i-1]
        yt1P    = Likelihood.data.data[i]  #view(Likelihood.data.data,2,:)
        yt      = Likelihood.data.data[i-1]

        ##

        InvMat_P[k] += appVar[k]
        Mean_P[k]   += appMean[k]*(yt1P-Dpsi[k]*yt)
    end

    for k in 1:kmax
        Covmat = inv(PDMat(Symmetric(InvMat_P[k])))
        #PDMat(Symmetric(InvMat_P[k])).mat-InvMat_P[k]
        #B = PDMat(Symmetric(InvMat_P[k])).mat
        Mean   = Covmat*Mean_P[k]

        mu0Par.parameteracc[k] = rand(MvNormal(Mean,Covmat))
        mu0Par.parameterprop[k] = deepcopy(mu0Par.parameteracc[k])
    end
    #isposdef(Likelihood.sigma.parameteraccinv[k])


    # delete objects
    InvMat_P    = nothing
    Dpsi        = nothing
    Mean_P      = nothing
    IdMatrix    = nothing
    appVar      = nothing
    appMean     = nothing
    yt1P        = nothing
    yt          = nothing
    Covmat      = nothing
    Mean        = nothing
end




###
function sample_mu0!(Likelihood::Likelihood_OU_CircLinmodel, mu0Par::VecParMvNormal)

    # Likelihood  = MCMCLikelihood
    # mu0Par      = MCMCLikelihood.mu0
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.psi.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.muC.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    InvMat_P      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(InvMat_P, deepcopy(mu0Par.prior_invmat.mat))
    end

    Mean_P      = Vector{Vector{Float64}}()
    for k in 1:kmax
        push!(Mean_P, deepcopy(mu0Par.prior.μ))
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
        app         =  Dpsirho[k]*MatR*sigmainv[k].mat*transpose(MatR)
        InvMat_P[k] += app*Dpsirho[k]
        Mean_P[k]   += app*(yt1P-yt+Dpsirho[k]*yt-Drho[k]*MatR2*muC[k]  )
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

        mu0Par.parameteracc[k] = rand(MvNormal(Mean,Covmat))
        mu0Par.parameterprop[k] = deepcopy(mu0Par.parameteracc[k])
    end
    #isposdef(Likelihood.sigma.parameteraccinv[k])

    return nothing

end
