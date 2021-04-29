
struct MCMCutils{T<:Integer,Td<:AbstractDict}
    iterations::T
    burnin::T
    thin::T
    nsamplesave::T
    mcmc_out::Td
    indexsave::Vector{T}

    function MCMCutils{T,Td}(iterations::T,burnin::T,thin::T,nsamplesave::T,mcmc_out::Td) where {T<:Integer,Td<:AbstractDict}
        if iterations <= burnin error("iterations must be greater than burnin") end
        if iterations < 1 error("iterations must greater than zero") end
        if burnin < 1 error("burnin must greater than zero") end
        if thin < 1 error("thin must greater than zero")  end
        indexsave = ones(T(1))
        new{T,Td}(iterations,burnin,thin,nsamplesave,mcmc_out,indexsave)
    end

end

# function MCMCutils(iterations::T,burnin::T,thin::T) where {T<:Integer}
#
#     MCMCutils{T}(iterations,burnin,thin,T((iterations-burnin)/thin))
# end


function OptionsMCMC(Likelihood::Likelihood_OUmodel, Clusterization::Clusterization_HMM;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2}} }(
    "mu0"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
    "psi"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
    "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
    "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
    "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax)
    )



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OUmodel, Clusterization::Clusterization_HMM)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt


    for k in 1:kmax
        MCMCout.mcmc_out["mu0"][indOut[1],:,k] = Likelihood.mu0.parameteracc[k]
        MCMCout.mcmc_out["psi"][indOut[1],:,k] = Likelihood.psi.parameteracc[k]
        MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
        MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
        MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
    end

    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])

    return nothing
end



#####

function OptionsMCMC(Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HMM;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2}} }(
    "mu0"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
    "psi"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
    "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
    "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
    "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
    "muC"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
    "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
    )



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HMM)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt


    for k in 1:kmax
        MCMCout.mcmc_out["mu0"][indOut[1],:,k] = Likelihood.mu0.parameteracc[k]
        MCMCout.mcmc_out["psi"][indOut[1],:,k] = Likelihood.psi.parameteracc[k]
        MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
        MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
        MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
        MCMCout.mcmc_out["muC"][indOut[1],:,k] = Likelihood.muC.parameteracc[k]
        MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]
    end

    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])

    return nothing
end





function OptionsMCMC(Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        mcmc_out = Dict{String,Union{Array{Float64,2},Array{Float64,3},Array{Int16,2}} }(
        "gammaDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "rhoDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "akDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "betaDP" => Array{Float64,3}(undef,nsamplesave,nanim,kmax),

        "mu0"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "psi"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "muC"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        )
    else
        mcmc_out = Dict{String,Union{Array{Float64,2},Array{Float64,3},Array{Int16,2},Array{Int32,1},  Array{Array{Int32,1},1}} }(
        "gammaDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "rhoDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "akDP" => Array{Float64,2}(undef,nsamplesave,nanim),
        "betaDP" => Array{Float64,3}(undef,nsamplesave,nanim,kmax),

        "mu0"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "psi"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "muC"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "missing" => Array{Float64,3}(undef,nsamplesave,size(Likelihood.miss.indexrow,1),Likelihood.data.ncol),
        "missingIndexRow" => Likelihood.miss.indexrow,
        "missingIndexCol" => Likelihood.miss.indexcol,
        )


    end
 #Likelihood = MCMCLikelihood



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt


#MCMCout.mcmc_out["gammaDP"][1,:] = MCMCClusterization.mcmc_gamma
# MCMCout.mcmc_out["betaDP"][1,:,:] = MCMCClusterization.mcmc_beta
# MCMCClusterization.pi.parameteracc
# MCMCout.mcmc_out["pi"][indOut[1],:,1]
    if Likelihood.SaveMissing == false
        MCMCout.mcmc_out["gammaDP"][indOut[1],:] = Clusterization.mcmc_gamma
        MCMCout.mcmc_out["rhoDP"][indOut[1],:] = Clusterization.mcmc_rho
        MCMCout.mcmc_out["akDP"][indOut[1],:] = Clusterization.mcmc_ak

        MCMCout.mcmc_out["betaDP"][indOut[1],:,:] = Clusterization.mcmc_beta
        for k in 1:kmax

            MCMCout.mcmc_out["mu0"][indOut[1],:,k] = Likelihood.mu0.parameteracc[k]
            MCMCout.mcmc_out["psi"][indOut[1],:,k] = Likelihood.psi.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["muC"][indOut[1],:,k] = Likelihood.muC.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]
        end
    else
        MCMCout.mcmc_out["gammaDP"][indOut[1],:] = Clusterization.mcmc_gamma
        MCMCout.mcmc_out["rhoDP"][indOut[1],:] = Clusterization.mcmc_rho
        MCMCout.mcmc_out["akDP"][indOut[1],:] = Clusterization.mcmc_ak

        MCMCout.mcmc_out["betaDP"][indOut[1],:,:] = Clusterization.mcmc_beta

        for k in 1:kmax
            
            MCMCout.mcmc_out["mu0"][indOut[1],:,k] = Likelihood.mu0.parameteracc[k]
            MCMCout.mcmc_out["psi"][indOut[1],:,k] = Likelihood.psi.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["muC"][indOut[1],:,k] = Likelihood.muC.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]

            MCMCout.mcmc_out["missing"][indOut[1],:,:]  =  Matrix(transpose(hcat(Likelihood.data.data[Likelihood.miss.indexrow,:]...)))
            # for irow in 1:size(Likelihood.miss.indexrow,1)
            #
            # end
        end
    end
    # "gammaDP" => Array{Float64,3}(undef,nsamplesave,nanim,1),
    # "rhoDP" => Array{Float64,3}(undef,nsamplesave,nanim,1),
    # "kappaDP" => Array{Float64,3}(undef,nsamplesave,nanim,1),
    # "kappaDP" => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
    #
    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])
#Likelihood.data.data[Likelihood.miss.indexrow,:]

#Matrix(Likelihood.data.data[Likelihood.miss.indexrow,:])
    return nothing
end
