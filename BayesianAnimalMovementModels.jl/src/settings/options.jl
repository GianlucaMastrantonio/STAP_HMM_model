
#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### LIKELIHOOD
#### #### #### #### #### #### #### #### #### #### #### ####
function OptionsLikelihood(;
    # numer of regimes
    kmax::Int16,
    # data
    data::Matrix{Union{Float64,Missing}},
    # model type
    likelihood_type::String,  # one of "OU"

    # Missing
    update_missing::Bool,
    # mu0
    update_mu0::Bool,
    inits_mu0::Vector{Vector{Float64}},
    prior_mu0::Dict,
    # psi
    update_psi::Bool,
    inits_psi::Vector{Vector{Float64}},
    prior_psi::Dict,
    # muC
    update_muC::Bool,
    inits_muC::Vector{Vector{Float64}},
    prior_muC::Dict,
    # rho
    update_rho::Bool,
    inits_rho::Vector{Vector{Float64}},
    prior_rho::Dict,
    # Sigma
    update_sigma::Bool,
    inits_sigma::Vector{Matrix{Float64}},
    prior_sigma::Dict,
    # zeta
    update_zeta::Bool,
    inits_zeta ::Vector{Int16}
)

    # missing
    datacopy        = deepcopy(data)
    app_missing = size(findall(ismissing, datacopy))[1] == 0

    if !update_missing
        mcmc_missing = MissingDoNotUpdate(datacopy)
        if app_missing

        else
            printstyled("ATTENTION: missing values have been interpolated and there will be not estimated \n", color=:blue)
        end
        #mcmc_missing = MissingDoNotUpdate(datacopy)
    else
        mcmc_missing = MissingDoUpdate(datacopy)
    end


    # dataset
    mcmc_data   = CoordinatesDataset(datacopy::Matrix{Union{Missing,Float64}})
    nanimales   = mcmc_data.nanimals
    nt          = mcmc_data.nt
    ncol        = mcmc_data.ncol

    datacopy = nothing

    # zeta
    if size(inits_zeta,1)!=nt error("length inits_zeta must be ", nt) end
    if update_zeta==false
        mcmc_zeta = ZetaDoNotUpdate()
    else
        mcmc_zeta = ZetaDoUpdate(inits_zeta,kmax)
    end


    # data model
    mcmc_muC    = VecParDoNotUpdate()
    mcmc_rhoC   = VecParDoNotUpdate()
    ### ### ### ### ### ### ###
    ### OU
    ### ### ### ### ### ### ###
    if likelihood_type=="OU"


        # mu0
        if update_mu0==false
            mcmc_mu0 = VecParDoNotUpdate()
        else
            checkPriorName  = size(findall(["MvNormal"] .== prior_mu0["name"]),1)== 0
            if checkPriorName error("mu0 prior must be one of MvNormal") end

            namecheck       = string("check_Mixture",prior_mu0["name"])
            nameprior       = string(prior_mu0["name"],"_mu0")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_mu0, prior_mu0, kmax, ncol)

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_mu0["name"])))
            mcmc_mu0        = structobject(p_aftercheck)

        end

        # psi
        if update_psi==false
            mcmc_psi = VecParDoNotUpdate()
        else
            checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_psi["name"]),1)== 0
            if checkPriorName error("psi prior must be one of MvBeta") end

            namecheck       = string("check_Mixture",prior_psi["name"])
            nameprior       = string(prior_psi["name"],"_psi")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_psi, prior_psi, kmax, Integer(ncol/2))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_psi["name"])))
            mcmc_psi        = structobject(p_aftercheck)

        end

        # Sigma
        if update_psi==false
            mcmc_sigma = PosDefMatParDoNotUpdate()
        else

            checkPriorName  = size(findall(["InverseWishart"] .== prior_sigma["name"]),1)== 0
            if checkPriorName error("sigma prior must be one of InverseWishart") end

            namecheck       = string("check_Mixture",prior_sigma["name"])
            nameprior       = string(prior_sigma["name"],"_sigma")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_sigma, prior_sigma, kmax, Integer(ncol))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("PosDefMat",  prior_sigma["name"])))
            mcmc_sigma        = structobject(p_aftercheck)

        end
        Likelihood = Likelihood_OUmodel(mcmc_data,mcmc_missing,mcmc_mu0,mcmc_psi,mcmc_sigma,mcmc_zeta,kmax)

    ### ### ### ### ### ### ###
    ### OU_CircLin
    ### ### ### ### ### ### ###
    elseif likelihood_type=="OU_CircLin"
        # mu0
        if update_mu0==false
            mcmc_mu0 = VecParDoNotUpdate()
        else
            checkPriorName  = size(findall(["MvNormal"] .== prior_mu0["name"]),1)== 0
            if checkPriorName error("mu0 prior must be one of MvNormal") end

            namecheck       = string("check_Mixture",prior_mu0["name"])
            nameprior       = string(prior_mu0["name"],"_mu0")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_mu0, prior_mu0, kmax, ncol)

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_mu0["name"])))
            mcmc_mu0        = structobject(p_aftercheck)

        end

        # psi
        if update_psi==false
            mcmc_psi = VecParDoNotUpdate()
        else
            checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_psi["name"]),1)== 0
            if checkPriorName error("psi prior must be one of MvBeta") end

            namecheck       = string("check_Mixture",prior_psi["name"])
            nameprior       = string(prior_psi["name"],"_psi")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_psi, prior_psi, kmax, Integer(ncol/2))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_psi["name"])))
            mcmc_psi        = structobject(p_aftercheck)

        end

        # muC
        if update_muC==false
            mcmc_muC = VecParDoNotUpdate()
        else
            checkPriorName  = size(findall(["MvNormal"] .== prior_muC["name"]),1)== 0
            if checkPriorName error("muC prior must be one of MvNormal") end

            namecheck       = string("check_Mixture",prior_muC["name"])
            nameprior       = string(prior_muC["name"],"_mu0")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_muC, prior_muC, kmax, ncol)

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_muC["name"])))
            mcmc_muC        = structobject(p_aftercheck)

        end

        # rho
        if update_rho==false

            namecheck       = string("check_Mixture",prior_rho["name"])
            nameprior       = string("MvUniform","_rho")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_rho, prior_rho, kmax, Integer(ncol/2))
            #print(p_aftercheck[3])

            mcmc_rho        = VecParMvNoPrior(p_aftercheck[3])

        else
            checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_rho["name"]),1)== 0
            if checkPriorName error("rho prior must be one of MvBeta") end

            namecheck       = string("check_Mixture",prior_rho["name"])
            nameprior       = string(prior_rho["name"],"_rho")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_rho, prior_rho, kmax, Integer(ncol/2))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_rho["name"])))
            mcmc_rho        = structobject(p_aftercheck)

        end

        # Sigma
        if update_psi==false
            mcmc_sigma = PosDefMatParDoNotUpdate()
        else

            checkPriorName  = size(findall(["InverseWishart"] .== prior_sigma["name"]),1)== 0
            if checkPriorName error("sigma prior must be one of InverseWishart") end

            namecheck       = string("check_Mixture",prior_sigma["name"])
            nameprior       = string(prior_sigma["name"],"_sigma")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_sigma, prior_sigma, kmax, Integer(ncol))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("PosDefMat",  prior_sigma["name"])))
            mcmc_sigma        = structobject(p_aftercheck)

        end

        # angle
        update_angle = true

        if update_angle==false
            mcmc_angle = VecParDoNotUpdate()
        else
            prior_angle = Dict("name"=>"MvUniform","a"=> -pi*ones(Float64,Integer(ncol/2)), "b"=> pi*ones(Float64,Integer(ncol/2)))
            inits_angle = [zeros(Integer(Integer(ncol/2)),kmax)[:,i] for i in 1:size(zeros(Integer(Integer(ncol/2)),kmax),2)]

            checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_angle["name"]),1)== 0
            if checkPriorName error("angle prior must be one of MvBeta") end

            namecheck       = string("check_Mixture",prior_angle["name"])
            nameprior       = string(prior_angle["name"],"_angle")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_angle, prior_angle, kmax, Integer(ncol/2))

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_angle["name"])))
            mcmc_angle        = structobject(p_aftercheck)

        end

        Likelihood = Likelihood_OU_CircLinmodel(mcmc_data,mcmc_missing,mcmc_mu0,mcmc_psi,mcmc_sigma,mcmc_zeta,kmax,mcmc_muC,mcmc_rho,mcmc_angle, true)
    else

        error("likelihood_type must be on of OU or OU_CircLin")
    end

    return Likelihood

end

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### LATENT CLASSIFICATION
#### #### #### #### #### #### #### #### #### #### #### ####
function OptionsClusterization(;
    Likelihood::Tl,
    clustering_type::String, # one of "HMM"
    # pi
    update_pi::Bool,
    inits_pi::Vector{Vector{Float64}},
    prior_pi::Dict,
    prior_ak = [1.0; 1.0]::Vector{Float64},
    prior_gamma = [1.0; 1.0]::Vector{Float64},
    prior_rho = [1.0 ;1.0]::Vector{Float64}
    ) where {Tl<:AbstractLikelihood}

    kmax = Likelihood.kmax
    if clustering_type == "HMM"

        # pi
        if update_pi==false
            mcmc_pi = PiDoNotUpdate()
        else

            checkPriorName  = size(findall(["Dirichlet"] .== prior_pi["name"]),1)== 0
            if checkPriorName error("pi prior must be one of Dirichlet") end

            namecheck       = string("check_Clusterization",prior_pi["name"])
            nameprior       = string(prior_pi["name"],"_pi")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_pi, prior_pi, kmax)

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_pi["name"])))
            mcmc_pi        = structobject(p_aftercheck)
        end

        mcmc_initpi              = 1.0/(Float64(kmax))*ones(Float64,kmax)
        mcmc_initz           = ones(Int16,1)
        Clusterization = Clusterization_HMM(Likelihood.clusterization,mcmc_pi,mcmc_initpi,mcmc_initz)
    elseif clustering_type == "HDP-HMM"
        # pi
        if update_pi==false
            mcmc_pi = PiDoNotUpdate()
        else

            checkPriorName  = size(findall(["Dirichlet"] .== prior_pi["name"]),1)== 0
            if checkPriorName error("pi prior must be one of Dirichlet") end

            namecheck       = string("check_Clusterization",prior_pi["name"])
            nameprior       = string(prior_pi["name"],"_pi")
            checkfunc       = getfield(BayesianAnimalMovementModels,Symbol(namecheck))
            p_aftercheck    = checkfunc(nameprior,inits_pi, prior_pi, kmax)

            structobject    = getfield(BayesianAnimalMovementModels,Symbol(string("VecPar",  prior_pi["name"])))
            mcmc_pi        = structobject(p_aftercheck)
        end

        mcmc_initpi              = 1.0/(Float64(kmax))*ones(Float64,kmax)
        mcmc_initz               = ones(Int16,1)
        mcmc_beta                = 1.0/kmax*ones(Float64,kmax)
        mcmc_rho                 = [0.5]
        mcmc_gamma               = [1.0]
        mcmc_ak                  = [1.0]



        Clusterization = Clusterization_HDPHMM(Likelihood.clusterization,mcmc_pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho)

    else

        error("clustering_type must be one of HMM or ...")
    end

    return Clusterization

end

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### MCMC interations
#### #### #### #### #### #### #### #### #### #### #### ####
