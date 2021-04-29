
#### GENERAL
function MCMCalgorithm(
    MCMCout::MCMCutils,
    MCMCLikelihood::AbstractLikelihood,
    MCMCClusterization::AbstractClusterization
)
    error("MCMC non implemented for Likelihood ",typeof(Likelihood), " and Clusterization ",typeof(Clusterization))
end

#### MODEL1

function MCMCalgorithm(
    MCMCout::MCMCutils,
    MCMCLikelihood::Likelihood_OUmodel,
    MCMCClusterization::Clusterization_HMM
)::Dict
    MCMCout.indexsave[1]    = 1
    appBurninOrThin         = MCMCout.burnin
    iterations = Int64(0)
    for saveMCMCindex in 1:MCMCout.nsamplesave

        for burnithinMCMCindex in 1:appBurninOrThin

            if mod(iterations,50) == 0 print(string("Iterations=",iterations,"\n"))  end
            #if mod(iterations,50) == 0 print(string("Clusters=",MCMCClusterization.clusterization.n_nonemptyC,"\n"))  end
            iterations += one(iterations)

            ### ### ### ### ### ### ###
            ### ### CLUSTERIZATION
            ### ### ### ### ### ### ###

            ### pi
            sample_pi!(MCMCClusterization, MCMCClusterization.pi)


            ### ### ### ### ### ### ###
            ### ### ZETA
            ### ### ### ### ### ### ###

            # zeta
            sample_zeta!(MCMCLikelihood, MCMCClusterization)


            ### ### ### ### ### ### ###
            ### ### LIKELIHOOD
            ### ### ### ### ### ### ###

            ### mu0
            sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)

            ### sigma
            sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)

            ### psi
            #sss = deepcopy(MCMCLikelihood.psi.parameteracc)
            sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)
#MCMCLikelihood.psi.parameteracc[1] == sss[1]
#MCMCLikelihood.psi.parameteracc[2] == sss[2]
            ### missing
            sample_missing!(MCMCLikelihood,MCMCLikelihood.miss)



        end # second for MCMC
        appBurninOrThin = MCMCout.thin

        save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)

    end # first for MCMC
    return Dict("PosteriorSamples"=>MCMCout.mcmc_out, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end




function MCMCalgorithm(
    MCMCout::MCMCutils,
    MCMCLikelihood::Likelihood_OU_CircLinmodel,
    MCMCClusterization::AbstractClusterization
)::Dict
    MCMCout.indexsave[1]    = 1
    appBurninOrThin         = MCMCout.burnin
    iterations = Int64(0)
    for saveMCMCindex in 1:MCMCout.nsamplesave

        for burnithinMCMCindex in 1:appBurninOrThin

            if mod(iterations,50) == 0 print(string("Iterations=",iterations,"\n"))  end
            if mod(iterations,50) == 0 print(string("Clusters=",MCMCClusterization.clusterization.n_nonemptyC[1],"\n"))  end
            iterations += one(iterations)

            ### ### ### ### ### ### ###
            ### ### ZETA
            ### ### ### ### ### ### ###

            # zeta
            sample_zeta!(MCMCLikelihood, MCMCClusterization)

            ### ### ### ### ### ### ###
            ### ### CLUSTERIZATION
            ### ### ### ### ### ### ###

            ### pi
            sample_pi!(MCMCClusterization, MCMCClusterization.pi)


            ### ### ### ### ### ### ###
            ### ### LIKELIHOOD
            ### ### ### ### ### ### ###

            ### mu0
            sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)

            ### sigma
            sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)

            ### psi
            sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)

            ### missing
            sample_missing!(MCMCLikelihood,MCMCLikelihood.miss)

            ### muC
            sample_muC!(MCMCLikelihood,MCMCLikelihood.muC)

            ### rho
            sample_rho!(MCMCLikelihood, MCMCLikelihood.rho)

        end # second for MCMC
        appBurninOrThin = MCMCout.thin

        save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)

    end # first for MCMC
    return Dict("PosteriorSamples"=>MCMCout.mcmc_out, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end
