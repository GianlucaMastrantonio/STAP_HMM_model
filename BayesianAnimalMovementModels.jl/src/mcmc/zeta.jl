

# function sample_zeta!(Likelihood::AbstractLikelihood, Clusterization::AbstractClusterization)
#     error(string("sample_zeta not defined for Likelihood ", typeof(Likelihood), " and clusterization ", typeof(Clusterization)) )
# end


function sample_zeta!(Likelihood::AbstractLikelihood, Clusterization::AbstractClusterization)

    # Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1

    LogPdfMatrix = zeros(Float64,nt,kmax)
    @inbounds for i  in 1:(nt-1)
        for k in 1:kmax
             LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
        end
    end

    @inbounds for i  in 1:(nt-2)
        # print(i)
        # print(" \n")
        probsVect       = zeros(Float64,kmax)

         zetasux         = Likelihood.clusterization.zeta[i+1]
        for k in 1:kmax
             probsVect[k] = LogPdfMatrix[i+1,k] #compute_logpdf(Likelihood, Int16(k), Int16(i+1))
             probsVect[k] += log( Clusterization.pi.parameteracc[zetaprex][k])
             probsVect[k] += log( Clusterization.pi.parameteracc[k][zetasux])
             # print(Clusterization.pi.parameteracc[zetaprex][k])
             # print("\n")
             # print(Clusterization.pi.parameteracc[k][zetasux])
             # print("\n")
        end

         Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)
         zetaprex = Likelihood.clusterization.zeta[i]

    end
    i = nt-1
     probsVect       = zeros(Float64,kmax)
    for k in 1:kmax
         probsVect[k] = LogPdfMatrix[i+1,k] # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
         probsVect[k] += log( Clusterization.pi.parameteracc[zetaprex][k])

    end
    Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)

    Update_ObsInClust!(Likelihood.clusterization)
    return nothing

end


#
# function sample_zeta!(Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)
#
#     # Likelihood      = MCMCLikelihood
#     # Clusterization  = MCMCClusterization
#
# #    Clusterization.pi.parameteracc
#     kmax        = Likelihood.kmax
#     nc          = Likelihood.data.ncol
#     nanim       = Likelihood.data.nanimals
#     nt          = Likelihood.data.nt
#
#     zetaprex = 1
#
#     # LogPdfMatrix = zeros(Float64,nt,kmax)
#     # @inbounds for i  in 1:(nt-1)
#     #     for k in 1:kmax
#     #          LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#     #     end
#     # end
#     veckacc = zeros(Int16,kmax)
#     probsVect       = zeros(Float64,kmax)
#     @inbounds for i  in 1:(nt-2)
#         nk              = 0
#
#
#          zetasux        = Likelihood.clusterization.zeta[i+1]
#          kold           = Likelihood.clusterization.zeta[i]
#          un             = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][kold]))
#         for k in 1:kmax
#
#             if Clusterization.pi.parameteracc[zetaprex][k]>un
#                 nk           += one(nk)
#                 veckacc[nk]   = k
#                 probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) #comp
#                 probsVect[nk] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
#                 probsVect[nk] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-300))
#             end
#         end
#
#         sampvec = probsVect[1:nk]
#         Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(sampvec)]
#         zetaprex = Likelihood.clusterization.zeta[i]
#
#     end
#     i       = nt-1
#     nk      = 0
#     kold           = Likelihood.clusterization.zeta[i]
#     un             = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][kold]))
#     for k in 1:kmax
#
#         if Clusterization.pi.parameteracc[zetaprex][k]>un
#             nk           += one(nk)
#             veckacc[nk]   = k
#             probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#             probsVect[nk] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
#         end
#
#     end
#     sampvec = probsVect[1:nk]
#     Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(sampvec)]
#
#     Update_ObsInClust!(Likelihood.clusterization)
#     return nothing
#
# end
#


#
#
# function sample_zeta!(Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)
#
#     # Likelihood      = MCMCLikelihood
#     # Clusterization  = MCMCClusterization
#
# #    Clusterization.pi.parameteracc
#     kmax        = Likelihood.kmax
#     nc          = Likelihood.data.ncol
#     nanim       = Likelihood.data.nanimals
#     nt          = Likelihood.data.nt
#
#     zetaprex = 1
#
#     LogPdfMatrix = zeros(Float64,nt,kmax)
#     @inbounds for i  in 1:(nt-1)
#         for k in 1:kmax
#              LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#         end
#     end
#
#     @inbounds for i  in 1:(nt-2)
#
#         probsVect       = zeros(Float64,kmax)
#
#          zetasux         = Likelihood.clusterization.zeta[i+1]
#         for k in 1:kmax
#              probsVect[k] = LogPdfMatrix[i+1,k] #compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#              probsVect[k] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-100))
#              probsVect[k] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-100))
#              # print(Clusterization.pi.parameteracc[zetaprex][k])
#              # print("\n")
#
#         end
#
#
#          Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)
#          zetaprex = Likelihood.clusterization.zeta[i]
#
#     end
#     i = nt-1
#      probsVect       = zeros(Float64,kmax)
#     for k in 1:kmax
#          probsVect[k] = LogPdfMatrix[i+1,k] # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#          probsVect[k] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-100))
#
#     end
#     Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)
#
#     Update_ObsInClust!(Likelihood.clusterization)
#     return nothing
#
# end






function sample_zeta!(Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)

    u = rand(Uniform(0.0,1.0))
    if u<0.95
        BayesianAnimalMovementModels.sample_zeta_beam!(Likelihood,Clusterization)
    else
        BayesianAnimalMovementModels.sample_zeta_marg!(Likelihood,Clusterization)
    end

end


###### BEAM SAMPLER

function sample_zeta_beam!(Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)

    # Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    Clust = Clusterization.clusterization
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    # LogPdfMatrix = zeros(Float64,nt,kmax)
    # @inbounds for i  in 1:(nt-1)
    #     for k in 1:kmax
    #          LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
    #     end
    # end


    un    = zeros(Float64,nt)
    zetaprex = 1
    for i  in 1:(nt-1)
        @inbounds k1          = Likelihood.clusterization.zeta[i]
        @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
        @inbounds zetaprex    = k1
    end

    zetaprex = 1
    @inbounds for i  in 1:(nt-2)
        nk              = 0
         @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
         @inbounds kold           = Likelihood.clusterization.zeta[i]
        for k in 1:kmax

            if (Clusterization.pi.parameteracc[zetaprex][k]>un[i]) && (Clusterization.pi.parameteracc[k][zetasux]>un[i+1])
                nk            += one(nk)
                @inbounds veckacc[nk]   = k
                #probsVect[nk] = LogPdfMatrix[i+1,k] #comp
                @inbounds probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) #comp
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-300))

                # if k == zetaprex
                #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
                # else
                #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
                # end
                #
                # if k == zetasux
                #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
                # else
                #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
                # end


            end
        end

        @inbounds sampvec = probsVect[1:nk]
        @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(sampvec)]
        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    for k in 1:kmax

        if Clusterization.pi.parameteracc[zetaprex][k]>un[i]
            nk           += one(nk)
            veckacc[nk]   = k
            @inbounds probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
            #@inbounds probsVect[nk] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))

            # if k == zetaprex
            #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
            # else
            #     @inbounds probsVect[nk] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
            # end


        end

    end
    sampvec = probsVect[1:nk]
    @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(sampvec)]

    Update_ObsInClust!(Likelihood.clusterization)
    return nothing

end

##### MARGINALIZATION

function sample_zeta_marg!(Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)

    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    # LogPdfMatrix = zeros(Float64,nt,kmax)
    # for i  in 1:(nt-1)
    #     for k in 1:kmax
    #         @inbounds LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
    #     end
    # end

    #
    # un    = zeros(Float64,nt)
    # zetaprex = 1
    # for i  in 1:(nt-1)
    #     @inbounds k1          = Likelihood.clusterization.zeta[i]
    #     @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
    #     @inbounds zetaprex    = k1
    # end
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1
    @inbounds for i  in 1:(nt-2)

        zetasux        = Likelihood.clusterization.zeta[i+1]
        kold           = Likelihood.clusterization.zeta[i]
        @inbounds Clust.n_itojC[zetaprex][kold] -= 1
        @inbounds Clust.n_itojC[kold][zetasux]  -= 1
        for k in 1:kmax

            @inbounds probsVect[k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))

            if k == zetaprex
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
            else
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
            end

            if k == zetasux
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
            else
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
            end
        end

        @inbounds Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)

        @inbounds Clust.n_itojC[zetaprex][Likelihood.clusterization.zeta[i]] += 1
        @inbounds Clust.n_itojC[Likelihood.clusterization.zeta[i]][zetasux]  += 1

        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    Clust.n_itojC[zetaprex][kold] -= 1
    for k in 1:kmax

        @inbounds probsVect[k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))

        if k == zetaprex
            @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
        else
            @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
        end


    end

    @inbounds Likelihood.clusterization.zeta[i] = sample_discretevar(probsVect)

    Update_ObsInClust!(Likelihood.clusterization)
    return nothing

end
