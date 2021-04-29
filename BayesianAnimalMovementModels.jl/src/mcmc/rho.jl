#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_rho!(Likelihood::AbstractLikelihood, rhoPar::AbstractVecPar)
    error(string("sample_rho not defined for Likelihood ", typeof(Likelihood), " and rho ", typeof(rhoPar)) )
end

function sample_rho!(Likelihood::Likelihood_OU_CircLinmodel, rhoPar::VecParMvNoPrior)

end

function sample_rho!(Likelihood::Likelihood_OU_CircLinmodel, rhoPar::VecParMvUniform)

    #Likelihood  = MCMCLikelihood
    # rhoPar      = MCMCLikelihood.rho

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.psi.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.muC.parameteracc
    mu0         = Likelihood.mu0.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv
    sigma       = Likelihood.sigma.parameteracc

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    # CovInv_AgivenB      = Vector{Matrix{Float64}}()
    # for k in 1:kmax
    #     push!(CovInv_AgivenB, zeros(Float64,2,2))
    # end
    #
    # CovABInvCovB      = Vector{Matrix{Float64}}()
    # for k in 1:kmax
    #     push!(CovABInvCovB, zeros(Float64,2,nc-2))
    # end


    for ianim in 1:nanim

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



        MatR = zeros(Float64,nc,nc)
        MatR_prop = zeros(Float64,nc,nc)
        MatR2 = zeros(Float64,nc,nc)
        MatR2_prop = zeros(Float64,nc,nc)

        Loglike = zeros(Float64,kmax)
        prop = zeros(Float64,kmax)
        lim1 = rhoPar.prior.v[ianim].a
        lim2 = rhoPar.prior.v[ianim].b

        lim1 = 0.0
        lim2 = 1.0


        prob_lim = Float64(1.0/3.0)*3.0
        for k in 1:kmax
            isinternal = Int16(0)
            uv   = rand(Uniform(0.0,3.0))
            if uv < prob_lim
                prop[k] = 0.0
            elseif uv > 3.0-prob_lim
                prop[k] = 1.0
            else
                u    = rand(Uniform(0.0,1.0))
                if u<0.25
                    prop[k] = rand(Uniform(rho[k][ianim]-0.5, rho[k][ianim]+0.5))
                elseif u<0.5
                    prop[k] = rand(Uniform(rho[k][ianim]-0.25, rho[k][ianim]+0.25))
                elseif u<0.75
                    prop[k] = rand(Uniform(rho[k][ianim]-0.1, rho[k][ianim]+0.1))
                else
                    prop[k] = rand(Uniform(rho[k][ianim]-0.01, rho[k][ianim]+0.01))
                end
                isinternal = 1
                #prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
            end


            if (prop[k]>0.0) && (prop[k]<1.0)
                Loglike[k] -= log(1.0-2.0*prob_lim/3.0)
                appP = 0.0
                if abs(prop[k]-rho[k][ianim])<0.5
                    appP += 0.25*1.0/(0.5*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.25
                    appP += 0.25*1.0/(0.25*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end
                Loglike[k] -= log(appP)
            else
                Loglike[k] -= log(prob_lim/3.0)
            end

            if (rho[k][ianim]>0.0) && (rho[k][ianim]<1.0)

                Loglike[k] += log(1.0-2.0*prob_lim/3.0)

                appP = 0.0
                if abs(prop[k]-rho[k][ianim])<0.5
                    appP += 0.25*1.0/(0.5*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.25
                    appP += 0.25*1.0/(0.25*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end

                Loglike[k] += log(appP)
            else
                Loglike[k] += log(prob_lim/3.0)
            end

            if isinternal == 1
                prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
            end


            #prop[k] = rho[k][ianim]
        end

        # for k in 1:kmax
        #     prop[k] = rho[k][ianim]+0.00001
        # end
        Dpsirho_prop     = deepcopy(Dpsirho)
        Drho_prop        = deepcopy(Drho)

        for k in 1:kmax
            WW = (ianim-1)*2 .+ [1;2]
            Dpsirho_prop[k][WW,WW]     = (1.0-prop[k]).*psi[k][ianim].*Matrix(I,2,2)
            Drho_prop[k][WW,WW]        = prop[k].*Matrix(I,2,2)
        end


        Cangle      = deepcopy(StartAngle[1])
        MatR        = zeros(Float64,nc,nc)
        MatRprop    = zeros(Float64,nc,nc)
        MatR2        = zeros(Float64,nc,nc)
        MatR2prop    = zeros(Float64,nc,nc)
        for i in 2:nt
            k       = zeta[i-1]
            yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
            yt      = Obs[i-1]

            for ianim2 in 1:nanim
                W = [1,2] .+(ianim2-1)*2
                MatR[W,W]       = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                MatR2[W,W]       = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                if ianim2!=ianim
                    MatR_prop[W,W]  = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                    MatR2_prop[W,W]  = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                else
                    MatR_prop[W,W]  = [ cos(Cangle[ianim2]*prop[k]) -sin(Cangle[ianim2]*prop[k]); sin(Cangle[ianim2]*prop[k])  cos(Cangle[ianim2]*prop[k])  ]
                    MatR2_prop[W,W]  = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                end

            end

            ##
            Mean_prop   = transpose(MatR_prop)*(yt+Dpsirho_prop[k]*(mu0[k]-yt) + Drho_prop[k]*MatR2_prop*muC[k])
            Mean_acc    = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k])


            Loglike[k] += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*yt1P)
            Loglike[k] -= logpdf(MvNormal(Mean_acc,  sigma[k]),   transpose(MatR)*yt1P)
            ##

            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
            end
            #print(Loglike)
            #print("\n\n")

        end

        for k in 1:kmax
            u = rand(Uniform(0.0,1.0))
            # print(exp(Loglike[k]))
            # print("\n")
            # print(prop[k])
            # print(" ")
            # print(rhoPar.parameteracc[k][ianim])
            # print("\n")
            if u<exp(Loglike[k])
                rhoPar.parameteracc[k][ianim]  = prop[k]
                rhoPar.parameterprop[k][ianim] = prop[k]
            end
        end

        # print("LogLike")
        # print(Loglike)
        # print("\n")
    end
    # print("RHO")
    # print(rhoPar.parameteracc)
    # print("\n")


    return nothing

end
