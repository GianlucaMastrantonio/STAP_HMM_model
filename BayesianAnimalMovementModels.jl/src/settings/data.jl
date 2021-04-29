
#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### DATASET
#### #### #### #### #### #### #### #### #### #### #### ####
abstract type AbstractDataset end

struct CoordinatesDataset <: AbstractDataset

    data::Vector{Vector{Union{Float64}}}
    nt::Int32
    ncol::Int32
    nanimals::Int16

    CoordinatesDataset(data::Vector{Vector{Union{Float64}}},nt::Int32,ncol::Int32,nanimals::Int16) = new(data,nt,ncol,nanimals)
end

function CoordinatesDataset(datacopy::Matrix{Union{Missing,Float64}})

    nt          = Int32(size(datacopy)[1])
    ncol        = Int32(size(datacopy)[2])
    nanimals    = Int16((ncol/2))
    data        = Vector{Vector{Union{Float64}}}()
    for i in 1:nt
        push!(data, zeros(Union{Float64},ncol))
        for j in 1:ncol
            data[i][j] = datacopy[i,j]
        end
    end


    CoordinatesDataset(data,nt,ncol,nanimals)
end

#### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### MODEL
#### #### #### #### #### #### #### #### #### #### ####


####
abstract type AbstractLikelihood end
struct Likelihood_OUmodel{
            Tmiss<:AbstractMissing,
            Tmu0<:AbstractVecPar,
            Tpsi<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar ,
            Tzeta<:AbstractZeta
                        }     <: AbstractLikelihood

    data::CoordinatesDataset
    miss::Tmiss
    mu0::Tmu0
    psi::Tpsi
    sigma::Tsigma
    clusterization::Tzeta
    kmax::Int16

    Likelihood_OUmodel{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta}(data::CoordinatesDataset,miss::Tmiss,mu0::Tmu0,psi::Tpsi,sigma::Tsigma,clusterization::Tzeta,kmax::Int16) where {
                Tmiss<:AbstractMissing,
                Tmu0<:AbstractVecPar,
                Tpsi<:AbstractVecPar,
                Tsigma<:AbstractPosDefMatPar,
                Tzeta<:AbstractZeta
                            } = new{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta}(data,miss,mu0,psi,sigma,clusterization,kmax)
end

function Likelihood_OUmodel(data::CoordinatesDataset, miss::Tmiss,mu0::Tmu0,psi::Tpsi,sigma::Tsigma,clusterization::Tzeta,kmax::Int16) where {
            Tmiss<:AbstractMissing,
            Tmu0<:AbstractVecPar,
            Tpsi<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar,
            Tzeta<:AbstractZeta
                        }
    Likelihood_OUmodel{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta}(data,miss,mu0,psi,sigma,clusterization,kmax)

end

####
struct Likelihood_OU_CircLinmodel{
            Tmiss<:AbstractMissing,
            Tmu0<:AbstractVecPar,
            Tpsi<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar ,
            Tzeta<:AbstractZeta,
            TmuC<:AbstractVecPar,
            Trho<:AbstractVecPar,
            Tangle<:AbstractVecPar
                        }     <: AbstractLikelihood

    data::CoordinatesDataset
    miss::Tmiss
    mu0::Tmu0
    psi::Tpsi
    sigma::Tsigma
    clusterization::Tzeta
    kmax::Int16
    muC::TmuC
    rho::Trho
    Angle::Tangle
    SaveMissing::Bool

    Likelihood_OU_CircLinmodel{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta,TmuC,Trho,Tangle}(data::CoordinatesDataset,miss::Tmiss,mu0::Tmu0,psi::Tpsi,sigma::Tsigma,clusterization::Tzeta,kmax::Int16,muC::TmuC,rho::Trho,Angle::Tangle,SaveMissing::Bool) where {
                Tmiss<:AbstractMissing,
                Tmu0<:AbstractVecPar,
                Tpsi<:AbstractVecPar,
                Tsigma<:AbstractPosDefMatPar,
                Tzeta<:AbstractZeta,
                TmuC<:AbstractVecPar,
                Trho<:AbstractVecPar,
                Tangle<:AbstractVecPar
                            } = new{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta,TmuC,Trho,Tangle}(data,miss,mu0,psi,sigma,clusterization,kmax,muC,rho,Angle,SaveMissing)
end

function Likelihood_OU_CircLinmodel(data::CoordinatesDataset, miss::Tmiss,mu0::Tmu0,psi::Tpsi,sigma::Tsigma,clusterization::Tzeta,kmax::Int16,muC::TmuC,rho::Trho,Angle::Tangle,SaveMissing::Bool) where {
            Tmiss<:AbstractMissing,
            Tmu0<:AbstractVecPar,
            Tpsi<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar,
            Tzeta<:AbstractZeta,
            TmuC<:AbstractVecPar,
            Trho<:AbstractVecPar,
            Tangle<:AbstractVecPar
                        }
    Likelihood_OU_CircLinmodel{Tmiss,Tmu0,Tpsi,Tsigma,Tzeta,TmuC,Trho,Tangle}(data,miss,mu0,psi,sigma,clusterization,kmax,muC,rho,Angle,SaveMissing)

end


#### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### Second level
#### #### #### #### #### #### #### #### #### #### ####
abstract type AbstractClusterization end
struct Clusterization_HMM{Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}   <: AbstractClusterization

    clusterization::Tzeta
    pi::Tpi
    mcmc_initpi::Vector{Float64}
    mcmc_initz::Vector{Int16}

    Clusterization_HMM{Tpi,Tzeta}(clusterization::Tzeta,pi::Tpi,mcmc_initpi,mcmc_initz) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta} = new{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz)
end

function Clusterization_HMM(clusterization::Tzeta,pi::Tpi,mcmc_initpi::Vector{Float64},mcmc_initz::Vector{Int16}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}
    Clusterization_HMM{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz)
end


struct Clusterization_HDPHMM{Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}   <: AbstractClusterization

    clusterization::Tzeta
    pi::Tpi
    mcmc_initpi::Vector{Float64}
    mcmc_initz::Vector{Int16}
    mcmc_beta::Vector{Float64}
    mcmc_rho::Vector{Float64}
    mcmc_gamma::Vector{Float64}
    mcmc_ak::Vector{Float64}

    prior_ak::Vector{Float64}
    prior_gamma::Vector{Float64}
    prior_rho::Vector{Float64}

    Clusterization_HDPHMM{Tpi,Tzeta}(clusterization::Tzeta,pi::Tpi,mcmc_initpi,mcmc_initz,mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta} = new{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho)
end

function Clusterization_HDPHMM(clusterization::Tzeta, pi::Tpi,mcmc_initpi::Vector{Float64},mcmc_initz::Vector{Int16},mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}
    Clusterization_HDPHMM{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho)
end


#Clusterization_HMM(Likelihood.clusterization,mcmc_pi,mcmc_initpi,mcmc_initz)
