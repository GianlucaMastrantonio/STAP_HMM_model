module BayesianAnimalMovementModels

##### Packages
using Distributions, Random
using LinearAlgebra, PDMats
using Impute
#import: Distributions:
##### Include
include(joinpath("checks.jl"))
include(joinpath("functions.jl"))
include(joinpath("settings/missing.jl"))

include(joinpath("settings/vectorparameters.jl"))
include(joinpath("settings/posdefmats.jl"))
include(joinpath("settings/zeta.jl"))
#
# include(joinpath("settings/mu0.jl"))
# include(joinpath("settings/psi.jl"))
# include(joinpath("settings/sigma.jl"))
# include(joinpath("settings/pi.jl"))

include(joinpath("settings/data.jl"))
include(joinpath("settings/options.jl"))
include(joinpath("settings/utils.jl"))

include(joinpath("mcmc/mu0.jl"))
include(joinpath("mcmc/muC.jl"))
include(joinpath("mcmc/sigma.jl"))
include(joinpath("mcmc/pi.jl"))
include(joinpath("mcmc/psi.jl"))
include(joinpath("mcmc/rho.jl"))
include(joinpath("mcmc/data.jl"))
include(joinpath("mcmc/zeta.jl"))
include(joinpath("mcmc/missing.jl"))
include(joinpath("mcmcalgorithms.jl"))
include(joinpath("informationalcriteria.jl"))
##### Functions
export
    OptionsLikelihood,
    OptionsClusterization,
    OptionsMCMC,
    MCMCalgorithm,
    InformationalCriteria


end # module BayesianAnimalMovementModels
