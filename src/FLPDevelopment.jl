module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, Dates, GLM, MixedModels, CSV, StatsBase, StatsPlots, HypothesisTests
import Distributions:Chisq
import Statistics: median

include("groups.jl")
include("utilities.jl")
include("DOFtests.jl")

export youngs, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export Likelyhood_Ratio_test, AIC_test, AICc_test
export median

end # module
