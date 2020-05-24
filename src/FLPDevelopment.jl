module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, Dates, MixedModels, CSV, StatsBase, StatsPlots, HypothesisTests
import Distributions:Chisq


include("groups.jl")
include("DOFtests.jl")

export youngs, age_exp, cnos_animals, wt_animals, below_thrs
export Likelyhood_Ratio_test, AIC_test, AICc_test

end # module
