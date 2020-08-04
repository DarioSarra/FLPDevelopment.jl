module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, CSV, Dates
@reexport using StatsBase, HypothesisTests, Bootstrap, StatsPlots
# @reexport using GLM, MixedModels
import Distributions:Chisq, Normal
import Statistics: median, std

include("groups.jl")
include("utilities.jl")
include("DOFtests.jl")
include("Analysis_type.jl")
include("pokes_plotting.jl")

export youngs, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export Likelyhood_Ratio_test, AIC_test, AICc_test
export median, std, dvAnalysis
export fraction_true, bootstrap_mean, mouse_summary,  group_summary, dvplot
export poke_plot!

end # module
