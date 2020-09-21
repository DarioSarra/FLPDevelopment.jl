module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, CSV, Dates
@reexport using StatsBase, HypothesisTests, Bootstrap, KernelDensity, StatsPlots, Optim
@reexport using Distributions
import Statistics: median, std

include("groups.jl")
include("DOFtests.jl")
include("Analysis_type.jl")
include("utilities.jl")
include("test_plotting.jl")
include("pokes_plotting.jl")
include("Mixture_fit.jl")

export youngs,dario_youngs, females, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export Likelyhood_Ratio_test, AIC_test, AICc_test
export median, std, kde, dvAnalysis
export fraction_true, bootstrap_mean, frequency
export mouse_summary,  group_summary, dvplot
export poke_plot!
export mixture_gamma, mixture_gamma_weighted, mixture_exp, mixture_exp_weighted

end # module
