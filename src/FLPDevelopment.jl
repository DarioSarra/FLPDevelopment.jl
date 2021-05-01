module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, CSV, Dates, CategoricalArrays
@reexport using StatsBase, HypothesisTests, Bootstrap, Optim, GLM, MixedModels, Random, Survival#GaussianMixtures
@reexport using Distributions, KernelDensity, StatsPlots, StatsPlots.PlotMeasures, LaTeXStrings
import Statistics: median, std

include("groups.jl")
include("DOFtests.jl")
include("Analysis_type.jl")
include("utilities.jl")
include("test_plotting.jl")
include("pokes_plotting.jl")
include("Mixture_fit.jl")
include("Figures_plotting.jl")

export youngs,dario_youngs, females, sixty_days_old, first_females_juveniles, second_females_juveniles, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export first_females_group, second_females_group
export Likelyhood_Ratio_test, AIC_test, AICc_test
export median, std, kde, dvAnalysis, DoubleAnalysis
export fraction_true, bootstrap_mean, frequency, dropnan, dropnan!, bin_axis, reallignpokes, rm_interpokes, reprocess_streaks, filter_pokestreak
export Heatmap_group, Heatmap_plot, Heatmap_difference, Leave_plots, bootstrapdf,leave_modelplt
export mixture_gamma, mixture_gamma_weighted, mixture_exp, mixture_exp_weighted, mixture_gamma_exp
export poke_plot!, maintitle!
export individual_summary, summary_xy, group_kde, group_cdf, group_frequency, group_summary, dvplot, check_distribution, check_distributions
export mediansurvival_analysis, survivalrate_algorythm, cumulative_algorythm, hazardrate_algorythm, function_analysis
export incorrect_fraction, median_ci_scatter
export add_bar!, add_info!, Difference

end # module
