module FLPDevelopment

using Reexport
@reexport using FLPprocess, DataFrames, CSV, Dates, CategoricalArrays
@reexport using StatsBase, HypothesisTests, Bootstrap, MixedModels, Random, Survival
@reexport using StatsPlots, StatsPlots.PlotMeasures, LaTeXStrings
import Statistics: median, std

include("groups.jl")
include("utilities.jl")
include("pokes_plotting.jl")
include("Analysis_fun.jl")
include("Plot_fun.jl")
# include("DOFtests.jl")
# include("Analysis_type.jl")
# include("test_plotting.jl")
# include("Mixture_fit.jl")

export youngs, dario_youngs, females, sixty_days_old, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export first_females_juveniles, second_females_juveniles
export first_females_group, second_females_group
export median, std#, kde, dvAnalysis, DoubleAnalysis
#export Likelyhood_Ratio_test, AIC_test, AICc_test
export frequency, fraction_true, bin_axis#, bootstrap_mean, dropnan, dropnan!, bin_axis, reallignpokes, rm_interpokes, reprocess_streaks, filter_pokestreak
#export Heatmap_group, Heatmap_plot, Heatmap_difference, Leave_plots, bootstrapdf,leave_modelplt
#export mixture_gamma, mixture_gamma_weighted, mixture_exp, mixture_exp_weighted, mixture_gamma_exp
export session_plot, maintitle!
export individual_summary, group_summary#, summary_xy, group_kde, group_cdf, group_frequency, dvplot, check_distribution, check_distributions
export function_analysis, survivalrate_algorythm, cumulative_algorythm, hazardrate_algorythm, mediansurvival_analysis
export incorrect_fraction, median_ci_scatter
export add_bar!, add_info!, Difference

end # module
