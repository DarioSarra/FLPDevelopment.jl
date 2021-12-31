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

export youngs, dario_youngs, females, sixty_days_old, age_exp, cnos_animals, wt_animals, below_thrs, VirusDict
export first_females_juveniles, second_females_juveniles
export first_females_group, second_females_group
export median, std
export frequency, fraction_true, bin_axis
export session_plot, maintitle!
export individual_summary, group_summary, bootstrapdf
export function_analysis, survivalrate_algorythm, cumulative_algorythm, hazardrate_algorythm, mediansurvival_analysis
export incorrect_fraction, median_ci_scatter
export add_bar!, add_info!, Difference

end # module
