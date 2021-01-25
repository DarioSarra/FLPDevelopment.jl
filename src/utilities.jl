function fraction_true(v::AbstractArray{Bool})
    sum(v)/length(v)
end

function bootstrap_mean(v; sampling = BalancedSampling, nsample = 1000)
    bs = bootstrap(mean, v, sampling(nsample))
    first(bs.t0)
end

function frequency(v)
    counting = countmap(v)
    freq = Dict()
    for (a,f) in counting
        freq[a] = round(f/length(v),digits = 5)
    end
    return [freq[a] for a in v]
end

"""
    _notnan!(res,col)
remove NaN from columns with element type <:Real
"""
function _notnan!(res, col)
    @inbounds for (i, el) in enumerate(col)
        if typeof(el) <:Real
            res[i] &= !isnan(el)
        end
    end
    return nothing
end

function complete_vals(df::AbstractDataFrame, col::Colon=:)
    if ncol(df) == 0
        throw(ArgumentError("Unable to compute complete vals of a data frame with no columns"))
    end
    res = trues(size(df, 1))
    for i in 1:size(df, 2)
        _notnan!(res, df[!, i])
    end
    res
end

function complete_vals(df::AbstractDataFrame, col::DataFrames.ColumnIndex)
    res = trues(size(df, 1))
    _notnan!(res, df[!, col])
    res
end

complete_vals(df::AbstractDataFrame, cols::Union{AbstractVector, Regex, Not, Between, All}) =
    complete_vals(df[!, cols])

"""
    dropnan(df::AbstractDataFrame, cols =:)
Makes a copy of df = newdf.
Remove NaN values from dataframe newdf. Proceeds on all column if cols is not specified
"""
function dropnan(df::AbstractDataFrame, cols =:)
    newdf = df[complete_vals(df, cols), :]
    # disallowmissing && disallowmissing!(newdf, cols)
    newdf
end

"""
    dropnan(df::AbstractDataFrame, cols =:)
Remove NaN values from dataframe df. Proceeds on all column if cols is not specified
"""
function dropnan!(df::AbstractDataFrame,
                      cols=:)
    delete!(df, (!).(complete_vals(df, cols)))
    df
end

"""
    bin_axis(v; length = 50)

bin a continuous vector in a number of cases defined by length
return a vector
"""
function bin_axis(v; length = 50, unit_step = nothing)
    if isnothing(unit_step)
        r = range(extrema(v)...;length = 50)
        return [round(r[findfirst(r .>= x )], digits=1) for x in v]
    else
        r = range(extrema(v)...;step = unit_step)
        return [isnothing(findfirst(r .>= x )) ? last(r) + step(r) : r[findfirst(r .>= x )] for x in v]
    end
end

"""
    process_filtered_streak(df::AbstractDataFrame, var::Symbol, val <:Number)

filter a pokes dataframe to chop trials according to the first encountered value val in column var
return a streak dataframe from the chopped pokes dataframe
"""
function process_filtered_streak(df::AbstractDataFrame, var::Symbol, val ::Number)
    gd = groupby(df,[:Streak,:MouseID])
    filt_df_p = combine(gd) do dd
        check = findfirst(dd[:, var] .>= val)
        if isnothing(check)
            limit = nrow(dd)
        elseif check == 1
            limit = 1
        elseif check > 1
            limit = check - 1
        end
        dd[1:limit,:]
    end
    gd = groupby(filt_df_p,:MouseID)
    filt_df_s = combine(gd) do dd
        prov = DataFrame(FLPprocess.process_streaks(dd))
        prov[!,:Gen] .= dd[1,:Gen]
        return prov
    end
    if in(:Age,propertynames(df))
        filt_df_s[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in filt_df_s.MouseID]
    elseif in(:Virus,propertynames(df))
        filt_df_s[!,:Virus] = [get(VirusDict,x,"Missing") for x in filt_df_s.MouseID]
    end
    filt_df_s[!,:Sex] = [x in females ? "F" : "M" for x in filt_df_s.MouseID]
    filt_df_s[!,:PreInterpoke] = [ismissing(x) ? 0.0 : x for x in filt_df_s.PreInterpoke]
    gd = groupby(filt_df_s,:MouseID)
    transform!(gd, :Streak => maximum => :Performance)
    Afreq = countmap(filt_df_s.AfterLast)
    Aprob = Dict()
    for (a,f) in Afreq
        Aprob[a] = round(f/nrow(filt_df_s),digits = 5)
    end
    filt_df_s[!,:IncorrectStart] = [!x for x in filt_df_s.CorrectStart]
    filt_df_s[!,:IncorrectLeave] = [!x for x in filt_df_s.CorrectLeave]
    filt_df_s[!,:P_AfterLast] = [Aprob[a] for a in filt_df_s.AfterLast]
    gd = groupby(filt_df_s,:MouseID)
    transform!(gd, :AfterLast => frequency)
    transform!(gd, :Num_Rewards => cumsum => :Cum_Rewards)
    filt_df_s[!,:RewRate] = filt_df_s.Cum_Rewards ./ filt_df_s.Stop
    return filt_df_s, filt_df_p

end

function joinfilter(df_s,df_p,var, value)
    # filter streak according to the value
    f_df_s = filter(var => x -> x <= value, df_s)
    # use column MouseID and Streak to select from the PokesDatafraim
    f_df_p = leftjoin(f_df_s[:,[:MouseID, :Streak]],df_p,on = [:MouseID, :Streak])
    dropmissing!(f_df_p,:PreInterpoke, disallowmissing = true)
    filter!(r -> 0 < r.PreInterpoke, f_df_p)
    return f_df_s, f_df_p
end
