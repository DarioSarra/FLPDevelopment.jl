function fraction_true(v::AbstractArray{Bool})
    sum(v)/length(v)
end
function incorrect_fraction(v::AbstractArray{Bool})
    sum(v)/(length(v)-sum(v))
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
        r = range(extrema(v)...;length = length)
        return [round(r[findfirst(r .>= x )], digits=1) for x in v]
    else
        r = range(extrema(v)...;step = unit_step)
        return [isnothing(findfirst(r .>= x )) ? last(r) + step(r) : r[findfirst(r .>= x )] for x in v]
    end
end

"""
    quantile_vec(vec, bins)

retunr the quantile position, from 1 to bins, of each element of a vectir vec
"""
function quantile_vec(vec, bins)
    step = 1/(bins)
    qbins = round.(collect(step .* (1:bins)), digits = 3)
    qs = quantile(vec,qbins)
    pre = [findfirst(x .<= qs) for x in vec]
    [isnothing(x) ? 1 : qbins[x] for x in pre]
end

function check_group(df)
    if in(:Age,propertynames(df))
        return grouping = :Age
    elseif in(:Virus,propertynames(df))
        return grouping = :Virus
    else
        error("Can't identify 2 groups for subtraction")
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
    # use column MouseID and Streak to select from the PokesDataframe
    f_df_p = leftjoin(f_df_s[:,[:MouseID, :Streak]],df_p,on = [:MouseID, :Streak])
    dropmissing!(f_df_p,:PreInterpoke, disallowmissing = true)
    filter!(r -> 0 < r.PreInterpoke, f_df_p)
    return f_df_s, f_df_p
end

"""
    `summary_df(streak_df, pokes_df)`
summarise per mouse info
"""
function summarydf(streak_df, pokes_df)
    grp = :Age in propertynames(streak_df) ? :Age : :Virus
    gds = groupby(streak_df,:MouseID)
    Summary_df = combine(gds,grp => first => grp,
        :Sex => first => :Sex,
        :Num_pokes => sum => :Pokes,
        :Streak => maximum => :Trials,
        :AfterLast => mean => :AvgAfterLast,
        :Num_Rewards => sum => :Rewards,
        :IncorrectLeave => sum => :ErrorTrials,
        :Travel_to => mean => :AvgTravel)
    Summary_df[!,:AvgAfterLast] = round.(Summary_df.AvgAfterLast, digits = 2)
    Summary_df[!,:AvgTravel] = round.(Summary_df.AvgTravel, digits = 2)
    Summary_df[!,:Perr] = Summary_df.ErrorTrials ./ Summary_df.Trials
    Summary_df[!,:Perr] = round.(Summary_df.Perr, digits = 2)
    Summary_df[!,:RewxPoke] = round.(Summary_df.Rewards ./ Summary_df.Pokes , digits = 2)
    Summary_df[!,:RewxTrial] = round.(Summary_df.Rewards ./ Summary_df.Trials , digits = 2)
    pokes_df = filter(r-> r.PreInterpoke > 0, pokes_df)
    gdp = groupby(pokes_df,:MouseID)
    Summary_df = leftjoin(Summary_df,combine(gdp, :PreInterpoke => mean => :AvgInterpoke), on = :MouseID)
    Summary_df[!,:AvgInterpoke] = round.(Summary_df.AvgInterpoke, digits = 2)
    sort!(Summary_df,grp)
    Summary_df = Summary_df[:,[:MouseID, grp, :Sex, :AvgAfterLast,
    :Pokes, :Trials, :ErrorTrials, :Rewards,
    :AvgInterpoke, :AvgTravel,
    :Perr,:RewxPoke,:RewxTrial]]
    return Summary_df
end

"""
    `reallignpokes(rv,pv)`

Using the reward vector, rv, and the poke_in vector of a trial
return the poke_in time realligned to the last reward
"""

function reallignpokes(rv,pv)
    s = findlast(rv)
    if isnothing(s)
        return pv
    else
        return pv .- pv[s]
    end
end

function reallignpokes(df::DataFrames.AbstractDataFrame)
    gd = groupby(df, [:MouseID,:Streak])
    pokes = transform(gd, [:Reward,:PokeOut] => ((rv,pv) -> reallignpokes(rv,pv)) => :PO_LR)
    pokes[!, :PI_LR] = pokes[:,:PO_LR] .- pokes[:,:PokeDur]
    # pokes[!,:PI_LR] = Int64.(round.(pokes[:,:PI_LR]))
    # pokes[!,:PO_LR] = Int64.(round.(pokes[:,:PO_LR]))
    pokes[!,:PI_LR] = pokes[:,:PI_LR]
    pokes[!,:PO_LR] = pokes[:,:PO_LR]
    return pokes
end

"""
    `rm_interpokes(session, shortlimit = 0.06, longlimit = 7)`

merge pokes separated by an interpoke interval shorter than shortimit
remove from the streak pokes that occures after an interpoke interval longer than longlimit
"""
function rm_shortinterpoke(trial::DataFrames.AbstractDataFrame; limit = 0.06)
    group_var = :Age in propertynames(trial) ? :Age : :Virus
    Same_vals = [:Session, :Streak, :Side, :Day, :MouseID, :Gen,
            :ExpSession, :ProtocolSession, group_var, :Sex, :Performance]
    PreviousPoke_vals = [:Poke, :SideHigh, :PokeHierarchy, :PokeIn, :PreInterpoke, :PokeNum, :PokeInStreak, :In]
    NextPoke_vals = [:Bout, :Leave, :PokeOut, :PostInterpoke, :Out]
    # find the index of interpokes shorter than 50ms.
    # It skips missing values because they are the last poke in a trial, therefore lacking a PostInterpoke value
    shortcases = findall((trial.PostInterpoke .< limit) .& (.!collect((ismissing.(trial.PostInterpoke)))))
    if length(shortcases) == 0
        return trial
    else
        for i in shortcases
            # println("Size = $(nrow(trial)), case = $i, Poke = $(trial[i,:Poke]), Session = $(trial[i,:Session])")
            check = mapcols(x -> x[1] == x[2], trial[i:i+1,Same_vals])
            if !all(check[1,:])
                # checks that the column that shouldn't change have the same value in the ith and ith+1 case
                error("some supposedely fixed values are dufferent")
            end
            # replace values in the ith case to values of the i+1 case in columns that need to be extended to the second poke
            trial[i,NextPoke_vals] = trial[i+1,NextPoke_vals]
            # assign reward true if any of the ith and ith+1 case was rewarded
            trial[i,:Reward] = any(trial[i:i+1,:Reward])
            # recalculate values
            trial[i,:PokeDur] = trial[i+1,:PokeOut] - trial[i,:PokeIn]
            trial[i+1,[:PokeIn,:PokeOut]] .= NaN
        end
        trial = filter(r -> !isnan(r.PokeIn), trial)
        trial[:,:PokeInStreak] = collect(1:nrow(trial))
        return trial
    end
end

function rm_longinterpoke(dd::DataFrames.AbstractDataFrame; limit = 7)
    longcase = findfirst((dd.PostInterpoke .> limit) .& (.!collect((ismissing.(dd.PostInterpoke)))))
    if isnothing(longcase)
        return dd
    else
        # println("Size = $(nrow(dd)), case = $longcase, Poke = $(dd[longcase,:Poke])")
        dd = dd[1:longcase,:]
        dd[end,:PostInterpoke] = missing
        return dd
    end
end

function rm_interpokes(session::DataFrames.AbstractDataFrame; shortlimit = 0.06, longlimit = 7)
    res3 = combine(groupby(session,:Streak)) do dd
        res1 = rm_shortinterpoke(dd; limit = shortlimit)
        res2 = rm_longinterpoke(res1; limit = longlimit)
    end
    transform!(res3,:Poke => (x -> collect(1:length(x))) => :Poke)
end


"""
    process_filtered_streak(df::AbstractDataFrame, var::Symbol, val <:Number)

filter a pokes dataframe to chop trials according to the first encountered value val in column var
return a streak dataframe from the chopped pokes dataframe
"""
function reprocess_streaks(pokes_df::AbstractDataFrame)
    gd = groupby(pokes_df,:MouseID)
    streak_df = combine(gd) do dd
        prov = DataFrame(FLPprocess.process_streaks(dd))
        prov[!,:Gen] .= dd[1,:Gen]
        return prov
    end
    if in(:Age,propertynames(pokes_df))
        streak_df[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in streak_df.MouseID]
    elseif in(:Virus,propertynames(pokes_df))
        streak_df[!,:Virus] = [get(VirusDict,x,"Missing") for x in streak_df.MouseID]
    end
    streak_df[!,:Sex] = [x in females ? "F" : "M" for x in streak_df.MouseID]
    streak_df[!,:PreInterpoke] = [ismissing(x) ? 0.0 : x for x in streak_df.PreInterpoke]
    gd = groupby(streak_df,:MouseID)
    transform!(gd, :Streak => maximum => :Performance)
    Afreq = countmap(streak_df.AfterLast)
    Aprob = Dict()
    for (a,f) in Afreq
        Aprob[a] = round(f/nrow(streak_df),digits = 5)
    end
    streak_df[!,:IncorrectStart] = [!x for x in streak_df.CorrectStart]
    streak_df[!,:IncorrectLeave] = [!x for x in streak_df.CorrectLeave]
    streak_df[!,:P_AfterLast] = [Aprob[a] for a in streak_df.AfterLast]
    gd = groupby(streak_df,:MouseID)
    transform!(gd, :AfterLast => frequency)
    transform!(gd, :Num_Rewards => cumsum => :Cum_Rewards)
    streak_df[!,:RewRate] = streak_df.Cum_Rewards ./ streak_df.Stop
    return streak_df
end

function filter_pokestreak(pokedf; min = 0.05, max = 50000)
    Remove_vals = [:Stim, :StimFreq, :Wall, :Block, :Delta, :Turn, :StreakInBlock, :ReverseStreak]
    N_pokes = pokedf[:,Not(Remove_vals)]
    F_pokes = combine(groupby(N_pokes, :Session)) do dd
        rm_interpokes(copy(dd); shortlimit = 0.1, longlimit = 50000)
    end
    F_streaks = reprocess_streaks(F_pokes)
    return F_pokes, F_streaks
end
