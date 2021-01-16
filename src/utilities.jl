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
