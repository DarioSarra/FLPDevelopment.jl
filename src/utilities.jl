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
