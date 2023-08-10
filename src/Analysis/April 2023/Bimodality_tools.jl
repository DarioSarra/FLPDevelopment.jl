using Clustering, Distributions, GaussianMixtures
import GaussianMixtures: nparams
#Use k means to cluster datas in nk means
function kcluster_df!(df,nk,columns; new_name = nothing)
    features = collect(Array(df[:,columns])') #turn vector into array for clustering
    result = kmeans(features, nk)
    isnothing(new_name) ? (col_name = :Assignment) : (col_name = newname)
    df[!, col_name] = result.assignments
end
#following kmeans replace cluster numbers with labels
function categorise_kclusters!(df, original_col,assignment_col,levels::AbstractVector{String}; new_name = nothing)
    @assert length(unique(df[:,assignment_col])) == length(levels)
    clust_df = combine(groupby(df, assignment_col), original_col=>mean=>:MeanValues)
    sort!(clust_df, :MeanValues)
    clust_df[!,:CatDur] = CategoricalArray(levels, levels = levels)
    clust_dict = Dict(r.Assignment => r.CatDur for r in eachrow(clust_df))
    isnothing(new_name) ? (col_name = :Cluster) : (col_name = new_name) 
    df[!, col_name] = [get(clust_dict, x, missing) for x in df[:,assignment_col]]
end
## using Akaike to distinguish unimodal and bimodal
#functions
params_counter(model) = sum([*(size(getproperty(model, k))...) for k in (:μ,:Σ,:w)])
function AICc(loglikelihood, n_params, n_obs)
    logL = loglikelihood
    k = n_params
    n = n_obs
    −2*logL + 2*k +2*k*(k−1)/(n−k−1)
end
function AIC(loglikelihood, n_params, n_obs)
    logL = loglikelihood
    k = n_params
    −2*logL + 2*k
end
function AkaikeGMM(x,modes)
    n = length(x)
    gausmix = GMM(modes, x; method=:split)
    L = loglikelihood(MixtureModel(gausmix),x)
    K = params_counter(gausmix)
    AIC(L, K, n)
end

function AkaikeBimodality(df, var; group = nothing)
    isnothing(group) ? (df1 = df) : (df1 = groupby(df,group))
    combine(df1, 
    var => (x -> AkaikeGMM(Vector(x), 1)) => :Unimodal_AICc,
    var => (x -> AkaikeGMM(Vector(x), 2)) => :Bimodal_AICc,
    )
end

## using LikelihoodRatioTest to distinguish unimodal and bimodal
#funtions
function LRT_GMM(simple, full,x)
    simpleL = loglikelihood(MixtureModel(simple),x)
    simpleK = params_counter(simple)
    fullL = loglikelihood(MixtureModel(full),x)
    fullK = params_counter(full)
    LRT = -2(simpleL - fullL)
    degrees = fullK-simpleK
    ccdf(Chisq(degrees), LRT)
end

function compareGMMs(x)
    s_GMM = GMM(1, x; method = :split)
    m_GMM = GMM(2, x; method = :split)
    LRT_GMM(s_GMM,m_GMM,x)
end

function LogLikeGMM(x,modes)
    n = length(x)
    gausmix = GMM(modes, x)
    L = loglikelihood(MixtureModel(gausmix),x)
end
