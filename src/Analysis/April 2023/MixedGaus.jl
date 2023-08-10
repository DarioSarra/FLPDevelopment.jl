using Distributions, GaussianMixtures
## using LikelihoodRatioTest to distinguish unimodal and bimodal
#funtions
function LRT_GMM(simple, full,x)
    simpleL = loglikelihood(MixtureModel(simple),x)
    simpleK = nparams(simple)
    fullL = loglikelihood(MixtureModel(full),x)
    fullK = nparams(full)
    LRT = -2(simpleL - fullL)
    degrees = fullK-simpleK
    ccdf(Chisq(degrees), LRT)
end

function compareGMMs(x)
    s_GMM = GMM(1, x)
    m_GMM = GMM(2, x)
    LRT_GMM(s_GMM,m_GMM,x)
end
#age
Age_Gauss = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => (x -> compareGMMs(Vector(x))) => :PGaussian)
transform!(Age_Gauss, :PGaussian => ByRow(x-> x<0.05) => :Bimodal)
sort!(Age_Gauss,:Age)
Age_Gaus_sum = combine(groupby(Age_Gauss, :Age), :Bimodal => sum => :Bimodal, nrow => :NAnimals)
Age_group = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => median => :LogDuration)
Age_Gaussgroup = combine(groupby(Age_group,:Age), :LogDuration => (x -> compareGMMs(Vector(x))) => :PGaussian)
transform!(Age_Gaussgroup, :PGaussian => ByRow(x-> x<0.05) => :Bimodal)

Cas_Gauss = combine(groupby(Cas_s,[:Virus,:MouseID]), :LogDuration => (x -> compareGMMs(Vector(x))) => :PGaussian)
transform!(Cas_Gauss, :PGaussian => ByRow(x-> x<0.05) => :Bimodal)
sort!(Cas_Gauss,:Virus)
Cas_Gaus_sum = combine(groupby(Cas_Gauss, :Virus), :Bimodal => sum => :Bimodal, nrow => :NAnimals)
Cas_group = combine(groupby(Cas_s,[:Virus,:MouseID]), :LogDuration => median => :LogDuration)
Cas_Gaussgroup = combine(groupby(Cas_group,:Virus), :LogDuration => (x -> compareGMMs(Vector(x))) => :PGaussian)
transform!(Cas_Gaussgroup, :PGaussian => ByRow(x-> x<0.05) => :Bimodal)
## using Akaike to distinguish unimodal and bimodal
#functions
nparams(model) = sum([*(size(getproperty(model, k))...) for k in (:μ,:Σ,:w)])
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
    gausmix = GMM(modes, x)
    L = loglikelihood(MixtureModel(gausmix),x)
    K = nparams(gausmix)
    AIC(L, K, n)
end
# Age
Age_AIC = combine(groupby(Age_s,[:Age,:MouseID]), 
    :LogDuration => (x -> AkaikeGMM(Vector(x), 1)) => :Unimodal_AICc,
    :LogDuration => (x -> AkaikeGMM(Vector(x), 2)) => :Bimodal_AICc,
    )
transform!(Age_AIC, [:Unimodal_AICc, :Bimodal_AICc] => ByRow((u,b) -> b < u ) => :Bimodal)
# open_html_table(Age_AIC)
Age_AIC_sum = combine(groupby(Age_AIC, :Age), :Bimodal => sum => :Bimodal, nrow => :NAnimals)

Age_group = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => median => :LogDuration)
Age_groupAIC = combine(groupby(Age_group,[:Age]), 
    :LogDuration => (x -> AkaikeGMM(Vector(x), 1)) => :Unimodal_AICc,
    :LogDuration => (x -> AkaikeGMM(Vector(x), 2)) => :Bimodal_AICc,
    )
transform!(Age_groupAIC, [:Unimodal_AICc, :Bimodal_AICc] => ByRow((u,b) -> b < u ) => :Bimodal)
#Virus
Cas_AIC = combine(groupby(Cas_s,[:Virus,:MouseID]), 
    :LogDuration => (x -> AkaikeGMM(Vector(x), 1)) => :Unimodal_AICc,
    :LogDuration => (x -> AkaikeGMM(Vector(x), 2)) => :Bimodal_AICc,
    )
transform!(Cas_AIC, [:Unimodal_AICc, :Bimodal_AICc] => ByRow((u,b) -> b < u ) => :Bimodal)
# open_html_table(Cas_AIC)
Cas_AIC_sum = combine(groupby(Cas_AIC, :Virus), :Bimodal => sum => :Bimodal, nrow => :NAnimals)

Cas_group = combine(groupby(Cas_s,[:Virus,:MouseID]), :LogDuration => mean => :LogDuration)
Cas_groupAIC = combine(groupby(Cas_group,[:Virus]), 
    :LogDuration => (x -> AkaikeGMM(Vector(x), 1)) => :Unimodal_AICc,
    :LogDuration => (x -> AkaikeGMM(Vector(x), 2)) => :Bimodal_AICc,
    )
transform!(Cas_groupAIC, [:Unimodal_AICc, :Bimodal_AICc] => ByRow((u,b) -> b < u ) => :Bimodal)


##
using UncertainData, KernelDensity
#transform Log duration data in uncertain value type
Age_JBT = combine(groupby(Age_s,[:Age,:MouseID]), 
    :LogDuration => (x -> UncertainValue(UnivariateKDE,Vector(x))) => :UV
)
#perform JarqueBeraTest on each animal uncertain value, extract pvalue, and outcome
transform!(Age_JBT, :UV => ByRow(JarqueBeraTest) => :JB) 
transform!(Age_JBT,:JB => ByRow(pvalue) => :P)
transform!(Age_JBT,:P => ByRow(x-> x<0.05) => :Bimodals)
# summarise results per group counting how many non-normally distributed animals there are 
Age_JBTsum = combine(groupby(Age_JBT, :Age), :Bimodals => sum => :Bimodals, nrow => :Total)
open_html_table(Age_JBT)
##
Age_dd = combine(groupby(Age_s,[:Age,:MouseID]), :LogDuration => mean => :LogDuration)
Age_JBTgroup = combine(groupby(Age_dd,:Age),
    :LogDuration => (x -> UncertainValue(UnivariateKDE,Vector(x))) => :UV
)
transform!(Age_JBTgroup, :UV => ByRow(JarqueBeraTest) => :JB)
transform!(Age_JBTgroup,:JB => ByRow(pvalue) => :P) 
transform!(Age_JBTgroup,:P => ByRow(x-> x<0.05) => :Bimodals)
Age_JBTgroup[:,[:Age,:P,:Bimodals]]
open_html_table(Age_JBTgroup)