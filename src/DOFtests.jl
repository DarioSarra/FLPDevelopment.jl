function Likelyhood_Ratio_test(simple,full)
    χdegrees = dof(full) - dof(simple)
    χ = deviance(simple) - deviance(full)
    Pχ = ccdf(Distributions.Chisq(χdegrees), χ)
    res = DataFrame(
        Formula = [simple.formula, full.formula],
        ModelDof = [dof(simple), dof(full)],
        Deviance = [deviance(simple), deviance(full)],
        Χ² = [nothing, χ],
        Χ²Dof = [nothing, χdegrees],
        PΧ = [nothing,Pχ])
    rename!(res, [:ModelDof => Symbol("Model-dof"), :PΧ=> Symbol("P(>Χ²)")])
    show(res)
end

function mixture_Likelyhood_Ratio_test(simple,full,vec)
    χdegrees = length(full.components)*2 - length(simple.components)*2
    χ = 2* loglikelihood(simple,vec)/loglikelihood(full,vec)
    Pχ = ccdf(Distributions.Chisq(χdegrees), χ)
end

function AIC_test(simple, full)
    exp((aic(full) - aic(simple))/2)
end

function AICc(model)
    aic(model) + ((2*(dof(model)^2) + 2*dof(model))/(nobs(model) - dof(model) - 1))
end
function AICc_test(simple, full)
    exp((AICc(full) - AICc(simple))/2)
end
