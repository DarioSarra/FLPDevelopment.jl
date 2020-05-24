function Likelyhood_Ratio_test(simple,full)
    degrees = dof(full) - dof(simple)
    ccdf(Distributions.Chisq(degrees), deviance(simple) - deviance(full))
end

function AIC_test(candidate, simpler)
    exp((aic(candidate) - aic(simpler))/2)
end
function AICc(model)
    aic(model) + ((2*(dof(model)^2) + 2*dof(model))/(nobs(model) - dof(model) - 1))
end
function AICc_test(candidate, simpler)
    exp((AICc(candidate) - AICc(simpler))/2)
end
