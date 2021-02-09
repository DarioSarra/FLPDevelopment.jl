function Likelyhood_Ratio_test(simple,full)
    degrees = dof(full) - dof(simple)
    ccdf(Distributions.Chisq(degrees), deviance(simple) - deviance(full))
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
