function table_coef(m)
    # main_ef = DataFrame(Variable = fixefnames(m),
    #     Coef = coef(m),
    #     Error = stderror(m))
    DataFrame(coeftable(m))
end

function table_ranef(m)
    main_ef = table_coef(m)
    ran_ef = DataFrame(only(raneftables(m)))
    for n in propertynames(ran_ef)[2:end]
        base_val = main_ef[findfirst(main_ef.Name .== string(n)), Symbol("Coef.")]
        ran_ef[!,Symbol("Net_"*string(n))] = ran_ef[:,Symbol(n)] .+ base_val
    end
    return ran_ef
end
