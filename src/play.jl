include("Young_to_run2.jl")
include("Caspase_to_run.jl")
for df in (Age_p, Age_b, Age_s, Cas_p, Cas_b, Cas_s)
    filter!(r -> r.Protocol == "90/90" &&
    # r.MouseID != "CD09" && # biting, see B1_CD09_2020-07-13 minute30
    !(r.MouseID in first_females_group) &&
    r.ProtocolSession == 1,
    df)
end
# for df in (Cas_p, Cas_b, Cas_s)
#     filter!(r -> r.Gen == "Rbp4-cre", df)
# end
##
check = combine(groupby(Cas_s, :MouseID), :Streak => maximum, :Performance => union)
##
fCas_p = filter(r->r.PokeDur > 0.3 &&
    (r.Reward || ismissing(r.PostInterpoke) || (r.PostInterpoke > 0.1)) &&
    (r.Reward || r.PreInterpoke == 0 || ismissing(r.PreInterpoke) || (r.PreInterpoke > 0.1)),
    Cas_p)
    gd = groupby(fCas_p,[:MouseID,:Session,:Virus,:Gen])
    casdf = combine(gd) do dd
        process_streaks(dd)
    end
categorical!(casdf,:Virus)
levels!(casdf.Virus,["tdTomato", "Caspase"])
categorical!(casdf,:Gen)
levels!(casdf.Gen,["Wild Type", "Rbp4-cre"])
ALCas0 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak  + (1|MouseID)),casdf))
ALCas1 = fit!(LinearMixedModel(@formula(AfterLast ~ 1 + Streak + Virus + Gen & Virus + (1|MouseID)),casdf))
