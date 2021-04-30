# using Revise, FLPDevelopment, BrowseTables

pathMac = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Pharmacology/CaspaseDevelopment"
pathLin = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Pharmacology/CaspaseDevelopment"
if ispath(pathMac)
    path = pathMac
elseif ispath(pathLin)
    path = pathLin
else
    error("unknown path")
end
##
Cas_d, Cas_p, Cas_b, Cas_s = process_dataset(path)

for df in (Cas_p, Cas_b, Cas_s)
    df[!,:Virus] = [get(VirusDict,x,"Missing") for x in df.MouseID]
    df[!,:Sex] = [x in females ? "F" : "M" for x in df.MouseID]
    # df[!,:PreInterpoke] = [ismissing(x) ? 0.0 : x for x in df.PreInterpoke]
    transform!(df, :Virus => categorical, renamecols=false)
    gd = groupby(df,:Session)
    transform!(gd, :Streak => maximum => :Performance)
end

levels!(Cas_p.Virus,["tdTomato", "Caspase"])
levels!(Cas_b.Virus,["tdTomato", "Caspase"])
levels!(Cas_s.Virus,["tdTomato", "Caspase"])

for dd in [Cas_s, Cas_p]
    dd[!,:Gen] = [g == "HET" ? "Rbp4-cre" : "Wild Type" for g in dd.Gen]
    dd[!,:Combo] = dd.Gen .* "\n" .* string.(dd.Virus)
    dd[!,:Group] = [r.Gen == "Wild Type" ? "Wild Type" : r.Gen .* "\n" .* string(r.Virus) for r in eachrow(dd)]
end

Afreq = countmap(Cas_s.AfterLast)
Aprob = Dict()
for (a,f) in Afreq
    Aprob[a] = round(f/nrow(Cas_s),digits = 5)
end
Cas_s[!,:IncorrectStart] = [!x for x in Cas_s.CorrectStart]
Cas_s[!,:IncorrectLeave] = [!x for x in Cas_s.CorrectLeave]
Cas_s[!,:P_AfterLast] = [Aprob[a] for a in Cas_s.AfterLast]
gd = groupby(Cas_s,:Session)
transform!(gd, :AfterLast => frequency)
transform!(gd, :AfterLast => (x -> x .<= quantile(x,0.95)) => :Limit)
transform!(gd, :Num_Rewards => cumsum => :Cum_Rewards)
Cas_s[!,:RewRate] = Cas_s.Cum_Rewards ./ Cas_s.Stop
Cas_s[!,:TrialRewRate] = Cas_s.Cum_Rewards ./ Cas_s.Streak
gdp = groupby(Cas_p,[:Session, :Streak])
transform!(gdp,:PokeIn => (x -> x .- x[1]) => :In)
transform!(gdp,[:PokeIn, :PokeOut] => ((x,y) -> y .- x[1]) => :Out)
# for (df,name) in zip([Cas_p, Cas_b, Cas_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end
