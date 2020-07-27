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
    gd = groupby(df,:Session)
    transform!(gd, :Streak => maximum => :Performance)
end

for df in (Cas_p, Cas_b, Cas_s)
    df[!,:Virus] = [get(VirusDict,x,"Missing") for x in df.MouseID]
end

Cas_s[!,:Gen] = [g == "HET" ? "Rbp4-cre" : "Wild Type" for g in Cas_s.Gen]
Cas_s[!,:Combo] = Cas_s.Gen .* "\n" .* Cas_s.Virus
Cas_s[!,:Group] = [r.Gen == "Wild Type" ? "Wild Type" : r.Gen .* "\n" .* r.Virus for r in eachrow(Cas_s)]
# categorical!(Cas_s,[:Gen,:Virus,:Group])
# levels!(Cas_s.Gen,["Wild Type","Rbp4-cre",])
# levels!(Cas_s.Virus,["tdTomato", "Caspase"])
# levels!(Cas_s.Group,["Rbp4-cre\ntdTomato", "Rbp4-cre\nCaspase", "Wild Type",])
Afreq = countmap(Cas_s.AfterLast)
Aprob = Dict()
for (a,f) in Afreq
    Aprob[a] = round(f/nrow(Cas_s),digits = 5)
end
Cas_s[!,:P_AfterLast] = [Aprob[a] for a in Cas_s.AfterLast]

# for (df,name) in zip([Cas_p, Cas_b, Cas_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end