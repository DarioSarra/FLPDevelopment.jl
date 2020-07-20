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

# for (df,name) in zip([Cas_p, Cas_b, Cas_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end
