# using Revise, FLPDevelopment, BrowseTables

# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Test"
# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data"
pathMac = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Pharmacology/DarioDevelopment"
pathLin = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Pharmacology/DarioDevelopment"
if ispath(pathMac)
    path = pathMac
elseif ispath(pathLin)
    path = pathLin
else
    error("unknown path")
end
##
Age_d, Age_p, Age_b, Age_s = process_dataset(path)
######### Filtering Datas #############
for df in (Age_p, Age_b, Age_s)
    df[!,:Age] = [x in dario_youngs ? "Juveniles" : "Adults" for x in df.MouseID]
    df[!,:Sex] = [x in females ? "F" : "M" for x in df.MouseID]
    # df[!,:PreInterpoke] = [ismissing(x) ? 0.0 : x for x in df.PreInterpoke]
    transform!(df, :Age => categorical, renamecols=false)
    gd = groupby(df,:Session)
    transform!(gd, :Streak => maximum => :Performance)
end

levels!(Age_p.Age,["Adults", "Juveniles"])
levels!(Age_b.Age,["Adults", "Juveniles"])
levels!(Age_s.Age,["Adults", "Juveniles"])

Afreq = countmap(Age_s.AfterLast)
Aprob = Dict()
for (a,f) in Afreq
    Aprob[a] = round(f/nrow(Age_s),digits = 5)
end
Age_s[!,:IncorrectStart] = [!x for x in Age_s.CorrectStart]
Age_s[!,:IncorrectLeave] = [!x for x in Age_s.CorrectLeave]
Age_s[!,:P_AfterLast] = [Aprob[a] for a in Age_s.AfterLast]
gd = groupby(Age_s,:Session)
transform!(gd, :AfterLast => frequency)
transform!(gd, :AfterLast => (x -> x .<= quantile(x,0.95)) => :Limit)
transform!(gd, :Num_Rewards => cumsum => :Cum_Rewards)
Age_s[!,:RewRate] = Age_s.Cum_Rewards ./ Age_s.Stop
Age_s[!,:TrialRewRate] = Age_s.Cum_Rewards ./ Age_s.Streak
gdp = groupby(Age_p,[:Session, :Streak])
transform!(gdp,:PokeIn => (x -> x .- x[1]) => :In)
transform!(gdp,[:PokeIn, :PokeOut] => ((x,y) -> y .- x[1]) => :Out)
# for (df,name) in zip([Age_p, Age_b, Age_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end
