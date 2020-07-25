# using Revise, FLPDevelopment, BrowseTables

# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Test"
# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data"
pathMac = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Pharmacology/BeatrizDevelopment"
pathLin = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/Datasets/Pharmacology/BeatrizDevelopment"
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
    df[!,:Exp_type] = [x in age_exp ? "Development" : "Projections" for x in df.MouseID]
    filter!(r->r.Exp_type == "Development",df)
    df[!,:Age] = [x in youngs ? "Young" : "Adult" for x in df.MouseID]
    gd = groupby(df,:Session)
    transform!(gd, :Streak => maximum => :Performance)
end

Afreq = countmap(Age_s.AfterLast)
Aprob = Dict()
for (a,f) in Afreq
    Aprob[a] = round(f/nrow(Age_s),digits = 5)
end
Age_s[!,:P_AfterLast] = [Aprob[a] for a in Age_s.AfterLast]

# for (df,name) in zip([Age_p, Age_b, Age_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end
