# using Revise, FLPDevelopment, BrowseTables

# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Test"
# Directory = "/home/beatriz/mainen.flipping.5ht@gmail.com/Flipping/run_task_2/Dev_raw_data"
path = "/Volumes/GoogleDrive/My Drive/Flipping/run_task_2/Dev_raw_data"
path = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Pharmacology/BeatrizDevelopment"
Age_d, Age_p, Age_b, Age_s = process_dataset(path)
######### Filtering Datas #############
for df in (Age_p, Age_b, Age_s)
    df[!,:Exp_type] = [x in age_exp ? "Development" : "Projections" for x in df.MouseID]
    filter!(r->r.Exp_type == "Development",df)
    df[!,:Age] = [x in youngs ? "Young" : "Adult" for x in df.MouseID]
    gd = groupby(df,:Session)
    transform!(gd, :Streak => maximum => :Performance)
end

# for (df,name) in zip([Age_p, Age_b, Age_s],["Pokes", "Bouts", "Streaks"])
#     filename = joinpath(path, "Results_" * string(today()),"FullInfo" * name * ".csv")
#     CSV.write(filename,df)
# end
