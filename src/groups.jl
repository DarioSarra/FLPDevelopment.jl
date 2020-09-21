# this file contains informations about animals grouping focusing on analysis
# about the first session. There are information about 2 experiments here,
# either Development or Projections. The possible categories spanning both
# experiments are Young or Adults, CNO or Saline, HET or WT. This file defines
# the animal belonging to one value of each pair.
# Animals either belong to a defined group or the opposite by exclusion

# CNO GROUP WITHOUT CNO THE FIRST DAY to be excluded from analysis
const exclude_mice = ["CN21", "CN22","CN23","CN24","CN41","CN42","CN43"]

#animals lists
const Ncs = ["NC"* string(i) for i in collect(1:14)]
const Ccs = ["CC"* string(i) for i in [3,4,5,9,10,14,15,16]]
const Rcs = ["RC"* string(i) for i in [1,2,3,4]]
const Nbs = ["NB"* string(i) for i in [21,22,23,24,25,41,42,43,44,45]]
const Bns = ["BN"* string(i) for i in [21,22,23,24,41,42,43,44]]
const Its = ["IT"* string(i) for i in [21,22,23,24,41,42,43,44]]

# Age information
# all animals of these groups were tested for age factor
const age_exp = vcat(Nbs,Bns,Its)

# specifies which animals in age_exp are young
const youngs = ["BN21", "BN22", "BN41", "BN42", "IT21", "IT22","IT41", "IT42", "NB21", "NB22","NB23","NB41","NB42"]
const dario_youngs = ["RJ01", "RJ02", "RJ03", "RJ04", "RJ05", "RJ06", "RJ07", "RJ08", "RJ09", "RJ10", "RJ11", "RJ12",
    "RJ25", "RJ26", "RJ27", "RJ28", "RJ29", "RJ30", "RJ31", "RJ32", "RJ33", "RJ34", "RJ35", "RJ36"]

# specifies sex
const females = ["RJ25", "RJ26", "RJ27", "RJ28", "RJ29", "RJ30", "RJ31", "RJ32", "RJ33", "RJ34", "RJ35", "RJ36",
    "RJ37", "RJ38", "RJ39", "RJ40", "RJ41", "RJ42", "RJ43", "RJ44", "RJ45", "RJ46", "RJ47", "RJ48",
    "CD03", "CD04", "CD07", "CD08", "CD11", "CD12", "CD14", "CD15", "CD16", "CD20", "CD21", "CD22"]

# Treatment information
# animals treated with CNO
const cnos_animals = vcat(Ncs,Ccs,Rcs)

# Genotype informations
const wt_animals = vcat(["NC"* string(i) for i in [1,2,3,4,9,10,11,12]],age_exp)

# Animals with less than 15 trials in the first session
const development_exp_below15_day1 = ["BN21","BN23","BN41","BN42","BN44","IT23","NB25","NB43","NB45"]
const projection_exp_below15_day1 = ["CC9","NC10","NC12","NC14","RC3","RC8"]
const below_thrs = vcat(development_exp_below15_day1,projection_exp_below15_day1)

## Caspase Virus

const VirusDict = Dict(
    "CD06" => "Caspase",
    "CD09" => "Caspase",
    "CD10" => "Caspase",
    "CD11" => "Caspase",
    "CD13" => "Caspase",
    "CD14" => "Caspase",
    "CD15" => "Caspase",
    "CD16" => "Caspase",
    "CD05" => "Caspase",
    "CD07" => "Caspase",
    "CD08" => "Caspase",
    "CD12" => "Caspase",
    "CD02" => "tdTomato",
    "CD04" => "tdTomato",
    "CD17" => "tdTomato",
    "CD18" => "tdTomato",
    "CD19" => "tdTomato",
    "CD20" => "tdTomato",
    "CD21" => "tdTomato",
    "CD22" => "tdTomato",
    "CD01" => "tdTomato",
    "CD03" => "tdTomato",
    "CD23" => "tdTomato",
    "CD24" => "tdTomato"
)
