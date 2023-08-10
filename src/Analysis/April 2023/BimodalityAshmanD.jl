BimCriteria(m1,s1,m2,s2) = abs(m1-m2) <= 2*minimum([s1,s2])
BimCriteria(d::GMM{Float64,Matrix{Float64}}) = BimCriteria(d.μ[1],d.Σ[1],d.μ[2],d.Σ[2])
AshmanD(m1,s1,m2,s2) =2^(1/2)((abs(m1-m2)/sqrt(s1^2+s2^2)))

Age_Bim_group = combine(groupby(Age_s,:Age),:LogDuration => (x->GMM(2, Vector(x); method=:split)) => :Bim)
transform!(Age_Bim_group, :Bim => ByRow(BimCriteria) => :BimCriteria)
Age_Bim_mouse = combine(groupby(Age_s,[:Age,:MouseID]),:LogDuration => (x->GMM(2, Vector(x); method=:split)) => :Bim)
transform!(Age_Bim_mouse, :Bim => ByRow(BimCriteria) => :BimCriteria)
Age_Bim_freq = combine(groupby(Age_Bim_mouse, :Age),:BimCriteria => sum => :Unimodal, nrow)
transform!(Age_Bim_freq, [:Unimodal,:nrow] => ByRow((b,n) -> b/n) => :Freq)

Age_Bim_group = combine(groupby(Age_k,:Age),:Diff => (x->GMM(2, Vector(x); method=:split)) => :Bim)
transform!(Age_Bim_group, :Bim => ByRow(BimCriteria) => :BimCriteria)

#=
We observe that the distribution of leaving times appears to be bimodal, devided in long and short leaving time.
This could reflect the fact that  animals are naive to the box environment and they all mightalternate trials of focused foraging
with trials mixed with investigating the box. Alternatively different animals might adopt different strategies, 
raising the possibility that differences across groups depends on the proportion of animals characterized by one particular 
strategy versus the other. To test within groups homogeneity, i.e. verify that individual animals within a group adopt both strategies,
we first use kmeans clustering to define long and short leaving time at the group level and categorize each trial accordingly.
Next we calculated the probability of performing a long trial at the group level and Perform Fisher's exact test of the null hypothesis 
that the probability of long leaving time at the individual and group levelels are equal. Using this approach we found that only one 
adult animal and one caspase animal had a significant lower chance to adopt long leaving times. Considering that our analysis showed
opposing tendencies in persistence between adults and caspase animals, these results indicate that, despite
the existence a potential distinction between short and long leaving time all animals perform both type of trials, resulting in 
an homogeneous behaviour within the groups.  


Next we calculated the difference in the frequency of short and long leaving time trials for each animals. Finally to test whether we can distinguish to separate strategies within group we tested whether 
the difference between the guasian mixture means was smaller than double the smaller standar deviation, a commonly used sufficient
criteria for unimodality (Behboodian J.,1970, doi:10.2307/1267357). 

we first fit a mixture of two gaussian on leaving time of each animals. Next, we tested 
whether the difference between means was smaller than double the smaller standar deviation, a commonly used sufficient 
criteria for bimodality (Behboodian J.,1970, doi:10.2307/1267357). Using this creiteria we found that
To confirm that all animals within
=#
