source("IS2/variants/standard.R")

#Make sure that any functions your LL calls are loaded here! (including packages)

result <- IS2(sampled, filter = "sample", n_cores = 25)
print(result$lw)
save(result, file = "IS2_1F.RData")