rm(list = ls())

interval = function(x, conf, ROUND) {
  # conf = 0.95
  # ROUND = 3
  # x = rnorm(100)
  al = 1 - conf
  int = quantile(x, c(al / 2, 1 - (al / 2)))
  int2 = paste0("(",
                round(int[1],ROUND),
                ",",
                round(int[2],ROUND),
                ")")
  return(int2)
}




#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
wd = getwd()
#Model 1--------
MODEL = c(1, 2,3,4,5,6)
for (M in MODEL) {

  file.path = paste0(wd,
                     "/Data_simulation/MC_estimates/mc_estimates_model_",
                     M,
                     ".rds")
  mc_estimates_model_1 <- readRDS(file.path)
  file.path.par = paste0(wd, "/Data_simulation/real_par_model_", M, ".txt")
  emv = read.table(file.path.par)

  mc_estimates_model_1$V3  = factor(mc_estimates_model_1$V3,
                                    levels = c(emv$par_names, "alpha"))

  estimates_model_1 <-
    mc_estimates_model_1 %>%
    group_by(V8, V6, V7, V2, V3,V4) %>%
    reframe(V5 = V5[!V5 %in% boxplot.stats(V5,1.40)$out]) %>%
    group_by(V8, V6, V7, V2, V3) %>%
    summarise(
      Par_real = mean(V4),
      Mean = round(mean(V5), 3),
      Median = round(median(V5), 3),
      #Interval = interval(V5,0.95,3),
      RB = round((Par_real - Mean) / Par_real, 3),
      MSE = round(mean(   (V5 - Par_real) ^ 2  ), 3),
      #Samples = length(V3)
    ) %>% data.frame()

  write.table(
    estimates_model_1,
    "Data_simulation/MC_estimates/results_statistics.txt",
    append = T,
    col.names = F,
    row.names = F,
  )

}
results_statistics = read.table("Data_simulation/MC_estimates/results_statistics.txt")
colnames(results_statistics) = c(
  "Model",
  "Method",
  "n",
  "alpha",
  "Parameter",
  "Par Real",
  "Mean",
  "Median",
 # "Interval",
  "RB",
  "MSE"
#  "Samples"
)
#write.csv(results_statistics,"Data_simulation/MC_estimates/results_statistics.csv")
write.csv(results_statistics,"Data_simulation/MC_estimates/results_statistics.csv")





