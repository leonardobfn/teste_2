METHOD = c("Method_1","Method_2")
MODEL = c(5,6)
MC = 3000
n. = c(84,144,168)
alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")

#alpha_value = alpha_value[c(1)]
wd <- getwd()
for (M in MODEL) {
  cat("Model=", M, "\n")
  estimates = NULL
  for (Method in METHOD) {
    #Method = METHOD[1]
    cat("Method=", Method, "\n")
    for (N. in n.) {
      #N. = n.[1]
      cat("N=", N., "\n")
      for (ALPHA in alpha_value) {
        #ALPHA = alpha_value[1]
        cat("ALPHA=", ALPHA, "\n")
        wd.estimates <-
          paste0(
            wd,
            "/Data_simulation/Model_",
            M,
            "/estimates/",
            Method,
            "/estimates/n",
            N.,
            "/",
            ALPHA
          )
        setwd(wd.estimates)
        for (i in 1:MC) {
          #i=1
          #cat("i=", i, "\n")
          names.files <-
            paste0(i, "_estimates_", ALPHA, "_", "n", N., ".txt")
          estimates.aux = read.table(names.files)
          na = length(which(is.na(estimates.aux$V5)) == T)
          if (na > 0) {
            next
          }
          estimates.aux$V3 <-
            stringr::str_replace(estimates.aux$V3, "_a", replacement = "_nu")
          estimates = rbind(estimates, estimates.aux)
        }
        setwd(wd)

      }
    }
  }

  setwd(wd)
  path = paste0(wd, "/Data_simulation/MC_estimates/")
  names.file = paste0("mc_estimates_model_", M, ".rds")
  setwd(path)
  saveRDS(estimates, names.file)
  setwd(wd)

}





dados = readRDS("estimates_model_1.rds")

require(dplyr)
require(tidyr)
head(dados)
est = dados %>% group_by(V8, V2, V6, V7, V3) %>%
  summarise(Par_real = mean(V4),
            Mean = mean(V5),
            Median = median(V5))
