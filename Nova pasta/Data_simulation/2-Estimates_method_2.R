

rm(list = ls())
devtools::load_all() # meu pacote
#devtools::install()
#require(RThesis)
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
n. = c(144,84,168)
#n. = c(144)
logs = NULL
MODEL = c(8)
MC = 3000
cpus <- 10
wd. = getwd()
alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")
#alpha_value = alpha_value[c(5)]
require(snowfall)
sfInit(cpus = cpus,
       type = "SOCK",
       parallel = TRUE)
sfExportAll()
#sfLibrary(RThesis)
sfLibrary(tidyr)
sfLibrary(dplyr)
sfLibrary(extraDistr)
sfLibrary(Formula)
for (M in MODEL) {
  if (M == 2 | M == 4) {
    FORMULA <- RH ~ sent + cost - 1 | semester
  }
  if (M == 3) {
    FORMULA <- RH ~ sent + cost | semester
  }
  if (M == 5) {
    FORMULA <- RH ~ sent + cost | semester-1
  }
  if (M == 6) {
    FORMULA <- RH ~ t | semester
  }
  if (M == 7|M==8) {
    FORMULA <- RH ~ log(t) | semester
  }
  for (k in n.) {
    for (alpha. in alpha_value) {
      tic <- tictoc::tic()
      sfLapply(
        1:MC,
        fun = snowfall.estimates_Method_2,
        model = M,
        alpha.value = alpha.,
        formula = FORMULA,
        wd = wd.,
        n = k
      )
      toc <- tictoc::toc()
      logs.aux = data.frame(
        MC = MC,
        cpus = cpus,
        n = k,
        alpha = paste0("0.", stringr::str_sub(alpha., 6, 7)),
        Model = M,
        Seconds = ((stringr::str_sub(
          toc$callback_msg,
          start = -stringr::str_length(toc$callback_msg),
          end = -13
        ) %>% as.numeric()))
      )
      #logs = rbind(logs, logs.aux)
      write.table(
        logs.aux,
        "Data_simulation/logs_estimates_method_2.txt",
        append = T,
        col.names = F,
        row.names = F
      )
    }

  }
}
sfStop()
