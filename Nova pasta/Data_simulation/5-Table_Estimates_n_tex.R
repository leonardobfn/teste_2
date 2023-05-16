rm(list = ls())
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(xtable)
require(kableExtra)

model = read.csv("Data_simulation/MC_estimates/results_statistics.csv", sep =
                   ";")
paste0("$\\beta_", 0:2, "$")
labels1_3 = c(
  paste0("$\\beta_",0:2,"$"),
  paste0("$\\gamma_",0:1,"$"),
  paste0("$\\alpha$")
)
labels2_4 = c(
  paste0("$\\beta_",1:2,"$"),
  paste0("$\\gamma_",0:1,"$"),
  paste0("$\\alpha$")
)

labels5 = c(
  paste0("$\\beta_",0:2,"$"),
  paste0("$\\gamma_",1:2,"$"),
  paste0("$\\alpha$")
)

labels6 = c(
  paste0("$\\beta_",0:1,"$"),
  paste0("$\\gamma_",0:1,"$"),
  paste0("$\\alpha$")
)

MODEL = c(5)
ALPHA = c(0.35, 0.50, 0.65, 0.80, 0.95)
n = c(84,144,168)
for (M in MODEL) {
  #M=4
  if (M == 1 | M == 3) {
    LAB = labels1_3
  }
  if (M == 2 | M == 4) {
    LAB = labels2_4
  }
  if (M == 5) {
    LAB = labels5
  }
  if (M == 6) {
    LAB = labels6
  }

  for (N in n) {
    #N=144
    path.model = paste0("Model ",M)

    # model_sub = model %>% filter(Model == path.model, alpha == AL) %>%
    #   select(-X, -Model, -alpha)
    #
    model_sub = model %>% filter(Model == path.model, n == N) %>%
      select(-X, -Model, -n)


    model_sub$Parameter = rep(LAB,nrow(model_sub)/length(LAB))
    col.names = c("Method","$\\alpha$","Parameter","Par Real","Mean","Median","Relative Bias","Mean Squared Error")
    caption_model = paste0("Estimates of parameters to Model ",
                           M,
                           " and ",
                           "$n=",
                           N,
                           "$")
    label_model = paste0("Est_model_", M, "_n", N)
    #model_6 = model %>% filter( Model=="Model 6",n=="144") %>% select(-X,-X.1,-Model,-n)
    table_model = kbl(
      model_sub,
      align = "c",
      format = "latex"  ,
      booktabs = T,
      col.names = col.names,
      escape = F,
      row.names = F,
      longtable = F,
      caption = caption_model,
      label = label_model,
      position = "h",
      table.envir = "table",

    ) %>%
      #kable_paper(full_width = F) %>%
      column_spec(1:2, bold = T) %>%
      #kable_styling(latex_options = c("repeat_header")) %>%
      #row_spec(c(1:5, 11:15)-1, extra_latex_after = "\\rowcolor{gray}") %>%
      collapse_rows(columns = 1:2, valign = c("top"))
    #al.names = stringr::str_sub(AL,start = 3,end = 4)
    path_model = paste0("Data_simulation/Tables/Table_model_",
                        M,
                        "_n",
                        N,
                        ".tex")
    save_kable(table_model, path_model)
  }
}
