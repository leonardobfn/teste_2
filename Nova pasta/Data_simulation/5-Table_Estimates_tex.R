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

MODEL = c(4, 5, 6)
ALPHA = c(0.35, 0.50, 0.65, 0.80, 0.95)
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

  for (AL in ALPHA) {
    #AL = 0.35
    path.model = paste0("Model ",M)

    model_sub = model %>% filter(Model == path.model, alpha == AL) %>%
      select(-X, -Model, -alpha)

    # model_sub = model %>% filter(Model == path.model, n == 144) %>%
    #   select(-X, -Model, -n)


    model_sub$Parameter = rep(LAB,nrow(model_sub)/length(LAB))
    col.names = c(colnames(model_sub)[1:3], "Par Real","Mean","Median","Relative Bias","Mean Squared Error")
    caption_model = paste0("Estimates of parameters to Model ",
                           M,
                           " and ",
                           "$\\alpha=",
                           AL,
                           "$")
    label_model = paste0("Est_model_", M, "_Alpha", AL)
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
      kable_classic(full_width=F) %>%
      #kable_paper(full_width = F) %>%
      #kable_styling(latex_options = c("scale_down")) %>%
      column_spec(1:2, bold = T) %>%
      #row_spec(c(1:5, 11:15)-1, extra_latex_after = "\\rowcolor{gray}") %>%
      collapse_rows(columns = 1:2,
                    latex_hline = 'major',
                    custom_latex_hline = c(2:8),
                    row_group_label_position = 'identity',
                    valign = "middle")
    al.names = stringr::str_sub(AL,start = 3,end = 4)
    path_model = paste0("Data_simulation/Tables/Table_model_",
                        M,
                        "_Alpha",
                        al.names,
                        ".tex")
    save_kable(table_model, path_model)
  }
}

\cmidrule{3-5}
