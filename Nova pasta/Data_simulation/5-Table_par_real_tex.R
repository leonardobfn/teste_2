rm(list = ls())
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(xtable)
require(kableExtra)

model = read.csv("Data_simulation/MC_estimates/results_statistics.csv", sep =
                   ";")
labels1_3 = c(
  paste0("$\\beta_",0:2,"$"),
  paste0("$\\gamma_",0:1,"$")
  #paste0("$\\alpha$")
)
labels2_4 = c(
  paste0("$\\beta_",1:2,"$"),
  paste0("$\\gamma_",0:1,"$")
  #paste0("$\\alpha$")
)

labels5 = c(
  paste0("$\\beta_",0:2,"$"),
  paste0("$\\gamma_",1:2,"$")
  #paste0("$\\alpha$")
)

labels6 = c(
  paste0("$\\beta_",0:1,"$"),
  paste0("$\\gamma_",0:1,"$")
  #paste0("$\\alpha$")
)
par.real = NULL

MODEL = c(4, 5, 6)
for(M in MODEL){
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
  path.par = paste0("Data_simulation/real_par_model_",M,".txt")
  PAR = read.table(path.par)
  variables = PAR$par_names
  PAR$par_names = LAB
  len.betas = sum(stringr::str_count(LAB,"beta"))
  len.gamma = sum(stringr::str_count(LAB,"gamma"))
  trend = rep("Trend",len.betas)
  disp = rep("Dispersion",len.gamma)
  predictor = factor(c(trend,disp),levels = c("Trend","Dispersion"))
  PAR = data.frame(Model=M,predictor,PAR,variables)
  par.real = rbind(par.real,PAR)
  rm(PAR)
}
par.real <- par.real  %>% dplyr::relocate(Model,predictor,variables,par_names,real_values)
col.names = c("Model","Predictor","Variable","Parameter","Real values")
caption_model = "Real values to Parameters of models 4, 5 and 6"
label_model = "real_value_par"
table_model = kbl(
  par.real,
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
  centering = T
) %>%
  kable_paper(full_width = F) %>%
  column_spec(1:2, bold = T
              ) %>%
  #row_spec(c(1:5, 11:15)-1, extra_latex_after = "\\rowcolor{gray}") %>%
  collapse_rows(columns = 1:2, valign = "middle",latex_hline = 'major')
path_model = paste0("Data_simulation/Tables/real_par",
                    ".tex")
save_kable(table_model, path_model)

