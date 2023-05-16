rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
devtools::load_all()
data(Manaus)
#st = read.table("dados.txt")/100
s = read.table("Data_simulation/Model_5/simulations/n144/alpha80/data100.txt")
data = Manaus
#data = "Data_simulation/"
data$RH = s$y#st[1:144,]
Models = readRDS("Data_real/Models.rds")
erro = 10 ^ (-4)
wd = getwd()
for(i in 1:length(Models)){
  #i=1
  formula = Models[[i]]
  results_1 = gkumar_method_1(formula,data,erro,n=132)
  results_2 = gkumar_method_2(formula,data,erro,n=132)
  path_1 = paste0(wd,"/Data_real/results_method_1_model_",i,".rds")
  path_2 = paste0(wd,"/Data_real/results_method_2_model_",i,".rds")
  saveRDS(results_1,path_1)
  saveRDS(results_2,path_2)
}

m12 =readRDS("Data_real/results_method_1_model_1.rds")
m22 =readRDS("Data_real/results_method_2_model_1.rds")

a = coef(m12)
b = coef(m22)
a
b

Method = rep( c("Methdo 1","Methodo 2"),each = nrow(coef(m12))  )
ab=rbind(a,b)
Variables = c(row.names(a),row.names(b))
Parameters = c(paste0("$\\beta_",0:4,"$") , paste0("$\\gamma_",0:1,"$"), paste0("$\\alpha","$"))

all = data.frame(Method,Variables,Parameters,ab)


require(kableExtra)
tab2 = all %>% kbl(
  escape = F,
  table.envir = T,
  align = "c",
  booktabs = T,
  format = "latex",
  row.names = F
) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "top",latex_hline = 'major')

save_kable(tab2,"Data_real/par_estimates.tex")

m12$AIC
m12$BIC

m22$AIC
m22$BIC
