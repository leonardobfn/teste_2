rm(list=ls())
require(tidyr)
require(dplyr)

dados1 = read.table("Data_simulation/Model_1/estimates/Method_2/estimates/estimates_alpha35_n144.txt")
dados1 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessian_Gamma_alpha95_n144.txt")
#dados3 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessianCovLambda_alpha95.txt")
#dados4 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessian_alpha95.txt")

nrow(dados1)
nrow(dados2)
#nrow(dados3)
#nrow(dados3)
#nrow(dados4)
conf = .05

id = dados1 %>% pull(V1) %>%sort()
length(table(id))
id.f = seq(1,1200) %in% id
which(id.f==F)
length(which(id.f==F) )
length(id)



est1 = dados1 %>% group_by(V2,V3,V6,V7,V8) %>% summarise(real_par = median(V4,na.rm=T)
                                                         ,m = median(V5,na.rm=T),
                                           l = quantile(V5,conf/2,na.rm=T),
                                           up = quantile(V5,1-(conf/2),na.rm=T))

est1
(dados1 %>% slice(which(is.na(V5)==T)) %>% nrow())/6
tt = dados1 %>% slice(which(V2>20)) %>% pull(V1)
est1 = dados1 %>% slice(tt) %>% group_by(V3) %>% summarise(m = median(V2,na.rm=T),
                                             l = quantile(V2,conf/2,na.rm=T),
                                             up = quantile(V2,1-(conf/2),na.rm=T))

est1
read.table("Data_simulation/real_par_model_1.txt")


files = list.files("Data_simulation/Model_1/estimates/Method_2/estimates/n168/alpha95/")
setwd("Data_simulation/Model_1/estimates/Method_2/estimates/n168/alpha95/")
estimates = NULL
for(FILES in files){
  estimates.aux = read.table(FILES)
  estimates = rbind(estimates,estimates.aux)
}

na = estimates$V5[is.na(estimates$V5)==T]
3000-(length(na)/6)
head(estimates)
estimates %>% group_by(V8,V2,V6,V7,V3) %>%
  summarise(Par_real = mean(V4),
            Mean = mean(V5,na.rm=T),
            Median = median(V5,na.rm=T))
