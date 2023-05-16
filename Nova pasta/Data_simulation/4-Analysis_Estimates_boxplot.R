


# rm(list = ls())
#
# interval = function(x, conf, ROUND) {
#   # conf = 0.95
#   # ROUND = 3
#   # x = rnorm(100)
#   al = 1 - conf
#   int = quantile(x, c(al / 2, 1 - (al / 2)))
#   int2 = paste0("(",
#                 round(int[1],ROUND),
#                 ",",
#                 round(int[2],ROUND),
#                 ")")
#   return(int2)
# }
#
# #packages--------------
# require(tidyr)
# require(dplyr)
# require(extraDistr)
# wd = getwd()
# #Model 1--------
# MODEL = c(1, 2,3,4,5,6)
# for (M in MODEL) {
#
#   file.path = paste0(wd,
#                      "/Data_simulation/MC_estimates/mc_estimates_model_",
#                      M,
#                      ".rds")
#   mc_estimates_model_1 <- readRDS(file.path)
#   file.path.par = paste0(wd, "/Data_simulation/real_par_model_", M, ".txt")
#   emv = read.table(file.path.par)
#
#   mc_estimates_model_1$V3  = factor(mc_estimates_model_1$V3,
#                                     levels = c(emv$par_names, "alpha"))
#
#
#   write.table(
#     mc_estimates_model_1,
#     "Data_simulation/MC_estimates/mc_estimates.txt",
#     append = T,
#     col.names = F,
#     row.names = F,
#   )
#
# }

#colnames(mc) = c("mc","alpha","Parameter","Par_real","Estimates","Method","n","Model","start")
#saveRDS(mc,"Data_simulation/MC_estimates/mc_estimates.rds")


rm(list = ls())
#packages--------------
require(ggplot2)
require(ggridges)
require(tidyr)
require(dplyr)
require(extraDistr)
devtools::load_all()
compiler::enableJIT(3)
alphas.plot = c("0.35","0.5" ,"0.65","0.8" ,"0.95")
sc = .80
width = 12 * sc
height = 7 * sc
mc = readRDS("Data_simulation/MC_estimates/mc_estimates.rds")
head(mc)
# width=NULL
# height=NULL
labels1_3 = c(
  expression(beta[0]),
  expression(beta[1]),
  expression(beta[2]),
  expression(gamma[0]),
  expression(gamma[1]),
  expression(alpha)
)
labels2_4 = c(
  expression(beta[1]),
  expression(beta[2]),
  expression(gamma[0]),
  expression(gamma[1]),
  expression(alpha)
)

labels5 = c(
  expression(beta[0]),
  expression(beta[1]),
  expression(beta[2]),
  expression(gamma[1]),
  expression(gamma[2]),
  expression(alpha)
)

labels6 = c(
  expression(beta[0]),
  expression(beta[1]),
  expression(gamma[0]),
  expression(gamma[1]),
  expression(alpha)
)
models = c(1,2,3,4,5,6)
wd = getwd()
theme.all = theme_bw() + theme(axis.text.x = element_text(size=5,angle = 90),
                               axis.text.y = element_text(size=5),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               strip.text = element_text(size=6))
# Model 1------------
for (M in models) {
  #M = 2
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

  file.path.par = paste0(wd, "/Data_simulation/real_par_model_", M, ".txt")
  emv = read.table(file.path.par)

  estimates_model_1 = mc %>% filter(Model == paste0("Model ", M))

  estimates_model_1$Parameter = factor(
    estimates_model_1$Parameter,
    levels = c(emv$par_names, "alpha"),
    labels = LAB
  )


  for (i in 1:length(alphas.plot)) {
    #i = 1
    fig_box_plot = estimates_model_1 %>% mutate(alpha = factor(alpha)) %>%
      filter(alpha == alphas.plot[i]) %>%
      dplyr::group_split(n, Parameter) %>% purrr::map(
        ~ ggplot(., aes(Method, Estimates)) +
          #geom_boxplot() +
          stat_summary(fun.data = calc_boxplot_stat , geom = "boxplot") +
          facet_grid(
            Parameter ~ n,
            shrink = T,
            scales = "free_y",
            labeller = label_parsed,
            space = "free"
          ) +
          geom_hline(data =  .,
                     aes(yintercept = Par_real)) +
          # xlab("Method")+
          # ylab("Estimates")+
          theme.all
      ) %>% cowplot::plot_grid(plotlist = ., nrow = 3)

    path1 = paste0("Data_simulation/Figures/Model_",
                   M,
                   "_",
                   alphas.plot[i],
                   "_box_plot.pdf")
    #path2 = paste0()
    pdf(path1, width = width, height = height)
    #postscript(path)
    print(fig_box_plot)
    dev.off()
  }

}
