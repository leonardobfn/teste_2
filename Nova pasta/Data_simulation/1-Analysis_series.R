require(ggplot2)
require(tidyverse)
require(dplyr)
alpha.values = c("alpha35", "alpha50", "alpha65", "alpha80" , "alpha95")
#alpha.values = alpha.values[1]
MODEL = c(4, 5, 6)
id=300
#MODEL = MODEL[1]
N = c(168)
sc = .80
width = 7 * sc
height = 7 * sc
y = NULL

for (M in MODEL) {
  for (n in N) {
    for (ALPHA in alpha.values) {
      path = paste0(wd,
                    "/Data_simulation/Model_",
                    M,
                    "/simulations/n",
                    n,
                    "/",
                    ALPHA,
                    "/data",id,".txt")
      y.aux = read.table(path)
      y = rbind(y, y.aux)
    }
  }
}
nrow(y)

#------------ séries------

fig_series = y %>%  mutate(t = rep(1:N, 15))  %>%
    ggplot(.) +
      geom_line(aes(t, y),linewidth=0.2) +
      facet_wrap(
        alpha ~ Model,
        scales = "free_y",
        labeller = label_bquote(alpha == .(alpha)),ncol = 3
      ) +
      theme_bw() +
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=5),
        axis.title.y = element_text(size=5),
        strip.text = element_text(size=5))


# fig_series = y %>%  mutate(t = rep(1:144, 30))  %>% dplyr::group_split(Model,alpha) %>%
#   purrr::map(
#     ~ggplot(.) +
#       geom_line(aes(t, y)) +
#       facet_grid(
#         alpha ~ Model,
#         scales = "free_y",
#         labeller = label_bquote(alpha == .(alpha))
#       ) +
#       theme_bw()
#   ) %>%
#   cowplot::plot_grid(plotlist = .,nrow=6)

#filter(Model == "Model 6") %>%

#fig_series

path1 = paste0("Data_simulation/Figures/series", ".pdf")
#path1 = paste0("Data_simulation/Figures/series",".eps")
#path2 = paste0()
pdf(path1,
    width = width,
    height = height,
    title = "y")
#cairo_ps(path1)
print(fig_series)
dev.off()
#------------ ACF series------

ic_alpha = function(alpha, acf_res) {
  return(qnorm((1 + (1 - alpha)) / 2) / sqrt(acf_res$n.used))
}

acff = y %>% group_by(Model, V3) %>%
  reframe(
    acff = acf(V1, plot = F)$acf[, , 1],
    lag = acf(V1, plot = F)$lag,
    lim1 = ic_alpha(0.05, acf(V1))
  )

fig_acf = acff %>%
  ggplot() +
  geom_segment(mapping = aes(
    x = lag,
    y = 0,
    xend = lag,
    yend = acff
  )) +
  geom_hline(aes(yintercept = 0)) +
  geom_hline(aes(yintercept = lim1), linetype = 2, color = 'blue') +
  geom_hline(aes(yintercept = -lim1), linetype = 2, color = 'blue') +
  facet_grid(V3 ~ Model  , labeller = label_bquote(alpha == .(V3))) +
  theme_bw()
#path1 = paste0("Data_simulation/Figures/series",".eps")
path1 = paste0("Data_simulation/Figures/acf", ".pdf")
#path2 = paste0()
pdf(path1,
    width = width,
    height = height,
    title = "acf")
#cairo_ps(path1)
print(fig_acf)
dev.off()
