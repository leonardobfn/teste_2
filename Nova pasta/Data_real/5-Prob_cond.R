rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
require(latex2exp)
devtools::load_all()
results = readRDS("Data_real/results_method_2_model_2.rds")
coef(results)
y = results$data.n$RH
yy = ts(y, frequency = 12, start = 2009)
n = results$n
prob_1 = gkumar_prob_cond_1(results, yt2 = 0.7)
prob_2 = prob_1
YLAB= TeX(sprintf(
  "$P( Y_t>0.70|Y_{t-1}=y_{t-1} )$"
))
sc = .80
width = 7 * sc
height = 7 * sc


month.complet = results$data.n$month_names

probs2 = data.frame(
  p = prob_2[c(-1)],
  t = paste0("t = ", month.complet[-c(1)]),
  t_1 = paste0("t-1 = ", month.complet[-c(132)])
#  t_2 = paste0("t-2 = ", month.complet[-c(131, 132)])
)
a = ts(c(NA,probs2$p),frequency = 12,start = 2009)

path1 = paste0("Data_real/Figures/p_ordem11", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
boxplot(a~results$data.n$month_names,xlab="Month (t) ",ylab=YLAB)
abline(h=0.50)
#abline(h=c(-2,2),lty = 2)
dev.off()


l = paste0("t = ", month.abb)
l1 = paste0("t-1 = ", month.abb)

lt = factor(paste0("t = ", month.complet),
            levels = l)

lt1 = factor(paste0("t-1 = ", month.complet),
             levels = l1)


probs2 = probs2 %>% mutate(t = lt[-c(1)],
                           t_1 = lt1[-c(132)])

teste = probs2 %>% filter(t == "Dec")
th = theme(
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 7),
  strip.text = element_text(size = 5)
)
r = function(y) {
  y = ifelse(y < 0.0001, "<0.0001", y)
  return(y)
}
plot. = probs2 %>% arrange(t) %>%
  ggplot(aes()) +
  geom_boxplot(aes(t_1, p), outlier.size = 0.2, varwidth = T) +
  facet_wrap(t ~ ., scales = "free_x", dir = "h") +
  scale_y_continuous(labels = r) +
  theme_bw() + th +
  ylab(TeX(sprintf(
    "$P( Y_t>0.70|Y_{t-1}=y_{t-1} )$"
  ))) +
  geom_hline(yintercept = 0.50)

path1 = paste0("Data_real/Figures/p_ordem1", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
plot.
#abline(h=c(-2,2),lty = 2)
dev.off()



prob_2 = NULL
for (t in 3:132) {
  t.vector = c(t - 1, t)
  prob_2[t] = gkumar_prob_cond_choose_2(
    results,
    h.before = 2,
    t.choose = t.vector,
    y.choose = c(0.70, 0.70)
  )$prob[t]
}
month.complet = results$data.n$month_names
b = ts(c(prob_2),start = 2009,frequency = 12)
YLAB = TeX(sprintf(
  "$P( Y_t>0.70,Y_{t-1}>0.70|Y_{t-2}=y_{t-2} )$"
))
path1 = paste0("Data_real/Figures/p_ordem22", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
boxplot(b~results$data.n$month_names,xlab="Month",ylab=YLAB)
abline(h=c(0.50),lty = 2)
dev.off()


probs2 = data.frame(
  p = prob_2[-c(1, 2)],
  t = paste0("t = ", month.complet[-c(1, 2)]),
  t_1 = paste0("t-1 = ", month.complet[-c(1, 132)]),
  t_2 = paste0("t-2 = ", month.complet[-c(131, 132)])
)

l = paste0("t = ", month.abb)
l1 = paste0("t-1 = ", month.abb)
l2 = paste0("t-2 = ", month.abb)

lt = factor(paste0("t = ", month.complet),
            levels = l)

lt1 = factor(paste0("t-1 = ", month.complet),
             levels = l1)
lt2 = factor(paste0("t-2 = ", month.complet),
             levels = l2)


probs2 = probs2 %>% mutate(t = lt[-c(1, 2)],
                           t_1 = lt1[-c(1, 132)],
                           t_2 = lt2[-c(131, 132)])

teste = probs2 %>% filter(t == "Dec")
th = theme(
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 7),
  strip.text = element_text(size = 5)
)
r = function(y) {
  y = ifelse(y < 0.0001, "<0.0001", y)
  return(y)
}
plot. = probs2 %>% arrange(t) %>%
  ggplot(aes()) +
  geom_boxplot(aes(t_2, p), outlier.size = 0.2, varwidth = T) +
  facet_wrap(t ~ t_1, scales = "free_x", dir = "h") +
  scale_y_continuous(labels = r) +
  theme_bw() + th +
  ylab(TeX(sprintf(
    "$P( Y_t>0.70,Y_{t-1}>0.70|Y_{t-2}=y_{t-2} )$"
  ))) +
  geom_hline(yintercept = 0.50)

path1 = paste0("Data_real/Figures/p_ordem2", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
plot.
#abline(h=c(-2,2),lty = 2)
dev.off()















prob_forecast_1 = gkumar_prob_cond_model_forecast_2(
  results,
  h.forecast = 12,
  h.before = 1,
  MC = T,
  SEED = 1234
)
prob_forecast_2 = gkumar_prob_cond_model_forecast_2(
  results,
  h.forecast = 12,
  h.before = 2,
  MC = T,
  SEED = 1234
)

res = round(
  data.frame(
    t = 133:144,
    y.real = results$data$RH[133:144],
    y.forecast = prob_forecast_2$y[133:144],
    P.1 = prob_forecast_1$prob[133:144],
    P.2 = prob_forecast_2$prob[133:144]
  ),
  3
)

mean((res$y.real - res$y.forecast) ^ 2)
mean((res$y.real - res$y.forecast) / res$y.real)


probs = matrix(0, 132, 6)
for (k in 1:6) {
  probs[, k] = gkumar_prob_cond_choose_2(results, h.before = k)$prob
}
dados = data.frame(
  ano = results$data.n$Year,
  mes = results$data.n$month_names,
  temp = results$data.n$TBS,
  y,
  probs
)
j = dados %>% filter(mes == "Jan")
boxplot(dados$X1 ~ dados$mes)


a = gkumar_forecast_mean(results, 12, MC = T)
plot.ts(a$y.forecast, ylim = c(0.40, 0.80))
L1 = a$y.forecast - 1.96 * a$SD
L2 = a$y.forecast + 1.96 * a$SD
lines(L1, col = 2)
lines(L2, col = 2)
lines(results$data$RH[133:144], col = 3)
