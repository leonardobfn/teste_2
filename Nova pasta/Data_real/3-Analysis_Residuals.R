rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
devtools::load_all()
results = readRDS("Data_real/results_method_1_model_2.rds")
coef(results)
resid = (gkumar_resid(results))
# par(mfrow = c(2, 2))
# plot(resid,main=" ",ylab="Residuls",ylim=c(-3,3))
# abline(h=c(-2,2),lty = 2)
# acf(resid,main=" ")
# test = Box.test(resid, lag = 20, type = "Ljung")
# test
# #mtext(paste0("Box-Ljung test - p.value: ", round(test$p.value, 4)))
# pacf(resid,main=" ")
# qqnorm((resid),main=" ")
# qqline((resid),main=" ")
#hist(resid)

sc = .80
width = 7 * sc
height = 7 * sc
#par(mfrow = c(2, 2))
path1 = paste0("Data_real/Figures/diag_resid_1", ".pdf")
path2 = paste0("Data_real/Figures/diag_acf_1", ".pdf")
path3 = paste0("Data_real/Figures/diag_pacf_1", ".pdf")
path4 = paste0("Data_real/Figures/diag_qqnorm_1", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
plot(resid,main=" ",ylab="Residuls",ylim=c(-3,3))
abline(h=c(-2,2),lty = 2)
dev.off()

pdf(path2,
    width = width,
    height = height,
    title = "y")

acf(resid,main=" ")
dev.off()

# test = Box.test(resid, lag = 20, type = "Ljung")
# test
#mtext(paste0("Box-Ljung test - p.value: ", round(test$p.value, 4)))
pdf(path3,
    width = width,
    height = height,
    title = "y")
pacf(resid,main=" ")
dev.off()

pdf(path4,
    width = width,
    height = height,
    title = "y")
qqnorm((resid),main=" ")
qqline((resid),main=" ")
dev.off()


#model 2

rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
devtools::load_all()
results = readRDS("Data_real/results_method_2_model_2.rds")
coef(results)
resid = (gkumar_resid(results))
# resid = qnorm(gkumar_resid(results))
par(mfrow = c(2, 2))
plot(resid,main=" ",ylab="Residuls",ylim=c(-3,3))
abline(h=c(-2,2),lty = 2)
acf(resid,main=" ")
test = Box.test(resid, lag = 20, type = "Ljung-Box")
test
# #mtext(paste0("Box-Ljung test - p.value: ", round(test$p.value, 4)))
# pacf(resid,main=" ")
# qqnorm((resid),main=" ")
# qqline((resid),main=" ")
#hist(resid)

sc = .80
width = 7 * sc
height = 7 * sc
#par(mfrow = c(2, 2))
path1 = paste0("Data_real/Figures/diag_resid_2", ".pdf")
path2 = paste0("Data_real/Figures/diag_acf_2", ".pdf")
path3 = paste0("Data_real/Figures/diag_pacf_2", ".pdf")
path4 = paste0("Data_real/Figures/diag_qqnorm_2", ".pdf")

pdf(path1,
    width = width,
    height = height,
    title = "y")
plot(resid,main=" ",ylab="Residuls",ylim=c(-3,3))
abline(h=c(-2,2),lty = 2)
dev.off()

pdf(path2,
    width = width,
    height = height,
    title = "y")

acf(resid,main=" ")
dev.off()

# test = Box.test(resid, lag = 20, type = "Ljung")
# test
#mtext(paste0("Box-Ljung test - p.value: ", round(test$p.value, 4)))
pdf(path3,
    width = width,
    height = height,
    title = "y")
pacf(resid,main=" ")
dev.off()

pdf(path4,
    width = width,
    height = height,
    title = "y")
qqnorm((resid),main=" ")
qqline((resid),main=" ")
dev.off()


