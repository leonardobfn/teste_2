rm(list=ls())
devtools::load_all()
results2 = readRDS("Data_real/results_method_2_model_2.rds")
gkumar_forecast_mean(results2,h=12,MC = F)
coef(results2)
  forecast_esp = gkumar_pred_quantile_forecast(results2,h=12,quant = 0.5)
  forecast_risco = gkumar_forecast_mean(results2,h=12,MC =F)
  y.real =  results2$data$RH[133:144]
  # MAPE_1 = mean((forecast1$y.forecast - y.real)/y.real)
  # MSE_1 = mean((forecast1$y.forecast - y.real))^2


  MAPE_2_esp = mean((forecast_esp$y.forecast - y.real)/y.real)
  MAPE_2_risco = mean((forecast_risco$y.forecast - y.real)/y.real)

  MSE_2_esp = mean((forecast_esp$y.forecast - y.real)^2)
  MSE_2_risco = mean((forecast_risco$y.forecast - y.real)^2)

  stat = data.frame(
                    Forecast = c("Median","Risco"),
                    MAPE = c(MAPE_2_esp,MAPE_2_risco),
                    MSE = c(MSE_2_esp,MSE_2_risco))

  require(kableExtra)
  kbl(stat,format = "latex",table.envir = "table",booktabs = T)

  sc = .80
  width = 7 * sc
  sc = .80
  width = 7 * sc
  height = 7 * sc
pred = gkumar_pred_quantile(results2,quant = 0.5)
pred2 = gkumar_pred_mean(results2,MC=F)
yy = pred[2:132]
y2 = pred2$y.pred
y = results2$data.n$RH[2:132]
path1 = paste0("Data_real/Figures/fitted", ".pdf")
pdf(path1,
    width = width,
    height = height,
    title = "y")
plot.ts(y,ylab="Relative Humidity")
lines(yy,col=1,lty=2,lwd=1.5)
lines(y2,col=1,lty=3,lwd=1.5)
legend("bottomleft", legend=c("Real", "Median","RR"),
       lty = c(1,2,3), cex=0.8)

dev.off()

