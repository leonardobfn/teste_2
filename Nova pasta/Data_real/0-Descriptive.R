
data(Manaus)
manaus = ts(Manaus$RH,frequency = 12,start = 2009)

sc = .80
width = 7 * sc
height = 7 * sc
path1 = paste0("Data_real/Figures/manaus", ".pdf")
path2 = paste0("Data_real/Figures/month", ".pdf")
path3 = paste0("Data_real/Figures/acf", ".pdf")
path4 = paste0("Data_real/Figures/pacf", ".pdf")
pdf(path1,
    width = width,
    height = height,
    title = "y")
plot.ts(manaus,ylab="Relative Humidity")
dev.off()

pdf(path2,
    width = width,
    height = height,
    title = "y")
#monthplot(manaus,ylab="Relative Humidity",xlab="Month")
#axis(1,labels = F,cex.axis=10)
boxplot(Manaus$RH~Manaus$month_names,ylab="Relative Humidity",xlab="Month")

dev.off()


pdf(path3,
    width = width,
    height = height,
    title = "y")
acf(manaus,main=" ")
#boxplot(Manaus$RH~Manaus$Month)
dev.off()

pdf(path4,
    width = width,
    height = height,
    title = "y")
pacf(manaus,main=" ",ylab="PACF")
#boxplot(Manaus$RH~Manaus$Month)
dev.off()


