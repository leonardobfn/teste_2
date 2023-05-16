devtools::load_all()
results = readRDS("Data_real/results_method_1_model_2.rds")
coef(results)
data  = results$data.n
formula = results$fomula
n = results$n
estimates = results$coefficients[, 1]
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 2)

ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)

beta <- estimates[c(1:ncx)]
lambda <- estimates[c((1 + ncx):(ncx + ncv))]
alpha <- estimates[-c((1):(ncx + ncv))]

nu = exp(cov_a %*% beta)
delta = exp(cov_delta %*% lambda)
b = 1 / delta

at = nu
bt = b

G.base.t = extraDistr::pkumar(y,
                              a = at,
                              b = bt ,
                              lower.tail = F)

LAMBDA.t = sort( -log( (c(G.base.t))  ))
F.barra = sort( exp( - (  (LAMBDA.t)^alpha  ) ))
LAMBDA.MAR = sort(-log(F.barra))

ReIns::ExpQQ(LAMBDA.t)
qqline(LAMBDA.t,distribution = function(p) qexp(p))

ReIns::WeibullQQ(LAMBDA.MAR)
qqline(LAMBDA.MAR,
       distribution = function(p)
         qweibull(p,shape = alpha,scale = 1))

ReIns::WeibullQQ(LAMBDA.t)
qqline(LAMBDA.t,
       distribution = function(p)
         qweibull(p,shape = alpha,scale = 1))





hist(LAMBDA.t ,prob=T)
curve(dexp(x,1),add=T)

hist(LAMBDA.MAR ,prob=T)
curve(dweibull(x,shape = alpha,scale = 1),add=T,col=2)


plot((1:132),LAMBDA.t)
plot(log(1:132),log(LAMBDA.MAR))


aa = qexp((gkumar_resid(results)))
bb = qweibull((gkumar_resid(results)),shape = alpha,scale = 1)
cc = qnorm((gkumar_resid(results)))


ReIns::ExpQQ(aa)
#abline(a=0,b=1)
qqline(aa,distribution = function(p) qexp(p))
#qqline
ReIns::WeibullQQ(bb)
qqline(bb,distribution = function(p) qweibull(p,shape = alpha,scale = 1))
#abline(a=0,b=1)
qqnorm(cc)
#abline(a=0,b=1)
qqline(cc)

