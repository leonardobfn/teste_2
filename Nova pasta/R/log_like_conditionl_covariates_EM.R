log.like_conditionl_covariates_EM = function(theta, y, x, w, Etil1,ncx,ncv) {
  # x: betas covariates
  # Etil1 = Esp.z
  # x = cov_a
  # w = cov_delta
  # theta = par.covariates.start

  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)] # real
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% Beta
  a = exp(xbeta)
  Gb = pkumar(y,
              a,
              b = 1 / delta,
              lower.tail = FALSE,
              log.p = FALSE)
  d = dkumar(y, a, b = 1 / delta, log = FALSE)
  r  = d / Gb
  l = Etil1 * log(Gb)
  ll = l+log(r)
  #ll = l - log(delta)+log(a) + a*log(y)-log(1-y^a) - log(y)
  #ll = l - log(delta)+log(a) - a*log(y)-log(1-y^a) - log(y)# como está no texto
  ll = ll[is.finite(ll)]
  L = sum(ll)
  #L = sum(l + log(r.), na.rm = T)

  return(L)
}
