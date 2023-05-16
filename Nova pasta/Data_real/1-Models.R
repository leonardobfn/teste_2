Models <- list(
  RH ~ sent + cost  | semester,
  #RH ~ sent + cost  |   semester-1
  RH ~ log(TBS) + log(t) + sent+cost | semester
  #okRH ~ log(TBS) + log(t) + sent | semester
  #RH ~   t+TBS + sent + cost  |  semester
  #RH ~ log(t)+TBS + sent + cost |log(t)
  #RH ~ I(log(t)) + TBS + sent + cost|semester
  #RH ~ t + log(TBS) + sent + cost|semester
  #RH ~ t + I(t^2)+TBS |semester
  #RH ~TBS+log(precp)|log(t)+log(precp)
  #RH ~TBS|log(t)

)
saveRDS(Models,"Data_real/Models.rds")
