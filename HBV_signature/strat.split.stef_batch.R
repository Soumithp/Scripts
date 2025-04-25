source("strat.split.stef.R")

strat.split.stef(
  "bMilan_tumor_percentP40_CubicSpline_symbolmed.gct",
  "bMilan_outcome_sampleOrderAsGCT.txt",

  output.name="Split",

  kfold=2,
  cls.var = "rec",
  strata.var="time.to.rec.days",
  strata.param=10,
  method=1,

  num.split=5,

  rnd.seed=2006
  )
