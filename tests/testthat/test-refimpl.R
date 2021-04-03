context('Reference implementation comparison')

nlsy97_1983_cohort <- function(omit.na=T) {
  data(nlsy97depression)
  df <- nlsy97depression[nlsy97depression$birth_y==1983,c('pid', 'occasion', 'nervous', 'calm', 'down', 'happy', 'depressed')]
  if (omit.na) {
    df_wide <- na.omit(reshape(df, timevar='occasion', idvar='pid', direction='wide',
                               v.names=c('nervous', 'calm', 'down', 'happy', 'depressed')))
    expect_equal(nrow(df_wide), 1397)
    reshape(df_wide)
  } else {
    df
  }
}

mxmmod_ref <- function(df, do_fiml=F) {
  require(OpenMx)
  require(mxmmod)
  structure <- list(
    F1 = c('nervous', 'down', 'depressed', 'calm', 'happy')
  )
  mmod_model <- mxMmodModel(data=df,
                            modelName='1 Factor MMOD',
                            idvar='pid', timevar='occasion', structure=structure, fiml=do_fiml)
  summary(mxRun(mmod_model))
}

estabrook_ref <- function(df) {
  # load required libraries
  require(OpenMx)

  # load raw and covariance data
  rawData <- reshape(df, timevar='occasion', idvar='pid', direction='wide', sep='0')[-1]
  covData <- cov(rawData, use='complete.obs')

  # make weight matrix with Deboeck ’s functions
  #source("GOLD.R")
  #weight <- t (solve (Contrasts (1:3, 2)))

  # make weight matrix without Deboeck ’s functions
  weight <- matrix(
    c(1/3, 1/3, 1/3,
      -0.5, 0, 0.5,
      1/2, -1, 1/2),
    nrow=3, ncol=3)

  # expand weight matrix for full data
  #fullWeight <- kronecker(weight, diag(5))

  # transform the data
  #rawDeriv <- as.matrix(rawData) %*% fullWeight
  #covDeriv <- t(fullWeight) %*% covData %*% fullWeight
  #dimnames(covDeriv) <- dimnames(covData)

  # one factor MMOD model, untransformed data
  obsVar <- colnames(covData)
  latDeriv <- c("nerv0", "down0", "dep0", "calm0", "happy0",
                 "nerv1", "down1", "dep1", "calm1", "happy1",
                 "nerv2", "down2", "dep2", "calm2", "happy2")
  longWeight <- as.vector(t(solve(weight)))

  mmod2 <- mxModel("One Factor MMOD 2",
                    mxData(covData, type="cov", numObs=1397),
                    type="RAM",
                    manifestVars=colnames(covData),
                    latentVars=c("F0", "F1", "F2", latDeriv),
                    # transformation
                    mxPath(from=c("nerv0", "nerv1", "nerv2"),
                           to=c("nervous00", "nervous01", "nervous02"),
                           connect="all.pairs",
                           free=FALSE, values=longWeight),
                    mxPath(from=c("down0", "down1", "down2"),
                           to=c("down00", "down01", "down02"),
                           connect="all.pairs",
                           free=FALSE, values=longWeight),
                    mxPath(from=c("dep0", "dep1", "dep2"),
                           to=c("depressed00", "depressed01", "depressed02"),
                           connect="all.pairs",
                           free=FALSE, values=longWeight),
                    mxPath(from=c("calm0", "calm1", "calm2"),
                           to=c("calm00", "calm01", "calm02"),
                           connect="all.pairs",
                           free=FALSE, values=longWeight),
                    mxPath(from=c("happy0", "happy1", "happy2"),
                           to=c("happy00", "happy01", "happy02"),
                           connect="all.pairs",
                           free=FALSE, values=longWeight),
                    # factor loadings
                    mxPath(from="F0",
                           to=c("nerv0", "down0", "dep0", "calm0", "happy0"),
                           values=0.5, free=TRUE),
                    mxPath(from="F1",
                           to=c("nerv1", "down1", "dep1", "calm1", "happy1"),
                           values =0.5, free=TRUE),
                    mxPath(from="F2",
                           to=c("nerv2", "down2", "dep2", "calm2", "happy2"),
                           values =0.5, free=TRUE),
                    # factor variances
                    mxPath(from=c("F0", "F1", "F2"),
                           arrows=2, values=1, free=FALSE),
                    # factor correlations
                    mxPath(from=c("F0", "F1", "F2"),
                           arrows=2, connect="unique.bivariate", free=TRUE),
                    # residual variances(only for latent derivatives !)
                    mxPath(from=latDeriv, arrows=2, values=1)
  )
  summary(mxRun(mmod2))
}

test_that('Same results as Estabrook 2015', {
  df <- nlsy97_1983_cohort()
  a <- mxmmod_ref(df)
  b <- estabrook_ref(df)
  expect_equal(a$parameters$Estimate, b$parameters$Estimate, tolerance = .001)
  expect_equal(a$parameters$Std.Error, b$parameters$Std.Error, tolerance = .001)
})

test_that('Floating point occasions', {
  df <- nlsy97_1983_cohort()

  scale_val <- 2
  scaled_nlsy <- df
  scaled_nlsy$occasion <- scaled_nlsy$occasion / scale_val

  trans <- rep(1, 33)
  trans[6:10] <- scale_val
  trans[11:15] <- scale_val^2
  trans[24:33] <- scale_val^2
  trans[29:33] <- scale_val^4

  a <- mxmmod_ref(scaled_nlsy)
  b <- estabrook_ref(df)

  a$parameters$Estimate <- a$parameters$Estimate / trans
  a$parameters$Std.Error <- a$parameters$Std.Error / trans

  expect_equal(a$parameters$Estimate, b$parameters$Estimate, tolerance = .001)
  expect_equal(a$parameters$Std.Error, b$parameters$Std.Error, tolerance = .001)
})

test_that('FIML', {
  df <- na.omit(nlsy97_1983_cohort(omit.na=F))
  a <- mxmmod_ref(df, do_fiml=T)
  b <- estabrook_ref(df)
  expect_equal(a$parameters$Estimate[1:33], b$parameters$Estimate, tolerance = .05)
  expect_equal(a$parameters$Std.Error[1:33], b$parameters$Std.Error, tolerance = .01)
})

test_that('Omit missing values warning', {
  df <- na.omit(nlsy97_1983_cohort(omit.na=F))
  expect_warning(mxmmod_ref(df, do_fiml=F), regexp="Missing values")
})

test_that('Argument check timevar numeric', {
  df <- nlsy97_1983_cohort()
  df$occasion <- as.factor(df$occasion)
  expect_error(mxmmod_ref(df), 'timevar must be a numeric type')
})

test_that('Argument check timevar', {
  df <- nlsy97_1983_cohort()
  df$occasion <- NULL
  expect_error(mxmmod_ref(df), "Column 'occasion' not found")
})

test_that('Argument check idvar', {
  df <- nlsy97_1983_cohort()
  df$pid <- NULL
  expect_error(mxmmod_ref(df), "Column 'pid' not found")
})
