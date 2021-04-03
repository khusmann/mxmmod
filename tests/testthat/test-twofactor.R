context('Two factor tests')

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

test_that("basic two factor model works", {
  require(OpenMx)
  require(mxmmod)

  df <- nlsy97_1983_cohort()

  structure <- list(
    F1 = c('nervous', 'down', 'depressed'),
    F2 = c('calm', 'happy')
  )
  mmod_model <- mxMmodModel(data=df,
                            modelName='2 Factor MMOD',
                            idvar='pid', timevar='occasion', structure=structure, fiml=F)
  summary(mxRun(mmod_model))

  dval1 <- c(1, 1, 1)
  dval2 <- c(-1, 0, 1)
  dval3 <- c(1/3, -2/3, 1/3)

  mmod_model_ref <- mxModel(
    type='RAM',
    manifestVars=c('nervous_1', 'nervous_2', 'nervous_3',
                   'down_1', 'down_2', 'down_3',
                   'depressed_1', 'depressed_2', 'depressed_3',
                   'calm_1', 'calm_2', 'calm_3',
                   'happy_1', 'happy_2', 'happy_3'
    ),
    latentVars=  c('F1_1', 'F2_1', 'F1_2', 'F2_2', 'F1_3', 'F2_3',
                   'dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                   'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                   'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ),

    # Factor structures
    mxPath(from='F1_1', to=c('dnervous_1', 'ddown_1', 'ddepressed_1'), values=0.5),
    mxPath(from='F1_2', to=c('dnervous_2', 'ddown_2', 'ddepressed_2'), values=0.5),
    mxPath(from='F1_3', to=c('dnervous_3', 'ddown_3', 'ddepressed_3'), values=0.5),

    mxPath(from='F2_1', to=c('dcalm_1', 'dhappy_1'), values=0.5),
    mxPath(from='F2_2', to=c('dcalm_2', 'dhappy_2'), values=0.5),
    mxPath(from='F2_3', to=c('dcalm_3', 'dhappy_3'), values=0.5),

    # Factor Covariances
    mxPath(from=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'),
           to=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'), 'unique.bivariate', arrows=2),

    # Factor variances fixed @ 1
    mxPath(from=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'),
           to=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'), 'single', arrows=2, values=1, free=F),

    # residual variances (for latent derivs)
    mxPath(from=c('dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                  'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                  'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ), arrows=2, values=1),

    # Build latent derivs
    mxPath(from='dnervous_1', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval1),
    mxPath(from='dnervous_2', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval2),
    mxPath(from='dnervous_3', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval3),

    mxPath(from='ddown_1', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval1),
    mxPath(from='ddown_2', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval2),
    mxPath(from='ddown_3', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval3),

    mxPath(from='ddepressed_1', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval1),
    mxPath(from='ddepressed_2', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval2),
    mxPath(from='ddepressed_3', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval3),

    mxPath(from='dcalm_1', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval1),
    mxPath(from='dcalm_2', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval2),
    mxPath(from='dcalm_3', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval3),

    mxPath(from='dhappy_1', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval1),
    mxPath(from='dhappy_2', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval2),
    mxPath(from='dhappy_3', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval3)
  )

  expect_equal(mmod_model$A, mmod_model_ref$A)
  expect_equal(mmod_model$S, mmod_model_ref$S)
  expect_equal(mmod_model$F, mmod_model_ref$F)
})

test_that("factors can have crossloadings", {
  require(OpenMx)
  require(mxmmod)

  df <- nlsy97_1983_cohort()

  structure <- list(
    F1 = c('nervous', 'down', 'depressed'),
    F2 = c('calm', 'happy', 'nervous')
  )
  mmod_model <- mxMmodModel(data=df,
                            modelName='2 Factor MMOD',
                            idvar='pid', timevar='occasion', structure=structure, fiml=F)
  summary(mxRun(mmod_model))

  dval1 <- c(1, 1, 1)
  dval2 <- c(-1, 0, 1)
  dval3 <- c(1/3, -2/3, 1/3)

  mmod_model_ref <- mxModel(
    type='RAM',
    manifestVars=c('nervous_1', 'nervous_2', 'nervous_3',
                   'down_1', 'down_2', 'down_3',
                   'depressed_1', 'depressed_2', 'depressed_3',
                   'calm_1', 'calm_2', 'calm_3',
                   'happy_1', 'happy_2', 'happy_3'
    ),
    latentVars=  c('F1_1', 'F2_1', 'F1_2', 'F2_2', 'F1_3', 'F2_3',
                   'dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                   'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                   'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ),

    # Factor structures
    mxPath(from='F1_1', to=c('dnervous_1', 'ddown_1', 'ddepressed_1'), values=0.5),
    mxPath(from='F1_2', to=c('dnervous_2', 'ddown_2', 'ddepressed_2'), values=0.5),
    mxPath(from='F1_3', to=c('dnervous_3', 'ddown_3', 'ddepressed_3'), values=0.5),

    mxPath(from='F2_1', to=c('dcalm_1', 'dhappy_1', 'dnervous_1'), values=0.5), # (crossloading!)
    mxPath(from='F2_2', to=c('dcalm_2', 'dhappy_2', 'dnervous_2'), values=0.5), # (crossloading!)
    mxPath(from='F2_3', to=c('dcalm_3', 'dhappy_3', 'dnervous_3'), values=0.5), # (crossloading!)

    # Factor Covariances
    mxPath(from=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'),
           to=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'), 'unique.bivariate', arrows=2),

    # Factor variances fixed @ 1
    mxPath(from=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'),
           to=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'), 'single', arrows=2, values=1, free=F),

    # residual variances (for latent derivs)
    mxPath(from=c('dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                  'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                  'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ), arrows=2, values=1),

    # Build latent derivs
    mxPath(from='dnervous_1', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval1),
    mxPath(from='dnervous_2', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval2),
    mxPath(from='dnervous_3', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval3),

    mxPath(from='ddown_1', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval1),
    mxPath(from='ddown_2', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval2),
    mxPath(from='ddown_3', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval3),

    mxPath(from='ddepressed_1', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval1),
    mxPath(from='ddepressed_2', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval2),
    mxPath(from='ddepressed_3', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval3),

    mxPath(from='dcalm_1', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval1),
    mxPath(from='dcalm_2', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval2),
    mxPath(from='dcalm_3', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval3),

    mxPath(from='dhappy_1', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval1),
    mxPath(from='dhappy_2', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval2),
    mxPath(from='dhappy_3', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval3)
  )

  expect_equal(mmod_model$A, mmod_model_ref$A)
  expect_equal(mmod_model$S, mmod_model_ref$S)
  expect_equal(mmod_model$F, mmod_model_ref$F)
})

test_that("orthogonal factors", {
  require(OpenMx)
  require(mxmmod)

  df <- nlsy97_1983_cohort()

  structure <- list(
    F1 = c('nervous', 'down', 'depressed'),
    F2 = c('calm', 'happy', 'nervous')
  )
  mmod_model <- mxMmodModel(data=df,
                            modelName='2 Factor orthogonal MMOD',
                            idvar='pid', timevar='occasion', structure=structure, fiml=F,
                            orthogonal=T)
  summary(mxRun(mmod_model))

  dval1 <- c(1, 1, 1)
  dval2 <- c(-1, 0, 1)
  dval3 <- c(1/3, -2/3, 1/3)

  mmod_model_ref <- mxModel(
    type='RAM',
    manifestVars=c('nervous_1', 'nervous_2', 'nervous_3',
                   'down_1', 'down_2', 'down_3',
                   'depressed_1', 'depressed_2', 'depressed_3',
                   'calm_1', 'calm_2', 'calm_3',
                   'happy_1', 'happy_2', 'happy_3'
    ),
    latentVars=  c('F1_1', 'F2_1', 'F1_2', 'F2_2', 'F1_3', 'F2_3',
                   'dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                   'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                   'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ),

    # Factor structures
    mxPath(from='F1_1', to=c('dnervous_1', 'ddown_1', 'ddepressed_1'), values=0.5),
    mxPath(from='F1_2', to=c('dnervous_2', 'ddown_2', 'ddepressed_2'), values=0.5),
    mxPath(from='F1_3', to=c('dnervous_3', 'ddown_3', 'ddepressed_3'), values=0.5),

    mxPath(from='F2_1', to=c('dcalm_1', 'dhappy_1', 'dnervous_1'), values=0.5), # (crossloading!)
    mxPath(from='F2_2', to=c('dcalm_2', 'dhappy_2', 'dnervous_2'), values=0.5), # (crossloading!)
    mxPath(from='F2_3', to=c('dcalm_3', 'dhappy_3', 'dnervous_3'), values=0.5), # (crossloading!)

    # Factor Covariances (orthogonal => within cov, NO BETWEEN COV)
    mxPath(from=c('F1_1', 'F1_2', 'F1_3'),
           to=c('F1_1', 'F1_2', 'F1_3'), 'unique.bivariate', arrows=2),
    mxPath(from=c('F2_1', 'F2_2', 'F2_3'),
           to=c('F2_1', 'F2_2', 'F2_3'), 'unique.bivariate', arrows=2),

    # Factor variances fixed @ 1
    mxPath(from=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'),
           to=c('F1_1', 'F1_2', 'F1_3', 'F2_1', 'F2_2', 'F2_3'), 'single', arrows=2, values=1, free=F),

    # residual variances (for latent derivs)
    mxPath(from=c('dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                  'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2',
                  'dnervous_3', 'ddown_3', 'ddepressed_3', 'dcalm_3', 'dhappy_3'
    ), arrows=2, values=1),

    # Build latent derivs
    mxPath(from='dnervous_1', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval1),
    mxPath(from='dnervous_2', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval2),
    mxPath(from='dnervous_3', to=c('nervous_1', 'nervous_2', 'nervous_3'), free=F, values=dval3),

    mxPath(from='ddown_1', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval1),
    mxPath(from='ddown_2', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval2),
    mxPath(from='ddown_3', to=c('down_1', 'down_2', 'down_3'), free=F, values=dval3),

    mxPath(from='ddepressed_1', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval1),
    mxPath(from='ddepressed_2', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval2),
    mxPath(from='ddepressed_3', to=c('depressed_1', 'depressed_2', 'depressed_3'), free=F, values=dval3),

    mxPath(from='dcalm_1', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval1),
    mxPath(from='dcalm_2', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval2),
    mxPath(from='dcalm_3', to=c('calm_1', 'calm_2', 'calm_3'), free=F, values=dval3),

    mxPath(from='dhappy_1', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval1),
    mxPath(from='dhappy_2', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval2),
    mxPath(from='dhappy_3', to=c('happy_1', 'happy_2', 'happy_3'), free=F, values=dval3)
  )

  expect_equal(mmod_model$A, mmod_model_ref$A)
  expect_equal(mmod_model$S, mmod_model_ref$S)
  expect_equal(mmod_model$F, mmod_model_ref$F)
})
