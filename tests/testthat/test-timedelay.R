context('Time Delay Embedding Tests')

set.seed(9000)

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

mock_data <- reshape(
  data.frame(
    occ1 = rnorm(10),
    occ2 = rnorm(10),
    occ3 = rnorm(10),
    occ4 = rnorm(10),
    occ5 = rnorm(10),
    occ6 = rnorm(10)
  ), varying=paste0("occ", 1:6), v.names="occ", direction="long")

lapply_rbind <- function(x, f) {
  Reduce(rbind, lapply(x, f))
}

time_delay_mapping <- function(ids, times, n_embed) {
  n <- length(times) - n_embed + 1
  lapply_rbind(ids, function(curr_id) {
    lapply_rbind(1:n, function(curr_n) {
      data.frame(time_embed = 1:n_embed,
                 embed_n = curr_n,
                 time = times[seq(curr_n, length.out = n_embed)],
                 id = curr_id
      )
    })
  })
}

time_delay_ref <- function(data, idvar, timevar, n_embed) {
  times <- sort(unique(data[[timevar]]))
  ids <- unique(data[[idvar]])
  mapping <- time_delay_mapping(ids, times, n_embed)
  mapping$embed_n <- paste0(mapping$id, "_", mapping$embed_n)
  names(mapping) <- c(paste0(timevar, "_embed"), paste0(idvar, "_embed"), timevar, idvar)
  merge(data, mapping)
}

shuffle <- function(data) {
  data[sample(nrow(data)), ]
}

test_that("embedding dimension mock data", {
  df <- mock_data
  for (i in 1:6) {
    expect_equal(
      time_delay_embed(df, "id", "time", i),
      time_delay_ref(df, "id", "time", i)
    )
  }
})

test_that("embedding dimension mock data shuffle", {
  df <- shuffle(mock_data)
  for (i in 1:6) {
    expect_equal(
      time_delay_embed(df, "id", "time", i),
      time_delay_ref(df, "id", "time", i)
    )
  }
})

test_that("embedding dimension real data", {
  df <- nlsy97_1983_cohort()
  pids <- sample(unique(df$pid), 100)
  df <- df[df$pid %in% pids, ]
  for (i in 3) {
    expect_equal(
      time_delay_embed(df, "pid", "occasion", i),
      time_delay_ref(df, "pid", "occasion", i)
    )
  }
})

test_that("basic two factor model works with embed_dim", {
  require(OpenMx)
  require(mxmmod)

  df <- nlsy97_1983_cohort()

  structure <- list(
    F1 = c('nervous', 'down', 'depressed'),
    F2 = c('calm', 'happy')
  )
  mmod_model <- mxMmodModel(data=df,
                            modelName='2 Factor MMOD',
                            idvar='pid', timevar='occasion', structure=structure, embed_dim=2, fiml=F)
  summary(mxRun(mmod_model))

  dval1 <- c(1, 1)
  dval2 <- c(-0.5, 0.5)

  mmod_model_ref <- mxModel(
    type='RAM',
    manifestVars=c('nervous_1', 'nervous_2',
                   'down_1', 'down_2',
                   'depressed_1', 'depressed_2',
                   'calm_1', 'calm_2',
                   'happy_1', 'happy_2'
    ),
    latentVars=  c('F1_1', 'F2_1', 'F1_2', 'F2_2',
                   'dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                   'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2'
    ),

    # Factor structures
    mxPath(from='F1_1', to=c('dnervous_1', 'ddown_1', 'ddepressed_1'), values=0.5),
    mxPath(from='F1_2', to=c('dnervous_2', 'ddown_2', 'ddepressed_2'), values=0.5),

    mxPath(from='F2_1', to=c('dcalm_1', 'dhappy_1'), values=0.5),
    mxPath(from='F2_2', to=c('dcalm_2', 'dhappy_2'), values=0.5),

    # Factor Covariances
    mxPath(from=c('F1_1', 'F1_2', 'F2_1', 'F2_2'),
           to=c('F1_1', 'F1_2', 'F2_1', 'F2_2'), 'unique.bivariate', arrows=2),

    # Factor variances fixed @ 1
    mxPath(from=c('F1_1', 'F1_2', 'F2_1', 'F2_2'),
           to=c('F1_1', 'F1_2', 'F2_1', 'F2_2'), 'single', arrows=2, values=1, free=F),

    # residual variances (for latent derivs)
    mxPath(from=c('dnervous_1', 'ddown_1', 'ddepressed_1', 'dcalm_1', 'dhappy_1',
                  'dnervous_2', 'ddown_2', 'ddepressed_2', 'dcalm_2', 'dhappy_2'
    ), arrows=2, values=1),

    # Build latent derivs
    mxPath(from='dnervous_1', to=c('nervous_1', 'nervous_2'), free=F, values=dval1),
    mxPath(from='dnervous_2', to=c('nervous_1', 'nervous_2'), free=F, values=dval2),

    mxPath(from='ddown_1', to=c('down_1', 'down_2'), free=F, values=dval1),
    mxPath(from='ddown_2', to=c('down_1', 'down_2'), free=F, values=dval2),

    mxPath(from='ddepressed_1', to=c('depressed_1', 'depressed_2'), free=F, values=dval1),
    mxPath(from='ddepressed_2', to=c('depressed_1', 'depressed_2'), free=F, values=dval2),

    mxPath(from='dcalm_1', to=c('calm_1', 'calm_2'), free=F, values=dval1),
    mxPath(from='dcalm_2', to=c('calm_1', 'calm_2'), free=F, values=dval2),

    mxPath(from='dhappy_1', to=c('happy_1', 'happy_2'), free=F, values=dval1),
    mxPath(from='dhappy_2', to=c('happy_1', 'happy_2'), free=F, values=dval2)
  )

  expect_equal(mmod_model$A, mmod_model_ref$A)
  expect_equal(mmod_model$S, mmod_model_ref$S)
  expect_equal(mmod_model$F, mmod_model_ref$F)
})
