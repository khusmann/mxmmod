# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

mxMmodModel <- function(data, modelName, idvar, timevar, structure, fiml=F) {
  derivName <- function(o, m) {paste0('d', m, '_', o)} # derivName(1, 'nervous') -> dnervous_1
  itemName <- function(o, m) {paste0(m, '_', o)} # itemName(1, 'nervous') -> nervous_1
  factorName <- function(o, f) {paste0(f, '_', o)} # factorName(1, 'F') -> F_1
  occasions <- unique(data[[timevar]])

  # factorStruct
  # input:
  #   (occasions = c(1, 2, 3))
  #   list (
  #     F = c('nervous', 'down', 'anxious')
  #   )
  # output:
  #   list (
  #     F1 = c('dnervous_1', 'ddown_1', 'danxious_1')
  #     F2 = c('dnervous_2', 'ddown_2', 'danxious_2')
  #     F3 = c('dnervous_3', 'ddown_3', 'danxious_3')
  #   )
  factorStruct <- unlist(lapply(occasions, function(o) {
    tmp <- lapply(structure, function(s) {derivName(o, s)}) # Create all derivs under factor
    names(tmp) <- factorName(o, names(tmp)) # Create factor names
    tmp
  }), recursive=F)

  # derivStruct
  # input:
  #   (occasions = c(1, 2, 3))
  #   list (
  #     F = c('nervous', 'down', 'anxious')
  #   )
  # output:
  #   list (
  #     d1nervious = c('nervous_1', 'nervous_2', 'nervous_3')
  #     d2nervious = c('nervous_1', 'nervous_2', 'nervous_3')
  #     d3nervious = c('nervous_1', 'nervous_2', 'nervous_3')
  #     d1down = c('down_1', 'down_2', 'down_3')
  #     d2down = c('down_1', 'down_2', 'down_3')
  #     d3down = c('down_1', 'down_2', 'down_3')
  #     d1anxious = c('anxious_1', 'anxious_2', 'anxious_3')
  #     d2anxious = c('anxious_1', 'anxious_2', 'anxious_3')
  #     d3anxious = c('anxious_1', 'anxious_2', 'anxious_3')
  #   )
  derivStruct <- lapply(occasions, function(o) {
    measures_flat <- unlist(structure, use.names=F)
    tmp <- lapply(measures_flat, function(m) {
      sapply(occasions, function(oo) {itemName(oo, m)})
    })
    names(tmp) <- derivName(o, measures_flat)
    tmp
  })

  factors <- names(factorStruct)
  derivatives <- unlist(factorStruct, use.names=F)
  manifests <- unique(unlist(derivStruct))

  data <- data[c(idvar, timevar, unlist(structure))]
  data <- reshape(as.data.frame(data), timevar=timevar, idvar=idvar, direction='wide', sep='_')[-1]
  stopifnot(setequal(manifests, names(data)))
  # TODO: Does the cols in the data have to be in the same order as the manifest list?
  data <- data[manifests]

  if (fiml) {
    mxd <- mxData(data, type="raw")
  } else {
    df_subset <- na.omit(data)
    # TODO: choose between polychoric / cor / cov
    #df_cors <- polychoric(df_subset)$rho
    df_cors <- cor(df_subset)
    mxd <- mxData(df_cors, type="cov", numObs=nrow(df_subset))
  }

  # TODO: CHECK THIS
  # Make weight matrix with Deboeck’s functions
  weight <- ContrastsGOLD(occasions, length(occasions) - 1)
  weightList <- as.list(as.data.frame(t(weight)))

  # Make weight matrix without Deboeck’s functions
  #weight <- matrix(c(1/3, 1/3, 1/3,
  #                   -1,   0, 1 ,
  #                   1/2,  -1, 1/2),
  #                 nrow=3, ncol=3)
  #weightList <- as.list(as.data.frame(t(solve(weight))))

  do.call('mxModel', c(list(
    modelName, mxd, type="RAM",
    manifestVars=manifests,
    latentVars=c(factors, derivatives),
    # factor loadings
    mapply(function(fct, drv) {
      mxPath(from=fct, to=drv, values=0.5, free=T)
    }, names(factorStruct), factorStruct),
    # factor variances
    mxPath(from=factors, arrows=2, values=1, free=F),
    # factor correlations
    mxPath(from=factors, arrows=2, connect="unique.bivariate", free=T),
    # residual variances(only for latent derivatives !)
    mxPath(from=derivatives, arrows=2, values=1)),
    # transformation
    mapply(function(dgrp, weight) {
       mapply(function(drv, mnf) {
        mxPath(from=drv, to=mnf, free=F, values=weight)
      }, names(dgrp), dgrp)
    }, derivStruct, weightList),
    # saturate model
    if (fiml) mxPath(from = 'one', to = manifests) else list()
  ))
}

# Code from Deboeck (2010)
ContrastsGOLD <- function(T, max) {
  Xi <- matrix(NA, length(T), length(T))
  for(r in 0:(length(T)-1)) {
    Xi[r+1,] <- (T^r)
    if(r>0) {
      for(p in 0:(r-1)) {
        Xi[r+1,] <- Xi[r+1,] - (Xi[p+1,]*(sum(Xi[p+1,]*(T^r)))/(sum(Xi[p+1,]*(T^p))))
      }
    }
  }
  return(Xi[1:(max+1),])
}
