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
  derivName <- function(o, m) {paste0('d', o, m)} # derivName(1, 'nervous') -> d1nervous
  itemName <- function(o, m) {paste0(m, '_', o)} # itemName(1, 'nervous') -> nervous_1
  occasions <- unique(data[[timevar]])

  data <- data[c(idvar, timevar, unlist(structure))]
  data <- reshape(as.data.frame(data), timevar=timevar, idvar=idvar, direction='wide', sep='_')

  # factorStruct
  # input:
  #   (occasions = c(1, 2, 3))
  #   list (
  #     F = c('nervous', 'down', 'anxious')
  #   )
  # output:
  #   list (
  #     F1 = c('d1nervous', 'd1down', 'd1anxious)
  #     F2 = c('d2nervous', 'd2down', 'd2anxious)
  #     F3 = c('d3nervous', 'd3down', 'd3anxious)
  #   )
  factorStruct <- unlist(lapply(occasions, function(o) {
    tmp <- lapply(structure, function(s) {derivName(o, s)}) # Create all derivs under factor
    names(tmp) <- paste0(names(tmp), o) # Append occasion to factor name
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
  #     d1nervious = c('nervous.1', 'nervous.2', 'nervous.3')
  #     d2nervious = c('nervous.1', 'nervous.2', 'nervous.3')
  #     d3nervious = c('nervous.1', 'nervous.2', 'nervous.3')
  #     d1down = c('down.1', 'down.2', 'down.3')
  #     d2down = c('down.1', 'down.2', 'down.3')
  #     d3down = c('down.1', 'down.2', 'down.3')
  #     d1anxious = c('anxious.1', 'anxious.2', 'anxious.3')
  #     d2anxious = c('anxious.1', 'anxious.2', 'anxious.3')
  #     d3anxious = c('anxious.1', 'anxious.2', 'anxious.3')
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

  if (fiml) {
    mxd <- mxData(data[manifests], type="raw")
  } else {
    df_subset <- na.omit(data[manifests])
    #df_cors <- polychoric(df_subset)$rho
    df_cors <- cor(df_subset)
    mxd <- mxData(df_cors, type="cov", numObs=nrow(df_subset))
  }

  # Make weight matrix without Deboeck’s functions
  weight <- matrix(c(1/3, 1/3, 1/3,
                     -1,   0, 1 ,
                     1/2,  -1, 1/2),
                   nrow=3, ncol=3)
  weightList <- as.list(as.data.frame(t(solve(weight))))

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
