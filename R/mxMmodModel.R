#' Create an MMOD
#'
#' This function builds a Measurement Model of Derivatives (MMOD; Estabrook 2015) with a given
#' factor structure.
#'
#' @param data a data frame with measurements in long format
#' @param modelName name for the resulting model
#' @param idvar name of column for subject IDs
#' @param timevar name of column for measurement occasion
#' @param structure factor structure, see 'Details'
#' @param orthogonal if true, fix correlations between factors to 0
#'                   (A factor and its derivatives will still intercorrelate)
#' @param embed_dim time delay embedding dimension
#' @param fiml if true, use raw data to fit model with FIML. Otherwise, fit using cov matrix
#'             (dropping missing values if necessary).
#' @return an MMOD as an mxModel object
#'
#' @details
#'
#' The \code{structure} argument is a list of latent factors and their mappings to manifest
#' variables. For example, a one factor structure would be:
#'
#' \code{list(F1 = c('m1', 'm2', 'm3', 'm4', 'm5', 'm6'))}
#'
#' And a two factor structure would be:
#'
#' \code{list(F1 = c('m1', 'm2', 'm3'), F2 = c('m4', 'm5', 'm6'))}
#'
#' @examples
#' data(nlsy97depression)
#' # Fit one factor MMOD
#' structure <- list(
#'   F1 = c('nervous', 'down', 'depressed', 'calm', 'happy')
#' )
#' mmod_model <- mxMmodModel(data=nlsy97depression,
#'                           modelName='1 Factor MMOD',
#'                           idvar='pid', timevar='occasion', structure=structure)
#' mmod_fit <- OpenMx::mxRun(mmod_model)
#' summary(mmod_fit)
#' @export

mxMmodModel <- function(data, modelName, idvar, timevar, structure, orthogonal=F, embed_dim=NULL, fiml=F) {
  if (!is.null(embed_dim)) {
    data <- time_delay_embed(data, idvar, timevar, embed_dim)
    idvar <- paste0(idvar, "_embed")
    timevar <- paste0(timevar, "_embed")
  }

  derivName <- function(o, m) {paste0('d', m, '_', o)} # derivName(1, 'nervous') -> dnervous_1
  itemName <- function(o, m) {paste0(m, '_', o)} # itemName(1, 'nervous') -> nervous_1
  factorName <- function(o, f) {paste0(f, '_', o)} # factorName(1, 'F') -> F_1

  if (!idvar %in% colnames(data)) {
    stop("Column '", idvar, "' not found")
  }

  if (!timevar %in% colnames(data)) {
    stop("Column '", timevar, "' not found")
  }

  if (!is.numeric(data[[timevar]])) {
    stop('timevar must be a numeric type')
  }

  # Convert data[[timevar]] to 1,2,3,etc.
  occasions_num <- sort(unique(data[[timevar]]))
  data[[timevar]] <- as.numeric(factor(data[[timevar]]))
  occasions <- 1:length(occasions_num)

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
    measures_flat <- unique(unlist(structure, use.names=F))
    tmp <- lapply(measures_flat, function(m) {
      sapply(occasions, function(oo) {itemName(oo, m)})
    })
    names(tmp) <- derivName(o, measures_flat)
    tmp
  })

  factors <- names(factorStruct)
  derivatives <- unique(unlist(factorStruct, use.names=F))
  manifests <- unique(unlist(derivStruct))

  data <- data[c(idvar, timevar, unique(unlist(structure)))]
  data <- stats::reshape(as.data.frame(data), timevar=timevar, idvar=idvar, direction='wide', sep='_')[-1]
  stopifnot(setequal(manifests, names(data))) # Sanity check

  # OpenMx manifest ordering bug: https://github.com/OpenMx/OpenMx/issues/247
  data <- data[manifests]

  if (fiml) {
    mxd <- OpenMx::mxData(data, type="raw")
  } else {
    if (any(is.na(data))) {
      warning('Missing values detected; omitting them.')
    }
    df_subset <- stats::na.omit(data)
    df_cov <- stats::cov(df_subset)
    mxd <- OpenMx::mxData(df_cov, type="cov", numObs=nrow(df_subset))
  }

  # Make weight matrix with Deboeck’s functions
  weight <- ContrastsGOLD(occasions_num, length(occasions_num) - 1)
  weightList <- as.list(as.data.frame(t(weight)))

  do.call(OpenMx::mxModel, c(list(
    modelName, mxd, type="RAM",
    manifestVars=manifests,
    latentVars=c(factors, derivatives),
    # factor loadings
    mapply(function(fct, drv) {
      OpenMx::mxPath(from=fct, to=drv, values=0.5, free=T)
    }, names(factorStruct), factorStruct),
    # factor variances
    OpenMx::mxPath(from=factors, arrows=2, values=1, free=F),
    # factor correlations
    if (orthogonal) {
      lapply(names(structure), function(f) {
        OpenMx::mxPath(from=factorName(occasions, f), arrows=2, connect='unique.bivariate', free=T)
      })
    } else {
      OpenMx::mxPath(from=factors, arrows=2, connect="unique.bivariate", free=T)
    },
    # residual variances(only for latent derivatives !)
    OpenMx::mxPath(from=derivatives, arrows=2, values=1)),
    # transformation
    mapply(function(dgrp, weight) {
       mapply(function(drv, mnf) {
        OpenMx::mxPath(from=drv, to=mnf, free=F, values=weight)
      }, names(dgrp), dgrp)
    }, derivStruct, weightList),
    # saturate model
    if (fiml) OpenMx::mxPath(from = 'one', to = manifests) else list()
  ))
}

#' Generate GOLD contrasts
#'
#' Code from Deboeck (2010)
#'
#' @param T timeseries
#' @param max max derivatives
#'
#' @keywords internal
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

#' Generate Time Delay Embeddings
#'
#' @param data a data frame with measurements in long format
#' @param idvar name of column for subject IDs
#' @param timevar name of column for measurement occasion
#' @param n_embed embedding dimension
#'
#' @keywords internal
time_delay_embed <- function(data, idvar, timevar, n_embed) {
  unique_occasions <- sort(unique(data[[timevar]]))
  unique_ids <- unique(data[[idvar]])

  n_copies <- length(unique_occasions) - n_embed + 1

  embedings <- unlist(lapply(1:n_copies, function (i) unique_occasions[i:(i+n_embed-1)]))

  embed_map <- data.frame(
    idvar = rep(unique_ids, each=length(embedings)),
    timevar = rep(embedings, length(unique_ids)),
    embedvar = rep(rep(1:n_embed, n_copies), length(unique_ids))
  )
  names(embed_map) <- c(idvar, timevar, paste0(timevar, "_embed"))
  embed_map[paste0(idvar, "_embed")] <-
    paste0(embed_map[[idvar]], "_", rep(rep(1:n_copies, each=n_embed), length(unique_ids)))

  merge(data, embed_map)
}
