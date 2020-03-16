# mxmmod

This package provides a convenient interface for building Measurement Model of Derivatives (MMOD; Estabrook, 2015) in OpenMx.

## Installing

Install in R with:

```
devtools::install_github('khusmann/mxmmod')
```

## Usage

```
data(nlsy97depression)
# Fit one factor MMOD
structure <- list(
  F1 = c('nervous', 'down', 'depressed', 'calm', 'happy')
)
mmod_model <- mxMmodModel(data=nlsy97depression,
                          modelName='1 Factor MMOD',
                          idvar='pid', timevar='occasion', structure=structure)
mmod_fit <- OpenMx::mxRun(mmod_model)
summary(mmod_fit)
```
