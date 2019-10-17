# Public domain data from nlsy97: https://www.bls.gov/nls/nlsy97.htm

var_list <- list(
  pid = 'R0000100',
  # sex = 'R0536300',
  # birth_m = 'R0536401',
  # birth_y = 'R0536402',
  # sample_type = 'R1235800',
  # eth = 'R1482600',
  nervous0 = 'R4893600',
  calm0 = 'R4893700',
  down0 = 'R4893800',
  happy0 = 'R4893900',
  depressed0 = 'R4894000',
  nervous1 = 'S0920800',
  calm1 = 'S0920900',
  down1 = 'S0921000',
  happy1 = 'S0921100',
  depressed1 = 'S0921200',
  nervous2 = 'S4681900',
  calm2 = 'S4682000',
  down2 = 'S4682100',
  happy2 = 'S4682200',
  depressed2 = 'S4682300'
)

nlsy97depression <- read.csv('data-raw/nlsy97_subset.csv')
nlsy97depression <- nlsy97depression[unlist(var_list)]
colnames(nlsy97depression) <- names(var_list)
nlsy97depression[nlsy97depression < 0] = NA
nlsy97depression <- reshape(nlsy97depression,
                            idvar='pid', timevar='occasion',
                            varying=names(var_list[-1]),
                            direction='long', sep='')

usethis::use_data(nlsy97depression, overwrite=T)
