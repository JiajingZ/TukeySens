## cleaning up raw data ##
library(tidyverse)
library(foreign)
NHANES <- read.dta(system.file("extdata", "NHANES3hbp_dbp.dta",
                               package = "TukeySensitivity")) %>% as_tibble

NHANES <- NHANES %>% dplyr::filter(ave_dbp > 30)

NHANES <- NHANES %>% select(-one_of(c("ave_sbp", "trt_sbp", "d_ctrl", "num_aht")))

## NHANES contains covariates , treatment "trt_dbp", outcome "ave_dbp"

## save the dataset in data/ ##
usethis::use_data(NHANES, overwrite = TRUE)
