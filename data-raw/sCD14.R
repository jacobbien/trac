## code to prepare `sCD14` dataset goes here

sCD14 <- readRDS("data-raw/sCD14.RDS")
usethis::use_data(sCD14, overwrite = TRUE)
