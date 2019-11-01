if (identical(Sys.getenv("BIGSNPR_CRAN"), "false") && bigsnpr:::get_os() == "Unix") {

  download.file("https://bit.ly/2TbSaEI", destfile = (liftOver <- tempfile()))

  # X    X chromosome                    -> 23
  # Y    Y chromosome                    -> 24
  # XY   Pseudo-autosomal region of X    -> 25
  # MT   Mitochondrial                   -> 26
  df <- dplyr::mutate(bigsnpr::LD.wiki34, pos = Start, chr = dplyr::case_when(
    Chr == 23 ~ "X",
    Chr == 24 ~ "Y",
    Chr == 25 ~ "XY",
    Chr == 26 ~ "M",
    TRUE ~ as.character(Chr)
  ))

  df.new <- snp_modifyBuild(df, liftOver, from = "hg18", to = "hg19")
  with(df.new, round(pos[Chr != 23], -5))

  # https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
  Start <- c(4.8e+07, 8.6e+07, 134500000, 1.83e+08, 47500000, 83500000,
             8.9e+07, 44500000, 9.8e+07, 1.29e+08, 135500000, 2.5e+07,
             5.7e+07, 1.4e+08, 5.5e+07, 7e+06, 4.3e+07, 1.12e+08, 3.7e+07,
             4.6e+07, 87500000, 3.3e+07, 109500000, 3.2e+07)
  Stop <- c(5.2e+07, 100500000, 1.38e+08, 1.9e+08, 5e+07, 8.7e+07,
            97500000, 50500000, 100500000, 1.32e+08, 138500000, 3.5e+07,
            6.4e+07, 142500000, 6.6e+07, 1.3e+07, 5e+07, 1.15e+08,
            4.3e+07, 5.7e+07, 90500000, 4e+07, 1.12e+08, 34500000)
  expect_equal(with(df.new, round(2 * pos[Chr != 23], -6) / 2), Start)
