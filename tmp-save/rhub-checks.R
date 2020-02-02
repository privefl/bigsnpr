ENV_VARS <-  c("_R_CHECK_FORCE_SUGGESTS_" = "FALSE",
               "_R_CHECK_CRAN_INCOMING_"  = "TRUE",
               "R_C_BOUNDS_CHECK"         = "yes")

# rhub::check_on_solaris(env_vars = ENV_VARS, show_status = FALSE) # NOK -> invalid permissions
rhub::check_on_windows(env_vars = ENV_VARS, show_status = FALSE) # OK
rhub::check_on_centos (env_vars = ENV_VARS, show_status = FALSE)

rhub::check_with_sanitizers(env_vars = ENV_VARS, show_status = FALSE) # TIMEOUT
rhub::check_with_valgrind  (env_vars = ENV_VARS, show_status = FALSE) # TIMEOUT

devtools::check_win_devel()  # OK
