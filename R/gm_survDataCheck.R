#' Checks if an object can be used to perform survival analysis
#'
#' The \code{gm_survDataCheck} function can be used to check if an object
#' containing survival data is formatted according to the expectations of the
#' \code{gm_survData} function.
#'
#'
#' @aliases gm_survDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, the function may produce diagnosis plots
#'
#' @return
#' @seealso \code{\link{survData}}
#'
#' @examples
#' # Run the check data function
#' data(zinc)
#' gm_survDataCheck(zinc)
#'
#' # Now we insert an error in the dataset, by artificially increasing the
#' # number of survivors at some time point, in such a way that the number
#' # of indivuals increases in some replicate
#' zinc[25, "Nsurv"] <- as.integer(20)
#' gm_survDataCheck(zinc)
#'
#'
#' @export
gm_survDataCheck <- function(data_surv, data_conc) {

  ##
  ## 0. check if data.frame
  ##
  if(!("data.frame" %in% class(data_surv) & "data.frame" %in% class(data_conc))) {
    stop("'data_surv' and/or 'data_conc' are not of class 'data.frame'.")
  }

  ##
  ## 1. check column names
  ##

    ### 1.1 data_surv
    ref.namesSurv <- c("profile","replicate","time","Nsurv")

    if ( ncol(data_surv) != 4  # check same number of column
         | !(all(colnames(data_surv) %in% ref.namesSurv))  ) { # check if the name are the same, order doesn't matter
      stop( "'data_surv' must have 4 columns with names: 'profile','replicate','time' and 'Nsurv' (order does not matter)." )
    }

    ### 1.2 data_conc
    ref.namesConc <- c("profile","time","conc")

    if ( ncol(data_conc) != 3  # check same number of column
         | !(all(colnames(data_conc) %in% ref.namesConc))  ) { # check if the name are the same, order doesn't matter
      stop( "'data_conc' must have 3 columns with names: 'profile', 'time' and 'conc' (order does not matter)." )
    }

  ##
  ## 2. assert the first time point is zero for each (profile, replicate)
  ##
  df_check0Surv = data_surv %>%
    group_by(profile,replicate) %>%
    summarize(check_0 = 0 %in% time)
  if(all(df_check0Surv$check_0) != TRUE ){
    stop( "In each '(profile, replicate)' of 'data_surv', the first time point must be '0'." )
  }

  df_check0Conc = data_conc %>%
    group_by(profile) %>%
    summarize(check_0 = 0 %in% time)

  if(all(df_check0Conc$check_0) !=TRUE ){
    stop( "In each 'profile' of 'data_conc', the first time point must be '0'." )
  }

  ##
  ## 3. time are numeric, Nsurv are interger, assert concentrations are numeric
  ##
  ### 3.1 data_surv
  if (!is.integer(data_surv$Nsurv) | any(data_surv$Nsurv < 0)
      | !is.numeric(data_surv$time) | any(data_surv$time < 0)      ) {
    stop( "In table 'data_surv', column 'time' must contain only positive numeric values and column 'Nsurv' must contain only positive integer values. " )
  }
  ### 3.2 data_conc
  if (!is.numeric(data_conc$conc) | any(data_conc$conc < 0)
      | !is.numeric(data_conc$time) | any(data_conc$time < 0)      ) {
    stop( "In table 'data_conc', column 'time' must contain only positive numeric values and column 'conc' must contain only positive integer values. " )
  }


  ##
  ## 4. assert Nsurv != 0 at time 0
  ##

  df_checkNeq0 = data_surv %>%
    filter(time == 0) %>%
    mutate( check_NsurvNeq0 = (Nsurv != 0))

  if(all(df_checkNeq0$check_NsurvNeq0) != TRUE){
    stop( "In each '(profile, replicate)', 'Nsurv' must be greater than '0'." )
  }

  # ##
  # ## 5. assert there is the same name and/or number of replicates between tables
  # ##
  #
  # df_checkReplicate_surv = data_surv %>%
  #   group_by(profile) %>%
  #   summarize( check_Replicate = length(unique(replicate)))
  #
  # df_checkReplicate_conc = data_conc %>%
  #   group_by(profile) %>%
  #   summarize( check_Replicate = length(unique(replicate)))
  #
  # if( all( df_checkReplicate_surv$check_Replicate == df_checkReplicate_conc$check_Replicate ) != TRUE){
  #   stop( "At least one 'replicate' is missing or different between tables 'data_surv' and 'data_conc'." )
  # }

  ##
  ## 6. assert there is the same name and/or number of profile between tables
  ##

  df_checkProfile_surv = data_surv %>%
    summarize( check_Profile = length(unique(profile)))

  df_checkProfile_conc = data_conc %>%
    summarize( check_Profile = length(unique(profile)))

  if( all( df_checkProfile_surv$check_Profile == df_checkProfile_conc$check_Profile ) != TRUE){
    stop( "At least one 'profile' is missing or different between tables 'data_surv' and 'data_conc'." )
  }

  ##
  ## 7. assert Nsurv never increases with time
  ##

  df_checkSurvIncrease = data_surv %>%
    group_by(profile, replicate) %>%
    arrange(time) %>%
    mutate(Nprec = ifelse(time == 0, Nsurv, lag(Nsurv))) %>%
    mutate( check_SurvIncrease = Nsurv <= Nprec)

  if(all(df_checkSurvIncrease$check_SurvIncrease) != TRUE){
    stop( "'Nsurv' increases at some time points." )
  }


  # ##
  # ## 8. check if the profile of concentration are the same for every `replicate` of a same `profile`
  # ##
  #
  # df_checkProfile = data_conc %>%
  #   group_by(profile, time) %>%
  #   summarize(checkProfile = length(unique(conc)))
  #
  # if(any(df_checkProfile$checkProfile) != 1){
  #   stop( "Elements in 'conc' must be the same between 'replicate' of a same '(profile, time)'." )
  # }

  ##
  ## 9. WARNING - assert max(time in data_conc) >= max(time in data_surv)
  ##

  df_checkMaxTimeSurv = data_surv %>%
    group_by(profile) %>%
    filter(time == max(time))

  df_checkMaxTimeConc = data_conc %>%
    group_by(profile) %>%
    filter(time == max(time))

  if(!all(df_checkMaxTimeConc$time >= df_checkMaxTimeSurv$time) ){
    warning( "In each 'profile', maximum time in 'data_conc' should be greater or equal to maximum time in 'data_surv'.
             Otherwise, last concentration is taken to fill concentration profile until the maximum time in 'data_surv'." )
  }

  ##
  ## 10. for every profile, 'time' must be the same in 'data_surv'
  ##

  df_checkTimeSurv = data_surv %>%
    mutate(concate = paste0(profile,replicate)) %>%
    mutate(time_ref = time)%>%
    select(concate, time, time_ref) %>%
    spread(key = concate, value = time)

  if(length(unique(as.list(df_checkTimeSurv))) != 1) {
    stop( "In 'data_surv' all 'time' vector must be the same." )
  }


}

