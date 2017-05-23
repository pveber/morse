#' Create a list of matrices from dataset to use in Bayesian modelling (JAGS or Stan)
#'
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return A list of matrix qith data for Bayesian modeling
#'
#' @export
#'
#'
#'


gm_matrix = function(gm_survData) {
  
  df.matrix = gm_survData %>%
    filter(time != 0) %>%
    arrange(profile, time) # ARRANGE To make proper data matrix
  
  count_data = count(df.matrix, profile)
  max.count_data = max(count_data$n)
  
  ls_out = split( df.matrix , f = df.matrix$profile )
  
  mat_time = ls_out%>%
    lapply("[[", "time")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_tprec = ls_out%>%
    lapply("[[", "tprec")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_Nsurv = ls_out%>%
    lapply("[[", "Nsurv")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_Nprec = ls_out%>%
    lapply("[[", "Nprec")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_conc = ls_out%>%
    lapply("[[", "conc")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  # mat_replicate = ls_out%>%
  #   lapply("[[", "replicate")%>%
  #   lapply(`length<-`, max.count_data)%>%
  #   unlist()%>%
  #   matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_profile = ls_out%>%
    lapply("[[", "profile")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  gm.listData = list(conc = mat_conc,
                     Nsurv = mat_Nsurv,
                     Nprec = mat_Nprec,
                     time = mat_time,
                     tprec = mat_tprec,
                     profile = mat_profile,
                     count_data = count_data$n)
  
  return(gm.listData)
  
}


gm_matrix_edo = function(gm_data_surv, gm_data_conc) {
  
  profile_Ncount = length(unique(gm_data_conc$profile)) ## Number of profile
  
  conc_Ncount = count(gm_data_conc, profile)$n ## Number of conc measure in each profile
  conc_Ncount_cumsum = c(1,cumsum(conc_Ncount))

  vector_conc = gm_data_conc$conc
  vector_conc_time = gm_data_conc$time

  
  df.matrix = gm_data_surv%>%
    filter(time !=0)%>%
    arrange(profile, time)
    
  surv_Ncount = count(df.matrix, profile)$n ## Number of Nsurv measure in each profile
  max.count_data = max(surv_Ncount)
  
  ls_out = split(df.matrix , f = df.matrix$profile )
  
  mat_time = ls_out%>%
    lapply("[[", "time")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  # mat_tprec = ls_out%>%
  #   lapply("[[", "tprec")%>%
  #   lapply(`length<-`, max.count_data)%>%
  #   unlist()%>%
  #   matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_Nsurv = ls_out%>%
    lapply("[[", "Nsurv")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_Nprec = ls_out%>%
    lapply("[[", "Nprec")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  # mat_conc = ls_out%>%
  #   lapply("[[", "conc")%>%
  #   lapply(`length<-`, max.count_data)%>%
  #   unlist()%>%
  #   matrix(ncol = max.count_data, byrow = TRUE)
  
  # mat_replicate = ls_out%>%
  #   lapply("[[", "replicate")%>%
  #   lapply(`length<-`, max.count_data)%>%
  #   unlist()%>%
  #   matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_profile = ls_out%>%
    lapply("[[", "profile")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  gm.listData = list(#conc = mat_conc,
                     conc = vector_conc,
                     conc_time = vector_conc_time,
                     Nsurv = mat_Nsurv,
                     Nprec = mat_Nprec,
                     time = mat_time,
                     #tprec = mat_tprec,
                     profile = mat_profile,
                     profile_Ncount = profile_Ncount,
                     conc_Ncount = conc_Ncount,
                     conc_Ncount_cumsum = conc_Ncount_cumsum,
                     surv_Ncount = surv_Ncount)
  
  return(gm.listData)
  
}

