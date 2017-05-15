#' Test if concentration of profile in 'gm_survData' is constant
#' 
#' @param gm_survData An object of class \code{gm_survData}
#' 
#' @return a boolean \code{TRUE} if concentration in profile is constant, or \code{FALSE} if the concentration in at least one of the profile is variable
#' 
#' @export

test_constant_conc = function(gm_survData){
  df.test = gm_survData %>%
    group_by(replicat)%>%
    summarise(count = n_distinct(conc))
  return(all(df.test$count==1))
}

#' Create a dataset for survival analysis when the profile of concentration is variable
#' 
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return A dataframe of class \code{gm_survData}.
#'
#' @export
#'
#'

gm_survData_interpolate = function (gm_survData){
  
  df.survData_interpolate = gm_survData %>%
    group_by(profile, replicate) %>%
    arrange(time) %>%
    mutate(conc.origin = conc) %>%
    mutate(conc = na.approx(conc,time, na.rm=FALSE)) %>%
    mutate(id_conc_interp = ifelse(is.na(Nsurv), NA, row_number()))%>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc)) %>%# FOR THE LAST VALUE
    arrange(profile, replicate, time)
  
  class(df.survData_interpolate) = c("gm_survData", "data.frame")
  return(df.survData_interpolate)    
}


gm_listData = function(gm_survData){
  
  df.listData = gm_survData %>%
    filter(time !=0) %>% # Remove time = 0
    group_by(profile, time) %>%
    summarise(Nsurv = sum(Nsurv, na.rm = TRUE),
              conc = unique(conc),)
  
  
    arrange(profile, replicate, time) # ARRANGE To make proper data matrix
  
  count_data = count(gm_survData, profile, replicate)
  max.count_data = max(count_data$n)
  
  ls_out = split( gm_survData , f = gm_survData$replicat )
  
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
  
  mat_replicat = ls_out%>%
    lapply("[[", "replicat")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  gm.listData = list(conc=mat_conc,
                     Nsurv=mat_Nsurv,
                     Nprec=mat_Nprec,
                     time=mat_time,
                     tprec=mat_tprec,
                     replicat = mat_replicat,
                     count_data = count_data$n)
  
  return(gm.listData)
  
}

#' Create a list of scalars giving priors to use in Bayesian modelling (JAGS or Stan)
#'
#' @param gm_survData An object of class 'gm_survData'
#'
#' @return A list of scalar for parameterization of priors for Bayesian modeling
#'
#' @export
#'


gm_bayesData = function(gm_survData){
  
  
  test_cst_conc = test_constant_conc(gm_survData = gm_survData)
  
  if(test_cst_conc){
    
    # return priors for model.
    gm.jagsPriors = gm_jagsPriors(gm_survData = gm_survData)
    # return list of Data used in JAGS model
    gm.listData = gm_listData(gm_survData = gm_survData)
    
  } else{
    
    df.long = df.long_interpolate(gm_survData,
                                  extend.time=extend_time)
    
    gm_survData = gm_survData%>%
      select(-tprec,-Nprec)%>%
      filter(!is.na(Nsurv))%>% # REMOVE NA
      arrange(replicat,time)%>%
      # 'lag' function copy values lagged by 1 (see 'dplyr' package)
      mutate(tprec = ifelse(time==0, time, lag(time))) %>%
      mutate(Nprec = ifelse(time==0, Nsurv, lag(Nsurv)))
    
    # return priors for model.
    gm.jagsPriors = gm_jagsPriors(gm_survData = gm_survData)
    # return list of Data used in JAGS model
    gm.listData = gm_listData(gm_survData = gm_survData)
    
    MYintCList = jags.intCList(gm_survData, df.long)
  }
  
  if(is.null(model_type)){
    stop("Please provide model.type SD or IT")
  }
  
  
  if(model_type=="IT"){
    if(test_cst_conc){ # return TRUE if conc are constant by replicat, FALSE otherwise
      jagsData.list = list(
        x = gm.listData$conc,
        y = gm.listData$Nsurv,
        
        Nprec = gm.listData$Nprec,
        time = gm.listData$time,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.alpha = gm.jagsPriors$meanlog10.alpha,
        taulog10.alpha = gm.jagsPriors$taulog10.alpha
      )
      jagsData.list_Null = list(
        x = gm.listData$conc,
        
        Nprec = gm.listData$Nprec,
        time = gm.listData$time,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.alpha = gm.jagsPriors$meanlog10.alpha,
        taulog10.alpha = gm.jagsPriors$taulog10.alpha
      )
    } else {
      
      jagsData.list = list(
        
        y = gm.listData$Nsurv,
        Nprec = gm.listData$Nprec,
        time = gm.listData$time,
        
        # Nécessaire lorsque la concentration est variable en fonction du temps
        intC = MYintCList$intC , 
        intC.time = MYintCList$intC.time,
        id.intCtime = MYintCList$id.intCtime,
        
        N.intC = MYintCList$N.intC,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        #n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        n.time =  gm.listData$count_data, # Au cas ou il n'y ai pas le même nombre d'observation
        
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.alpha = gm.jagsPriors$meanlog10.alpha,
        taulog10.alpha = gm.jagsPriors$taulog10.alpha
        
      )
      jagsData.list_Null = list(
        
        Nprec = gm.listData$Nprec,
        time = gm.listData$time,
        
        # Nécessaire lorsque la concentration est variable en fonction du temps
        intC = MYintCList$intC , 
        intC.time = MYintCList$intC.time,
        id.intCtime = MYintCList$id.intCtime,
        
        N.intC = MYintCList$N.intC,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        #n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        n.time =  gm.listData$count_data, # Au cas ou il n'y ai pas le même nombre d'observation
        
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.alpha = gm.jagsPriors$meanlog10.alpha,
        taulog10.alpha = gm.jagsPriors$taulog10.alpha
        
      )}
  }
  else if(model_type =="SD"){
    if(test_cst_conc){ # return TRUE if conc are constant by replicat, FALSE otherwise
      jagsData.list = list(
        x = gm.listData$conc,
        y = gm.listData$Nsurv,
        
        Nprec = gm.listData$Nprec,
        tprec = gm.listData$tprec,
        time = gm.listData$time,
        
        # POURQUOI + 10 ? si on mesure à la minute sur 10 jours: 14400 minutes, on devrait plutôt MULTIPLIER non ?   
        bigtime = max(gm.listData$time)+10,
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        
        meanlog10.kk = gm.jagsPriors$meanlog10.kk,
        taulog10.kk = gm.jagsPriors$taulog10.kk,
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.z = gm.jagsPriors$meanlog10.z,
        taulog10.z = gm.jagsPriors$taulog10.z
      )
      jagsData.list_Null = list(
        x = gm.listData$conc,
        
        Nprec = gm.listData$Nprec,
        tprec = gm.listData$tprec,
        time = gm.listData$time,
        
        # POURQUOI + 10 ? si on mesure à la minute sur 10 jours: 14400 minutes, on devrait plutôt MULTIPLIER non ?   
        bigtime = max(gm.listData$time)+10,
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        
        meanlog10.kk = gm.jagsPriors$meanlog10.kk,
        taulog10.kk = gm.jagsPriors$taulog10.kk,
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.z = gm.jagsPriors$meanlog10.z,
        taulog10.z = gm.jagsPriors$taulog10.z
      )
      
    } else{
      
      jagsData.list = list(
        
        y = gm.listData$Nsurv,
        
        Nprec = gm.listData$Nprec,
        
        # Nécessaire lorsque la concentration est variable en fonction du temps
        intC = MYintCList$intC , 
        intC.time = MYintCList$intC.time,
        id.intCtime = MYintCList$id.intCtime,
        
        N.intC = MYintCList$N.intC,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        #n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        n.time =  gm.listData$count_data, # Au cas ou il n'y ai pas le même nombre d'observation
        
        meanlog10.kk = gm.jagsPriors$meanlog10.kk,
        taulog10.kk = gm.jagsPriors$taulog10.kk,
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.z = gm.jagsPriors$meanlog10.z,
        taulog10.z = gm.jagsPriors$taulog10.z
      )
      jagsData.list_Null = list(
        
        Nprec = gm.listData$Nprec,
        
        # Nécessaire lorsque la concentration est variable en fonction du temps
        intC = MYintCList$intC , 
        intC.time = MYintCList$intC.time,
        id.intCtime = MYintCList$id.intCtime,
        
        N.intC = MYintCList$N.intC,
        
        n.group = nrow(gm.listData$conc), # MUST BE IMPROVED
        #n.time = ncol(gm.listData$conc), # MUST BE IMPROVED
        n.time =  gm.listData$count_data, # Au cas ou il n'y ai pas le même nombre d'observation
        
        meanlog10.kk = gm.jagsPriors$meanlog10.kk,
        taulog10.kk = gm.jagsPriors$taulog10.kk,
        meanlog10.ke = gm.jagsPriors$meanlog10.ke,
        taulog10.ke = gm.jagsPriors$taulog10.ke,
        meanlog10.hb = gm.jagsPriors$meanlog10.hb,
        taulog10.hb = gm.jagsPriors$taulog10.hb,
        meanlog10.z = gm.jagsPriors$meanlog10.z,
        taulog10.z = gm.jagsPriors$taulog10.z
      )
      
    }}
  
  else{stop('Please provide the model.type ("IT", "SD",...)')}
  
  if(test_cst_conc){
    return(list(jagsData=jagsData.list,
                jagsData_Null = jagsData.list_Null,
                replicat = gm.listData$replicat,
                model_type = model_type,
                cst_conc = test_cst_conc))
    
  } else{
    return(list(jagsData=jagsData.list,
                jagsData_Null = jagsData.list_Null,
                replicat =  MYintCList$replicat, # C'est pour le replicat !!!
                model_type = model_type,
                cst_conc = test_cst_conc))
  }
}