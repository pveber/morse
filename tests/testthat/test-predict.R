# test predict
library(morse)

setwd("~/Documents/morse_test")

# -----------------------------------------------------------------------------
# --- Propiconazole

data_varPRZ <- survData(propiconazole_pulse_exposure)
# survFit_varSD_PRZ_new <- survFit(data_varPRZ, model_type = "SD", ode_solver = "new")
# save(survFit_varSD_PRZ_new, file = "fit/survFit_varSD_PRZ_new.rda")
load(file = "fit/survFit_varSD_PRZ_new.rda")
plot(survFit_varSD_PRZ_new)
predict_PRZ_var_to_var <- predict(survFit_varSD_PRZ_new, data_predict = propiconazole_pulse_exposure, check_to_Nsurv = TRUE)
predict_PRZ_var_to_cst <- predict(survFit_varSD_PRZ_new, data_predict = propiconazole, check_to_Nsurv = TRUE)
save(predict_PRZ_var_to_var, file = "fit/predict_PRZ_var_to_var.rda")
save(predict_PRZ_var_to_cst, file = "fit/predict_PRZ_var_to_cst.rda")
data_cstPRZ <- survData(propiconazole)
# survFit_cstSD_PRZ <- survFit(data_cstPRZ, model_type = "SD")
# save(survFit_cstSD_PRZ, file = "fit/survFit_cstSD_PRZ.rda")
load(file = "fit/survFit_cstSD_PRZ.rda")
plot(survFit_cstSD_PRZ, adddata = TRUE)

predict_PRZ_cst_to_var <- predict(survFit_cstSD_PRZ, data_predict = propiconazole_pulse_exposure, check_to_Nsurv = TRUE)
predict_PRZ_cst_to_cst <- predict(survFit_cstSD_PRZ, data_predict = propiconazole, check_to_Nsurv = TRUE)
save(predict_PRZ_cst_to_var, file = "predict_PRZ_cst_to_var.rda")
save(predict_PRZ_cst_to_cst, file = "predict_PRZ_cst_to_cst.rda")

plot(predict_PRZ_var_to_var)
plot(predict_PRZ_cst_to_var)

plot(predict_PRZ_var_to_cst)
plot(predict_PRZ_cst_to_cst)

# -----------------------------------------------------------------------------
# --- Diazinon

load(file = "data/data_Diaz_var.rda")
data_varDIAZ <- survData(data_Diazinon)
survFit_varSD_DIAZ_new <- survFit(data_varDIAZ, model_type = "SD", ode_solver = "new")
save(survFit_varSD_DIAZ_new, file = "survFit_varSD_DIAZ_new.rda")
load(file = "fit/survFit_varSD_DIAZ_new.rda")
predict_DIAZ_var_to_var <- predict(survFit_varSD_DIAZ_new, data_predict = data_Diazinon, check_to_Nsurv = TRUE)
save(predict_DIAZ_var_to_var, file = "predict_DIAZ_var_to_var.rda")



# -----------------------------------------------------------------------------
# --- TEST AS PAPER: Car ; Cyp ; Dim ; Mal

load(file = "data/data_Car_cst.rda")
load(file = "data/data_Car_var.rda")

load(file = "data/data_Cyp_cst.rda")
load(file = "data/data_Cyp_var.rda")

load(file = "data/data_Dim_cst.rda")
load(file = "data/data_Dim_var.rda")

load(file = "data/data_Mal_cst.rda")
load(file = "data/data_Mal_var.rda")


load(file = "fit/survFit_cstSD_Car.rda")
load(file = "fit/survFit_cstIT_Car.rda")
load(file = "fit/survFit_varSD_Car.rda")
load(file = "fit/survFit_varIT_Car.rda")

load(file = "fit/survFit_cstSD_Cyp.rda")
load(file = "fit/survFit_cstIT_Cyp.rda")
load(file = "fit/survFit_varSD_Cyp.rda")
load(file = "fit/survFit_varIT_Cyp.rda")

load(file = "fit/survFit_cstSD_Dim.rda")
load(file = "fit/survFit_cstIT_Dim.rda")
load(file = "fit/survFit_varSD_Dim.rda")
load(file = "fit/survFit_varIT_Dim.rda")

load(file = "fit/survFit_cstSD_Mal.rda")
load(file = "fit/survFit_cstIT_Mal.rda")
load(file = "fit/survFit_varSD_Mal.rda")
load(file = "fit/survFit_varIT_Mal.rda")
### Car
# predict_cstSD_to_cst_Car <- predict(survFit_cstSD_Car, data_predict = Car_cst, check_to_Nsurv = TRUE)
# predict_varSD_to_var_Car <- predict(survFit_varSD_Car, data_predict = Car_var, check_to_Nsurv = TRUE)
# predict_cstSD_to_var_Car <- predict(survFit_cstSD_Car, data_predict = Car_var, check_to_Nsurv = TRUE)
# predict_varSD_to_cst_Car <- predict(survFit_varSD_Car, data_predict = Car_cst, check_to_Nsurv = TRUE)
# 
# predict_cstIT_to_cst_Car <- predict(survFit_cstIT_Car, data_predict = Car_cst, check_to_Nsurv = TRUE)
# predict_varIT_to_var_Car <- predict(survFit_varIT_Car, data_predict = Car_var, check_to_Nsurv = TRUE)
# predict_cstIT_to_var_Car <- predict(survFit_cstIT_Car, data_predict = Car_var, check_to_Nsurv = TRUE)
# predict_varIT_to_cst_Car <- predict(survFit_varIT_Car, data_predict = Car_cst, check_to_Nsurv = TRUE)
# 
# save(predict_cstSD_to_cst_Car, file = "predict_cstSD_to_cst_Car.rda")
# save(predict_varSD_to_var_Car, file = "predict_varSD_to_var_Car.rda")
# save(predict_cstSD_to_var_Car, file = "predict_cstSD_to_var_Car.rda")
# save(predict_varSD_to_cst_Car, file = "predict_varSD_to_cst_Car.rda")
# 
# save(predict_cstIT_to_cst_Car, file = "predict_cstIT_to_cst_Car.rda")
# save(predict_varIT_to_var_Car, file = "predict_varIT_to_var_Car.rda")
# save(predict_cstIT_to_var_Car, file = "predict_cstIT_to_var_Car.rda")
# save(predict_varIT_to_cst_Car, file = "predict_varIT_to_cst_Car.rda")
# ### Cyp
# predict_cstSD_to_cst_Cyp <- predict(survFit_cstSD_Cyp, data_predict = Cyp_cst, check_to_Nsurv = TRUE)
# predict_varSD_to_var_Cyp <- predict(survFit_varSD_Cyp, data_predict = Cyp_var, check_to_Nsurv = TRUE)
# predict_cstSD_to_var_Cyp <- predict(survFit_cstSD_Cyp, data_predict = Cyp_var, check_to_Nsurv = TRUE)
# predict_varSD_to_cst_Cyp <- predict(survFit_varSD_Cyp, data_predict = Cyp_cst, check_to_Nsurv = TRUE)
# 
# predict_cstIT_to_cst_Cyp <- predict(survFit_cstIT_Cyp, data_predict = Cyp_cst, check_to_Nsurv = TRUE)
# predict_varIT_to_var_Cyp <- predict(survFit_varIT_Cyp, data_predict = Cyp_var, check_to_Nsurv = TRUE)
# predict_cstIT_to_var_Cyp <- predict(survFit_cstIT_Cyp, data_predict = Cyp_var, check_to_Nsurv = TRUE)
# predict_varIT_to_cst_Cyp <- predict(survFit_varIT_Cyp, data_predict = Cyp_cst, check_to_Nsurv = TRUE)
# 
# save(predict_cstSD_to_cst_Cyp, file = "predict_cstSD_to_cst_Cyp.rda")
# save(predict_varSD_to_var_Cyp, file = "predict_varSD_to_var_Cyp.rda")
# save(predict_cstSD_to_var_Cyp, file = "predict_cstSD_to_var_Cyp.rda")
# save(predict_varSD_to_cst_Cyp, file = "predict_varSD_to_cst_Cyp.rda")
# 
# save(predict_cstIT_to_cst_Cyp, file = "predict_cstIT_to_cst_Cyp.rda")
# save(predict_varIT_to_var_Cyp, file = "predict_varIT_to_var_Cyp.rda")
# save(predict_cstIT_to_var_Cyp, file = "predict_cstIT_to_var_Cyp.rda")
# save(predict_varIT_to_cst_Cyp, file = "predict_varIT_to_cst_Cyp.rda")
# 
# ### Mal
# predict_cstSD_to_cst_Mal <- predict(survFit_cstSD_Mal, data_predict = Mal_cst, check_to_Nsurv = TRUE)
# predict_varSD_to_var_Mal <- predict(survFit_varSD_Mal, data_predict = Mal_var, check_to_Nsurv = TRUE)
# predict_cstSD_to_var_Mal <- predict(survFit_cstSD_Mal, data_predict = Mal_var, check_to_Nsurv = TRUE)
# predict_varSD_to_cst_Mal <- predict(survFit_varSD_Mal, data_predict = Mal_cst, check_to_Nsurv = TRUE)
# 
# predict_cstIT_to_cst_Mal <- predict(survFit_cstIT_Mal, data_predict = Mal_cst, check_to_Nsurv = TRUE)
# predict_varIT_to_var_Mal <- predict(survFit_varIT_Mal, data_predict = Mal_var, check_to_Nsurv = TRUE)
# predict_cstIT_to_var_Mal <- predict(survFit_cstIT_Mal, data_predict = Mal_var, check_to_Nsurv = TRUE)
# predict_varIT_to_cst_Mal <- predict(survFit_varIT_Mal, data_predict = Mal_cst, check_to_Nsurv = TRUE)
# 
# save(predict_cstSD_to_cst_Mal, file = "predict_cstSD_to_cst_Mal.rda")
# save(predict_varSD_to_var_Mal, file = "predict_varSD_to_var_Mal.rda")
# save(predict_cstSD_to_var_Mal, file = "predict_cstSD_to_var_Mal.rda")
# save(predict_varSD_to_cst_Mal, file = "predict_varSD_to_cst_Mal.rda")
# 
# save(predict_cstIT_to_cst_Mal, file = "predict_cstIT_to_cst_Mal.rda")
# save(predict_varIT_to_var_Mal, file = "predict_varIT_to_var_Mal.rda")
# save(predict_cstIT_to_var_Mal, file = "predict_cstIT_to_var_Mal.rda")
# save(predict_varIT_to_cst_Mal, file = "predict_varIT_to_cst_Mal.rda")
# 
# ### Dim
# predict_cstSD_to_cst_Dim <- predict(survFit_cstSD_Dim, data_predict = Dim_cst, check_to_Nsurv = TRUE)
# predict_varSD_to_var_Dim <- predict(survFit_varSD_Dim, data_predict = Dim_var, check_to_Nsurv = TRUE)
# predict_cstSD_to_var_Dim <- predict(survFit_cstSD_Dim, data_predict = Dim_var, check_to_Nsurv = TRUE)
# predict_varSD_to_cst_Dim <- predict(survFit_varSD_Dim, data_predict = Dim_cst, check_to_Nsurv = TRUE)
# 
# predict_cstIT_to_cst_Dim <- predict(survFit_cstIT_Dim, data_predict = Dim_cst, check_to_Nsurv = TRUE)
# predict_varIT_to_var_Dim <- predict(survFit_varIT_Dim, data_predict = Dim_var, check_to_Nsurv = TRUE)
# predict_cstIT_to_var_Dim <- predict(survFit_cstIT_Dim, data_predict = Dim_var, check_to_Nsurv = TRUE)
# predict_varIT_to_cst_Dim <- predict(survFit_varIT_Dim, data_predict = Dim_cst, check_to_Nsurv = TRUE)
# 
# save(predict_cstSD_to_cst_Dim, file = "predict_cstSD_to_cst_Dim.rda")
# save(predict_varSD_to_var_Dim, file = "predict_varSD_to_var_Dim.rda")
# save(predict_cstSD_to_var_Dim, file = "predict_cstSD_to_var_Dim.rda")
# save(predict_varSD_to_cst_Dim, file = "predict_varSD_to_cst_Dim.rda")
# 
# save(predict_cstIT_to_cst_Dim, file = "predict_cstIT_to_cst_Dim.rda")
# save(predict_varIT_to_var_Dim, file = "predict_varIT_to_var_Dim.rda")
# save(predict_cstIT_to_var_Dim, file = "predict_cstIT_to_var_Dim.rda")
# save(predict_varIT_to_cst_Dim, file = "predict_varIT_to_cst_Dim.rda")

load(file = "predict_cstSD_to_cst_Car.rda")
load(file = "predict_varSD_to_var_Car.rda")
load(file = "predict_cstSD_to_var_Car.rda")
load(file = "predict_varSD_to_cst_Car.rda")
load(file = "predict_cstIT_to_cst_Car.rda")
load(file = "predict_varIT_to_var_Car.rda")
load(file = "predict_cstIT_to_var_Car.rda")
load(file = "predict_varIT_to_cst_Car.rda")

load(file = "predict_cstSD_to_cst_Mal.rda")
load(file = "predict_varSD_to_var_Mal.rda")
load(file = "predict_cstSD_to_var_Mal.rda")
load(file = "predict_varSD_to_cst_Mal.rda")
load(file = "predict_cstIT_to_cst_Mal.rda")
load(file = "predict_varIT_to_var_Mal.rda")
load(file = "predict_cstIT_to_var_Mal.rda")
load(file = "predict_varIT_to_cst_Mal.rda")

load(file = "predict_cstSD_to_cst_Cyp.rda")
load(file = "predict_varSD_to_var_Cyp.rda")
load(file = "predict_cstSD_to_var_Cyp.rda")
load(file = "predict_varSD_to_cst_Cyp.rda")
load(file = "predict_cstIT_to_cst_Cyp.rda")
load(file = "predict_varIT_to_var_Cyp.rda")
load(file = "predict_cstIT_to_var_Cyp.rda")
load(file = "predict_varIT_to_cst_Cyp.rda")

load(file = "predict_cstSD_to_cst_Dim.rda")
load(file = "predict_varSD_to_var_Dim.rda")
load(file = "predict_cstSD_to_var_Dim.rda")
load(file = "predict_varSD_to_cst_Dim.rda")
load(file = "predict_cstIT_to_cst_Dim.rda")
load(file = "predict_varIT_to_var_Dim.rda")
load(file = "predict_cstIT_to_var_Dim.rda")
load(file = "predict_varIT_to_cst_Dim.rda")

# --- NRMSE SD
predict_cstSD_to_cst_Car$ls_check_on_Nsurv$nrmse
predict_cstSD_to_cst_Cyp$ls_check_on_Nsurv$nrmse
predict_cstSD_to_cst_Dim$ls_check_on_Nsurv$nrmse
predict_cstSD_to_cst_Mal$ls_check_on_Nsurv$nrmse

predict_varSD_to_var_Car$ls_check_on_Nsurv$nrmse
predict_varSD_to_var_Cyp$ls_check_on_Nsurv$nrmse
predict_varSD_to_var_Dim$ls_check_on_Nsurv$nrmse
predict_varSD_to_var_Mal$ls_check_on_Nsurv$nrmse

predict_cstSD_to_var_Car$ls_check_on_Nsurv$nrmse
predict_cstSD_to_var_Cyp$ls_check_on_Nsurv$nrmse
predict_cstSD_to_var_Dim$ls_check_on_Nsurv$nrmse
predict_cstSD_to_var_Mal$ls_check_on_Nsurv$nrmse

predict_varSD_to_cst_Car$ls_check_on_Nsurv$nrmse
predict_varSD_to_cst_Cyp$ls_check_on_Nsurv$nrmse
predict_varSD_to_cst_Dim$ls_check_on_Nsurv$nrmse
predict_varSD_to_cst_Mal$ls_check_on_Nsurv$nrmse

# --- NRMSE IT
predict_cstIT_to_cst_Car$ls_check_on_Nsurv$nrmse
predict_cstIT_to_cst_Cyp$ls_check_on_Nsurv$nrmse
predict_cstIT_to_cst_Dim$ls_check_on_Nsurv$nrmse
predict_cstIT_to_cst_Malt$ls_check_on_Nsurv$nrmse

predict_varIT_to_var_Car$ls_check_on_Nsurv$nrmse
predict_varIT_to_var_Cyp$ls_check_on_Nsurv$nrmse
predict_varIT_to_var_Dim$ls_check_on_Nsurv$nrmse
predict_varIT_to_var_Mal$ls_check_on_Nsurv$nrmse

predict_cstIT_to_var_Car$ls_check_on_Nsurv$nrmse
predict_cstIT_to_var_Cyp$ls_check_on_Nsurv$nrmse
predict_cstIT_to_var_Dim$ls_check_on_Nsurv$nrmse
predict_cstIT_to_var_Mal$ls_check_on_Nsurv$nrmse

predict_varIT_to_cst_Car$ls_check_on_Nsurv$nrmse
predict_varIT_to_cst_Cyp$ls_check_on_Nsurv$nrmse
predict_varIT_to_cst_Dim$ls_check_on_Nsurv$nrmse
predict_varIT_to_cst_Mal$ls_check_on_Nsurv$nrmse

# --- PPC SD
predict_cstSD_to_cst_Car$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_cst_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_cst_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_cst_Malt$ls_check_on_Nsurv$percent_ppc_check

predict_varSD_to_var_Car$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_var_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_var_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_var_Mal$ls_check_on_Nsurv$percent_ppc_check

predict_cstSD_to_var_Car$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_var_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_var_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_cstSD_to_var_Mal$ls_check_on_Nsurv$percent_ppc_check

predict_varSD_to_cst_Car$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_cst_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_cst_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_varSD_to_cst_Mal$ls_check_on_Nsurv$percent_ppc_check

# --- PPC IT
predict_cstIT_to_cst_Car$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_cst_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_cst_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_cst_Malt$ls_check_on_Nsurv$percent_ppc_check

predict_varIT_to_var_Car$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_var_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_var_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_var_Mal$ls_check_on_Nsurv$percent_ppc_check

predict_cstIT_to_var_Car$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_var_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_var_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_cstIT_to_var_Mal$ls_check_on_Nsurv$percent_ppc_check

predict_varIT_to_cst_Car$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_cst_Cyp$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_cst_Dim$ls_check_on_Nsurv$percent_ppc_check
predict_varIT_to_cst_Mal$ls_check_on_Nsurv$percent_ppc_check

