library(morse)
library(dplyr)
library(tidyr)

Data.Set.C_Cal.Conc_modif <- read.csv("~/Documents/PostDoc-LBBE/GUTS_RING_TEST/Data-Set_GUTS-ring-test/DATA-SET_C/Data-Set-C_Cal-Conc_modif.csv") %>%
  gather(replicate, conc, -time)
Data.Set.C_Survival <- read.csv("~/Documents/PostDoc-LBBE/GUTS_RING_TEST/Data-Set_GUTS-ring-test/DATA-SET_C/Data-Set-C_Survival.csv") %>%
  mutate(time = Day) %>%
  select(- Day) %>%
  gather(replicate, Nsurv,-time)

data_diazinon_cst <- full_join(Data.Set.C_Cal.Conc_modif, Data.Set.C_Survival, by = c("replicate", "time"))
fit_diaz_cstSD <- survFit(survData(data_diazinon_cst), model_type = "SD")

plot(fit_diaz_cstSD)

FOCUS.SW <- read.csv("~/Documents/GUTS_shinyapp/FOCUS-SW.csv")
FOCUS.SW$replicate <- rep("A", nrow(FOCUS.SW))
pred_cstDiaz_FocusSW <- predict(fit_diaz_cstSD, data_predict = FOCUS.SW)
plot(pred_cstDiaz_FocusSW)

library(ggplot2)
ggplot(data = pred_cstDiaz_FocusSW$df_quantile )+
  geom_ribbon(aes(x = time, ymin = qinf95, ymax = qsup95)) +
  geom_line(aes(x = time, y = q50, color = "orange"))

# FOCUS.SW$conc <- FOCUS.SW$conc * 1000
# pred_cstDiaz_FocusSW_1000 <- predict(fit_diaz_cstSD, data_predict = FOCUS.SW)
# plot(pred_cstDiaz_FocusSW_10)

# --- deSolve function
library(deSolve)
predODE_cstDiaz_FocusSW <- predict_ode(fit_diaz_cstSD, data_predict = FOCUS.SW, mcmc_size = 1000)
plot(predODE_cstDiaz_FocusSW)

# ---
setwd("~/Documents/morse_test")
load(file = "data/data_PRZ_cst.rda")
load(file = "data/data_PRZ_var.rda")
load(file = "fit/survFit_cstSD_PRZ.rda")
load(file = "fit/survFit_varSD_PRZ.rda")
load(file = "fit/survFit_cstIT_PRZ.rda")
load(file = "fit/survFit_varIT_PRZ.rda")
fit_cstIT <- survFit_cstIT_PRZ
fit_varIT <- survFit_varIT_PRZ
fit_cstSD <- survFit_cstSD_PRZ
fit_varSD <- survFit_varSD_PRZ

predODE_SD_cstTOcst <- predict(survFit_cstSD_PRZ, data_predict = PRZ_cst)
predODE_SD_cstTOvar <- predict(survFit_cstSD_PRZ, data_predict = PRZ_var)
predODE_SD_varTOcst <- predict(survFit_varSD_PRZ, data_predict = PRZ_cst)
predODE_SD_varTOvar <- predict(survFit_varSD_PRZ, data_predict = PRZ_var)
save(predODE_SD_cstTOcst, file = "predict/predODE_SD_cstTOcst.rda")
save(predODE_SD_cstTOvar, file = "predict/predODE_SD_cstTOvar.rda")
save(predODE_SD_varTOcst, file = "predict/predODE_SD_varTOcst.rda")
save(predODE_SD_varTOvar, file = "predict/predODE_SD_varTOvar.rda")

predODE_SD_cstTOcst_ode <- predict_ode(survFit_cstSD_PRZ, data_predict = PRZ_cst, mcmc_size = 1000)
predODE_SD_cstTOvar_ode <- predict_ode(survFit_cstSD_PRZ, data_predict = PRZ_var, mcmc_size = 1000)
predODE_SD_varTOcst_ode <- predict_ode(survFit_varSD_PRZ, data_predict = PRZ_cst, mcmc_size = 1000)
predODE_SD_varTOvar_ode <- predict_ode(survFit_varSD_PRZ, data_predict = PRZ_var, mcmc_size = 1000)
save(predODE_SD_cstTOcst_ode, file = "predict/predODE_SD_cstTOcst_ode.rda")
save(predODE_SD_cstTOvar_ode, file = "predict/predODE_SD_cstTOvar_ode.rda")
save(predODE_SD_varTOcst_ode, file = "predict/predODE_SD_varTOcst_ode.rda")
save(predODE_SD_varTOvar_ode, file = "predict/predODE_SD_varTOvar_ode.rda")

plot(predODE_SD_cstTOcst_ode)
plot(predODE_SD_cstTOcst)
plot(predODE_SD_varTOcst_ode)
plot(predODE_SD_varTOcst)
plot(predODE_SD_cstTOvar_ode)
plot(predODE_SD_cstTOvar)
plot(predODE_SD_varTOvar_ode)
plot(predODE_SD_varTOvar)


predODE_IT_cstTOcst <- predict(survFit_cstIT_PRZ, data_predict = PRZ_cst)
predODE_IT_cstTOvar <- predict(survFit_cstIT_PRZ, data_predict = PRZ_var)
predODE_IT_varTOcst <- predict(survFit_varIT_PRZ, data_predict = PRZ_cst)
predODE_IT_varTOvar <- predict(survFit_varIT_PRZ, data_predict = PRZ_var)
save(predODE_IT_cstTOcst, file = "predict/predODE_IT_cstTOcst.rda")
save(predODE_IT_cstTOvar, file = "predict/predODE_IT_cstTOvar.rda")
save(predODE_IT_varTOcst, file = "predict/predODE_IT_varTOcst.rda")
save(predODE_IT_varTOvar, file = "predict/predODE_IT_varTOvar.rda")
predODE_IT_cstTOcst_ode <- predict_ode(survFit_cstIT_PRZ, data_predict = PRZ_cst, mcmc_size = 1000)
predODE_IT_cstTOvar_ode <- predict_ode(survFit_cstIT_PRZ, data_predict = PRZ_var, mcmc_size = 1000)
predODE_IT_varTOcst_ode <- predict_ode(survFit_varIT_PRZ, data_predict = PRZ_cst, mcmc_size = 1000)
predODE_IT_varTOvar_ode <- predict_ode(survFit_varIT_PRZ, data_predict = PRZ_var, mcmc_size = 1000)
save(predODE_IT_cstTOcst_ode, file = "predict/predODE_IT_cstTOcst_ode.rda")
save(predODE_IT_cstTOvar_ode, file = "predict/predODE_IT_cstTOvar_ode.rda")
save(predODE_IT_varTOcst_ode, file = "predict/predODE_IT_varTOcst_ode.rda")
save(predODE_IT_varTOvar_ode, file = "predict/predODE_IT_varTOvar_ode.rda")

plot(predODE_IT_cstTOcst_ode)
plot(predODE_IT_cstTOcst)

plot(predODE_IT_varTOcst_ode)
plot(predODE_IT_varTOcst)

plot(predODE_IT_cstTOvar_ode)
plot(predODE_IT_cstTOvar)

plot(predODE_IT_varTOvar_ode)
plot(predODE_IT_varTOvar)


