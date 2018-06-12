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
  geom_line(aes(x = time, y = q50, color = "orange"))

FOCUS.SW$conc <- FOCUS.SW$conc * 1000
pred_cstDiaz_FocusSW_1000 <- predict(fit_diaz_cstSD, data_predict = FOCUS.SW)
plot(pred_cstDiaz_FocusSW_10)

# --- deSolve function
predODE_cstDiaz_FocusSW <- predict_ode.survFit(fit_diaz_cstSD, data_predict = FOCUS.SW)
plot(predODE_cstDiaz_FocusSW)

