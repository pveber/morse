data("cadmium1")
cadmium1 
cadmium1 <- cadmium1 %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(cadmium1, file = "data/cadmium1_2.rda")
cadmium1


data("cadmium2")
cadmium2
cadmium2 <- cadmium2 %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(cadmium2, file = "data/cadmium2_2.rda")
cadmium2


data("chlordan")
chlordan
chlordan <- chlordan %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(chlordan, file = "data/chlordan_2.rda")
chlordan

data("copper")
copper
copper <- copper %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(copper, file = "data/copper_2.rda")
copper

data("dichromate")
dichromate
dichromate <- dichromate %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(dichromate, file = "data/dichromate_2.rda")
dichromate

data("propiconazole")
propiconazole
propiconazole <- propiconazole %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(propiconazole, file = "data/propiconazole_2.rda")
propiconazole

data("zinc")
zinc
zinc <- zinc %>%
  arrange(conc,replicate,time) %>%
  rename(replicate_ID = replicate) %>%
  mutate(replicate = as.numeric(interaction(replicate_ID,conc)))
save(zinc, file = "data/zinc_2.rda")
zinc
