############################################################################
# Datos  Cebada - ceiboxcarumbe                                           ##
# Datos tomados por Nicolas Mastandrea                                    ##
# Gaston Quero  - Nicolas Mastandrea                                      ##
# marzo 2019                                                            ##
############################################################################

getwd()
setwd ("R:/Tesis_Nicolas_Mastandrea")

# Paquetes 
library(lme4)
library(emmeans)
library("car")
library("nlmrt")
library("easynls")
library("plotrix")  
library("lattice")
library("latticeExtra")
library(multcompView)
library(effects)
library("dplyr")
library("corrplot")
library ("ggplot2")
library ("FactoMineR")
library (ggjoy)
#library (hrbrthemes)
library(tidyverse)
#library(forcats)
library("viridis")
library("lmerTest")
library("lubridate")
library("stringr")

library (forcats)
today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos
# cargar datos 

E1_ciclo <- read.table ("./Data/procdata/ciclo.gd.RIL.E1.txt" ,
                        header = TRUE, sep = ",",dec = ".",
                        na.strings = "NA" )
#View(E1_ciclo)
#E1_ciclo.1 <- E1_ciclo %>%
#             rename (env = ambiente)

E1_ciclo.1 <- E1_ciclo %>%
              dplyr::mutate (env = ambiente ) %>%
              dplyr::select (-ambiente)

#View(E1_ciclo.1.1)

E1_ciclo.1.1 <- E1_ciclo.1 %>%
                dplyr::mutate (pot =str_c (E1_ciclo.1$env,E1_ciclo.1$bloque,E1_ciclo.1$parcela)) %>%
                dplyr::mutate (ambiente = "E1")


ciclo.gd.E1 <- E1_ciclo.1.1 %>%
               dplyr::select (genotipo, anio,ambiente, env, pot, bloque,parcela,e.z21.gd, 
                               z21.z31.gd,z31.a.gd,e.a.gd) %>% 
               dplyr::arrange (pot)

dim(ciclo.gd.E1)

hist(x = ciclo.gd.E1$e.a.gd)
hist(x = ciclo.gd.E1$e.a.gd, main = "Ciclo E1", 
     xlab = "GD", ylab = "Frecuencia", col = "grey")



####  E2 #############

E2_ciclo <- read.table ("./Data/procdata/outpout.HAUN/ciclo.gd.RIL.E2.txt" ,
                        header = TRUE, sep = ",",dec = ".",
                        na.strings = "NA" )


#View(E2_ciclo)
#E2_ciclo.1 <- E2_ciclo %>%
#             rename (env = ambiente)


E2_ciclo.1 <- E2_ciclo %>%
  mutate (env = ambiente ) %>%
  select (-ambiente)
#View(E2_ciclo.1.1)



E2_ciclo.1.1 <- E2_ciclo.1 %>%
  mutate (pot =str_c (E2_ciclo.1$env,E2_ciclo.1$bloque,E2_ciclo.1$parcela)) %>%
  mutate (ambiente = "E2")


ciclo.gd.E2 <- E2_ciclo.1.1 %>%
  select(genotipo, anio,ambiente, env, pot, bloque,parcela,e.z21.gd, 
         z21.z31.gd,z31.a.gd,e.a.gd) %>% 
  arrange (pot)


dim(ciclo.gd.E2)

#### ESTO LO PRECISAMOS?? ANTES LO USABAMOS APARENTEMENTE #####
pot.E2 <- as.character (unique(NFH.completo.E2.1$pot))

ciclo.gd.E2.1 <- ciclo.gd.E2 %>%
  filter  (pot %in% pot.E2)

dim (ciclo.gd.E2.1)


hist(x = ciclo.gd.E2$e.a.gd)
hist(x = ciclo.gd.E2$e.a.gd, main = "Ciclo E2", 
     xlab = "GD", ylab = "Frecuencia", col = "grey")


#### E3 ###################


E3_ciclo <- read.table ("./Data/procdata/outpout.HAUN/ciclo.gd.RIL.E3.txt" ,
                        header = TRUE, sep = ",",dec = ".",
                        na.strings = "NA" )

#View(E3_ciclo)
#E3_ciclo.1 <- E3_ciclo %>%
#             rename (env = ambiente)

E3_ciclo.1 <- E3_ciclo %>%
  select (-ambiente) %>% 
  mutate (env=3)


#E3_ciclo.1 <- E3_ciclo %>%
#  mutate (env = ambiente ) %>%
#  select (-ambiente)
#View(E3_ciclo.1.1)


E3_ciclo.1.1 <- E3_ciclo.1 %>%
  mutate (pot =str_c (E3_ciclo.1$env,E3_ciclo.1$bloque,E3_ciclo.1$parcela)) %>%
  mutate (ambiente = "E3")


ciclo.gd.E3 <- E3_ciclo.1.1 %>%
  select(genotipo, anio,ambiente, env, pot, bloque,parcela,e.z21.gd, 
         z21.z31.gd,z31.a.gd,e.a.gd) %>% 
  arrange (pot)

dim(ciclo.gd.E3)

#### ESTO LO PRECISAMOS?? ANTES LO USABAMOS APARENTEMENTE #####

pot.E3 <- as.character (unique(NFH.completo.E3$pot))

ciclo.gd.E3$pot

ciclo.gd.E3.1 <- ciclo.gd.E3 %>%
  dplyr::filter (pot%in%pot.E3)

hist(x = ciclo.gd.E3$e.a.gd)
hist(x = ciclo.gd.E3$e.a.gd, main = "Ciclo E3", 
     xlab = "GD", ylab = "Frecuencia", col = "grey")

#### E4 ##############

E4_ciclo <- read.table ("./Data/procdata/outpout.HAUN/ciclo.gd.RIL.E4.txt" ,
                        header = TRUE, sep = ",",dec = ".",
                        na.strings = "NA" )

#View(E4_ciclo)
#E4_ciclo.1 <- E4_ciclo %>%
#             rename (env = ambiente)

E4_ciclo.1 <- E4_ciclo %>%
  select (-ambiente) %>% 
  mutate (env=4)


#E4_ciclo.1 <- E4_ciclo %>%
#  mutate (env = ambiente ) %>%
#  select (-ambiente)
#View(E4_ciclo.1.1)


E4_ciclo.1.1 <- E4_ciclo.1 %>%
  mutate (pot =str_c (E4_ciclo.1$env,E4_ciclo.1$bloque,E4_ciclo.1$parcela)) %>%
  mutate (ambiente = "E4")


ciclo.gd.E4 <- E4_ciclo.1.1 %>%
  select(genotipo, anio,ambiente, env, pot, bloque,parcela,e.z21.gd, 
         z21.z31.gd,z31.a.gd,e.a.gd) %>% 
  arrange (pot)

dim(ciclo.gd.E4)

#### ESTO LO PRECISAMOS?? ANTES LO USABAMOS APARENTEMENTE #####

pot.E4 <- as.character (unique(NFH.completo.E4$pot))

ciclo.gd.E4$pot

ciclo.gd.E4.1 <- ciclo.gd.E4 %>%
  dplyr::filter (pot%in%pot.E4)


dim(ciclo.gd.E4.1)

hist(x = ciclo.gd.E4$e.a.gd)
hist(x = ciclo.gd.E4$e.a.gd, main = "Ciclo E4", 
     xlab = "GD", ylab = "Frecuencia", col = "grey")
###########################################

gd.E1  <- ciclo.gd.E1 %>%
  mutate (year = 2016)%>%
  mutate (SD = "early")%>%
  mutate (temp =12.28)

dim (NFH.completo.E1)

gd.E2 <- ciclo.gd.E2 %>%
  mutate (year = 2016)%>%
  mutate (SD = "late")%>%
  mutate (temp =17.44)

gd.E3 <- ciclo.gd.E3 %>%
  mutate (year = 2017)%>%
  mutate (SD = "early")%>%
  mutate (temp =14.66)

gd.E4 <-  ciclo.gd.E4 %>%
  mutate (year = 2017)%>%
  mutate (SD = "late")%>%
  mutate (temp =19.23)

gd.env <- rbind (gd.E1,
                 gd.E2,
                 gd.E3,
                 gd.E4)

gd.env$pot <- as.character(gd.env$pot)

gd.env <- gd.env %>%
  arrange (pot)

pot.NFH.env.1 <- as.character(unique (NFH.env.1$pot))




gd.env.1 <- gd.env %>%
  filter (pot %in% pot.NFH.env.1) 

gd.env.1 <- gd.env.1 %>%
  arrange (pot)
dim(gd.env.1)

NFH.env.1$pot <- as.character(NFH.env.1$pot)
NFH.env.1 <- NFH.env.1 %>%
  arrange (pot) 

dim(NFH.env.1)
NFH.env.1$pot == gd.env.1$pot 

head(gd.env.1)
summary(gd.env.1)

gd.env.1 %>%
  mutate (YearFct = fct_rev(as.factor(year))) %>%
  ggplot(aes(y = YearFct)) +
  geom_density_ridges(
    aes(x = e.a.gd, fill = paste(YearFct, SD)), 
    alpha = .8, color = "white", from = 500, to = 1700
  ) +
  labs(
    x = "°Cd",
    y = "Year"
    #title = "Indy vs Unionist vote in Catalan elections",
    #subtitle = "Analysis unit: municipalities (n = 949)",
    #caption = "Marc Belzunces (@marcbeldata) | Source: Idescat"
  ) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_fill_cyclical(
    breaks = c("2016 early", "2016 late"),
    labels = c(`2016 early` = "early", `2016 late` = "late"),
    values = c("#ff0000", "#0000ff"),
    name = "SD", guide = "legend"
  ) +
  theme_ridges(grid = FALSE)




write.table (gd.env.1, file = "./Data/procdata/gd.env.1.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)


######## #################### aca termina el codigo para duracion de ciclo ###########

