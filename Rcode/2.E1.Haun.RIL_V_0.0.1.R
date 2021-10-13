############################################################################
# Datos  Cebada - ceiboxcarumbe                                           ##
# Datos tomados por Nicolas Mastandrea                                    ##
# Gaston Quero  - Nicolas Mastandrea                                      ##
# 9-04-2018                                                             ##
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


#########  se carga los datos crudos
# cargar datos 
#### ambiente 1 
E1_RIL_haun <- read.table ("./Data/rawdata/E1_RIL_haun.txt" ,
                           header = TRUE, sep = "\t",dec = ".",
                           na.strings = "NA" )

summary (E1_RIL_haun)
str (E1_RIL_haun)

E1_hojayalt <- read.table ("./Data/rawdata/E1_hojayalt.txt" ,
                           header = TRUE, sep = "\t",dec = ".",
                           na.strings = "NA" )

str (E1_hojayalt)

# funcion para calcular el numero de hoja 
# en cada momento y los deltas de hojas

####
E1_RIL_NH <- NH (E1_RIL_haun, E1_hojayalt) 
                                       

### vuelvo a ingresar los datos

#E1_RIL_NH <- read.table ("./Data/procdata/outpout.HAUN/datos.NHE1.csv" ,
 #                        header = TRUE, sep = ",",dec = ".",
  #                       na.strings = "NA" )


E1_RIL_NH <- E1_RIL_NH %>%
             dplyr::mutate (cod.pl = str_c (E1_RIL_NH$planta,E1_RIL_NH$pot))


E1_RIL_NH.z21  <- E1_RIL_NH  %>% 
                  dplyr::filter  (estado.feno =="z2.1") %>%
                  dplyr::arrange (bloque) %>% 
                  dplyr::arrange (parcela) %>%
                  dplyr::mutate (haun.Z21 = haun) %>%
                  dplyr::select (-haun)%>%
                  dplyr::arrange (pot)

dim(E1_RIL_NH.z21)

E1_RIL_NH.z31  <- E1_RIL_NH  %>% 
                  dplyr::filter  (estado.feno=="z3.1") %>%
                  dplyr::arrange (bloque) %>% 
                  dplyr::arrange (parcela) %>%
                  dplyr::mutate (haun.Z31 = haun) %>%
                  dplyr::select (-haun)%>%
                  dplyr::arrange (pot)


dim (E1_RIL_NH.z31)

E1_RIL_NH.NFH  <- E1_RIL_NH  %>% 
                  dplyr::filter  (estado.feno=="flag.leaf") %>%
                  dplyr::arrange (bloque) %>% 
                  dplyr::arrange (parcela) %>%
                  dplyr::mutate (haun.NFH = haun) %>%
                  dplyr::select (-haun) %>%
                  dplyr::arrange (pot)

dim (E1_RIL_NH.NFH)

hojas.RIL.E1 <- E1_RIL_NH.z21 %>% 
                dplyr::select(-estado.feno) %>% 
                dplyr::mutate (haun.Z31 = E1_RIL_NH.z31$haun.Z31)%>%
                dplyr::mutate (haun.NFH = E1_RIL_NH.NFH$haun.NFH) %>%
                dplyr::mutate (haun.NFH.1 = haun.NFH -1)    %>%
                dplyr::mutate (deltaH.Z31.Z21 =  haun.Z31 - haun.Z21) %>%
                dplyr::mutate (deltaH.NFH.Z31 =  haun.NFH - haun.Z31)   %>%
                dplyr::mutate (deltaH.NFH.1.Z31 =  haun.NFH.1 - haun.Z31)

write.table (hojas.RIL.E1, file = "./Data/procdata/outpout.HAUN/delta.NH.RIL.E1.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)

# aca generamos el data.frame de numero de hojas en cada subfase
# y la diferencias entre entre el numero de hoja de una fase menos el
# numero de hojas de la fase anterior

