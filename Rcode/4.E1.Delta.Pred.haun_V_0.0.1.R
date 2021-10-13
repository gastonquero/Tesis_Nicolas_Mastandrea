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
today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos
# cargar datos 

NH.estado.pred.E1 <- read.table ("./Data/procdata/NH.estado.pred.E1.testigos.txt" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA")

summary (NH.estado.pred.E1)

str (NH.estado.pred.E1)


NH.pred.E1.z21  <- NH.estado.pred.E1  %>% 
                   dplyr::filter (estado.feno =="z2.1") %>%
                   dplyr::arrange (parcela.1) %>%
                   dplyr::mutate (haun.Z21.pred = haun.pred) %>%
                   dplyr::select (-haun.pred)

dim (NH.pred.E1.z21)

NH.pred.E1.z31  <- NH.estado.pred.E1 %>%
                   dplyr::filter (estado.feno=="z3.1")%>%
                   dplyr::arrange (parcela.1) %>% 
                   dplyr::mutate (haun.Z31 = haun.pred) %>%
                   select (-haun.pred)

dim (NH.pred.E1.z31) 

NH.pred.E1.NFH  <- NH.estado.pred.E1 %>% 
                   dplyr::filter (estado.feno=="hb")%>%
                   dplyr::arrange (parcela.1) %>%
                   dplyr::mutate (haun.NFH = haun.pred) %>%
                   dplyr::select (-haun.pred)

dim (NH.pred.E1.NFH)

hojas.pred.E1.testigos <- NH.pred.E1.z21 %>% 
                 dplyr::select(-estado.feno) %>% 
                 dplyr::mutate (haun.Z31.pred = NH.pred.E1.z31$haun.Z31)%>%
                 dplyr::mutate (haun.NFH.pred = NH.pred.E1.NFH$haun.NFH) %>%
                 dplyr::mutate (haun.NFH.1.pred = haun.NFH.pred -1)    %>%
                 dplyr::mutate (deltaH.Z31.Z21.pred =  haun.Z31.pred - haun.Z21.pred) %>%
                 dplyr::mutate (deltaH.NFH.Z31.pred =  haun.NFH.pred - haun.Z31.pred)   %>%
                 dplyr::mutate (deltaH.NFH.1.Z31.pred =  haun.NFH.1.pred - haun.Z31.pred)

write.table (hojas.pred.E1.testigos, 
             file = "./Data/procdata/delta.NH.pred.estado.E1.testigos.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)


############################ aca termina el codigo para el delta del    ####
############################  numero de hojas predicho en cada subfase ###########
