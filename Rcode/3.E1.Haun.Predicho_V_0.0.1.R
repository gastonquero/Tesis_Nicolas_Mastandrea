###############################################################
# predicho de NH segun los filocrones calculados             #
# Gaston Quero - Nicolas Mastandrea                          #
# 21-06-2018                                                 #
############################################################# 

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
####aca podríamos cargar desde el outpout pero no tenemos s-e ###
### entonces usamos desde el rawdata como estaba antes #####

E1_ciclo <- read.table ("./Data/rawdata/E1_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )
summary(E1_ciclo)
str(E1_ciclo)

testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 

E1_ciclo.testigos <- E1_ciclo %>% 
                     dplyr::filter (genotipo %in% tes.geno)%>%
                     dplyr::arrange (parcela)
str(E1_ciclo.testigos)
head(E1_ciclo.testigos)

E1_ciclo.testigos <- E1_ciclo.testigos %>%
                     dplyr::mutate (s.z21.gd = s.e.gd + e.z21.gd)%>%
                     dplyr::mutate (s.z31.gd = s.e.gd + e.z21.gd + z21.z31.gd) %>%
                     dplyr::mutate (s.a.gd = s.e.gd + e.z21.gd + z21.z31.gd + z31.a.gd)
 
E1_ciclo.testigos$parcela


#E1_ciclo.testigos.1 <- E1_ciclo.testigos %>% filter(parcela !=4)%>% 
                       #filter(parcela !=70)

#dim (E1_ciclo.testigos.1)

######### corremos el model GD ########
## para el ciclo completo

E1_RIL_haun <- read.table ("./Data/rawdata/E1_RIL_haun.txt" ,
                           header = TRUE, sep = "\t",dec = ".",
                           na.strings = "NA" )

## asegurse de correr antes las funcion pre.fil

E1_RIL_datos.nfh.ciclo <- pre.fil (E1_RIL_haun, subfase = "ciclo")

### Vamos a calcular el filocron 

filo.ciclo.E1 <- GD.model.3 (E1_RIL_datos.nfh.ciclo, subfase = "ciclo" )


filocron.E1a.testigos <- filo.ciclo.E1 %>% 
                         dplyr::filter (genotipo %in% tes.geno)%>%
                         dplyr::arrange (parcela) %>% 
                         dplyr::mutate (parcela.1 = parcela) %>%
                         dplyr::select (-parcela)

unique(filocron.E1a.testigos$subfase)
str(filocron.E1a.testigos)

filocron.E1b.testigos <-  filocron.E1a.testigos %>% 
                         dplyr::arrange (parcela.1) 
                         #%>%dplyr::filter (subfase == "ciclo")


filocron.E1b.testigos$parcela.1 %in% E1_ciclo.testigos$parcela
E1_ciclo.testigos$parcela %in% filocron.E1b.testigos$parcela.1
filocron.E1b.testigos$parcela.1  == E1_ciclo.testigos$parcela



filocron.E1.pred.testigos <- filocron.E1b.testigos %>% 
                             dplyr::mutate (s.z21.gd = E1_ciclo.testigos$s.z21.gd)%>% 
                             dplyr::mutate (s.z31.gd = E1_ciclo.testigos$s.z31.gd)%>%
                             dplyr::mutate (s.a.gd = E1_ciclo.testigos$s.a.gd)

NH.estado.pred.E1.testigos <- filocron.E1.pred.testigos %>% 
                              dplyr::mutate (haun.z21.pred = a + b*s.z21.gd)%>%
                              dplyr::mutate (haun.z31.pred = a + b*s.z31.gd)%>%
                              dplyr::mutate (haun.NFH.pred = a + b*s.a.gd)
                      
    

########## hasta aca calculamos NH predicho en cada fase #########

NH.E1.z21.testigos  <- NH.estado.pred.E1.testigos %>% 
                       dplyr::arrange (parcela.1) %>%
                       dplyr::mutate (haun.pred = NH.estado.pred.E1.testigos$haun.z21.pred) %>%
                       dplyr::select (-haun.z21.pred) %>%
                       dplyr::mutate (estado.feno ="z2.1") %>%
                       dplyr::select (genotipo, id.diseno, bloque, parcela.1,
                             estado.feno,haun.pred)


NH.E1.z31.testigos  <- NH.estado.pred.E1.testigos %>% 
                       dplyr::arrange (parcela.1) %>%
                       dplyr::mutate (haun.pred = haun.z31.pred)%>%
                       dplyr::select (-haun.z31.pred) %>%
                       dplyr::mutate (estado.feno ="z3.1") %>%
                       dplyr::select (genotipo, id.diseno, bloque, parcela.1,
                            estado.feno,haun.pred)


NH.E1.NFH.testigos  <- NH.estado.pred.E1.testigos %>% 
                       dplyr::arrange (parcela.1) %>%
                       dplyr::mutate (haun.pred = haun.NFH.pred) %>%
                       dplyr::select (-haun.NFH.pred)%>%
                       dplyr::mutate (estado.feno ="hb") %>%
                       dplyr::select (genotipo, id.diseno, bloque, parcela.1,
                            estado.feno,haun.pred)


hojas.E1.testigos.pred  <- rbind (NH.E1.z21.testigos, 
                                  NH.E1.z31.testigos, 
                                  NH.E1.NFH.testigos)

write.table (hojas.E1.testigos.pred, 
             file = "./Data/procdata/NH.estado.pred.E1.testigos.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)               
                      
                      
####### saca el numero de hojas predicho en cada fase #
