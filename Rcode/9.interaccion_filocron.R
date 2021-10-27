############################################################################
# Datos  Cebada - ceiboxcarumbe                                           ##
# Datos tomados por Nicolas Mastandrea                                    ##
# Gaston Quero  - Nicolas Mastandrea                                      ##
# 13-03-2019                                                            ##
############################################################################

getwd()
setwd ("E:/tesis_Copia_Nueva_Mastandrea")

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
library(multcomp)
library(forcats)
library(ggpubr)

today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos

filocron.E1a <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE1ciclo.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


filocron.E2a <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE2ciclo.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


filocron.E3a <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE3ciclo.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


filocron.E4a <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE4ciclo.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


### Los cuatro ambientes #####

filocron.env <-  bind_rows (filocron.E1a, 
                            filocron.E2a,
                            filocron.E3a, 
                            filocron.E4a )

unique (filocron.env$ambiente )

write.table (filocron.env, file = "./Data/procdata/filocron.env.RIL.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)

testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 

filocron.testigos <- filocron.env  %>% 
                     dplyr::filter (genotipo %in% tes.geno)

filocron.testigos <- filocron.testigos %>%
                     dplyr::mutate (env = ambiente) %>%
                     dplyr::filter (filocron < 150) %>%
                     dplyr::filter (genotipo != "Berolina")




list.geno <- c("Kenia","Prior","Logan",  
               "Baronesse",
               "Bowman","Quebracho",   
               "Danuta", "Carumbe","Ceibo") 


head (filocron.testigos)

filocron.testigos$bloque <- as.factor (filocron.testigos$bloque)



##### filocron #####
filocron.mod.1 <- lmer (filocron ~ (1|bloque) + genotipo + env + genotipo * env , data = filocron.testigos)

anova (filocron.mod.1)

### funcion que corre los contrastes 
run_contrastes <- function (data.model = NULL, trait = NULL){
  
  dt.c <- bind_rows (lapply( list.geno, function (filt.geno) {
    
    #filt.geno ="Prior"
    print (str_c (trait,"_",filt.geno))
    em.geno <- emmeans (data.model, "env",
                        at = list (genotipo = filt.geno))
    
    cr.geno <- contrast (em.geno ,
                         method = "pairwise")
    
    cr.geno.1 <- as.data.frame (cr.geno )
    
    contrastes <- tibble( genotipo = filt.geno, trait =trait, cr.geno.1)
    
    write_csv2 (contrastes, file= str_c ("./Data/procdata/contrastes","_", trait,"_", filt.geno, ".csv"))
    
    em.geno <- cbind (genotipo = filt.geno, trait =trait,
                      cld (em.geno, sort=FALSE))
    
    
  }))
  
  
}

contrastes_filocron <- run_contrastes (data.model = filocron.mod.1, trait = "filocron" )


write_csv2 (contrastes_filocron , file= "./Data/procdata/contrastes_filocron.csv")

