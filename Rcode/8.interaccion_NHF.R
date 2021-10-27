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
today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos

E1_hoja <- read.table ("./Data/procdata/outpout.HAUN/delta.NH.RIL.E1.txt" ,
                           header = TRUE, sep = ",",dec = ".",
                           na.strings = "NA" )

E2_hoja <- read.table ("./Data/procdata/outpout.HAUN/delta.NH.RIL.E2.txt" ,
                           header = TRUE, sep = ",",dec = ".",
                           na.strings = "NA" )

E3_hoja <- read.table ("./Data/procdata/outpout.HAUN/delta.NH.RIL.E3.txt" ,
                           header = TRUE, sep = ",",dec = ".",
                           na.strings = "NA" )

E4_hoja <- read.table ("./Data/procdata/outpout.HAUN/delta.NH.RIL.E4.txt",
                           header = TRUE, sep = ",",dec = ".",
                           na.strings = "NA" )


### Los cuatro ambientes #####

hoja.env <-  bind_rows (E1_hoja , 
                        E2_hoja  ,
                        E3_hoja , 
                        E4_hoja  )


write.table (hoja.env, file = "./Data/procdata/hoja.env.RIL.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)

testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 

hoja.testigos <- hoja.env %>% 
                 dplyr::filter (genotipo %in% tes.geno)

hoja.testigos <-hoja.testigos %>%
                dplyr::rename (env = ambiente)

list.geno <- c("Kenia","Prior","Logan",  
               "Baronesse",
               "Bowman","Quebracho",   
               "Danuta", "Carumbe","Ceibo") 


head (hoja.testigos)

hoja.testigos$bloque <- as.factor (hoja.testigos$bloque)

head (hoja.testigos )

##### e.a.gd #####
haun.NFH.mod.1 <- lmer (haun.NFH ~ (1|bloque) + genotipo + env + genotipo * env , data = hoja.testigos)

anova(haun.NFH.mod.1)

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


contrastes_haun.NFH <- run_contrastes (data.model = haun.NFH.mod.1, trait = "haun.NFH" )


write_csv2 (contrastes_haun.NFH , file= "./Data/procdata/contrastes_hojas.csv")

###### haun.Z31 #####
haun.Z31.mod.1 <- lmer (haun.Z31 ~ (1|bloque) + genotipo + env + genotipo * env , data = hoja.testigos)

anova (haun.Z31.mod.1)

contrastes_haun.Z31 <- run_contrastes (data.model = haun.Z31.mod.1, trait = "haun.Z31" )


write_csv2 (contrastes_haun.Z31 , file= "./Data/procdata/contrastes_haun.Z31.csv")


###### deltaH.NFH.Z31 #####
deltaH.NFH.Z31.mod.1 <- lmer (deltaH.NFH.Z31 ~ (1|bloque) + genotipo + env + genotipo * env , data = hoja.testigos)

anova (deltaH.NFH.Z31.mod.1)

contrastes_deltaH.NFH.Z31 <- run_contrastes (data.model = deltaH.NFH.Z31.mod.1, trait = "deltaH.NFH.Z31" )


write_csv2 (contrastes_deltaH.NFH.Z31 , file= "./Data/procdata/contrastes_deltaH.NFH.Z31.csv")





## obtener una media ajusta por modelo fijo por ambiente  ###
str (E1_hojayalt)


#### NFH ajustado E1

E1_hojayalt$pot <- as.factor (E1_hojayalt$pot)

nfh.E1.mod.1 <- lm (haun.NFH ~ pot, data=E1_hojayalt)

em.nfh.E1.mod.1 <- emmeans (nfh.E1.mod.1, ~ pot)

em.nfh.E1.mod1 <- summary (em.nfh.E1.mod.1)



em.nfh.E1.mod1 <- em.nfh.E1.mod1 %>%
  arrange (pot)

em.nfh.E1.mod1 <- em.nfh.E1.mod1 %>%
  mutate (em.nfh = emmean) %>%
  select (-emmean)


View(em.nfh.E1.mod1)

hist(x = em.nfh.E1.mod1$em.nfh)
hist(x = em.nfh.E1.mod1$em.nfh, main = "NFH E1", 
     xlab = "leaf number", ylab = "Frecuencia", col = "grey")

#### NFH ajustado E2

E2_hojayalt$pot <- as.factor (E2_hojayalt$pot)

nfh.E2.mod.1 <- lm (haun.NFH ~ pot, data=E2_hojayalt)

em.nfh.E2.mod.1 <- emmeans (nfh.E2.mod.1, ~ pot)

em.nfh.E2.mod1 <- summary (em.nfh.E2.mod.1)


em.nfh.E2.mod1 <- em.nfh.E2.mod1 %>%
  arrange (pot)

em.nfh.E2.mod1 <- em.nfh.E2.mod1 %>%
  mutate (em.nfh = emmean) %>%
  select (-emmean)

View(em.nfh.E2.mod1)

hist(x = em.nfh.E2.mod1$em.nfh)
hist(x = em.nfh.E2.mod1$em.nfh, main = "NFH E2", 
     xlab = "leaf number", ylab = "Frecuencia", col = "grey")

#### NFH ajustado E3

E3_hojayalt$pot <- as.factor (E3_hojayalt$pot)

nfh.E3.mod.1 <- lm (haun.NFH ~ pot, data=E3_hojayalt)

em.nfh.E3.mod.1 <- emmeans (nfh.E3.mod.1, ~ pot)

em.nfh.E3.mod1 <- summary (em.nfh.E3.mod.1)

em.nfh.E3.mod1 <- em.nfh.E3.mod1 %>%
  arrange (pot)

em.nfh.E3.mod1 <- em.nfh.E3.mod1 %>%
  mutate (em.nfh = emmean) %>%
  select (-emmean)

View(em.nfh.E3.mod1)

hist(x = em.nfh.E3.mod1$em.nfh)
hist(x = em.nfh.E3.mod1$em.nfh, main = "NFH E3", 
     xlab = "leaf number", ylab = "Frecuencia", col = "grey")

#### NFH ajustado E4

E4_hojayalt$pot <- as.factor (E4_hojayalt$pot)

nfh.E4.mod.1 <- lm (haun.NFH ~ pot, data=E4_hojayalt)

em.nfh.E4.mod.1 <- emmeans (nfh.E4.mod.1, ~ pot)

em.nfh.E4.mod1 <- summary (em.nfh.E4.mod.1)

em.nfh.E4.mod1 <- em.nfh.E4.mod1 %>%
  arrange (pot)

em.nfh.E4.mod1 <- em.nfh.E4.mod1 %>%
  mutate (em.nfh = emmean) %>%
  select (-emmean)

View(em.nfh.E4.mod1)

hist(x = em.nfh.E4.mod1$em.nfh)
hist(x = em.nfh.E4.mod1$em.nfh, main = "NFH E4", 
     xlab = "leaf number", ylab = "Frecuencia", col = "grey")

#### 
str(E1_hojayalt)

hist (em.nfh.E1.mod1$em.nfh)
hist (em.nfh.E2.mod1$em.nfh)
hist (em.nfh.E3.mod1$em.nfh)
hist (em.nfh.E4.mod1$em.nfh)

########

class(em.nfh.E1.mod1)
dim(em.nfh.E1.mod1)

NFH.completo.E1  <- em.nfh.E1.mod1 %>%
  mutate (year = 2016)%>%
  mutate (SD = "early")%>%
  mutate (temp =12.28)

dim (NFH.completo.E1)

NFH.completo.E2 <- em.nfh.E2.mod1 %>%
  mutate (year = 2016)%>%
  mutate (SD = "late")%>%
  mutate (temp =17.44)

#NFH.completo.E2.1 <- NFH.completo.E2 %>%
#  filter (pot!=2114)%>%
#  filter (pot!=2115)%>%
#  filter (pot!=2117)%>%
#  filter (pot!=2118)%>%
#  filter (pot!=2121)%>%
#  filter (pot!=2133)%>%
#  filter (pot!=2255)%>%
#  filter (pot!=2266)%>%
#  filter (pot!=2275)%>%
#  filter (pot!=2391)

pot.E2 <- unique(E2_RIL_ciclo$pot)


NFH.completo.E2.1 <- NFH.completo.E2 %>%
  filter  (pot %in% pot.E2)

dim (NFH.completo.E2.1)

setdiff (NFH.completo.E2.1$pot , E2_RIL_ciclo$pot)
setdiff (E2_RIL_ciclo$pot, NFH.completo.E2.1$pot)


dim (NFH.completo.E2.1 )

NFH.completo.E3 <- em.nfh.E3.mod1 %>%
  mutate (year = 2017)%>%
  mutate (SD = "early")%>%
  mutate (temp =14.66)

 #### despues ver si esto queda así, saque un pot de filocron ###
E3_RIL_ciclo.1 <- E3_RIL_ciclo  %>%
 filter (pot!=3382)
#####################################

dim(NFH.completo.E3)

setdiff (NFH.completo.E3$pot , E3_RIL_ciclo.1$pot)
setdiff (E3_RIL_ciclo.1$pot, NFH.completo.E3$pot)


NFH.completo.E4 <-  em.nfh.E4.mod1 %>%
  mutate (year = 2017)%>%
  mutate (SD = "late")%>%
  mutate (temp =19.23)

E4_RIL_ciclo.1 <- E4_RIL_ciclo %>%
  filter (pot!=4237)%>%
  filter (pot!=4232)%>%
  filter (pot!=4126)%>%
  filter (pot!=4118)



dim (NFH.completo.E4)

setdiff (NFH.completo.E4$pot , E4_RIL_ciclo.1$pot)
setdiff (E4_RIL_ciclo.1$pot, NFH.completo.E4$pot)



NFH.env <- rbind (NFH.completo.E1,
                  NFH.completo.E2.1,
                  NFH.completo.E3,
                  NFH.completo.E4)


pot.filocron.env.2  <- unique(filocron.env.2$pot)

NFH.env.1 <- NFH.env %>%
  filter (pot %in% pot.filocron.env.2) 

NFH.env.1$pot <- as.character(NFH.env.1$pot)
NFH.env.1 <- NFH.env.1 %>%
  arrange (pot) 

summary(NFH.env.1)

#### sobrelamarcha ####

filocron.env.2.1 <- filocron.env.2 %>%
  filter (pot!=4237)%>%
  filter (pot!=4232)%>%
  filter (pot!=4126)%>%
  filter (pot!=4118)%>%
  filter (pot!=3382)%>%
  filter (pot!=3249)

## linea 251 la cambie en feb 2020 ## no estaba filtrada ##
NFH.env.1$pot == filocron.env.2.1$pot 

setdiff (NFH.env.1$pot , filocron.env.2.1$pot)
setdiff (filocron.env.2.1$pot, NFH.env.1$pot)



NFH.env.1 %>%
  mutate (YearFct = fct_rev(as.factor(year))) %>%
  ggplot(aes(y = YearFct)) +
  geom_density_ridges(
    aes(x = em.nfh, fill = paste(YearFct, SD)), 
    alpha = .8, color = "white", from = 6, to = 16
  ) +
  labs(
    x = "Leaf (FN)",
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

write.table (NFH.env.1, file = "./Data/procdata/NFH.env.1.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)


###### hasta aca parece que viene bien #####
### comparar con figura.NFH.N.R lo que falta #######

