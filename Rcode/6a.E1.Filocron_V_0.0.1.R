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

######### cargamos los datos desde outpout.HAUN xa filo ####
### esta es la salidad del GDmodel.3

filocron.E1a <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE1ciclo.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 


filocron.E1b <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE1e_z21.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 



filocron.E1c <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE1z21_z31.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )

testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 


filocron.E1d <- read.table ("./Data/procdata/outpout.HAUN/datos.filocronE1z31_f.csv" ,
                            header = TRUE, sep = ",",dec = ".",
                            na.strings = "NA" )


testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 


filocron.E1a.testigos <- filocron.E1a %>% 
                         dplyr::filter (genotipo %in% tes.geno)%>%
                         dplyr::arrange (parcela)


filocron.E1a.testigos$bloque <- as.factor (filocron.E1a.testigos$bloque)

filocron.E1a.testigos <- filocron.E1a.testigos %>% 
                         dplyr::mutate (parcela.1 = parcela) %>%
                         dplyr::select (-parcela)


filocron.E1b.testigos <- filocron.E1b %>%
                         dplyr::filter (genotipo %in% tes.geno)%>%
                         dplyr::arrange (parcela)

filocron.E1b.testigos$bloque <- as.factor (filocron.E1b.testigos$bloque)

filocron.E1b.testigos <- filocron.E1b.testigos %>% 
                         dplyr::mutate (parcela.1 = parcela) %>%
                         dplyr::select (-parcela)


filocron.E1c.testigos <- filocron.E1c %>% 
                         dplyr::filter (genotipo %in% tes.geno)%>%
                         dplyr::arrange (parcela)

filocron.E1c.testigos$bloque <- as.factor (filocron.E1c.testigos$bloque)

filocron.E1c.testigos <- filocron.E1c.testigos %>% 
                        dplyr::mutate (parcela.1 = parcela) %>%
                        dplyr::select (-parcela)


filocron.E1d.testigos <- filocron.E1d %>% filter (genotipo %in% tes.geno)%>%
                         dplyr::arrange (parcela)

filocron.E1d.testigos$bloque <- as.factor (filocron.E1d.testigos$bloque)

filocron.E1d.testigos <- filocron.E1d.testigos %>% 
                         dplyr::mutate (parcela.1 = parcela) %>%
                         dplyr::select (-parcela)

str(filocron.E1b.testigos)

#filocron.E1b <-  filocron.E1a.testigos %>% arrange (parcela.1) %>%
#  filter (subfase== "ciclo")
#filocron.E1b$parcela.1 %in% E1_ciclo.testigos$parcela
#E1_ciclo.testigos$parcela %in% filocron.E1b$parcela.1
#filocron.E1b$parcela.1  == E1_ciclo.testigos$parcela

#summary (filocron.E1b.testigos)
#View(filocron.E1b.testigos)

filocron.E1a.testigos.1 <- filocron.E1a.testigos %>%
                           dplyr::filter (rsq > 0.51)

filocron.E1a.testigos.2 <- filocron.E1a.testigos.1 %>%
                           dplyr::filter (filocron < 151)


filocron.E1b.testigos.1 <- filocron.E1b.testigos %>%
                           dplyr::filter (filocron < 151)

filocron.E1c.testigos.1 <- filocron.E1c.testigos %>%
                           dplyr::filter (rsq > 0.51)

filocron.E1c.testigos.2 <- filocron.E1c.testigos.1 %>%
                           dplyr::filter (filocron < 151)


filocron.E1d.testigos.1 <- filocron.E1d.testigos %>%
                           dplyr::filter (rsq > 0.51)

filocron.E1d.testigos.2 <- filocron.E1d.testigos.1 %>%
                           dplyr::filter (filocron < 151)


############ graficos exploratorios  ###################3


dotplot (filocron  ~ genotipo, data = filocron.E1a.testigos.2)
bwplot (filocron  ~ genotipo, data = filocron.E1a.testigos.2)
dotplot (filocron  ~ genotipo, data = filocron.E1b.testigos.1)
bwplot (filocron  ~ genotipo, data = filocron.E1b.testigos.1)
dotplot (filocron  ~ genotipo, data = filocron.E1c.testigos.2)
bwplot (filocron  ~ genotipo, data = filocron.E1c.testigos.2)
dotplot (filocron ~ genotipo, data = filocron.E1d.testigos.2)
bwplot (filocron  ~ genotipo, data = filocron.E1d.testigos.2)

###############  Prueba ##################
## Hay que ver si esto esta bien ### el modelo , etc #####
##### obtener un delta filo por parcela ciclo ## para los abline lo precisaría

filocron.E1a.testigos.mod <- lmer (filocron ~  genotipo + (1|bloque), data = filocron.E1a.testigos.2)
anova(filocron.E1a.testigos.mod)
em.filocron.E1a.testigos <- emmeans (filocron.E1a.testigos.mod, ~ genotipo)

em.filocron.E1a.testigos.sum <- summary (em.filocron.E1a.testigos  , infer = c(TRUE,TRUE),
                                         level = .90, adjust = "bon", 
                                         by = c("genotipo"))

em.filocron.E1a.testigos.sum <- em.filocron.E1a.testigos.sum %>% 
  mutate (em.filocron = em.filocron.E1a.testigos.sum$emmean) %>%
  select (genotipo, em.filocron)


plot (em.filocron.E1a.testigos, comparisons = TRUE, alpha = .05)

(cr.filocron <- contrast (em.filocron.E1a.testigos, method = "pairwise"))

write.table (em.filocron.E1a.testigos.sum, file = "./Data/procdata/em.filocron.E1a.testigos.sum.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)


#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.filocron, file = "./Data/procdata/cr.filocron.ciclo.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.filocron <- read.table ("./Data/procdata/cr.filocron.ciclo.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.ciclo.filo<- df.cr.filocron %>%
    dplyr::filter (p.value < 0.05))

####### para e-z21 #############

filocron.E1b.testigos.mod <- lmer (filocron ~  genotipo + (1|bloque), data = filocron.E1b.testigos.1)
anova(filocron.E1b.testigos.mod)
em.filocron.E1b.testigos <- emmeans (filocron.E1b.testigos.mod, ~ genotipo)

em.filocron.E1b.testigos.sum <- summary (em.filocron.E1b.testigos  , infer = c(TRUE,TRUE),
                                         level = .90, adjust = "bon", 
                                         by = c("genotipo"))

em.filocron.E1b.testigos.sum <- em.filocron.E1b.testigos.sum %>% 
  mutate (em.filocron = em.filocron.E1b.testigos.sum$emmean) %>%
  select (genotipo, em.filocron)


plot (em.filocron.E1b.testigos, comparisons = TRUE, alpha = .05)

(cr.filocron <- contrast (em.filocron.E1b.testigos, method = "pairwise"))

write.table (em.filocron.E1b.testigos.sum , file = "./Data/procdata/em.filocron.E1a.testigos.e_z21.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)
#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.filocron, file = "./Data/procdata/cr.filocron.e_21.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.filocron <- read.table ("./Data/procdata/cr.filocron.e_21.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.ciclo.filo<- df.cr.filocron %>%
    dplyr::filter (p.value < 0.05))

###############  z21_z31  ###################

filocron.E1c.testigos.mod <- lmer (filocron ~  genotipo + (1|bloque), data = filocron.E1c.testigos.2)
anova(filocron.E1c.testigos.mod)
em.filocron.E1c.testigos <- emmeans (filocron.E1c.testigos.mod, ~ genotipo)

em.filocron.E1c.testigos.sum <- summary (em.filocron.E1c.testigos  , infer = c(TRUE,TRUE),
                                         level = .90, adjust = "bon", 
                                         by = c("genotipo"))

em.filocron.E1c.testigos.sum <- em.filocron.E1c.testigos.sum %>% 
  mutate (em.filocron = em.filocron.E1c.testigos.sum$emmean) %>%
  select (genotipo, em.filocron)


plot (em.filocron.E1c.testigos, comparisons = TRUE, alpha = .05)

(cr.filocron <- contrast (em.filocron.E1c.testigos, method = "pairwise"))

write.table (em.filocron.E1c.testigos.sum, file = "./Data/procdata/em.filocron.E1a.testigos.z21_z31.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)
#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.filocron, file = "./Data/procdata/cr.filocron.z21_z31.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.filocron <- read.table ("./Data/procdata/cr.filocron.z21_z31.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.ciclo.filo<- df.cr.filocron %>%
    dplyr::filter (p.value < 0.05))

###############  z31_f  ###################

filocron.E1d.testigos.mod <- lmer (filocron ~  genotipo + (1|bloque), data = filocron.E1d.testigos.2)
anova(filocron.E1d.testigos.mod)
em.filocron.E1d.testigos <- emmeans (filocron.E1d.testigos.mod, ~ genotipo)

em.filocron.E1d.testigos.sum <- summary (em.filocron.E1d.testigos  , infer = c(TRUE,TRUE),
                                         level = .90, adjust = "bon", 
                                         by = c("genotipo"))

em.filocron.E1d.testigos.sum <- em.filocron.E1d.testigos.sum %>% 
  mutate (em.filocron = em.filocron.E1d.testigos.sum$emmean) %>%
  select (genotipo, em.filocron)


plot (em.filocron.E1d.testigos, comparisons = TRUE, alpha = .05)

(cr.filocron <- contrast (em.filocron.E1d.testigos, method = "pairwise"))

write.table (em.filocron.E1d.testigos.sum, file = "./Data/procdata/em.filocron.E1a.testigos.z31_f.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)
#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.filocron, file = "./Data/procdata/cr.filocron.z31_f.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.filocron <- read.table ("./Data/procdata/cr.filocron.z31_f.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.ciclo.filo<- df.cr.filocron %>%
    dplyr::filter (p.value < 0.05))



# Armmo la matriz para el plot 

delta.filo.E1 <- cbind ( em.filocron.E1b.testigos.sum [,1:2], 
                       delta.Z21_31.gd = em.filocron.E1c.testigos.sum [,2],
                       delta.Z31_a.gd = em.filocron.E1d.testigos.sum [,2])

names (delta.filo.E1 )
colnames(delta.filo.E1 )
rownames(delta.filo.E1)

delta.filo.E1.t <- t (delta.filo.E1)

#View (delta.em.pred.haun.E1.t)
colnames (delta.filo.E1.t) <- delta.filo.E1.t[1,]
#dim (delta.em.pred.haun.E1.t)

delta.filo.E1.t.1 <- delta.filo.E1.t [-1,]

#View (delta.em.pred.haun.E1.t.1)

delta.filo.E1.t.1a <- cbind (
  Berolina = delta.filo.E1.t.1 [,2],
  Prior = delta.filo.E1.t.1  [,9],
  Quebracho = delta.filo.E1.t.1  [,10],
  Ceibo = delta.filo.E1.t.1  [,5],
  Carumbe = delta.filo.E1.t.1 [,4],
  Kenia = delta.filo.E1.t.1  [,7],
  Bowman = delta.filo.E1.t.1  [,3],
  Logan = delta.filo.E1.t.1  [,8],
  Baronesse = delta.filo.E1.t.1  [,1],
  Danuta = delta.filo.E1.t.1  [,6]
)

#View (delta.em.pred.haun.E1.t.1a)
class(delta.filo.E1.t.1a )
# haun.Z21 promedio 
num.filo.Z21 <- as.numeric (delta.filo.E1.t.1a[1,] )
prom.filo.Z21 <- mean (num.filo.Z21)

## gd.Z21.z31 promedio
num.filo.Z21_31 <- as.numeric (delta.filo.E1.t.1a[2,] )
prom.filo.Z21_31 <- mean (num.filo.Z21_31)


## delta.haun.Z31_NFH.1 promedio
num.filo.Z31.a <- as.numeric (delta.filo.E1.t.1a[3,] )
prom.filo.Z31.a <- mean (num.filo.Z31.a)



svg (filename="./Figures/Fig.2/filo.mayo.E1.svg", 
     width=7, 
     height=5, 
     pointsize=12)

barCenters <- barplot (delta.filo.E1.t.1a , 
                       xlim =c(0,500),
                       beside = FALSE,
                       las = 1,
                       cex.names = 0.75,
                       ylab = NULL,
                       xlab = "filo",
                       col=c("gold","orange1","red3"),
                       border = "black", 
                       axes = TRUE,
                       horiz = TRUE,
                       legend.text = TRUE,
                       args.legend = list(title = "filo.E1", 
                                          x = "topleft",
                                          cex = .7))

abline (v=prom.filo.Z21, col="black", lty=3, lwd=2)
abline (v=prom.filo.Z21 + prom.filo.Z21_31 , col="black", lty=3, lwd=2)
abline (v=prom.filo.Z21 + prom.filo.Z21_31 + prom.filo.Z31.a  , col="black", lty=3, lwd=2)
box()
dev.off()




#########hasta aca #####mayo 2019 ##########
##################################################