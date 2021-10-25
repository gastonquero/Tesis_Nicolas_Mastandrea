############################################################################
# Datos  Cebada - ceiboxcarumbe                                           ##
# Datos tomados por Nicolas Mastandrea                                    ##
# Gaston Quero  - Nicolas Mastandrea                                      ##
# 27-06-2018                                                             ##
############################################################################

getwd()
setwd ("R:/tesis_Copia_Nueva_Mastandrea")

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
E1_ciclo <- read.table ("./Data/rawdata/E1_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )

E1_ciclo.1 <- E1_ciclo %>%
              dplyr::mutate (env = E1_ciclo$ambiente)%>%
              dplyr::select (-ambiente) %>%
              dplyr::mutate (pot =str_c (bloque,parcela))



head(E1_ciclo.1)
summary(E1_ciclo)
str(E1_ciclo.1)


ciclo.gd.E1 <- E1_ciclo %>%
               dplyr::select (genotipo, anio,ambiente,bloque,parcela,e.z21.gd, 
                              z21.z31.gd,z31.a.gd,e.a.gd) %>% 
               dplyr::arrange (parcela)

write.table (ciclo.gd.E1, file = "./Data/procdata/ciclo.gd.RIL.E1.txt",
              append = FALSE, quote = TRUE, sep = ",",
              eol = "\n", na = "NA", dec = ".", row.names =FALSE,
             col.names = TRUE)


#E1.ciclo.gd <- read.table ("./Data/procdata/outpout.HAUN/ciclo.gd.RIL.E1.txt" ,
 #                          header = TRUE, sep = ",",dec = ".",
  #                         na.strings = "NA" )

testigos <- c("testigo1","testigo2","testigo3","testigo4",
              "testigo5","testigo6","testigo7","testigo8",
              "testigo9","testigo10")

tes.geno <- c("Kenia","Prior","Logan",  
              "Baronesse","Berolina",
              "Bowman","Quebracho",   
              "Danuta", "Carumbe","Ceibo") 

E1_ciclo.testigos <- ciclo.gd.E1 %>% 
                     dplyr::filter (genotipo %in% tes.geno)%>%
                     dplyr::arrange (parcela)



#View (E1.ciclo.gd)
unique(E1_ciclo.testigos$genotipo)
head(E1_ciclo.testigos )

#plot (z31.a.gd ~  z21.z31.gd, ylim=c(300,500), data = E1.ciclo.gd)

#xyplot(z31.a.gd ~  z21.z31.gd | genotipo,data = E1.ciclo.gd )
  
E1_ciclo.testigos <-  E1_ciclo.testigos %>% 
                      dplyr::arrange (parcela)

#str(fenologia.E1)
#dim(fenologia.E1)


E1_ciclo.testigos$bloque <- as.factor (E1_ciclo.testigos$bloque)

##E1.ciclo.gd


############ graficos exploratorios  ###################3


dotplot (e.z21.gd  ~ genotipo, data = E1_ciclo.testigos)
bwplot (e.z21.gd  ~ genotipo, data = E1_ciclo.testigos)
dotplot (z21.z31.gd  ~ genotipo, data = E1_ciclo.testigos)
bwplot (z21.z31.gd  ~ genotipo, data = E1_ciclo.testigos)
dotplot (z31.a.gd  ~ genotipo, data = E1_ciclo.testigos)
bwplot (z31.a.gd   ~ genotipo, data = E1_ciclo.testigos)


#dotplot (s.a.gd   ~ genotipo, data = fenologia.E1)


##### e.z21.gd #####
e.z21.gd.mod.1 <- lmer (e.z21.gd ~ genotipo + (1|bloque), data = E1_ciclo.testigos)
anova(e.z21.gd.mod.1)
em.e.z21.gd.E1 <- emmeans (e.z21.gd.mod.1, ~ genotipo)

em.e.z21.gd.E1.sum <- summary (em.e.z21.gd.E1  , infer = c(TRUE,TRUE),
                                    level = .90, adjust = "bon", 
                                    by = c("genotipo"))

em.e.z21.gd.E1.sum <- em.e.z21.gd.E1.sum %>% 
                      dplyr::mutate (e.z21.gd = em.e.z21.gd.E1.sum$emmean) %>%
                      dplyr::select (genotipo, e.z21.gd)


plot (em.e.z21.gd.E1 , comparisons = TRUE, alpha = .05)

(cr.e.z21.gd.E1 <- contrast (em.e.z21.gd.E1 , method = "pairwise"))

#(cld (cr.NFH.Z31.pred, sort=FALSE))

write.table (cr.e.z21.gd.E1, file = "./Data/procdata/cr.e.z21.gd.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.e.z21.gd.E1 <- read.table ("./Data/procdata/cr.e.z21.gd.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.e.z21.gd <- df.cr.e.z21.gd.E1  %>%
    dplyr::filter (p.value < 0.05))


################ 
##### z21.z31.gd #####
z21.z31.gd.mod.1 <- lmer (z21.z31.gd ~ genotipo + (1|bloque), data = E1_ciclo.testigos)
anova(z21.z31.gd.mod.1)
em.z21.z31.gd.E1 <- emmeans (z21.z31.gd.mod.1, ~ genotipo)

em.z21.z31.gd.E1.sum <- summary (em.z21.z31.gd.E1  , infer = c(TRUE,TRUE),
                               level = .90, adjust = "bon", 
                               by = c("genotipo"))

em.z21.z31.gd.E1.sum <- em.z21.z31.gd.E1.sum %>% 
                        dplyr::mutate (z21.z31.gd = em.z21.z31.gd.E1.sum$emmean) %>%
                        dplyr::select (genotipo, z21.z31.gd)


plot (em.z21.z31.gd.E1 , comparisons = TRUE, alpha = .05)

(cr.z21.z31.gd.E1 <- contrast (em.z21.z31.gd.E1 , method = "pairwise"))


write.table (cr.z21.z31.gd.E1, file = "./Data/procdata/cr.z21.z31.gd.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.z21.z31.gd.E1 <- read.table ("./Data/procdata/cr.z21.z31.gd.E1.txt" ,
                                 header = TRUE, sep = ",",dec = ".",
                                 na.strings = "NA" )

(sig.z21.z31.gd <- df.cr.z21.z31.gd.E1  %>%
    dplyr::filter (p.value < 0.05))

################ 
##### z31.a.gd #####

z31.a.gd.mod.1 <- lmer (z31.a.gd ~ genotipo + (1|bloque), data = E1_ciclo.testigos)
anova(z31.a.gd.mod.1)
em.z31.a.gd.E1 <- emmeans (z31.a.gd.mod.1, ~ genotipo)

em.z31.a.gd.E1.sum <- summary (em.z31.a.gd.E1  , infer = c(TRUE,TRUE),
                                 level = .90, adjust = "bon", 
                                 by = c("genotipo"))

em.z31.a.gd.E1.sum <- em.z31.a.gd.E1.sum %>% 
                      dplyr::mutate (z31.a.gd = em.z31.a.gd.E1.sum$emmean) %>%
                      dplyr::select (genotipo, z31.a.gd)


plot (em.z31.a.gd.E1 , comparisons = TRUE, alpha = .05)

(cr.z31.a.gd.E1 <- contrast (em.z31.a.gd.E1 , method = "pairwise"))

#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.z31.a.gd.E1, file = "./Data/procdata/cr.z31.a.gd.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.z31.a.gd.E1 <- read.table ("./Data/procdata/cr.z31.a.gd.E1.txt" ,
                                   header = TRUE, sep = ",",dec = ".",
                                   na.strings = "NA" )

(sig.z31.a.gd <- df.cr.z31.a.gd.E1  %>%
    dplyr::filter (p.value < 0.05))


### emergencia antesis #### 

E1_ciclo.testigos <- E1_ciclo.testigos %>%
                     dplyr::mutate (e.a.gd = E1_ciclo.testigos$e.z21.gd +
                         E1_ciclo.testigos$z21.z31.gd +E1_ciclo.testigos$z31.a.gd)


##### e.a.gd #####
e.a.gd.mod.1 <- lmer (e.a.gd ~ genotipo + (1|bloque), data = E1_ciclo.testigos)

anova(e.a.gd.mod.1)
em.e.a.gd.E1 <- emmeans (e.a.gd.mod.1, ~ genotipo)

em.e.a.gd.E1.sum <- summary (em.e.a.gd.E1  , infer = c(TRUE,TRUE),
                               level = .90, adjust = "bon", 
                               by = c("genotipo"))

em.e.a.gd.E1.sum <- em.e.a.gd.E1.sum %>% 
                    dplyr::mutate (e.a.gd = em.e.a.gd.E1.sum$emmean) %>%
                    dplyr::select (genotipo, e.a.gd)


plot (em.e.a.gd.E1 , comparisons = TRUE, alpha = .05)

(cr.e.a.gd.E1 <- contrast (em.e.a.gd.E1 , method = "pairwise"))

write.table (em.e.a.gd.E1.sum, file = "./Data/procdata/em.e.a.gd.E1.sum.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)
#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.e.a.gd.E1, file = "./Data/procdata/cr.e.a.gd.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.e.a.gd.E1 <- read.table ("./Data/procdata/cr.e.a.gd.E1.txt" ,
                                 header = TRUE, sep = ",",dec = ".",
                                 na.strings = "NA" )

(sig.e.a.gd <- df.cr.e.a.gd.E1  %>%
    dplyr::filter (p.value < 0.05))


########## 

# Armmo la matriz para el plot 

delta.gd.E1 <- cbind ( em.e.z21.gd.E1.sum  [,1:2], 
                       delta.Z21_31.gd = em.z21.z31.gd.E1.sum [,2],
                       delta.Z31_a.gd = em.z31.a.gd.E1.sum [,2])

names (delta.gd.E1 )
colnames(delta.gd.E1 )
rownames(delta.gd.E1)

delta.gd.E1.t <- t (delta.gd.E1)

#View (delta.em.pred.haun.E1.t)
colnames (delta.gd.E1.t) <- delta.gd.E1.t[1,]
#dim (delta.em.pred.haun.E1.t)

delta.gd.E1.t.1 <- delta.gd.E1.t [-1,]

#View (delta.em.pred.haun.E1.t.1)

delta.gd.E1.t.1a <- cbind (
  Berolina = delta.gd.E1.t.1 [,2],
  Prior = delta.gd.E1.t.1  [,9],
  Quebracho = delta.gd.E1.t.1  [,10],
  Ceibo = delta.gd.E1.t.1  [,5],
  Carumbe = delta.gd.E1.t.1 [,4],
  Kenia = delta.gd.E1.t.1  [,7],
  Bowman = delta.gd.E1.t.1  [,3],
  Logan = delta.gd.E1.t.1  [,8],
  Baronesse =delta.gd.E1.t.1  [,1],
  Danuta = delta.gd.E1.t.1  [,6]
)

#View (delta.em.pred.haun.E1.t.1a)
class(delta.gd.E1.t.1a )
# haun.Z21 promedio 
num.gd.Z21 <- as.numeric (delta.gd.E1.t.1a[1,] )
prom.gd.Z21 <- mean (num.gd.Z21)

## gd.Z21.z31 promedio
num.gd.Z21_31 <- as.numeric (delta.gd.E1.t.1a[2,] )
prom.gd.Z21_31 <- mean (num.gd.Z21_31)


## delta.haun.Z31_NFH.1 promedio
num.gd.Z31.a <- as.numeric (delta.gd.E1.t.1a[3,] )
prom.gd.Z31.a <- mean (num.gd.Z31.a)



svg (filename="./Figures/Fig.1/GD.testigos.E1.svg", 
     width=7, 
     height=5, 
     pointsize=12)

barCenters <- barplot (delta.gd.E1.t.1a , 
                       xlim =c(0,1800),
                       beside = FALSE,
                       las = 1,
                       cex.names = 0.75,
                       ylab = NULL,
                       xlab = "grados día",
                       col=c("gold","orange1","red3"),
                       border = "black", 
                       axes = TRUE,
                       horiz = TRUE,
                       legend.text = TRUE,
                       args.legend = list(title = "gd.E1", 
                                          x = "topleft",
                                          cex = .7))

abline (v=prom.gd.Z21, col="black", lty=3, lwd=2)
abline (v=prom.gd.Z21 + prom.gd.Z21_31 , col="black", lty=3, lwd=2)
abline (v=prom.gd.Z21 + prom.gd.Z21_31 + prom.gd.Z31.a  , col="black", lty=3, lwd=2)
box()
dev.off()




#########hasta aca #####mayo 2019 ##########
##################################################
