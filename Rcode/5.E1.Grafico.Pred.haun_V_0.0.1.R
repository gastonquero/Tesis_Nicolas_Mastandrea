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
library("ggcorrplot")
library(multcomp)
library(ggpubr)

today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos
# cargar datos 
#delta.NH.estado.E1.txt
# delta de numero de hojas

delta.pred.haun.E1 <- read.table ("./Data/procdata/delta.NH.pred.estado.E1.testigos.txt" ,
                                  header = TRUE, sep = ",",dec = ".",
                                  na.strings = "NA" )


parcela.1 <- unique (delta.pred.haun.E1$parcela)

delta.pred.haun.E1$parcela <-  as.factor(delta.pred.haun.E1$parcela)
delta.pred.haun.E1$bloque  <-  as.factor(delta.pred.haun.E1$bloque)


str(delta.pred.haun.E1)
head (delta.pred.haun.E1)
unique (delta.pred.haun.E1$genotipo)


##### obtener el haun ajustado a Z21 
#dotplot(haun.Z21.pred ~ genotipo , data = delta.pred.haun.E1 )

Z21.pred.haun.E1.mod.2 <- lmer (haun.Z21.pred ~ genotipo + (1|bloque), data = delta.pred.haun.E1)

anova (Z21.pred.haun.E1.mod.2)

em.Z21.pred.haun.E1 <- emmeans (Z21.pred.haun.E1.mod.2, ~ genotipo)

plot (em.Z21.pred.haun.E1, comparisons = TRUE, alpha = .05)
(cr.haun.Z21 <- contrast (em.Z21.pred.haun.E1, method = "pairwise"))


write.table (cr.haun.Z21, file = "./Data/procdata/cr.haun.Z21.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

df.cr.haun.Z21 <- read.table ("./Data/procdata/cr.haun.Z21.E1.txt" ,
                                  header = TRUE, sep = ",",dec = ".",
                                  na.strings = "NA" )
df.cr.haun.Z21 <- df.cr.haun.Z21 %>%
                  dplyr::arrange (p.value)

### este filtro te dice quienes son significativamente diferentes       
(sig.leaf.num.z21 <- df.cr.haun.Z21 %>%
                    dplyr::filter (p.value < 0.05))
  

em.Z21.pred.haun.E1.sum <- summary (em.Z21.pred.haun.E1  , infer = c(TRUE,TRUE),
                                    level = .90, adjust = "bon", 
                                    by = c("genotipo"))
head(em.Z21.pred.haun.E1.sum)

em.Z21.pred.haun.E1.sum <- em.Z21.pred.haun.E1.sum %>% 
                           dplyr::mutate (haun.Z21 = em.Z21.pred.haun.E1.sum$emmean) %>%
                           dplyr::select (genotipo, haun.Z21)

### haun a Z31 

Z31.pred.haun.E1.mod.2 <- lmer (haun.Z31.pred ~ genotipo + (1|bloque), data = delta.pred.haun.E1)

anova (Z31.pred.haun.E1.mod.2)

em.Z31.pred.haun.E1 <- emmeans (Z31.pred.haun.E1.mod.2, ~ genotipo)

plot (em.Z31.pred.haun.E1, comparisons = TRUE, alpha = .05)
(cr.haun.Z31 <- contrast (em.Z31.pred.haun.E1, method = "pairwise"))


write.table (cr.haun.Z31, file = "./Data/procdata/cr.haun.Z31.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.haun.Z31 <- read.table ("./Data/procdata/cr.haun.Z31.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

### este filtro te dice quienes son significativamente diferentes       
(sig.leaf.num.z31 <- df.cr.haun.Z31 %>%
    dplyr::filter (p.value < 0.05))


em.Z31.pred.haun.E1.sum <- summary (em.Z31.pred.haun.E1  , infer = c(TRUE,TRUE),
                                    level = .90, adjust = "bon", 
                                    by = c("genotipo"))
head(em.Z31.pred.haun.E1.sum)

em.Z31.pred.haun.E1.sum <- em.Z31.pred.haun.E1.sum %>% 
                           dplyr::mutate (haun.Z31 = em.Z31.pred.haun.E1.sum$emmean) %>%
                           dplyr::select (genotipo, haun.Z31)


##### obterer un delta haun por parcela Z21-Z31

Z31.Z21.pred.haun.E1.mod.2 <- lmer (deltaH.Z31.Z21.pred ~ genotipo + (1|bloque), data = delta.pred.haun.E1)
anova (Z31.Z21.pred.haun.E1.mod.2)

em.Z31.Z21.pred.haun.E1 <- emmeans (Z31.Z21.pred.haun.E1.mod.2, ~ genotipo)


plot (em.Z31.Z21.pred.haun.E1, comparisons = TRUE, alpha = .05)

(cr.Z31.Z21.pred <- contrast (em.Z31.Z21.pred.haun.E1, method = "pairwise"))

#(cld (cr.Z31.Z21.pred, sort=FALSE))


write.table (cr.Z31.Z21.pred, file = "./Data/procdata/cr.Z31.Z21.pred.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.Z31.Z21.pred <- read.table ("./Data/procdata/cr.Z31.Z21.pred.E1.txt" ,
                              header = TRUE, sep = ",",dec = ".",
                              na.strings = "NA" )

(sig.leaf.num.z21.Z31 <- df.cr.Z31.Z21.pred %>%
    dplyr::filter (p.value < 0.05))


em.Z31.Z21.pred.haun.E1.sum <- summary (em.Z31.Z21.pred.haun.E1 , infer = c(TRUE,TRUE),
                                        level = .90, adjust = "bon", 
                                        by = c("genotipo"))


em.Z31.Z21.pred.haun.E1.sum <- em.Z31.Z21.pred.haun.E1.sum %>% 
                               dplyr::mutate (delta.haun.Z21_31 = em.Z31.Z21.pred.haun.E1.sum$emmean) %>%
                               dplyr::select (genotipo, delta.haun.Z21_31)


##### obterer un delta haun por parcela Z31-a

NFH.Z31.pred.haun.E1.mod.2 <- lmer (deltaH.NFH.Z31.pred ~ genotipo + (1|bloque), data = delta.pred.haun.E1)
anova (NFH.Z31.pred.haun.E1.mod.2)
em.NFH.Z31.pred.haun.E1 <- emmeans (NFH.Z31.pred.haun.E1.mod.2, ~ genotipo)


em.NFH.Z31.pred.haun.E1.sum <- summary (em.NFH.Z31.pred.haun.E1 , infer = c(TRUE,TRUE),
                                          level = .90, adjust = "bon", 
                                          by = c("genotipo"))


plot (em.NFH.Z31.pred.haun.E1, comparisons = TRUE, alpha = .05)

(cr.NFH.Z31.pred <- contrast (em.NFH.Z31.pred.haun.E1, method = "pairwise"))

#(cld (cr.NFH.Z31.pred, sort=FALSE))


write.table (cr.NFH.Z31.pred, file = "./Data/procdata/cr.NFH.Z31.pred.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.NFH.Z31.pred <- read.table ("./Data/procdata/cr.NFH.Z31.pred.E1.txt" ,
                                  header = TRUE, sep = ",",dec = ".",
                                  na.strings = "NA" )

(sig.leaf.num.NFH.Z31.pred<- df.cr.NFH.Z31.pred %>%
    dplyr::filter (p.value < 0.05))



em.NFH.Z31.pred.haun.E1.sum <- em.NFH.Z31.pred.haun.E1.sum %>% 
                                dplyr::mutate (delta.haun.Z31_NFH = em.NFH.Z31.pred.haun.E1.sum$emmean) %>%
                                dplyr::select (genotipo, delta.haun.Z31_NFH)


##### obterer un delta haun por parcela antesis ## para los abline lo precisaría

NFH.pred.haun.E1.mod <- lmer (haun.NFH.pred ~  genotipo + (1|bloque), data = delta.pred.haun.E1)
anova (NFH.pred.haun.E1.mod)
em.NFH.pred.haun.E1 <- emmeans (NFH.pred.haun.E1.mod, ~ genotipo)

em.NFH.pred.haun.E1.sum <- summary (em.NFH.pred.haun.E1  , infer = c(TRUE,TRUE),
                                      level = .90, adjust = "bon", 
                                      by = c("genotipo"))

em.NFH.pred.haun.E1.sum <- em.NFH.pred.haun.E1.sum %>% 
                          dplyr::mutate (em.NFH = em.NFH.pred.haun.E1.sum$emmean) %>%
                          dplyr::select (genotipo, em.NFH)


plot (em.NFH.pred.haun.E1, comparisons = TRUE, alpha = .05)

(cr.NFH.pred <- contrast (em.NFH.pred.haun.E1, method = "pairwise"))

write.table (em.NFH.pred.haun.E1.sum, file = "./Data/procdata/em.NFH.pred.haun.E1.sum.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)


write.table (cr.NFH.pred, file = "./Data/procdata/cr.NFH.pred.E1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

df.cr.NFH.pred <- read.table ("./Data/procdata/cr.NFH.pred.E1.txt" ,
                                  header = TRUE, sep = ",",dec = ".",
                                  na.strings = "NA" )

(sig.leaf.num.NFH.pred<- df.cr.NFH.pred %>%
    dplyr::filter (p.value < 0.05))


#########################################
# Armmo la matriz para el plot 

delta.em.pred.haun.E1 <- cbind ( em.Z21.pred.haun.E1.sum [,1:2], 
                                 delta.haun.Z21_31 = em.Z31.Z21.pred.haun.E1.sum[,2],
                                 delta.haun.Z31_NFH = em.NFH.Z31.pred.haun.E1.sum[,2])

names (delta.em.pred.haun.E1 )
colnames(delta.em.pred.haun.E1)
rownames(delta.em.pred.haun.E1)

delta.em.pred.haun.E1.t <- t (delta.em.pred.haun.E1)

#View (delta.em.pred.haun.E1.t)
colnames (delta.em.pred.haun.E1.t) <- delta.em.pred.haun.E1.t [1,]
#dim (delta.em.pred.haun.E1.t)

delta.em.pred.haun.E1.t.1 <- delta.em.pred.haun.E1.t [-1,]

#View (delta.em.pred.haun.E1.t.1)

delta.em.pred.haun.E1.t.1a <- cbind (
  Berolina = delta.em.pred.haun.E1.t.1 [,2],
  Prior = delta.em.pred.haun.E1.t.1 [,9],
  Quebracho = delta.em.pred.haun.E1.t.1 [,10],
  Ceibo = delta.em.pred.haun.E1.t.1 [,5],
  Carumbe = delta.em.pred.haun.E1.t.1 [,4],
  Kenia = delta.em.pred.haun.E1.t.1 [,7],
  Bowman = delta.em.pred.haun.E1.t.1 [,3],
  Logan = delta.em.pred.haun.E1.t.1 [,8],
  Baronesse = delta.em.pred.haun.E1.t.1 [,1],
  Danuta = delta.em.pred.haun.E1.t.1 [,6]
  )

#View (delta.em.pred.haun.E1.t.1a)
class(delta.em.pred.haun.E1.t.1a)
# haun.Z21 promedio 
num.haun.Z21 <- as.numeric (delta.em.pred.haun.E1.t.1a[1,] )
prom.haun.Z21 <- mean (num.haun.Z21)

## haun.Z21 promedio
num.delta.haun.Z21_31 <- as.numeric (delta.em.pred.haun.E1.t.1a[2,] )
prom.delta.haun.Z21_31 <- mean (num.delta.haun.Z21_31 )


## delta.haun.Z31_NFH.1 promedio
num.delta.haun.Z31_NFH <- as.numeric (delta.em.pred.haun.E1.t.1a[3,] )
prom.delta.haun.Z31_NFH <- mean (num.delta.haun.Z31_NFH)

svg (filename="./Figures/Fig.2/haun.pred.E1.testigos.svg", 
     width=7, 
     height=5, 
     pointsize=12)

barCenters <- barplot (delta.em.pred.haun.E1.t.1a , 
                       xlim =c(0,15),
                       beside = FALSE,
                       las = 1,
                       cex.names = 0.75,
                       ylab = NULL,
                       xlab = "haun",
                       col=c("gold","orange1","red3"),
                       border = "black", 
                       axes = TRUE,
                       horiz = TRUE,
                       legend.text = TRUE,
                       args.legend = list(title = "haun.pred.E1a", 
                                          x = "topleft",
                                          cex = .7))

abline (v=prom.haun.Z21, col="black", lty=3, lwd=2)
abline (v=prom.haun.Z21 + prom.delta.haun.Z21_31 , col="black", lty=3, lwd=2)
abline (v=prom.haun.Z21 + prom.delta.haun.Z21_31+prom.delta.haun.Z31_NFH  , col="black", lty=3, lwd=2)
box()
dev.off()

#####################################################












