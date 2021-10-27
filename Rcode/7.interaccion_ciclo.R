############################################################################
# Datos  Cebada - ceiboxcarumbe                                           ##
# Datos tomados por Nicolas Mastandrea                                    ##
# Gaston Quero  - Nicolas Mastandrea                                      ##
# 27-06-2018                                                             ##
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
library(ggpubr)
library(lmerTest)
library(nlme)
library(emmeans)
library("car")
#library("gplots")    
library("ggplot2")       
library("plotrix")  
library("lattice")
library("latticeExtra")
library(multcompView)
#library(effects)
library(dplyr)
library(xtable)
library(multcomp)
??cld.emmGrid
today ()
col.l2 <- colorRampPalette (c("orange","blue"))(50)

#########  se carga los datos crudos
# cargar datos 
#E1
E1_ciclo <- read.table ("./Data/rawdata/E1_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )

E1_ciclo.1 <- E1_ciclo %>%
              dplyr::mutate ( env= "E1") %>%
              dplyr::mutate (bloque = str_c (env, bloque)) %>%
              dplyr::mutate (pot =str_c (bloque,parcela))

unique (E1_ciclo.1$env)

#E2
E2_ciclo <- read.table ("./Data/rawdata/E2_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )

E2_ciclo.1 <- E2_ciclo %>%
              dplyr::mutate (env= "E2") %>%
              dplyr::mutate (bloque = str_c (env, bloque)) %>%
              dplyr::mutate (pot =str_c (bloque,parcela))

unique (E2_ciclo.1$env)

#E3
E3_ciclo <- read.table ("./Data/rawdata/E3_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )

E3_ciclo.1 <- E3_ciclo %>%
              dplyr::mutate (env = "E3")%>%
              dplyr::mutate (bloque = str_c (env, bloque)) %>%
              dplyr::mutate (pot =str_c (bloque,parcela))

unique (E3_ciclo.1$env)

#E4
E4_ciclo <- read.table ("./Data/rawdata/E4_ciclo.txt" ,
                        header = TRUE, sep = "\t",dec = ".",
                        na.strings = "NA" )

E4_ciclo.1 <- E4_ciclo %>%
              dplyr::mutate (env="E4")%>%
              dplyr::mutate (bloque = str_c (env, bloque)) %>%
              dplyr::mutate (pot =str_c (bloque,parcela))

unique (E4_ciclo.1$env)

### Los cuatro ambientes #####
 
ciclo.env <-  bind_rows (E1_ciclo.1, 
                         E2_ciclo.1 ,
                         E3_ciclo.1, 
                         E4_ciclo.1 )

unique (ciclo.env$env)

ciclo.gd <- ciclo.env %>%
            dplyr::select (genotipo, env, anio,ambiente,bloque,parcela,e.z21.gd, 
                              z21.z31.gd,z31.a.gd,e.a.gd) 

write.table (ciclo.gd, file = "./Data/procdata/ciclo.gd.RIL.txt",
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

ciclo.testigos <- ciclo.gd %>% 
                  dplyr::filter (genotipo %in% tes.geno) %>%
                  dplyr::filter (genotipo != "Berolina")

max(ciclo.testigos$e.a.gd, na.rm = TRUE)
min(ciclo.testigos$e.a.gd, na.rm = TRUE)
#View (E1.ciclo.gd)
unique(ciclo.testigos$genotipo)
unique(ciclo.testigos$env)
head(ciclo.testigos )

# 
############ graficos exploratorios  ###################3
amb.12 <- c("E1", "E2")
E1.E2 <- ciclo.testigos %>%
         dplyr::filter (env %in% amb.12)

plot.1 <- ggstripchart(E1.E2, "env", "e.a.gd", 
              ylim=c(500, 1800),
             #shape = "supp",
             size	=1,
             color = "genotipo", 
             add = c("mean_sd"),
             palette = c("chartreuse4", "gray48","darkorange","navyblue", "azure3",  "slategray4", "darkolivegreen3", "royalblue", "darkolivegreen4", "red1"),
             position = position_dodge(0.8)
             )

amb.34 <- c("E3", "E4")
E3.E4 <- ciclo.testigos %>%
        dplyr::filter (env %in% amb.34)

plot.2 <- ggstripchart(E3.E4, "env", "e.a.gd", 
                       ylim=c(500, 1800),
                       #shape = "supp",
                       size	=1,
                       color = "genotipo", 
                       add = c("mean_sd"),
                       palette = c("chartreuse4", "gray48","darkorange","navyblue", "azure3",  "slategray4", "darkolivegreen3", "royalblue", "darkolivegreen4", "red1"),
                       position = position_dodge(0.8)
)
ggarrange(plot.1, plot.2, ncol=2, nrow=1, common.legend = TRUE )

ciclo.testigos$bloque <- as.factor (ciclo.testigos$bloque)
ciclo.testigos$env <- as.factor (ciclo.testigos$env)

list.geno <- c("Kenia","Prior","Logan",  
               "Baronesse",
               "Bowman","Quebracho",   
               "Danuta", "Carumbe","Ceibo") 


head (ciclo.testigos )

##### e.a.gd #####
e.a.gd.mod.1 <- lmer (e.a.gd ~ (1|bloque) + genotipo + env + genotipo * env , data = ciclo.testigos)

anova(e.a.gd.mod.1)


##e.z21.gd 
e.z21.gd.mod.1 <- lmer (e.z21.gd ~ (1|bloque) + genotipo + env + genotipo * env , data = ciclo.testigos)

anova(e.z21.gd.mod.1)

## z21.z31.gd 
z21.z31.gd.mod.1 <- lmer ( z21.z31.gd  ~ (1|bloque) + genotipo + env + genotipo * env , data = ciclo.testigos)

anova (z21.z31.gd.mod.1)


##z31.a.gd 
z31.a.gd.mod.1 <- lmer ( z31.a.gd   ~ (1|bloque) + genotipo + env + genotipo * env , data = ciclo.testigos)
anova (z31.a.gd.mod.1)


## funcion que calcula los contrastes 
#data.model = e.a.gd.mod.1
#trait = "e.a.gd"
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

####### 
contrastes_e.a.gd <- run_contrastes (data.model = e.a.gd.mod.1, trait = "e.a.gd" )

contrastes_z21.z31.gd <- run_contrastes (data.model = z21.z31.gd.mod.1, trait = "z21.z31.gd" )

contrastes_z31.a.gd <- run_contrastes (data.model = z31.a.gd.mod.1, trait = "z31.a.gd" )

 contrastes_suma_termica <- bind_rows (contrastes_e.a.gd,contrastes_z21.z31.gd, contrastes_z31.a.gd  )

 write_csv2 (contrastes_suma_termica, file= "./Data/procdata/contrastes_suma_termica.csv")

 
 ###############3 hasat cas el 27 /10/2021
em.e.a.gd.Kenia <- emmeans (e.a.gd.mod.1, "env",
                                     at = list (genotipo = "Kenia"))





em.e.a.gd.Kenia <- cbind (genotipo = "Kenia", parametro ="e.a.gd",
                          cld (em.e.a.gd.Kenia, sort=FALSE))











e.a.gd.mod.1.em <- emmeans (e.a.gd.mod.1, ~ genotipo*env)

e.a.gd.mod.1.sum <- summary (e.a.gd.mod.1.em  , infer = c(TRUE,TRUE),
                             level = .90, adjust = "bon", 
                             by = c("genotipo", "env"))

e.a.gd.mod.1.sum.1 <- e.a.gd.mod.1.sum %>% 
                     dplyr::mutate (e.a.gd = e.a.gd.mod.1.sum$emmean) %>%
                     dplyr::select (genotipo, env, e.a.gd,SE, lower.CL, upper.CL)




write.table (e.a.gd.mod.1.sum.1, file = "./Data/procdata/e.a.gd.mod.1.sum.1.txt", 
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)

unique (e.a.gd.mod.1.em$)








e.a.gd.mod.1.sum.2 <- read.table ("./Data/procdata/e.a.gd.mod.1.sum.1.txt" ,
                                 header = TRUE, sep = ",",dec = ".",
                                 na.strings = "NA" )

anio1 <- c("E1", "E2")

e.a.gd.anio1 <- e.a.gd.mod.1.sum.2 %>%
                dplyr::filter (env %in% anio1)

plot.3 <- ggpaired (e.a.gd.anio1, 
          ylim=c(500, 1800),
          x = "env", y = "e.a.gd",
          id = "genotipo",
          shape="geotipo",
          color = "genotipo", 
          line.color = "genotipo", 
          line.size = 1,
          width =0,
          point.size = 3,
          xlab = "zafra",
          ylab = "gd",
          label="genotipo",
          repel = TRUE,
          palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"))


anio2 <- c("E3", "E4")

e.a.gd.anio2 <- e.a.gd.mod.1.sum.2 %>%
                dplyr::filter (env %in% anio2)

plot.4 <- ggpaired (e.a.gd.anio2, 
                    ylim=c(500, 1800),
                    x = "env", y = "e.a.gd",
                    id = "genotipo",
                    shape="geotipo",
                    color = "genotipo", 
                    line.color = "genotipo", 
                    line.size = 1,
                    width =0,
                    point.size = 3,
                    xlab = "zafra",
                    ylab = "gd",
                    label="genotipo",
                    repel = TRUE,
                    palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"))

ggarrange(plot.3, plot.4, ncol=2, nrow=1, common.legend = TRUE )

########### barplot ############

e.a.gd.anio1.1 <- e.a.gd.mod.1.sum.2 %>%
                  dplyr::filter (env == "E1")

e.a.gd.anio1.2 <- e.a.gd.mod.1.sum.2 %>%
                  dplyr::filter (env == "E2")
                 
delta.gd.anio1 <- e.a.gd.anio1.1 %>%
                  dplyr::inner_join(e.a.gd.anio1.2, by="genotipo" )%>%
                  dplyr::mutate (esta =e.a.gd.y  -  e.a.gd.x)


plot.5 <- ggbarplot (delta.gd.anio1 , x = "genotipo", y = "esta",
                 fill = "genotipo",           # change fill color by mpg_level
                 #color = "white",            # Set bar border colors to white
                 palette = c("orchid", "palegreen1", "yellow4", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"),            # jco journal color palett. see ?ggpar
                 sort.val = "asc",           # Sort the value in ascending order
                 sort.by.groups = FALSE,     # Don't sort inside each group
                 x.text.angle = 90,          # Rotate vertically x axis texts
                 ylab = "paravos",
                 xlab = FALSE,
                 rotate = TRUE
                 ,legend.title = "genotipo")

ggarrange(plot.5, plot.6, ncol=2, nrow=1, common.legend = TRUE )




e.a.gd.anio2.3 <- e.a.gd.mod.1.sum.2 %>%
                  dplyr::filter (env == "E3")

e.a.gd.anio2.4 <- e.a.gd.mod.1.sum.2 %>%
                  dplyr::filter (env == "E4")

delta.gd.anio2 <- e.a.gd.anio2.3 %>%
                  dplyr::inner_join(e.a.gd.anio2.4, by="genotipo" )%>%
                  dplyr::mutate (esta =e.a.gd.y  -  e.a.gd.x)

x2 <- mean (delta.gd.anio2$esta )
plot.6 <- ggbarplot (delta.gd.anio2 , x = "genotipo", y = "esta",
                     fill = "genotipo",           # change fill color by mpg_level
                     #color = "white",            # Set bar border colors to white
                     palette = c("steelblue3", "tan1", "violetred3", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"),            # jco journal color palett. see ?ggpar
                     sort.val = "asc",           # Sort the value in ascending order
                     sort.by.groups = FALSE,     # Don't sort inside each group
                     x.text.angle = 90,          # Rotate vertically x axis texts
                     ylab = "paravos",
                     xlab = FALSE,
                     rotate = TRUE
                     ,legend.title = "genotipo") +
                geom_vline (xintercept = 400, linetype = 2)
                 



##E1.ciclo.gd


############ graficos exploratorios  ###################3

plot.3 <- ggline (E1.E2,
                   x="env", y="e.a.gd", 
                   #size = 0.9,
                   #shape = 19,
                  label = "genotipo",
                   #repel = TRUE,
                   plot_type="b",
                   color = "genotipo",
                   add = c("mean"),
                   palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"))

x2 <- ggbarplot (xx.1, x = "genotipo", y = "esta",
                 fill = "genotipo",           # change fill color by mpg_level
                 #color = "white",            # Set bar border colors to white
                 palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"),            # jco journal color palett. see ?ggpar
                 sort.val = "asc",           # Sort the value in ascending order
                 sort.by.groups = FALSE,     # Don't sort inside each group
                 x.text.angle = 90,          # Rotate vertically x axis texts
                 ylab = "paravos",
                 xlab = FALSE,
                 rotate = TRUE
                 ,legend.title = "genotipo")



ggline (E3.E4,
        x="env", y="e.a.gd", 
        #size = 0.9,
        #shape = 19,
        repel = TRUE,
        plot_type="b",
        color = "genotipo",
        add = c("mean"),
        palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"))



+
  geom_hline (yintercept = 1,  linetype = 2)







dotplot (e.a.gd  ~ genotipo + env, data = ciclo.testigos)
bwplot (e.z21.gd  ~ genotipo, data = E1_ciclo.testigos)
dotplot (z21.z31.gd  ~ genotipo, data = E1_ciclo.testigos)
bwplot (z21.z31.gd  ~ genotipo, data = E1_ciclo.testigos)
dotplot (z31.a.gd  ~ genotipo, data = E1_ciclo.testigos)
bwplot (z31.a.gd   ~ genotipo, data = E1_ciclo.testigos)

ggp

#dotplot (s.a.gd   ~ genotipo, data = fenologia.E1)




e.a.gd.mod.1.sum.1$env <- as.character(e.a.gd.mod.1.sum.1$env )



amb.12 <- c("E1", "E2")
E1.E2a <- e.a.gd.mod.1.sum.1 %>%
          dplyr::filter (env %in% amb.12)



unique ( E1.E2a$env)


E1a  <-  e.a.gd.mod.1.sum.1 %>%
         dplyr::filter (env== "E1")

E2a  <-  e.a.gd.mod.1.sum.1 %>%
        dplyr::filter (env== "E2")


xx <- E1a %>%
      dplyr::inner_join(E2a, by="genotipo" ) 

head(xx )
xx.1 <- xx %>%
        dplyr::mutate (esta =e.a.gd.y  -  e.a.gd.x)



E3a  <-  e.a.gd.mod.1.sum.1 %>%
  dplyr::filter (env== "E3")

E4a  <-  e.a.gd.mod.1.sum.1 %>%
  dplyr::filter (env== "E4")


xx2 <- E3a %>%
  dplyr::inner_join(E4a, by="genotipo" ) 

head(xx )
xx.2 <- xx2 %>%
  dplyr::mutate (esta =e.a.gd.y  -  e.a.gd.x)
ggbarplot (xx.2, x = "genotipo", y = "esta",
           fill = "genotipo",           # change fill color by mpg_level
           #color = "white",            # Set bar border colors to white
           palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"),            # jco journal color palett. see ?ggpar
           sort.val = "asc",           # Sort the value in ascending order
           sort.by.groups = FALSE,     # Don't sort inside each group
           x.text.angle = 90,          # Rotate vertically x axis texts
           ylab = "paravos",
           xlab = FALSE,
           rotate = TRUE
           ,legend.title = "genotipo")


levels( E1.E2a$env)
         
x2 <- ggpaired (E1.E2, 
          x = "env", y = "e.a.gd",
          id = "genotipo",
          shape="geotipo",
          color = "genotipo", 
          line.color = "genotipo", 
          line.size = 0.5,
          width =0,
          point.size = 3,
          xlab = "zafra",
          ylab = "gd",
          #label="genotipo",
          repel = TRUE,
          palette = c("black", "red", "blue", "orange", "pink", "green", "gray", "darkgreen", "darkred", "yellow"))



ggarrange(x1, x2, ncol=2, nrow=1, common.legend = TRUE )




+
  geom_hline (yintercept = 1,  linetype = 2)



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
