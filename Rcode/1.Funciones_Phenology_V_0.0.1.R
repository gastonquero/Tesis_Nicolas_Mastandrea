################################################################
#   Funciones individuales para el paquete de Fenologia       ##
# Gaston Quero - Nicolas Mastandrea                           ##
# 24-5- 2019                                                  ##
################################################################

getwd()
#setwd ("E:/tesis_Copia_Nueva_Mastandrea")
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
library("stringr")


# funcion para calcular el numero de hoja 
# en cada momento y los deltas de hojas
NH <- function (dt.1 = NULL, dt.2 = NULL) {
  dir.create (file.path ("Data"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata", "outpout.HAUN"), showWarnings = FALSE)
  
  listapot.1 <- unique(dt.1$pot)
  
  dt.haun <- lapply (listapot.1, function(filtro){
    datos.1 <- filter (dt.1, pot== filtro)
    datos.1 <- select (datos.1, genotipo,id.diseno, anio,ambiente,bloque,parcela, pot,planta, estado.feno, haun)
    datos.1a <- datos.1 %>% filter (estado.feno == "z2.1")
    datos.1b <- datos.1 %>% filter (estado.feno == "z3.1")
    datos.1c <- rbind (datos.1a, datos.1b)
    print (datos.1c) 
  })
  
  datos.haun <- do.call(rbind.data.frame,dt.haun)
  
  listapot.2 <- unique(dt.2$pot)
  
  
  dt.nfh <- lapply (listapot.2, function(filtro){
    datos.2 <- filter (dt.2, pot== filtro)
    datos.2b <- select(datos.2,genotipo,id.diseno, anio,ambiente,bloque,parcela,pot,planta)
    datos.2c <- datos.2b %>% 
      mutate (estado.feno = "flag.leaf")
    
    datos.2d <- datos.2c %>%
      mutate (haun = datos.2$nfh) 
    print (datos.2d) 
  })
  
  datos.nfh <- do.call(rbind.data.frame,dt.nfh)
  
  datos.NH <- rbind (datos.haun, datos.nfh)
  
  datos.NH <- datos.NH %>%  
    arrange (bloque) %>%
    arrange (parcela) %>%
    arrange (estado.feno)
  
  ambiente <- datos.NH$ambiente[1] 
  
  write.table (datos.NH, file =paste("./Data/procdata/outpout.HAUN/datos.NH",ambiente,".csv",sep=""),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               col.names = TRUE)
  return (datos.NH)
}


# esta es la funcion para el previo de la funcion GDmodel
# en un futuro FF tiene que ser mas general

FF <- function (dt = NULL, estado.feno = NULL) {

  e.f <-  estado.feno
  listapl.1 <- unique(dt$cod.pl)
  
  dt.nh.x <- lapply (listapl.1, function(filtro){
    datos.x1 <- filter (dt, cod.pl==filtro)
    
    x.gd <- datos.x1 %>%
      dplyr::filter (estado.feno == e.f) %>%
      dplyr::select (grados.dias)
    x.gd [1,]
  if (is.na (x.gd [1,]) == FALSE ){
    if(e.f == "z3.1") {
      datos.x2 <- datos.x1 %>%
        dplyr::filter(grados.dias <=  x.gd [1,])
    }
    if(e.f == "antesis") {
      datos.x1$estado.feno <- as.character(datos.x1$estado.feno)
      etq <- unique (datos.x1$estado.feno)
      etq.1 <- etq [!is.na(etq )]
      etqx <- which(etq.1 == "antesis")
      
      etqx.1 <- etq.1 [etqx -1]
      x.gd.etqx.1 <- datos.x1 %>%
        dplyr::filter (estado.feno == etqx.1) %>%
        dplyr::select (grados.dias)
      x.gd.etqx.1 [1,]
      
      datos.x2 <- datos.x1 %>%
        dplyr::filter (   grados.dias >= x.gd.etqx.1 [1,]    )%>%
        dplyr::filter (   grados.dias <=  x.gd [1,])
    }

    print (datos.x2) 
  }
  })
  
  
  datos.nh.ef <- do.call(rbind.data.frame,  dt.nh.x)

 # unique (datos.nh.ef$estado.feno)
  #datos.nh.ef$estado.feno <- as.character(datos.nh.ef$estado.feno)
 # etq <- unique (datos.nh.ef$estado.feno)
  #etq.1 <- etq [!is.na(etq )]
 # etq.2 <- str_c(etq.1, collapse = "_")
  
  #datos.nh.ef <- datos.nh.ef %>%
                  #dplyr::mutate (sub.est = etq.2) %>%
                  #dplyr::select (genotipo, id.diseno, anio, 
                   #ambiente,sub.est, everything())
    
  return (datos.nh.ef)
  #write.table (datos.nh.ef, file =paste("./Data/procdata/outpout.HAUN/datos.nfh",ambiente, e.f,".csv",sep=""),
              # append = FALSE, quote = TRUE, sep = ",",
               #eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               #col.names = TRUE)
}


pre.fil <- function (dt = NULL, subfase = NULL) {
  dir.create (file.path ("Data"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata", "outpout.HAUN"), showWarnings = FALSE)
  
  listapl.1 <- unique(dt$cod.pl)
  
  dt.nfh.1 <- lapply (listapl.1, function(filtro){
    datos.x1 <- filter (dt, cod.pl==  filtro)
    datos.x2 <- select (datos.x1, genotipo,id.diseno, anio,ambiente,bloque,parcela, pot,cod.pl,planta, grados.dias,estado.feno, haun)

    if (is.na (datos.x2[1,"haun"]) == FALSE ) {
    if (is.na (datos.x2[3,"haun"]) == FALSE) {
      if (length(na.omit(datos.x2$haun)) > 2) {
        datos.x2.a <- datos.x2
        #datos.x2.a
        if (length(na.omit(datos.x2.a$haun)) > 3) {
          datos.x3 <-   datos.x2.a %>% filter (haun < max(na.omit(haun)))
        } else {
          datos.x3 <- datos.x2.a
        }
      }
      print (datos.x3) 
     }
     }
  
  })
  
  datos.nfh.1 <- do.call(rbind.data.frame, dt.nfh.1)
  datos.nfh.1 <-  datos.nfh.1 %>%
                  dplyr::mutate (sub.est = subfase) %>%
                  dplyr::select (genotipo, id.diseno, anio, 
                          ambiente,sub.est, everything())
  
  ambiente <-  datos.nfh.1$ambiente[1]
  
  write.table (datos.nfh.1, file =paste("./Data/procdata/outpout.HAUN/datos.nfh",ambiente,subfase,".csv",sep=""),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               col.names = TRUE)
  
  return(datos.nfh.1)
}


# esta es la funcion que calcula el filocron

GD.model.3 <- function(x = NULL, subfase = NULL) {
  
  dir.create (file.path ("Data"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata"), showWarnings = FALSE)
  dir.create (file.path ("Data", "procdata", "outpout.HAUN"), showWarnings = FALSE)
  
  listapot.1 <- unique(x$pot)
  dt.nfh.1 <- lapply (listapot.1, function(filtro){
    datos.x1 <- filter (x, pot== filtro)
    NH  <- datos.x1  %>% select(haun)
    NH  <- NH[,1]
    GD  <- datos.x1 %>% select(grados.dias)
    GD  <- GD[,1]
    mod1 <- lm ( NH ~  GD, na.action=na.exclude)
    if (is.na  (coef (mod1)[2]) == FALSE) {
      if (summary (mod1)$r.squared < 0.9999 ) {
        if (is.nan (summary(mod1)$adj.r.square) == FALSE) {
          if (summary(mod1)$adj.r.square > 0) {
            
            ss <- rbind (coef (mod1))
            colnames (ss) <- c("a","b")
            
            
            plot (GD, NH, pch=16, col="black",
                  main =paste("pot_",unique(datos.x1$pot)))
            abline (mod1,  col="red", lwd=2)
            res <- signif (residuals (mod1), 5)
            pre <- predict (mod1) # plot distances between points and the regression line
            segments ( GD, NH,GD, pre, col="red")
            #etiq (GD, NH, res, cex=0.7)
            rsq <- summary(mod1)$r.squared
            rsq.adj <- summary(mod1)$adj.r.square
            
            ss <- data.frame (ss) %>%
              mutate (filocron = 1/b)
            bloque  <- unique (datos.x1$bloque)
            parcela <- unique(datos.x1$parcela)
            genotipo <- unique(datos.x1$genotipo)
            id.diseno <- unique(datos.x1$id.diseno)
            pot <- unique(datos.x1$pot)
            ambiente <- unique(datos.x1$ambiente)
            ss <- cbind(ss,rsq, rsq.adj)
            #class(ss[,1])
            ss.1 <- data.frame (ambiente,genotipo,id.diseno,bloque,parcela,pot,ss)
            print(ss.1)
          }
        }
      }
    }
  })
  
  
  datos.filocron <- do.call(rbind.data.frame, dt.nfh.1)
  envx <-  datos.filocron$ambiente[1]
  env.1 <- str_c (envx)
  env.1a <- str_split (envx, "E")
  #class( env.1a)
  etv <- env.1a [[1]][2]
  subfase <- subfase
  datos.filocron <-  datos.filocron %>%
    dplyr::mutate (env =  etv)%>%
    dplyr::mutate (subfase = subfase)%>% 
    dplyr::select (ambiente, env,	genotipo, subfase	,id.diseno,
                   bloque,	parcela,	pot, everything())
  
  write.table (datos.filocron , file =paste("./Data/procdata/outpout.HAUN/datos.filocron",envx,subfase,".csv",sep=""),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = FALSE,
               col.names = TRUE)
  return  (datos.filocron)
}





