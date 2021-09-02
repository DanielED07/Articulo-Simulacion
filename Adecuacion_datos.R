setwd("./Articulo")
rm(list=ls())
#-----------------------#
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(sp)
library(sf)
library(geoR)
library(gstat)
library(rgdal)
library(raster)
library(ape)
#-----------------------#
#-------------------------------------------#
# Leamos los datos y el contorno de brazil.
#-------------------------------------------#
load("viento.RData")

plot(viento[,c(2,3)])

cont_bzl <- shapefile("brazil_administrative_boundaries_national_polygon.shp")
cont <- cont_bzl@polygons[[2]]@Polygons[[1]]@coords

brazil_contorno <-SpatialPolygons(list(Polygons(list(Polygon( cont )), "x")))
brazil_contorno <- st_as_sf(brazil_contorno)
st_crs(brazil_contorno) = 4326

brazil_espacial <- st_as_sf(viento[,-1],coords = c("Longitude","Latitude"), crs = 4326, agr = "constant")
brazil_espacial$Value <- as.numeric(brazil_espacial$Value)

ggplot() + 
  geom_sf(data = brazil_contorno) + 
  geom_sf(data = brazil_espacial,aes(color = Value), size = 2)+
  scale_color_viridis_c(option = "C")+ 
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank())

#---------------------------------------------------------------------#
# Proyectemos los datos y el contorno de brazil a coordenadas planas.
#---------------------------------------------------------------------#
# 29101

brazil_contorno <- st_transform(brazil_contorno,29101)
brazil_espacial <- st_transform(brazil_espacial,29101)

ggplot() + 
  geom_sf(data = brazil_contorno) + 
  geom_sf(data = brazil_espacial,aes(color = Value), size = 2)+
  scale_color_viridis_c(option = "C")+coord_sf(datum=st_crs(29101))+ 
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank())

#--------#
# E.D.A
#--------#
#---------------------------------#
# Datos con coordenadas separadas
#---------------------------------#

datos <- cbind("viento"=brazil_espacial$Value,st_coordinates(brazil_espacial)[,c(1,2)])
datos <- na.omit(datos)
datos <- as.data.frame(datos)
coordinates(datos) <- ~X+Y
contorno <- st_coordinates(brazil_contorno)[,c(1,2)]

#----------------------------------#
# Transformar datos a formato geoR
#----------------------------------#

viento_geo <- as.geodata(datos,coords=2:3,var=1,borders=T)
viento_geo$borders <- contorno

plot.geodata(viento_geo)

#----------------------------------#
# Distancimxima y semivariograma
#----------------------------------#
max_dist <-max(dist(datos[,-1]))
# 4374773


vbin <- variogram(viento ~ 1, locations = coordinates(datos),datos)
plot(vbin)
# Fit Bin
vbin.fit <- fit.variogram(vbin,vgm(c("Gau", "Sph", "Mat", "Exp")), fit.kappa = TRUE)
plot(vbin,vbin.fit,main="Ajuste Semivariograma Bin")

#----------------------------------#
library(moments)
library(forecast)
library(rcompanion)
library(HDInterval)
library(nortest)
library(normtest)

test <- na.omit(datos$viento)
P <- ecdf(test) # función de distribución acumulada
skewness(test)


n_t <- blom(test) # Transformación Normal Scores
xs_nt <- cbind(test,n_t) # Real y transformado
xs <- inverseCDF(pnorm(n_t),P)
View(cbind(xs_nt,xs))
sum(abs(test-xs)) # Sesgo < 1% muy pequeño

ggdensity(n_t)
skewness(n_t)
ggqqplot(n_t)

shapiro.test(n_t) # Shapiro test - NO Rechaza a un nivel alpha = 0.03
ad.test(n_t) # Anderson Darling  - NO Rechaza a un nivel alpha = 0.05
jb.norm.test(n_t) #J arque Bera - NO Rechaza a un nivel alpha = 0.05
pearson.test(n_t) # Pearson      - Rechaza
cvm.test(n_t) # Cramer-Von Mises - NO Rechaza a un nivel alpha = 0.05
sf.test(n_t) # Shapiro-Francia - NO Rechaza a un nivel alpha = 0.05
lillie.test(n_t) # Kolmogorov-Smirnov - Rechaza
agostino.test(n_t) # D´Agostino - NO Rechaza a un nivel alpha = 0.05

# 6 de las 8 pruebas no rechazan la  hipotesis de normalidad
# La transformación funciona y da fe de normalidad.


