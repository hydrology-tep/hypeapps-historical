#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-additional-utils.R
#
# Copyright 2016-2017 SMHI
#
# This file is part of H-TEP Hydrological Modelling Application, which is open source 
# and distributed under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version. The Hydrology TEP Hydrological Modelling Application is distributed 
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU 
# General Public License for more details. You should have received a copy of the Lesser 
# GNU General Public License along with the Hydrology TEP Hydrological Modelling Application. 
# If not, see <http://www.gnu.org/licenses/>.
#
# hypeapps-additional-utils.R: Additional supporting utility functions for offline preparation of model data, etc.
# Author:                      David Gustafsson, SMHI
# Version:                     2017-05-15





# small script to make the hype2csv from standard WHIST generated HYPE model shapefiles
library(sp)
library(rgdal)
library(rgeos)
shape2csvkey<-function(shpdf){
  csvdf=data.frame("SUBID"=shpdf@data$SUBID)
  csvdf$CENTERX=shpdf@data$CENTERX
  csvdf$CENTERY=shpdf@data$CENTERY
  csvdf$POURX=shpdf@data$POURX
  csvdf$POURY=shpdf@data$POURY
  csvdf$LON=round(coordinates(shpdf)[,1],digits=6)
  csvdf$LAT=round(coordinates(shpdf)[,2],digits=6)
  csvdf$ELEV=shpdf@data$elevation
  return(csvdf)
}
writeCsvKey<-function(csvdf,filename="hype2csv.txt"){
  write.table(x=csvdf,file = filename,quote = F,sep = "\t",row.names = F,col.names = T)
  return(0)
}

# Example usage:

# read shapefile
shapefile = readOGR(dsn="/hypeapps-historical/src/main/app-resources/model/my-hype/shapefiles",layer="my-hype")

# make csv key
csvdf = shape2csvkey(shapefile)

# write hype2csv.txt
writeCsvKey(csvdf = csvdf,filename = "/hypeapps-historical/src/main/app-resources/model/my-hype/hype2csv/my-hype2csv.txt")
