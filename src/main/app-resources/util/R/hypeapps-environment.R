#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-environment.R
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
# hypeapps-environment.R: Paths and other settings for h-TEP hydrological modelling applications
# Author:                 David Gustafsson, SMHI
# Version:                2017-05-16
#

## set system flag, if not set
if(!exists("app.sys")){
  app.sys ="tep"
}

## load rciop library
if(app.sys=="tep"){
  library("rciop")
}

if(app.sys=="tep"){rciop.log ("DEBUG", paste("rciop package loaded"), "/util/R/hypeapps-environment.R")}


## Model settings 
model.name = "niger-hype"
if(app.sys=="win"){
  model.bin  = "HYPE_assimilation.exe"
}else{
#  model.bin  = "hype_assimilation"
  model.bin  = "hype-4.12.0"
}
if(app.sys=="tep"){rciop.log ("DEBUG", paste("hype binary set to:", model.bin,sep="" ), "/util/R/hypeapps-environment.R")}

## Forcing data source
if(app.sys=="win"){forcing.data.source =  "local"}
if(app.sys=="tep"){forcing.data.source =  "hydro-smhi"}
if(app.sys=="tep"){rciop.log ("DEBUG", paste(" forcing data soruce set to:", forcing.data.source,sep="" ), "/util/R/hypeapps-environment.R")}


## folder settings
if(app.sys=="tep"){
  app.app_path = Sys.getenv("_CIOP_APPLICATION_PATH")
  app.tmp_path = TMPDIR
}else if(app.sys=="win"){
  app.app_path = paste(getwd(),'application',sep="/")
  app.tmp_path = paste(getwd(),'tmp',sep="/")
}else{
  app.app_path = ""
  app.tmp_path = ""
}
if(app.sys=="tep"){rciop.log ("DEBUG", paste("application paths set"), "/util/R/hypeapps-environment.R")}

## Flag that this file has been sourced
app.envset=T
