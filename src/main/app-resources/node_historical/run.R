#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-historical/src/main/app-resources/node_historical/run.R

# Copyright 2017 SMHI
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

# Application 1: "Niger-HYPE historical period" (hypeapps-historial)
# Author:         David Gustafsson, SMHI
# Version:        2017-05-16

# Remaining issues: 
#
# PLOTTING is disabled due to (a) lack of X11 server on sandbox, or 
#   (b) conflict on libjeg-version between r-cairo and python pil (needed for gdal and r-rgdal) 
#
# solution: waiting for Terradues response to support issue #5720
#
# XOBS file integration still to be developed


# Workflow overview:
# ------------------
# 1 Initialization          (load environmental variables, libraries, utility functions, etc)
# 2 Application inputs      (read all user inputs)
# 3 Application setup       (setup working folders, copy model parameters, generate forcing data, etc)
# 4 Run application         (run the model)
# 5 Output                  (pepare and publish output data)
# 6 End of workflow

#################################################################################
## 1 - Initialization
## ------------------------------------------------------------------------------
## set application name
app.name = "historical"
## ------------------------------------------------------------------------------
## flag which environment is used, if not set
if(!exists("app.sys")){
  app.sys ="tep"
}
## ------------------------------------------------------------------------------
## load rciop package and set working directory to TMPDIR when running on TEP 
if(app.sys=="tep"){
  library("rciop")
  
  rciop.log ("DEBUG", " *** hypeapps-historical *** TEP hydrological modelling applications ***", "/node_historical/run.R")
  rciop.log ("DEBUG", " rciop library loaded", "/node_historical/run.R")
  
  setwd(TMPDIR)
  rciop.log("DEBUG", paste(" R session working directory set to ",TMPDIR,sep=""), "/node_historical/run.R")
}
## ------------------------------------------------------------------------------
## load hypeapps environment and additional R utility functions
if(app.sys=="tep"){
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/hypeapps-environment.R",sep="/"))
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/hypeapps-utils.R", sep="/"))
  rciop.log ("DEBUG", paste(" libraries loaded and utilities sourced"), "/node_historical/run.R")
}else if(app.sys=="win"){
  source("application/util/R/hypeapps-environment.R")  
  source("application/util/R/hypeapps-utils.R")
}
#################################################################################
## 2 - Application user inputs
## ------------------------------------------------------------------------------
## application input parameters
app.inputs <- getHypeAppInputs(appName = app.name)

if(app.sys=="tep"){rciop.log ("DEBUG", paste(" hypeapps inputs and parameters read"), "/node_historical/run.R")}

#################################################################################
## 3 - Application setup
## ------------------------------------------------------------------------------
## Prepare basic model setup (static input files and hype model executable copied to working folder)
app.setup <- getHypeAppSetup(modelName = model.name,
                             modelBin  = model.bin,
                             tmpDir    = app.tmp_path,
                             appDir    = app.app_path,
                             appName   = app.name)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("HypeApp setup read"), "/node_historical/run.R")}

## ------------------------------------------------------------------------------
## forcing data
model.forcing <- getModelForcing(appSetup   = app.setup,
                                 appInputs  = app.inputs,
                                 dataSource = forcing.data.source)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("model forcing set"), "/node_historical/run.R")}

## ------------------------------------------------------------------------------
## modify some model files based on input parameters
model.input <- updateModelInput(appSetup = app.setup, appInputs = app.inputs)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("model inputs modified"), "/node_historical/run.R")}

#################################################################################
## 4 - Run application
## ------------------------------------------------------------------------------
##  run model
if(model.input==0){
  if(app.sys=="tep"){rciop.log ("DEBUG", " starting model run ...", "/node_historical/run.R")}
  model.run = system(command = app.setup$runCommand,intern = T)
}
if(app.sys=="tep"){rciop.log ("DEBUG", " ... model run ready", "/node_historical/run.R")}

#################################################################################
## 5 - Output
## ------------------------------------------------------------------------------
## post-process output data
app.outdir <- prepareHypeAppsOutput(appSetup  = app.setup,appInput = app.inputs, runRes = model.run)

## ------------------------------------------------------------------------------
## publish postprocessed results
if(app.sys=="tep"){
  rciop.publish(path=app.outdir, recursive=TRUE, metalink=TRUE)
}
#################################################################################
## 6 - End of workflow
## ------------------------------------------------------------------------------
## exit with appropriate status code
q(save="yes", status = 0)
