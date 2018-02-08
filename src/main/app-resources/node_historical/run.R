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
# Version:        2018-02-06

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
## create a date tag to include in output filenames
app.date = format(Sys.time(), "%Y%m%d_%H%M")

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
  
  rciop.log("DEBUG", " *** hypeapps-historical *** TEP hydrological modelling applications ***", "/node_historical/run.R")
  rciop.log("DEBUG", " rciop library loaded", "/node_historical/run.R")
  
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

## open application logfile
logFile=appLogOpen(appName = app.name,tmpDir = getwd(),appDate = app.date,prefix="000")

#################################################################################
## 2 - Application user inputs
## ------------------------------------------------------------------------------
## application input parameters
app.input <- getHypeAppInput(appName = app.name)

if(app.sys=="tep"){rciop.log("DEBUG", paste(" hypeapps inputs and parameters read"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "Inputs and parameters read",fileConn = logFile$fileConn)

#################################################################################
## 3 - Application setup
## ------------------------------------------------------------------------------
## Prepare basic model setup (static input files and hype model executable copied to working folder)
app.setup <- getHypeAppSetup(modelName = model.name,
                             modelBin  = model.bin,
                             tmpDir    = app.tmp_path,
                             appDir    = app.app_path,
                             appName   = app.name,
                             appInput  = app.input,
                             modelFilesURL = model.files.url,
                             forcingArchiveURL = forcing.archive.url,
                             stateFilesURL = state.files.url,
                             stateFilesIN = state.files)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("HypeApp setup read"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "HypeApp setup read",fileConn = logFile$fileConn)

## ------------------------------------------------------------------------------
## forcing data
model.forcing <- getModelForcing(appSetup   = app.setup,
                                 appInput   = app.input,
                                 dataSource = forcing.data.source)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("model forcing set"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "Model forcing prepared",fileConn = logFile$fileConn)

## ------------------------------------------------------------------------------
## get Xobs input file(s) from open catalogue
xobs.data <- getXobsData(appInput = app.input,
                         appSetup = app.setup)
if(app.sys=="tep"){rciop.log ("DEBUG", paste("xobs data downloaded from catalogue"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "xobs data (if any) downloaded from catalogue",fileConn = logFile$fileConn)

## ------------------------------------------------------------------------------
## read downloaded Xobs input file(s) - merge into one Xobs.txt in the model run folder
xobs.input <- readXobsData(appSetup = app.setup,
                         xobsData = xobs.data)
if(app.sys=="tep"){rciop.log ("DEBUG", paste("xobs data merged to model rundir"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "Xobs data (if any) merged into model directory",fileConn = logFile$fileConn)

## ------------------------------------------------------------------------------
## modify some model files based on input parameters
model.input <- updateModelInput(appSetup = app.setup, appInput = app.input, modelForcing = model.forcing, xobsInput = xobs.input)

if(app.sys=="tep"){rciop.log ("DEBUG", paste("model inputs modified"), "/node_historical/run.R")}
log.res=appLogWrite(logText = "model inputs modified",fileConn = logFile$fileConn)

#################################################################################
## 4 - Run application
## ------------------------------------------------------------------------------
##  run model
if(app.sys=="tep"){rciop.log ("DEBUG", " starting model run ...", "/node_historical/run.R")}
log.res=appLogWrite(logText = "starting model run ...",fileConn = logFile$fileConn)

model.run = system(command = app.setup$runCommand,intern = T,)

log.res=appLogWrite(logText = "..model run ready",fileConn = logFile$fileConn)
if(app.sys=="tep"){rciop.log ("DEBUG", paste(" ... model run ready, exit code: ",attr(model.run,"status"),sep=""), "/node_historical/run.R")}



#################################################################################
## 5 - Output
## ------------------------------------------------------------------------------
## post-process output data
#app.outdir <- prepareHypeAppsOutput(appSetup  = app.setup,
#                                    appInput = app.input,
#                                    modelInput = model.input,
#                                    modelForcing = model.forcing,
#                                    runRes = attr(model.run,"status"))
app.outfiles <- prepareHypeAppsOutput(appSetup  = app.setup,
                                    appInput = app.input,
                                    modelInput = model.input,
                                    modelForcing = model.forcing,
                                    runRes = attr(model.run,"status"),
                                    appDate = app.date)

if(length(app.outfiles)>1){
  app.outfiles=sort(app.outfiles,decreasing = F)
}

log.res=appLogWrite(logText = "HypeApp outputs prepared",fileConn = logFile$fileConn)

## ------------------------------------------------------------------------------
## publish postprocessed results (adding /* to avoid duplicate outputs to the user)
#if(app.sys=="tep"){
#  rciop.publish(path=paste(app.outdir,"/*",sep=""), recursive=FALSE, metalink=TRUE)
#}
for(k in 1:length(app.outfiles)){
  rciop.publish(path=app.outfiles[k], recursive=FALSE, metalink=TRUE)
}
log.res=appLogWrite(logText = "HypeApp outputs published",fileConn = logFile$fileConn)

## close and publish the logfile
log.file=appLogClose(appName = app.name,fileConn = logFile$fileConn)
if(app.sys=="tep"){
  rciop.publish(path=logFile$fileName, recursive=FALSE, metalink=TRUE)
}

#################################################################################
## 6 - End of workflow
## ------------------------------------------------------------------------------
## exit with appropriate status code

q(save="no", status = 0)
