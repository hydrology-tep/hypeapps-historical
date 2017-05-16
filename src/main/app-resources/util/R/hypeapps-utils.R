#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-utils.R
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
# hypeapps-utils.R: R tools for the HTEP hydrological modelling application 
# Author:           David Gustafsson, SMHI
# Version:          2017-05-16
#

# Remaining issues: 
#
# PLOTTING is disabled due to (a) lack of X11 server on sandbox, or 
#   (b) conflict on libjeg-version between r-cairo and python pil (needed for gdal and r-rgdal) 
#
# solution: waiting for Terradues response to support issue #5720

## --------------------------------------------------------------------------------
## initial settings
## --------------------------------------------------------------------------------
# set system flag if not set
if(!exists("app.sys")){
  app.sys ="tep"
}

# load rciop library
if(app.sys=="tep"){
  library("rciop")
}

# source hypeapps environment file, if needed
if(!exists("app.envset")){
  if(app.sys=="tep"){
    source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-environment.R",sep="/"))
  }else if(app.sys=="win"){
    source("application/util/R/hypeapps-environment.R")  
  }
  if(app.sys=="tep"){rciop.log ("DEBUG", paste("hypeapps-environment.R sourced"), "/util/R/hypeapps-utils.R")}
}

# source csv and hype-utils.R
if(app.sys=="tep"){
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-hype-utils.R",sep="/"))
  rciop.log ("DEBUG", paste("hypeapps-hype-utils.R sourced"), "/util/R/hypeapps-utils.R")
  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-csv-utils.R",sep="/"))
  rciop.log ("DEBUG", paste("hypeapps-csv-utils.R sourced"), "/util/R/hypeapps-utils.R")

# removing plotting until Support #5720 is solved by Terradue
#  source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-plotting-utils.R",sep="/"))
#  rciop.log ("DEBUG", paste("hypeapps-plotting-utils.R sourced"), "/util/R/hypeapps-utils.R")
  if(app.name=="eodata"){
    source(paste(Sys.getenv("_CIOP_APPLICATION_PATH"),"util/R/hypeapps-eo-utils.R",sep="/"))
    rciop.log ("DEBUG", paste("hypeapps-eo-utils.R sourced"), "/util/R/hypeapps-utils.R")
    }
}else if(app.sys=="win"){
  source(paste("application","util/R/hypeapps-hype-utils.R",sep="/"))
  source(paste("application","util/R/hypeapps-csv-utils.R",sep="/"))
# removing plotting until Support #5720 is solved by Terradue
#  source(paste("application","util/R/hypeapps-plotting-utils.R",sep="/"))
  if(app.name=="eodata"){source(paste("application","util/R/hypeapps-eo-utils.R",sep="/"))}
}else{
}

## --------------------------------------------------------------------------------
## functions
## -------------------------------------------------------------------------------
## getHypeAppInputs - function to load user parameter inputs, depending on application
getHypeAppInputs<-function(appName){
  
  # get input parameters for the different hypeapps applications
  if(appName=="historical"){
    
    if(app.sys=="tep"){
      # get parameters with the rciop function when running on the TEP system
      bdate     <- rciop.getparam("bdate")
      cdate     <- rciop.getparam("cdate")
      edate     <- rciop.getparam("edate")
      outvars   <- rciop.getparam("variables")
      outbasins <- rciop.getparam("basins")
    }else if(app.sys=="win"){
      # set to default values, if not set
      # set some test values for development on windows:
      if(!exists("win.bdate")|
         !exists("win.cdate")|
         !exists("win.edate")|
         !exists("win.outvars")|
         !exists("win.outbasins")){
        bdate     <- "2006-01-01"
        cdate     <- "2016-01-01"
        edate     <- "2016-12-01"
        outvars   <- "cout,cprc,ctmp,evap,epot"
        outbasins <- "37"
      }else{
        bdate     <- win.bdate
        cdate     <- win.cdate
        edate     <- win.edate
        outvars   <- win.outvars
        outbasins <-win.outbasins
      }
    }else{
      bdate     <- NULL
      cdate     <- NULL
      edate     <- NULL
      outvars   <- NULL
      outbasins <- NULL
      print("WARNING: hypeapps.sys not set, allowed values are 'tep' or 'win' ")
    }
    appInput=list("bdate"     = bdate,
                  "cdate"     = cdate,
                  "edate"     = edate,
                  "outvars"   = outvars,
                  "outbasins" = outbasins)
    
  }else if(appName=="forecast"){
    
  }else if(appName=="returnperiod"){
    
  }else if(appName=="eodata"){
    
    if(app.sys=="tep"){
      # get parameters and data files
      wlData  <- as.character(rciop.getparam("wlData"))   # Water level data
      wlSubid <- as.character(rciop.getparam("wlSubid"))  # water level data subbasin identifiers
      
      flData  <- as.character(rciop.getparam("flData"))   # flood map data
      flSubid <- as.character(rciop.getparam("flSubid"))  # flood map data subbasin identifiers
      
      # parse the water level inputs
      wlDataNum=0
      for(i in 1:length(wlData)){
        if(nchar(wlData[i])>1){
          # update the number of wldata inputs
          wlDataNum=wlDataNum+1
          
          # extract the URL from wlData input string:
          if(wlDataNum==1){
            wlDataURL=strsplit(wlData[1],split = "&")[[1]][1]
          }else{
            wlDataURL=c(wlDataURL,strsplit(wlData[1],split = "&")[[1]][1])
          }
        }
      }
      if(wlDataNum>0){
        wlDataSubid=as.integer(strsplit(wlSubid,","))
        if(length(wlDataSubid)==wlDataNum & length(wlDataURL)==wlDataNum){
          wlDataInput=T
        }else{
          wlDataNum=0
          wlDataInput=F
          wlDataURL=NULL
          wlDataSubid=NULL
        }
      }else{
        wlDataNum=0
        wlDataInput=F
        wlDataURL=NULL
        wlDataSubid=NULL
      }
      
      #      # parse the flood mapping inputs
      #      if(nchar(flData)>1){
      #        flDataInput=T
      #        
      #        # extract the URL from wlData input string:
      #        flDataURL=strsplit(flData[1],split = "&")[[1]][1]
      #        
      #        # subid list
      #        flDataSubid = flSubid
      #        
      #      }else{
      flDataNum=0
      flDataInput=F
      flDataURL=NULL
      flDataSubid=NULL
      #      }
    }else{
      wlDataNum    <- 0
      wlDataInput  <- F
      wlDataURL    <- NULL
      wlDataSubid  <- NULL
      flDataNum    <- 0
      flDataInput  <- F
      flDataURL    <- NULL
      flDataSubid  <- NULL
    }
    appInput=list("wlDataNum"=wlDataNum,
                  "wlDataInput"=wlDataInput,
                  "wlDataURL"=wlDataURL,
                  "wlDataSubid"=wlDataSubid,
                  "flDataNum"=flDataNum,
                  "flDataInput"=flDataInput,
                  "flDataURL"=flDataURL,
                  "flDataSubid"=flDataSubid)
    
  }else{
    appInputs=list("appName"=NULL)
  }
  return(appInput)
}

## -------------------------------------------------------------------------------
## prepare work directories and copy basic model files
getHypeAppSetup<-function(modelName,modelBin,tmpDir,appDir,appName){
  
  ## model input files (common for all applications)
  # model files source directory
  modelFilesAppDir=paste(appDir,'model',modelName,sep="/")
  
  # model files run directory
  modelFilesRunDir=paste(tmpDir,'model',modelName,sep="/")
  dir.create(modelFilesRunDir,recursive = T,showWarnings = F)
  
  # model files results directory (not necessary for "eodata" and "stats")
  if(appName=="historical"){
    modelResDir=paste(modelFilesRunDir,'results',sep="/")
    dir.create(modelResDir,recursive = T,showWarnings = F)
    
  }else if(appName=="forecast"){
    modelResDir=paste(modelFilesRunDir,'results/hindcast',sep="/")
    modelResDir=c(paste(modelFilesRunDir,'results/forecast',sep="/"))
    for(i in 1:2){
      dir.create(modelResDir[i],recursive = T,showWarnings = F)
    }
  }else{
    modelResDir=NULL
  }
  
  # eodata run directory
  if(appName=="eodata"){
    eodataResDir=paste(modelFilesRunDir,"eodata",sep="/")
    dir.create(eodataResDir,recursive = T,showWarnings = F)
  }else{
    eodataResDir=NULL
  }
  
  # copy model files to working directory
  fileNames=c("par.txt",
              "GeoData.txt",
              "GeoClass.txt",
              "BranchData.txt",
              "FloodData.txt",
              "LakeData.txt")
  
  if(appName=="historical"|appName=="eodata"){
    fileNames=c(fileNames,"info-historical.txt")
  }else if(appName=="forecast"){
    fileNames=c(fileNames,"info-hindcast.txt","info-forecast.txt")
  }
  
  for(i in 1:length(fileNames)){
    if(app.sys=="tep"){
      res <- rciop.copy(paste(modelFilesAppDir,fileNames[i],sep="/"), modelFilesRunDir, uncompress=TRUE)
    }else{
      file.copy(from=paste(modelFilesAppDir,fileNames[i],sep="/"),
                to =paste(modelFilesRunDir,fileNames[i],sep="/"),
                overwrite = T)
    }
  }
  
  ## model binary file (in application folder)
  if(appName=="historical"|appName=="forecast"){
    # model binary source file
    modelBinaryFile=paste(appDir,'util/bin',modelBin,sep="/")
    
    # command line to run the model "hype rundir"
    sysCommand = paste(paste(modelBinaryFile,modelFilesRunDir,sep=" "),"/",sep="")
    
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("HYPE RUN system command=",sysCommand,sep=""), "/util/R/hypeapps-utils.R")}
  }else{
    sysCommand=NULL
  }

  ## return list with application setup
  appSetup = list("runDir"=modelFilesRunDir,
                  "resDir"=modelResDir,
                  "runCommand"=sysCommand,
                  "eodataResDir"=eodataResDir,
                  "tmpDir"=tmpDir,
                  "appDir"=appDir,
                  "appName"=appName,
                  "modelName"=modelName,
                  "modelBin"=modelBin,
                  "hype2csv"=paste(modelFilesAppDir,"hype2csv",paste(modelName,"2csv.txt",sep=""),sep="/"))
  return(appSetup)
}

## -------------------------------------------------------------------------------
## get eo data from open catalogue
getEodata<-function(appInput,appSetup){
  
  #  appInput=list("wlDataNum"=wlDataNum
  #                "wlDataInput"=wlDataInput,
  #                "wlDataURL"=wlDataURL,
  #                "wlDataSubid"=wlDataSubid,
  #                "flDataNum"=flDataNum  
  #                "flDataInput"=flDataInput,
  #                "flDataURL"=flDataURL,
  #                "flDataSubid"=flDataSubid)
  
  #  appSetup = list("runDir"=modelFilesRunDir,
  #                  "resDir"=modelResDir,
  #                  "runCommand"=sysCommand,
  #                  "eodataResDir"=eodataResDir,
  #                  "tmpdir"=runDir)
  
  # loop over wlDataInput
  nDownLoad=0
  if(appInput$wlDataNum>0){
    for(i in 1:appInput$wlDataNum){
      # download data
      sysCmd=paste("opensearch-client",appInput$wlDataURL[i],"enclosure | ciop-copy -s -U -O /tmp/ -", sep=" ")
      wlFile=system(command = sysCmd,intern = T)
      if(file.exists(wlFile)){
        nDownLoad=nDownLoad+1
        if(nDownLoad==1){
          wlFiles=wlFile
          wlSubid=appInput$wlDataSubid[i]
        }else{
          wlFiles=c(wlFiles,wlFile)
          wlSubid=c(wlSubid,appInput$wlDataSubid[i])
        }
      }
    }
    eoData=list("wlDataNum"=nDownLoad,
                "wlDataFile"=wlFiles,
                "wlDataSubid"=wlSubid)
  }else{
    eoData=list("wlDataNum"=0,
                "wlDataFile"=NULL,
                "wlDataSubid"=NULL)
  }
  #      #wlData="http://sb-10-15-36-31.hydro.terradue.int/sbws/production/run/water-levels/0000081-170411124531167-oozie-oozi-W/products/search?uid=0000081-170411124531167-oozie-oozi-W/outputs/lakes_summary_multi_ATK_Kainji_L2.csv&format=atom"
  #      if(nchar(wlData)>1){
  #        # download data
  #        sysCmd=paste("opensearch-client",wlDataURL,"enclosure | ciop-copy -s -U -O /tmp/ -", sep=" ")
  #        wlFile=system(command = sysCmd,intern = T)
  
  return(eoData)
}

## -------------------------------------------------------------------------------
## Read eodata and prepare data on the Xobs format
readEodata<-function(appSetup,eoData){
  for(i in 1:eoData$wlDataNum){
    xobs = readWaterLevelCSV(fname=eoData$wlDataFile[i],
                             subid=eoData$wlDataSubid[i],
                             vars="wstr",
                             lwl0=0,
                             t1=NULL,
                             t2=NULL,
                             modelSHP=NULL,
                             limData=T,
                             lowLim=-0.9999,
                             highLim=9999)
    if(i==1){
      xobsData=xobs
    }else{
      #      xobsData= mergeXobs
    }
  }
  return(xobsData)
}

## -------------------------------------------------------------------------------
## write eo data in Xobs format
writeEodata<-function(appSetup,xobsData){
  # output directory
  outputDir=paste(appSetup$tmpdir,"output",sep="/")
  dir.create(outputDir,recursive = T,showWarnings = F)
  
  # output filename
  xobsFile=paste(outputDir,"Xobs-eodata.txt",sep="/")
  
  # write to file
  outres = WriteXobs(xobsData, filename = xobsFile)
  
  # return output list
  appOutput = list("outDir"=outputDir,"files"=xobsFile)
  return(appOutput)
}


## -------------------------------------------------------------------------------
## prepare model forcing data
getModelForcing<-function(appSetup,appInputs,dataSource="local",hindcast=T){
  
  # Forcing data for historical simulation
  if(appSetup$appName=="historical"){
    
    # A. get data from local file
    if(dataSource=="local"){
      forcingFile=paste(appSetup$appDir,"gfd",appSetup$modelName,"hindcast/archive.zip",sep="/")
      unzip(zipfile = forcingFile,overwrite = T,exdir = appSetup$runDir)
      if(app.sys=="tep"){rciop.log ("DEBUG", paste("forcingFile from local file = ",forcingFile,sep=""), "/util/R/hypeapps-utils.R")}
      return(0)
    }
    
    # B. get archive data from catalogue.terradue.com/hydro-smhi/
    if(dataSource=="hydro-smhi"){
      # download archive file to tmpDir
      sysCmd=paste("curl -o ", appSetup$tmpDir,"/archive.zip https://store.terradue.com//smhi/gfd/niger-hype/hindcast/files/v1/archive.zip",sep="")
      if(app.sys=="tep"){rciop.log ("DEBUG", paste(" trying command >> ",sysCmd,sep=""), "/util/R/hypeapps-utils.R")}
      a=system(sysCmd,intern=T)
      #unzip to runDir
      forcingFile=paste(appSetup$tmpDir,"archive.zip",sep="/")
      if(app.sys=="tep"){rciop.log ("DEBUG", paste("forcingFile from https://catalogue.terradue.com/hydro-smhi/ = ",forcingFile,sep=""), "/util/R/hypeapps-utils.R")}
      unzip(zipfile = forcingFile,overwrite = T,exdir = appSetup$runDir)
      return(0)
    }
    
    # C. get data from ftp.smhi.se
    if(dataSource=="ftp-smhi"){
      
    }
  }
  
}
## -------------------------------------------------------------------------------
## modify some model input files
updateModelInput<-function(appSetup=NULL,appInputs=NULL,hindcast=NULL){
  
  if(appSetup$appName=="historical"){
    
    ## INFO.TXT ##
    # read info.txt template
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("appSetup$runDir = ",appSetup$runDir,sep=""), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("TMPDIR=",appSetup$tmpDir,sep=""), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("dir(TMPDIR)=",dir(appSetup$tmpDir),sep=" "), "/util/R/hypeapps-utils.R")}
    if(app.sys=="tep"){rciop.log ("DEBUG", paste("dir(appSetup$runDir)=",dir(appSetup$runDir),sep=" "), "/util/R/hypeapps-utils.R")}
    
    info=readInfo(paste(appSetup$runDir,"info-historical.txt",sep="/"))
    
    # modify info according to inputs
    # resultdir
    if(info$isResDir){
      info$info.lines[info$resultdir.lineNr]=paste(paste('resultdir', appSetup$resDir, sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('resultdir', appSetup$resDir, sep=" "),"/",sep=""))
    }
    # modeldir
    if(info$isModDir){
      info$info.lines[info$modeldir.lineNr]=paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep="")
    }else{
      info$info.lines=c(info$info.lines,paste(paste('modeldir', appSetup$runDir, sep=" "),"/",sep=""))
    }
    # bdate,cdate,edate
    info$info.lines[info$bdate.lineNr]=paste('bdate',as.character(appInputs$bdate),sep=" ")
    info$info.lines[info$cdate.lineNr]=paste('cdate',as.character(appInputs$cdate),sep=" ")
    info$info.lines[info$edate.lineNr]=paste('edate',as.character(appInputs$edate),sep=" ")
    # basinoutput
    info$info.lines[info$basinoutput_variable.lineNr]=paste('basinoutput variable',gsub(pattern=",",replacement = " ",appInputs$outvars),sep=" ")
    info$info.lines[info$basinoutput_subbasin.lineNr]=paste('basinoutput subbasin',gsub(pattern=",",replacement = " ",appInputs$outbasins),sep=" ")
    # timeoutput
    info$info.lines[info$timeoutput_variable.lineNr]=paste('timeoutput variable',gsub(pattern=",",replacement = " ",appInputs$outvars),sep=" ")
    
    # write info file
    return(writeInfo(info$info.lines,filenm = paste(appSetup$runDir,"info.txt",sep="/")))
    
  }else{
    return(-1)
  }
}

## -------------------------------------------------------------------------------
## prepare application outputs
prepareHypeAppsOutput<-function(appSetup=NULL,appInputs=NULL,runRes=NULL){
  
  ## Create folder for data to be published
  outDir = paste(appSetup$tmpDir,'output',sep="/")
  dir.create(outDir,recursive = T,showWarnings = F)
  
  ## Post-process requested outputs (copy some files...)
  if(appSetup$appName=="historical"){
    
    # copy log-file in model working folder
    hyssLogFile = dir(path = appSetup$runDir , pattern =".log")
    if(length(hyssLogFile)>0){
      if(app.sys=="tep"){
        res <- rciop.copy(paste(appSetup$runDir,hyssLogFile[1],sep="/"), outDir, uncompress=TRUE)
      }
    }
    
    # list files in model results folder
    timeFiles   = dir(path = appSetup$resDir , pattern ="time")
    mapFiles    = dir(path = appSetup$resDir , pattern ="map")
    subassFiles = dir(path = appSetup$resDir , pattern ="subass")
    simassFile  = dir(path = appSetup$resDir , pattern ="simass")
    allFiles    = dir(path = appSetup$resDir , pattern =".txt")
    
    # copy time files
    if(length(timeFiles)>0){
      for(i in 1:length(timeFiles)){
        if(app.sys=="tep"){rciop.copy(paste(appSetup$resDir,timeFiles[i],sep="/"), outDir, uncompress=TRUE)}
      }
    }
    # copy map files
    if(length(mapFiles)>0){
      for(i in 1:length(mapFiles)){
        if(app.sys=="tep"){rciop.copy(paste(appSetup$resDir,mapFiles[i],sep="/"), outDir, uncompress=TRUE)}
      }
    }
    # copy subass files
    if(length(subassFiles)>0){
      for(i in 1:length(subassFiles)){
        if(app.sys=="tep"){rciop.copy(paste(appSetup$resDir,subassFiles[i],sep="/"), outDir, uncompress=TRUE)}
      }
    }
    # copy simass files
    if(length(simassFile)>0){
      if(app.sys=="tep"){rciop.copy(paste(appSetup$resDir,simassFile[1],sep="/"), outDir, uncompress=TRUE)}
    }
    
    # basin outputfiles
    outbasins = strsplit(appInputs$outbasins,split = ",")[[1]]
    zeroString="0000000000000000000000000000000"
    basinFiles=NULL
    if(length(outbasins)>0){
      for(i in 1:length(outbasins)){
        outFile=paste(outbasins[i],".txt",sep="")
        ni=nchar(outFile)
        for(j in 1:length(allFiles)){
          nj=nchar(allFiles[j])
          if(nj>=ni){
            if(substr(allFiles[j],nj-ni+1,nj)==outFile){
              if(nj>ni){
                if(substr(allFiles[j],1,nj-ni)==substr(zeroString,1,nj-ni)){
                  if(app.sys=="tep"){
                    rciop.copy(paste(appSetup$resDir,allFiles[j],sep="/"), outDir, uncompress=TRUE)
                    basinFiles=c(basinFiles,allFiles[j])
                    }
                }
              }else{
                if(app.sys=="tep"){
                  rciop.copy(paste(appSetup$resDir,allFiles[j],sep="/"), outDir, uncompress=TRUE)
                  basinFiles=c(basinFiles,allFiles[j])
                }
              }
            }
          }
        }
      }
    }
    # transform basinoutput files to csv format
    if(length(basinFiles)>0){
      for(i in 1:length(basinFiles)){
        resCsv = basinfiles2csv(hypeFile=paste(appSetup$resDir,basinFiles[i],sep="/"),
                                csvFile=paste(outDir,paste(substr(basinFiles[i],1,nchar(basinFiles[j])-3),"csv",sep=""),sep="/"),
                                hype2csvFile=appSetup$hype2csv)
      }
    }
    
# removing plotting until Support #5720 is solved by Terradue
    
#     # finally, plot content of basinoutput files
#     if(length(basinFiles)>0){
#       for(i in 1:length(basinFiles)){
#         resPlot = plotBasinOutputFile(fileDir=appSetup$resDir,
#                                       plotDir=outDir,
#                                       basinFile=basinFiles[i],
#                                       modelName=appSetup$modelName,
#                                       graphScale=1.6, lgndScale = 0.98,lineWidth=2,
#                                       maxmarg=0.1)
#       }
#     }
    
    
    
    
  }
  
  
  
  ## return outDir
  return(outDir)
}


if(app.sys=="tep"){rciop.log ("DEBUG", paste("all functions sourced"), "/util/R/hypeapps-utils.R")}
