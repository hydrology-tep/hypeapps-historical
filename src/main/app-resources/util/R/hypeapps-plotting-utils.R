#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet
#
# /hypeapps-[appName]/src/main/app-resources/util/R/hypeapps-plotting-utils.R
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
# hypeapps-plotting-utils.R: Utilities to plot HYPE model dala for TEP Hydrology.
# Author:                    David Gustafsson, Jafet Andersson SMHI
# Version:                   2017-05-15

# -------------------------------------------------------------------
# dependancies
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# color schemes
# -------------------------------------------------------------------
## makeTransparent plotting color
.makeTransparent <- function(someColor, alpha=60) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, maxColorValue = 255)})
}

## function to make character vector with days from a posix vector
.dates2days<-function(dates){
  days=as.character(dates)
  days=substr(days,6,10)
  return(days)
}

# ---------------------------------------------------------------------
# (internal) plot hindcast+forecast river discgarhe with warning levels
# ---------------------------------------------------------------------
.plotHindForWL<-function(hindcast=NULL,
                         forecast=NULL,
                         obsdata=NULL,
                         annsim=NULL,
                         annobs=NULL,
                         wl.sim=NULL,
                         wl.obs=NULL,
                         wl.col=c("yellow","orange","red"),
                         maxsim=NULL,
                         maxobs=NULL,
                         maintitle=NULL,
                         hinddays=NULL,
                         maxmarg=0.10,
                         wl.alpha=120,
                         wl.alpha.obs=120,
                         hind.col="blue",
                         fore.col="red",
                         obs.col="black",
                         lineWidth=3,
                         graphScale=2,
                         lgndScale=1,
                         lgndDY=-0.2){
  
  ## initiate plot
  if(!is.null(hinddays)){
    nh=nrow(hindcast)
    dateAxis=c(hindcast$DATE[(nh-hinddays+1):nh],forecast$DATE)
    iniData=c(hindcast[(nh-hinddays+1):nh,2],forecast[,2])
  }else{
    hinddays=nrow(hindcast)
    nh=nrow(hindcast)
    dateAxis=c(hindcast$DATE,forecast$DATE)
    iniData=c(hindcast[,2],forecast[,2])
  }
  maxy=max(c(maxobs,maxsim,max(hindcast[,2]),max(forecast[,2]),max(obsdata[,2],na.rm = T)),na.rm = T)*(1+maxmarg)
  miny=0
  
  plot(dateAxis,iniData,type="n",ylim=c(miny,maxy),main=maintitle,xlab="",ylab="Discharge (m3/s)",xaxs="i",cex.lab=graphScale,cex.axis=graphScale,cex.main=graphScale,xpd=T)
  axis(side = 1,at=dateAxis-60*60*2,labels=F,tcl=-0.2)  # add daily tickmarks for simplicity
  
  ## warning levels as polygons with transparent colors
  if(!is.null(wl.sim)){
    polcol <- .makeTransparent(wl.col[1], wl.alpha)
    if(length(wl.col)>1){
      for(j in 2:length(wl.col)){
        polcol <- c(polcol,.makeTransparent(wl.col[j], wl.alpha))
      }
    }
    for(i in 1:length(wl.sim)){
      x=c(dateAxis[1],dateAxis[length(dateAxis)],dateAxis[length(dateAxis)],dateAxis[1])
      y=c(wl.sim[i],wl.sim[i])
      if(i==length(wl.sim)){
        y=c(y,maxy*2,maxy*2)
      }else{
        y=c(y,wl.sim[i+1],wl.sim[i+1])
      }
      polygon(x, y, col = polcol[i], border = NA)
    }
  }
  ## warning levels from observations
  if(!is.null(wl.obs)){
    polcol <- .makeTransparent(wl.col[1], wl.alpha.obs)
    if(length(wl.col)>1){
      for(j in 2:length(wl.col)){
        polcol <- c(polcol,.makeTransparent(wl.col[j], wl.alpha))
      }
    }
    for(i in 1:length(wl.obs)){
      lines(c(min(dateAxis),max(dateAxis)),c(1,1)*wl.obs[i],col=polcol[i],lty=1,lwd=lineWidth)
    }
  }    
  
  ## maximum ever observed and simulated
  if(!is.null(maxsim)){
    lines(c(min(dateAxis),max(dateAxis)),c(1,1)*maxsim,col=hind.col,lty=2,lwd=lineWidth)
  }
  if(!is.null(maxobs)){
    lines(c(min(dateAxis),max(dateAxis)),c(1,1)*maxobs,col=obs.col,lty=2,lwd=lineWidth)
  }
  
  ## regime ribbons
  days2plot = .dates2days(dateAxis)
  ind2plot  = match(days2plot,annsim$mean$day)
  
  # simualtion regime
  polcol <- .makeTransparent(hind.col, 30)
  mindata=cbind(dateAxis,annsim$minimum[ind2plot, 3])
  maxdata=cbind(dateAxis,annsim$maximum[ind2plot,3])
  polcoor.minmax <- rbind(mindata, maxdata[nrow(maxdata):1,])
  mindata=cbind(dateAxis,annsim$p25[ind2plot, 3])
  maxdata=cbind(dateAxis,annsim$p75[ind2plot,3])
  polcoor.p25p75 <- rbind(mindata, maxdata[nrow(maxdata):1,])
  polygon(polcoor.minmax[, 1], polcoor.minmax[, 2], col = polcol, border = NA)
  
  # hindcast
  lines(rbind(hindcast[(nh-hinddays+1):nh,],forecast[1,]),col=.makeTransparent(hind.col, 200),lwd=lineWidth,lty=1)
  # forecast
  lines(forecast,col=.makeTransparent(fore.col, 200),lwd=lineWidth,lty=1,)
  # observation
  if(!is.null(obsdata)){
    lines(obsdata,col=.makeTransparent(obs.col, 200),lwd=lineWidth,lty=1)
    #  points(obsdata,col=.makeTransparent("black", 200),pch=20)
  }
  
  # Legend 1: forecast, hindcast, maxima and observations (lines)
  legend(x=par("usr")[2],y=par("usr")[3],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth, yjust=0,pt.cex=graphScale*2,
         legend = c("10 day forecast","Current hindcast","Hindcast max (1979-2015)","Hindcast daily range (1979-2015)",if(!is.null(obsdata)) {"Observed"},if(!is.null(maxobs)) {"Observed max (1979-2015)"}), 
         col = c(.makeTransparent(fore.col, alpha = 200),.makeTransparent(hind.col, alpha = 200),hind.col,.makeTransparent(hind.col, 30),.makeTransparent(obs.col, 200),obs.col),
         lty = c(1,1,2,NA,1,2),
         pch = c(NA,NA,NA,15,NA,NA))
  
  # Legend 2: warning levels and hindcast range for this period (polygons)
  legend(x=par("usr")[2],y=par("usr")[4],xpd=T, bty = "n", cex=graphScale*lgndScale, lwd = lineWidth,lty=NA,pch=15, pt.cex=graphScale*2,
         legend = c(paste("Warning 3 (",wl.rp[3]," yr RP)",sep=""),paste("Warning 2 (",wl.rp[2]," yr RP)",sep=""),paste("Warning 1 (",wl.rp[1]," yr RP)",sep="")),
         col=c(.makeTransparent(wl.col[3], wl.alpha),.makeTransparent(wl.col[2], wl.alpha),.makeTransparent(wl.col[1], wl.alpha)))
  return(0)
}

# --------------------------------------------------------------------------
# (external) function to plot hindcast-forecast river discharga data
# --------------------------------------------------------------------------
plotBasinDischargeHindcastForecast<-function(){
  
}

# --------------------------------------------------------------------------
# (external) function to plot content of a basin output file to output dir
# --------------------------------------------------------------------------
plotBasinOutputFile<-function(fileDir,plotDir=NULL,basinFile,modelName="*any-hype*",
                              graphScale=1.6, lgndScale = 0.98,lineWidth=2,
                              maxmarg=0.1){
  
  # set plotDir=resDir if plotDir not given on input
  if(is.null(plotDir)){plotDir=resDir}
  
  # create plotDir if missing
  if(!file.exists(plotDir)) {dir.create(plotDir)}
  
  # Read basin output file
  basinData = ReadBasinOutput(paste(fileDir,basinFile,sep="/"))
  
  # List variables and their units
  variables = attr(basinData,"variables")
  units     = attr(basinData,"units")
  
  # jpeg filename base
  jpgFileBase = paste(plotDir,substr(basinFile,1,nchar(basinFile)-4),sep="/")
  
  # main title base
  mainBase = paste(substr(basinFile,1,nchar(basinFile)-4)," (",modelName,")",sep="")
  
  # loop over variables
  for(i in 1:length(variables)){
    
    # initiate jpg file for plotting
    jpgFileName = paste(jpgFileBase,"-",variables[i],".jpg",sep="")
    jpeg(filename = jpgFileName, res = 300, width = 31, height = 20, units = "cm")
    
    # plotting parameters (Jafet's settings)
    oldpar<-par()
    par(mgp=c(2.5,0.8,0),mar=c(3.1,4.1,3.1,23))
    
    # simple plot
    plot(basinData[,1],basinData[,variables[i]],type="l",col="black",
         main = paste(variables[i], "(",units[i],")"," - ",mainBase,sep=""),
         xlab="",ylab=paste(variables[i], "(",units[i],")",sep=""),
         cex.lab=graphScale,cex.axis=graphScale,cex.main=graphScale,
         lwd=lineWidth)
    
    # reset plotting parameters
    suppressWarnings(par(oldpar)) # warns for some parameters that cannot be modified... 
    dev.off()
  }
  return(0)
}

