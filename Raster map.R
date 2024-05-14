library(gdistance); library(dplyr); library(lubridate); library(reshape); library(ggplot2); library(ggrepel)
SL=read.csv(file='00_raw_data/industry data.csv'); SL=na.omit(SL, cols='Average.caligus.per.fish')
SL=subset(SL,Finfish.Aquaculture.Management.Unit=='Discovery Islands')
SL$Incident.Date=as.Date(SL$Incident.Date,"%d-%b-%y")
SL$isoweek=isoweek(SL$Incident.Date)
SL=subset(SL,isoweek>17 & isoweek<27 & Longitude> -125.8)
SL=subset(SL, Year==2018)

SL_unique=SL%>% distinct(Site.Common.Name, .keep_all = TRUE)

#SHP=readOGR(dsn="00_raw_data/Water_Background.shp"); save(SHP,file='01_tidy_data/SHP.RData')
load(file='01_tidy_data/SHP.RData')
LongMin=-126; LongMax=-124.7911; LatiMin=49.9; LatiMax=50.51
ext=extent(LongMin,LongMax,LatiMin,LatiMax)
PixelLong=1.411973/(10^3); PixelLati=0.8983153/(10^3)
rr <- raster(ext, res=c(PixelLong,PixelLati)) # length of border of a pixel = 100 meters
#r<-rasterize(SHP, rr); save(r, file='01_tidy_data/r.RData')
load(file='01_tidy_data/r.RData')
tr=transition(r,mean,directions=4)
tr_Discovery=geoCorrection(tr, type="c")

LongDivision=seq(LongMin,LongMax,PixelLong); LatiDivision=seq(LatiMin,LatiMax,PixelLati)
LongID=1:length(LongDivision); LatiID=1:length(LatiDivision)
FullDF=data.frame(Lon=rep(LongDivision,length(LatiDivision)),Lat=rep(LatiDivision,each=length(LongDivision)),LongID=rep(LongID,length(LatiDivision)), LatiID=rep(LatiID,each=length(LongDivision)))
#mm=matrix(NA, nrow = length(LongDivision),ncol = length(LatiDivision)); for(i in 1:length(LongDivision)){for(j in 1:length(LatiDivision)){mm[i,j]=r[length(LatiDivision)-j+1,i,1]}}; 
#save(mm, file='01_tidy_data/mm.RData')
load(file='01_tidy_data/mm.RData')

DF=cbind(FullDF,Ocean=melt(mm)[,3]) # DF = all pixels
df=subset(DF,Ocean==1)
# Distance matrix; cost_Discovery=matrix(NA,nrow(df),nrow(SL_unique)); for(j in 1:nrow(SL_unique)){for(i in 1:40000){cost_Discovery[i,j]=costDistance(tr_Discovery,rbind(as.numeric(df[i,c(1,2)]),as.numeric(SL_unique[j,c(7,6)])))}; print(j)}; 
#load(file='01_tidy_data/cost_Discovery1.RData'); cost_Discovery1=cost_Discovery; load(file='01_tidy_data/cost_Discovery2.RData'); cost_Discovery2=cost_Discovery; load(file='01_tidy_data/cost_Discovery3.RData'); cost_Discovery3=cost_Discovery; load(file='01_tidy_data/cost_Discovery4.RData'); cost_Discovery4=cost_Discovery
#cost_Discovery[1:40000,]=cost_Discovery1[1:40000,]; cost_Discovery[40001:80000,]=cost_Discovery2[40001:80000,]; cost_Discovery[80001:120000,]=cost_Discovery3[80001:120000,]; cost_Discovery[120001:nrow(cost_Discovery),]=cost_Discovery4[120001:nrow(cost_Discovery),]
#save(cost_Discovery, file='01_tidy_data/cost_Discovery.RData')
load(file='01_tidy_data/cost_Discovery.RData')

InfCost=which(cost_Discovery[,1]==Inf) # Inf in distance matrix
DistanceMatrix=cost_Discovery[-InfCost,] # to remove rows with Inf
PixelInfo=df[-InfCost,] # to remove rows with Inf
