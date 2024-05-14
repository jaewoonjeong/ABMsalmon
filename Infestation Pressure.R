rm(list=ls())
source(file='02_r_scripts/Raster map.R')

PixelPoint=function(Point){c(LongDivision[which.min(abs(LongDivision-Point[1]))], LatiDivision[which.min(abs(LatiDivision-Point[2]))])}
StartPoint=c(-125.020121, 49.939181); Destination=c(-125.813433, 50.392261); 
SP=PixelPoint(StartPoint) # Start Point
DP=PixelPoint(Destination) # Destination Point

W=function(d){(1/(sqrt(2*3.14)*(Bandwidth/4)))*exp((-d^2)/(2*(Bandwidth/4)^2))} # Gaussian kernel density
d=DistanceMatrix/1000 # to change m to km
Bandwidth=30
dSum=matrix(NA,nrow(PixelInfo),nrow(SL_unique)); dim(dSum)=c(nrow(PixelInfo),nrow(SL_unique))
ip=matrix(NA,nrow(PixelInfo),nrow(SL_unique)); dim(ip)=c(nrow(PixelInfo),nrow(SL_unique))
for(i in 1:nrow(PixelInfo)){for(j in 1:nrow(SL_unique)){dSum[i,j]=W(d[i,j])/1000}}
for(i in 1:nrow(PixelInfo)){for(j in 1:nrow(SL_unique)){ip[i,j]=dSum[i,j]/sum(dSum[,j])};print(i)}
DistanceToEnd=numeric(); for(n in 1:nrow(PixelInfo)){DistanceToEnd[n]=costDistance(tr_Discovery,DP,as.numeric(PixelInfo[n,c(1,2)]))/1000; print(n)}
#save(dSum,ip,Bandwidth, DistanceToEnd, R3DP1_DM, R3DP2_DM, file='01_tidy_data/ip.RData')
load(file='01_tidy_data/ip.RData')

LiceCountsPerFarm=62500000
LiceCounts=LiceCountsPerFarm*8 # Lice counts from all farms
iipp=matrix(NA,nrow(PixelInfo),nrow(SL_unique)*length(LiceCounts)); dim(iipp)=c(nrow(PixelInfo),nrow(SL_unique),length(LiceCounts))
IP=matrix(NA,nrow(PixelInfo),length(LiceCounts))
for(i in 1:nrow(PixelInfo)){for(j in 1:nrow(SL_unique)){for(s in 1:length(LiceCounts)){iipp[i,j,s]=ip[i,j]*LiceCounts[s]}};print(i)}
for(i in 1:nrow(PixelInfo)){for(s in 1:length(LiceCounts)){IP[i,s]=sum(iipp[i,,s])};print(i)}

Trawl=function(IP){exp(-5.49+log(0.0032141)+0.425*log(IP))}
TP=Trawl(IP)

PI=cbind(PixelInfo,IP=IP,TP=TP,DistanceToEnd=DistanceToEnd,R3DP1_DM=R3DP1_DM,R3DP2_DM=R3DP2_DM)

#save(dSum,ip,iipp,IP,PI,LiceCounts,Bandwidth, TP, file='01_tidy_data/PI.RData')
load(file='01_tidy_data/PI.RData')
library(ggplot2)
IP=ggplot(PI,mapping=aes(Lon,Lat))+geom_raster(aes(fill=IP))+scale_fill_gradientn(colors = c("skyblue", "yellow", "red"))+coord_fixed()+
  geom_point(data=SL_unique,mapping=aes(x=Longitude,y=Latitude))+theme(legend.position='bottom', legend.text=element_text(angle=90))
TP=ggplot(PI,mapping=aes(Lon,Lat))+geom_raster(aes(fill=TP))+scale_fill_gradientn(colors = c("skyblue", "yellow", "red"))+coord_fixed()+
  geom_point(data=SL_unique,mapping=aes(x=Longitude,y=Latitude))+theme(legend.position='bottom', legend.text=element_text(angle=90))
library(gridExtra)
grid.arrange(IP,TP)

xx=seq(min(PI$IP),max(PI$IP),length.out=100)
plot(xx,Trawl(xx),type='l',xlab='Infestation Pressure',ylab='Transmission Pressure')
