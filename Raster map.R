############# Raster 

SL=read.csv(file='') # lice count data from salmon farms in British Columbia, Canada 
SL=na.omit(SL, cols='Average.caligus.per.fish')
SL=subset(SL,Finfish.Aquaculture.Management.Unit=='Discovery Islands')
SL$Incident.Date=as.Date(SL$Incident.Date,"%d-%b-%y")
SL$isoweek=isoweek(SL$Incident.Date)
SL=subset(SL,Year==2018 & isoweek>17 & isoweek<27 & Longitude> -125.8)
SL_unique=SL%>% distinct(Site.Common.Name, .keep_all = TRUE)

LongMin=-126; LongMax=-124.7911; LatiMin=49.9; LatiMax=50.51
PixelLong=1.411973/(10^3); PixelLati=0.8983153/(10^3)
LongDivision=seq(LongMin,LongMax,PixelLong); LatiDivision=seq(LatiMin,LatiMax,PixelLati)

PixelPoint=function(Point){c(LongDivision[which.min(abs(LongDivision-Point[1]))], LatiDivision[which.min(abs(LatiDivision-Point[2]))])}
StartPoint=c(-125.020121, 49.939181); Destination=c(-125.813433, 50.392261); 
SP=PixelPoint(StartPoint) # Start Point
DP=PixelPoint(Destination) # Destination Point

W=function(d){(1/(sqrt(2*3.14)*(Bandwidth/4)))*exp((-d^2)/(2*(Bandwidth/4)^2))} # Gaussian kernel density

SHP=readOGR(dsn="") # coastline geographic file data in shp file
LongMin=-126; LongMax=-124.7911; LatiMin=49.9; LatiMax=50.51
ext=extent(LongMin,LongMax,LatiMin,LatiMax)
PixelLong=1.411973/(10^3); PixelLati=0.8983153/(10^3)
rr <- raster(ext, res=c(PixelLong,PixelLati)) # length of border of a pixel = 100 meters
r<-rasterize(SHP, rr)

LongID=1:length(LongDivision); LatiID=1:length(LatiDivision)
FullDF=data.frame(Lon=rep(LongDivision,length(LatiDivision)),Lat=rep(LatiDivision,each=length(LongDivision)),LongID=rep(LongID,length(LatiDivision)), LatiID=rep(LatiID,each=length(LongDivision)))
mm=matrix(NA, nrow = length(LongDivision),ncol = length(LatiDivision)); for(i in 1:length(LongDivision)){for(j in 1:length(LatiDivision)){mm[i,j]=r[length(LatiDivision)-j+1,i,1]}}; 
DF=cbind(FullDF,Ocean=melt(mm)[,3]) # DF = all pixels
df=subset(DF,Ocean==1)

#Distance matrix
cost_Discovery=matrix(NA,nrow(df),nrow(SL_unique)); for(j in 1:nrow(SL_unique)){for(i in 1:40000){cost_Discovery[i,j]=costDistance(tr_Discovery,rbind(as.numeric(df[i,c(1,2)]),as.numeric(SL_unique[j,c(7,6)])))}; print(j)}; 
load(file='01_tidy_data/cost_Discovery1.RData'); cost_Discovery1=cost_Discovery; load(file='01_tidy_data/cost_Discovery2.RData'); cost_Discovery2=cost_Discovery; load(file='01_tidy_data/cost_Discovery3.RData'); cost_Discovery3=cost_Discovery; load(file='01_tidy_data/cost_Discovery4.RData'); cost_Discovery4=cost_Discovery
cost_Discovery[1:40000,]=cost_Discovery1[1:40000,]; cost_Discovery[40001:80000,]=cost_Discovery2[40001:80000,]; cost_Discovery[80001:120000,]=cost_Discovery3[80001:120000,]; cost_Discovery[120001:nrow(cost_Discovery),]=cost_Discovery4[120001:nrow(cost_Discovery),]

InfCost=which(cost_Discovery[,1]==Inf) # Inf in distance matrix
DistanceMatrix=cost_Discovery[-InfCost,] # to remove rows with Inf

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

Trawl=function(IP){exp(-5.49+log(0.0032141)+0.425*log(IP))} # To change IP to TP
TP=Trawl(IP)

PI=cbind(PixelInfo,IP=IP,TP=TP,DistanceToEnd=DistanceToEnd,R3DP1_DM=R3DP1_DM,R3DP2_DM=R3DP2_DM)
