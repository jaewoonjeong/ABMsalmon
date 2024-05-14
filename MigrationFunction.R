PixelPoint=function(Point){c(LongDivision[which.min(abs(LongDivision-Point[1]))], LatiDivision[which.min(abs(LatiDivision-Point[2]))])}
StartPoint=c(-125.020121, 49.939181); Destination=c(-125.813433, 50.392261); 
SP=PixelPoint(StartPoint) # Start Point
DP=PixelPoint(Destination) # Destination Point
R3DP1=PixelPoint(c(-125.023030, 50.223122)); R3DP2=PixelPoint(c(-125.329973,50.408208)); R4DP1=PixelPoint(c(-124.8, 50.1)); R4DP2=c(-125.191415, 50.306418); R5DP1=PixelPoint(c(-124.957136,50.210867)); R5DP2=PixelPoint(c(-125.415109,50.454108)) #  intermediate Destination Point
SPDP=t(data.frame(SP,DP)); colnames(SPDP)=c("Long","Lati"); SPDP=data.frame(SPDP)

#########################################################################################################################################
ClemensiFunction=function(MaxTime, Route){
  R1_linger=0; R2_linger=0.29; R3_linger=0.40 
  
  time=0
  salmon=PI[which.min(abs(PI$Lon-SP[1])+abs(PI$Lat-SP[2])),]
  salmon$Time=time; salmon$Distance=time*0.1
  salmon$RealTime=ymd_hms("2018-06-01 00:00:00")
  UpdateSalmon=salmon
  
  TimeDiff=15 # time gap to detect congestion
  MinDistance=5 # minimum distance to acknowledge progress
  TimePerPixel=277.6982 # time in seconds per pixel
  
  while(UpdateSalmon[1]>DP[1] & time<=MaxTime ){
    time=time+1
    
    ddf=PI[which.min(abs(salmon[time,1]-PI[,1])+abs(salmon[time,2]-PI[,2])),]
    North=subset(PI,LongID==ddf$LongID & LatiID==ddf$LatiID+1)
    East=subset(PI,LongID==ddf$LongID+1 & LatiID==ddf$LatiID)
    South=subset(PI,LongID==ddf$LongID & LatiID==ddf$LatiID-1)
    West=subset(PI,LongID==ddf$LongID-1 & LatiID==ddf$LatiID)
    dir=rbind(North,East,South,West)
    dir$Time=time; dir$Distance=time*0.1
    
    # decided point to move  c(N,S,E,W)
    if(Route==1){ # Route 1: one sub-section
      DistOrder=dir$DistanceToEnd
      
      if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
        NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
      } else {
        if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger,R1_linger/2,R1_linger/2,0))
        } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger*3/4,R1_linger*3/4,0))
        } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger/2,R1_linger/2))
        } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}}
    } 
    
    if(Route==2){ # Route 2: two sub-sections
      if(UpdateSalmon[1,2]<50.03){
        DistOrder=dir$DistanceToEnd
        NewLocation=sample(1:4,1,prob=c(.7,0,0,.3))
      } else {
        DistOrder=dir$DistanceToEnd
        
        if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
          NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
        } else {
          if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger,R2_linger/2,R2_linger/2,0))
          } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger*3/4,R2_linger*3/4,0))
          } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger/2,R2_linger/2))
          } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
        }
      }
    }
    
    if(Route==3){ # Route 3 : three sub-sections
      if(UpdateSalmon[1,2]<50.2){
        DistOrder=dir$R5DP1_DM
        if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
          NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
        } else {
          if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
          } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
          } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
          } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
        }
      } else {
        if(UpdateSalmon[1,1]> -125.4 & UpdateSalmon[1,2]>=50.2){
          DistOrder=dir$R5DP2_DM
          if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
            NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
          } else {
            if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
            } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
            } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
            } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}    
          } 
        } else {
          DistOrder=dir$DistanceToEnd
          
          if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
            NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
          } else {
            if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
            } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
            } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
            } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
          }
        }
      }
    }
    
    dir$RealTime=tail(salmon,1)$RealTime+seconds(TimePerPixel)
    UpdateSalmon=dir[NewLocation,] # Updated salmon with a new location
    
    PixelID=which(UpdateSalmon[1,3]==PI[,3] & UpdateSalmon[1,4]==PI[,4]) # number of row when updatesalmon is in PI
    
    salmon=rbind(salmon,UpdateSalmon)
    #print(UpdateSalmon)
  }
  return(salmon)
}

##################################################################
# Sensitivity Analysis: Progression Rate
SA_function=function(MaxTime, Route, beta){
  R1_linger=ifelse(beta<1,0,(beta-1)*1.12) ; R2_linger=0.29*beta; R3_linger=0.40*beta
  
  time=0
  salmon=PI[which.min(abs(PI$Lon-SP[1])+abs(PI$Lat-SP[2])),]
  salmon$Time=time; salmon$Distance=time*0.1
  salmon$RealTime=ymd_hms("2018-06-01 00:00:00")
  UpdateSalmon=salmon
  
  TimeDiff=15 # time gap to detect congestion
  MinDistance=5 # minimum distance to acknowledge progress
  TimePerPixel=277.6982 # time in seconds per pixel
  
  while(UpdateSalmon[1]>DP[1] & time<=MaxTime ){
    # & ifelse(time>(TimeDiff+30),sum(abs(salmon[time-(TimeDiff+30),]$LongID-UpdateSalmon$LongID),abs(salmon[time-(TimeDiff+30),]$LatiID-UpdateSalmon$LatiID))>MinDistance,TRUE)  
    time=time+1
    
    ddf=PI[which.min(abs(salmon[time,1]-PI[,1])+abs(salmon[time,2]-PI[,2])),]
    North=subset(PI,LongID==ddf$LongID & LatiID==ddf$LatiID+1)
    East=subset(PI,LongID==ddf$LongID+1 & LatiID==ddf$LatiID)
    South=subset(PI,LongID==ddf$LongID & LatiID==ddf$LatiID-1)
    West=subset(PI,LongID==ddf$LongID-1 & LatiID==ddf$LatiID)
    dir=rbind(North,East,South,West)
    dir$Time=time; dir$Distance=time*0.1
    #dir=subset(dir, SurroundingValue != 0)
    
    # decided point to move  c(N,S,E,W)
    if(Route==1){ # Route 1: one sub-section
      DistOrder=dir$DistanceToEnd
      if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
        NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
      } else {
        if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger,R1_linger/2,R1_linger/2,0))
        } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger*3/4,R1_linger*3/4,0))
        } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R1_linger/2,R1_linger/2))
        } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}}
    } 
    
    if(Route==2){ # Route 2: two sub-sections
      if(UpdateSalmon[1,2]<50.03){
        DistOrder=dir$DistanceToEnd
        NewLocation=sample(1:4,1,prob=c(.7,0,0,.3))
      } else {
        DistOrder=dir$DistanceToEnd
        if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
          NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
        } else {
          if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger,R2_linger/2,R2_linger/2,0))
          } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger*3/4,R2_linger*3/4,0))
          } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R2_linger/2,R2_linger/2))
          } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
        }
      }
    }
    
    if(Route==3){ # Route 3 : three sub-sections
      if(UpdateSalmon[1,2]<50.2){
        DistOrder=dir$R5DP1_DM
        if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
          NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
        } else {
          if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
          } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
          } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
          } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
        }
      } else {
        if(UpdateSalmon[1,1]> -125.4 & UpdateSalmon[1,2]>=50.2){
          DistOrder=dir$R5DP2_DM
          if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
            NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
          } else {
            if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
            } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
            } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
            } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}    
          } 
        } else {
          DistOrder=dir$DistanceToEnd
          if(time>TimeDiff & sum(abs(salmon[time-TimeDiff,]$LongID-UpdateSalmon$LongID),abs(salmon[time-TimeDiff,]$LatiID-UpdateSalmon$LatiID))<MinDistance){
            NewLocation=sample(1:nrow(dir),1,prob=c(rep(1,nrow(dir))))
          } else {
            if(nrow(dir)==4){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger,R3_linger/2,R3_linger/2,0))
            } else {if(nrow(dir)==3){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger*3/4,R3_linger*3/4,0))
            } else {if(nrow(dir)==2){NewLocation=sample(order(DistOrder),size=1,prob=c(1-R3_linger/2,R3_linger/2))
            } else {NewLocation=sample(order(DistOrder),size=1,prob=c(1))}}}
          }}}}
    
    dir$RealTime=tail(salmon,1)$RealTime+seconds(TimePerPixel)
    UpdateSalmon=dir[NewLocation,] # Updated salmon with a new location
    
    PixelID=which(UpdateSalmon[1,3]==PI[,3] & UpdateSalmon[1,4]==PI[,4]) # number of row when updatesalmon is in PI
    
    salmon=rbind(salmon,UpdateSalmon)
    #print(UpdateSalmon)
  }
  return(salmon)
}
