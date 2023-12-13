############################ Migration Model & Infestation Model
source(file='Raster map.R')
source(file='MigrationFunction.R')

MaxTime=12000
iteration=1
RR1=list()
for(r in 1){II=list()
for(i in 1:iteration){II[[i]]=ClemensiFunction(MaxTime, r)}
RR1[[r]]=II
print(r)}

iteration=3000
RR=list()
for(r in 2:3){II=list()
for(i in 1:iteration){II[[i]]=ClemensiFunction(MaxTime, r);print(i)}
RR[[r]]=II; print(r)}

for(i in 1:iteration){RR[[1]][[i]]=RR1[[1]][[1]]}

# To accumulate IP
for(r in 1:3){for(i in 1:iteration){for(p in 1:nrow(RR[[r]][[i]])){RR[[r]][[i]]$AccumulatedIP[p]=sum(RR[[r]][[i]]$IP[1:p])}};print(r)}

# How Infestation pressure has changed during migration depending on Route
r=1; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_1 <- do.call(rbind, df_a)
r=2; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_2 <- do.call(rbind, df_a)
r=3; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_3 <- do.call(rbind, df_a)
combined_df=cbind(rbind(combined_df_1,combined_df_2,combined_df_3),Route=c(rep('1',nrow(combined_df_1)),rep('2',nrow(combined_df_2)),rep('3',nrow(combined_df_3))))
combined_df$Iteration=as.factor(combined_df$Iteration)
ss=subset(combined_df,Iteration==1|Iteration==2|Iteration==3|Iteration==4|Iteration==5)
ss$Route=paste('Route',ss$Route)
ss$AccumulatedIP=ss$AccumulatedIP/1000000

# Figure 5
F5=ggplot(data=ss,mapping=aes(x=Distance,y=AccumulatedIP,fill=Iteration,color=factor(Iteration)))+geom_line(linewidth=.8)+
  guides(color = guide_legend(title = "Iteration"))+scale_color_manual(values=c("1"=2,"2"=3,"3"=4,"4"=5,"5"=1))+
  facet_wrap(.~Route,nrow=1)+xlab("Distance (km)")+ylab("Cumulative Infestation Pressure per Pixel in millions")+theme_bw()+theme(legend.position='none',axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))


sim=1

Infestation=matrix(NA,sim*3*iteration); dim(Infestation)=c(sim,iteration,3)
for(r in 1:3){for(i in 1:iteration){NewInf=matrix(NA,nrow((RR[[r]][[i]])),sim)
for(t in 1:sim){for(n in 1:nrow((RR[[r]][[i]]))){NewInf[n,t]=rnbinom(1,size=2.733,mu=RR[[r]][[i]]$TP[n])}
  Infestation[t,i,r]=sum(NewInf[,t])};print(i)};print(r)}

############################## Modeling outcome
# Migration period
TotalPeriod=matrix(NA,iteration,3); for(i in 1:iteration){for(r in 1:3){TotalPeriod[i,r]=nrow(RR[[r]][[i]])/10}}
summary(TotalPeriod)
md=melt(TotalPeriod); colnames(md)=c('Iteration','Route','Distance')
md$Route <- factor(md$Route, levels = c(3,2,1))
a=ggplot(data=subset(md,as.factor(Route)!=1),aes(y=as.factor(Route),x=Distance,fill=as.factor(Route)))+geom_density_ridges2(bandwidth=2.2)+scale_fill_manual(values=c("2"=3,"3"=4))+
  theme_bw()+theme(legend.position="none",plot.title.position = "plot",axis.title.x = element_text(size = 8))+xlab("Migration Distance (km)")+ylab("Route")+ggtitle("")+xlim(TotalPeriod[1,1],max(TotalPeriod))+geom_vline(xintercept=TotalPeriod[1,1], linetype = "dashed", color = 2,lwd=1)
# Progression rate
RouteDistance=c(106.7, 113.7, 143.7)
ProgressionRate=matrix(NA,iteration,3); for(i in 1:iteration){for(r in 1:3){ProgressionRate[i,r]=RouteDistance[r]/(nrow(RR[[r]][[i]])*0.0032141)}}
summary(ProgressionRate)
pr=melt(ProgressionRate); colnames(pr)=c('Iteration','Route','PR')
pr$Route <- factor(pr$Route, levels = c(3,2,1))
b=ggplot(data=subset(pr,as.factor(Route)!=1),aes(y=as.factor(Route),x=PR,fill=as.factor(Route)))+geom_density_ridges2(bandwidth=0.2)+scale_fill_manual(values=c("2"=3,"3"=4))+
  theme_bw()+theme(legend.position="none",plot.title.position = "plot",axis.title.x = element_text(size = 8))+xlab("Progression Rate (km/day)")+ylab("Route")+ggtitle("")+xlim(min(ProgressionRate),ProgressionRate[1,1])+geom_vline(xintercept=ProgressionRate[1,1], linetype = "dashed", color = 2,lwd=1)
# Total Infestation Pressure
TotalIP=matrix(NA,iteration,3); for(i in 1:iteration){for(r in 1:3){TotalIP[i,r]=sum(RR[[r]][[i]]$IP)/1000000}}
summary(TotalIP)
tip=melt(TotalIP); colnames(tip)=c('Iteration','Route','TIP')
tip$Route <- factor(tip$Route, levels = c(3,2,1))
c=ggplot(data=subset(tip,as.factor(Route)!=1),aes(y=as.factor(Route),x=TIP,fill=as.factor(Route)))+geom_density_ridges2(bandwidth=2)+scale_fill_manual(values=c("2"=3,"3"=4))+
  theme_bw()+theme(legend.position="none",plot.title.position = "plot",axis.title.x = element_text(size = 8))+xlab("Total Infestation Pressure (millions)")+ylab("Route")+ggtitle("")+xlim(TotalIP[1,1],max(TotalIP))+geom_vline(xintercept=TotalIP[1,1], linetype = "dashed", color = 2,lwd=1)

# Figure 4
F4=a+b+c+plot_annotation(tag_levels='A') 

# Abundance by Routes
Abun=numeric(); for(r in 1:3){Abun[r]=mean(Infestation[,,r])}

Minf=melt(Infestation)
colnames(Minf)=c('Simulation','Iteration','Route','Abundance')
Minf$Route=factor(Minf$Route,levels=c('1','2','3'))
Minf$Route=ifelse(Minf$Route=='1', 'Route 1', ifelse(Minf$Route=='2','Route 2','Route 3'))
Minf$Iteration=as.factor(Minf$Iteration); Minf$Route=as.factor(Minf$Route)

# Figure 6
F6=ggplot(data=subset(Minf, Abundance<10),mapping=aes(x=Abundance))+
  geom_histogram(binwidth=0.5)+facet_wrap(Route~.,nrow=1)+theme(legend.position='none')+
  xlab('Count of sea lice per fish')+ylab('Proportion')+theme_bw()+
  scale_x_continuous(breaks=0:9,labels=c(0:9))+
  scale_y_continuous(breaks=c(0,300,600,900,1200),labels=c(0,300,600,900,1200)/3000)

