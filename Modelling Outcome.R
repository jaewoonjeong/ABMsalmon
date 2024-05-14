library(reshape); library(ggplot2); library(ggridges); library(gridExtra)
#save(RR, Infestation, sim, iteration, file='01_tidy_data/InfestationModel.RData')
load(file='01_tidy_data/InfestationModel.RData'); iteration; length(RR[[1]])
#dim(Infestation)=c(sim,iteration,5); Infestation[t,i,r]

######## Range of migration distance  ##################### IN Manuscript
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

library(patchwork)
# Abundance, Total IP, and Migration distance
tiff(filename="03_plots/Fig. 4.tif", width=19.05, height=7, units='cm', res=600)
a+b+c+plot_annotation(tag_levels='A')
dev.off()

# Abundance by Routes
Abun=numeric(); for(r in 1:3){Abun[r]=mean(Infestation[,,r])}
Minf=melt(Infestation)
colnames(Minf)=c('Simulation','Iteration','Route','Abundance')
Minf$Route=factor(Minf$Route,levels=c('1','2','3'))
Minf$Route=ifelse(Minf$Route=='1', 'Route 1', ifelse(Minf$Route=='2','Route 2','Route 3'))
Minf$Iteration=as.factor(Minf$Iteration); Minf$Route=as.factor(Minf$Route)

Abun_Fixed=numeric(); for(r in 1:3){Abun_Fixed[r]=mean(Infestation_Fixed[,,r])}
Minf_Fixed=melt(Infestation_Fixed)
colnames(Minf_Fixed)=c('Simulation','Iteration','Route','Abundance')
Minf_Fixed$Route=factor(Minf_Fixed$Route,levels=c('1','2','3'))
Minf_Fixed$Route=ifelse(Minf_Fixed$Route=='1', 'Route 1', ifelse(Minf_Fixed$Route=='2','Route 2','Route 3'))
Minf_Fixed$Iteration=as.factor(Minf_Fixed$Iteration); Minf_Fixed$Route=as.factor(Minf_Fixed$Route)

Mchoice=cbind(rbind(Minf, Minf_Fixed), Choice=c(rep('Log',nrow(Minf)),rep('Fixed',nrow(Minf_Fixed))))
Mchoice$Choice=factor(Mchoice$Choice, levels=c('Log','Fixed'))
custom_labels <- as_labeller(c("Route 1"="Route 1","Route 2"="Route 2","Route 3"="Route 3", Fixed = "Proportional model", Log = "Logarithmic model"))

F6=ggplot(data=subset(Mchoice, Abundance<10),mapping=aes(x=Abundance))+
  geom_histogram(binwidth=0.5)+facet_grid(Choice~Route, labeller = custom_labels)+theme(legend.position='none')+
  xlab('Count of sea lice per fish')+ylab('Proportion')+theme_bw()+
  scale_x_continuous(breaks=0:9,labels=c(0:9))+
  scale_y_continuous(breaks=c(0,300,600,900,1200),labels=c(0,300,600,900,1200)/3000)

tiff(filename="03_plots/Fig. 6.tif", width=19.05, height=11, units='cm', res=600)
F6
dev.off()

aggregate(Abundance~Route, data=Minf, mean)
aggregate(Abundance~Route, data=Minf_Fixed, mean)

####################### Variation within Route
Minf_c=cbind(Minf,Value=rep(melt(TotalIP)[,3],each=sim))
Minf_d=cbind(Minf,Value=rep(melt(TotalPeriod)[,3],each=sim))
Minf_cd=cbind(rbind(Minf_c,Minf_d),Type=c(rep('Total Infestation Pressure',iteration*3),rep('Migration Distance',iteration*3)))
ggplot(subset(Minf_cd,Simulation==1&Route!=1),aes(Value,Abundance))+geom_point()+facet_grid(Route~Type,scales = "free_x")+xlab('Migration Distance (km)                   Total Infestation Pressure')
ggplot(subset(Minf_cd,Simulation==1&Route!=1),aes(Value,as.factor(Abundance)))+geom_violin()+facet_grid(Route~Type,scales = "free_x")+xlab('Migration Distance (km)                   Total Infestation Pressure')

# correlation between abundance and total IP or migration distance
par(mfrow=c(2,2))
df=subset(Minf_cd,Route==2 & Type=="Total Infestation Pressure"); cor(df$Value,df$Abundance); plot(df$Value,df$Abundance)
df=subset(Minf_cd,Route==2 & Type=="Migration Distance"); cor(df$Value,df$Abundance); plot(df$Value,df$Abundance)
df=subset(Minf_cd,Route==3 & Type=="Total Infestation Pressure"); cor(df$Value,df$Abundance); plot(df$Value,df$Abundance)
df=subset(Minf_cd,Route==3 & Type=="Migration Distance"); cor(df$Value,df$Abundance); plot(df$Value,df$Abundance)

par(mfrow=c(2,1))
TIP=colMeans(TotalIP)
matplot(TIP,Abun,type='p',col = rgb(1, 0, 0, alpha = 0),las=1,xlab='Total Infestation Pressure',ylab='Abundance')
for(r in 1:3){text(TIP[r],Abun[r],r,col=1)}

TIP=colMeans(TotalPeriod)
matplot(TIP,Abun,type='p',col = rgb(1, 0, 0, alpha = 0),las=1,xlab='Migration distance',ylab='Abundance')
for(r in 1:3){text(TIP[r],Abun[r],r,col=1)}
