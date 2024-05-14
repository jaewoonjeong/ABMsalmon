rm(list=ls())
source(file='02_r_scripts/Raster map.R')

################################################################
# Migration model Simulation
source(file='02_r_scripts/MigrationFunction.R')
load(file='01_tidy_data/PI.RData') # PI: IP, TP
# ClemensiFunction=function(MaxTime, Route)
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

# To add fixed proportional relation of IP and TP for sensitivity analysis
Fixed=function(IP){exp(-11.6+log(0.0032141)+1*log(IP))}
for(r in 1:3){for(i in 1:3000){RR[[r]][[i]]$TP_Fixed=Fixed(RR[[r]][[i]]$IP)};print(r)}

#save(RR1, RR, MaxTime, iteration, file='01_tidy_data/RR.RData')
load(file='01_tidy_data/RR.RData')
length(RR); length(RR[[2]]); head(RR[[1]][[1]]); dim((RR1[[1]][[1]]))

# How Infestation pressure has changed during migration depending on Route
r=1; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_1 <- do.call(rbind, df_a)
r=2; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_2 <- do.call(rbind, df_a)
r=3; df_a=list(); for(i in 1:iteration){df=data.frame(Distance=RR[[r]][[i]]$Distance,AccumulatedIP=RR[[r]][[i]]$AccumulatedIP); df_a[[i]]=cbind(df,Iteration=rep(i,nrow(df)))}; combined_df_3 <- do.call(rbind, df_a)
combined_df=cbind(rbind(combined_df_1,combined_df_2,combined_df_3),Route=c(rep('1',nrow(combined_df_1)),rep('2',nrow(combined_df_2)),rep('3',nrow(combined_df_3))))
combined_df$Iteration=as.factor(combined_df$Iteration)
ss=subset(combined_df,Iteration==1|Iteration==2|Iteration==3|Iteration==4|Iteration==5)
ss$Route=paste('Route',ss$Route)
ss$AccumulatedIP=ss$AccumulatedIP/1000000
F5=ggplot(data=ss,mapping=aes(x=Distance,y=AccumulatedIP,fill=Iteration,color=factor(Iteration)))+geom_line(linewidth=.8)+
  guides(color = guide_legend(title = "Iteration"))+scale_color_manual(values=c("1"=2,"2"=3,"3"=4,"4"=5,"5"=1))+
  facet_wrap(.~Route,nrow=1)+xlab("Distance (km)")+ylab("Cumulative Infestation Pressure per Pixel in millions")+theme_bw()+theme(legend.position='none',axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))

tiff(filename="03_plots/Fig. 5.tif", width=19.05, height=7, units='cm', res=600)
F5
dev.off()
# Non-accumulated IP
ggplot(data=ss,mapping=aes(x=Distance,y=IP,fill=Iteration,color=factor(Iteration)))+geom_line(linewidth=.01)+guides(color = guide_legend(title = "Iteration"))+facet_grid(Route~.)+xlab("Distance (km)")+ylab("Infestation Pressure per Pixel in millions")+theme_bw()+theme(legend.position='none')
######################################################################################################
# Infestation model
load(file='01_tidy_data/PI.RData')
load(file='01_tidy_data/RR.RData'); iteration; length(RR); length(RR[[1]]); head(RR[[1]][[1]])
iteration; length(RR[[1]]); length(RR[[2]]); length(RR[[3]])

sim=1
# Logarithmic relationship between IP and TP
Infestation=matrix(NA,sim*3*iteration); dim(Infestation)=c(sim,iteration,3)
for(r in 1:3){for(i in 1:iteration){NewInf=matrix(NA,nrow((RR[[r]][[i]])),sim)
for(t in 1:sim){for(n in 1:nrow((RR[[r]][[i]]))){NewInf[n,t]=rnbinom(1,size=0.921,mu=RR[[r]][[i]]$TP[n])}
  Infestation[t,i,r]=sum(NewInf[,t])}};print(r)}

# Proportional relationship between IP and TP
Infestation_Fixed=matrix(NA,sim*3*iteration); dim(Infestation_Fixed)=c(sim,iteration,3)
for(r in 1:3){for(i in 1:iteration){NewInf_Fixed=matrix(NA,nrow((RR[[r]][[i]])),sim)
for(t in 1:sim){for(n in 1:nrow((RR[[r]][[i]]))){NewInf_Fixed[n,t]=rnbinom(1,size=0.921,mu=RR[[r]][[i]]$TP_Fixed[n])}
  Infestation_Fixed[t,i,r]=sum(NewInf_Fixed[,t])}};print(r)}

#save(RR, Infestation, Infestation_Fixed, NewInf, sim, iteration, file='01_tidy_data/InfestationModel.RData')
load(file='01_tidy_data/InfestationModel.RData')
iteration; sim

################## Sensitivity analysis of dispersion parameter
# Different beta values in negative binomial distribution
beta_para=c(1/100, 1/10, 1/2, 1, 2, 10, 100)
Infestation_beta=matrix(NA,sim*3*iteration*length(beta_para)); dim(Infestation_beta)=c(sim,iteration,3,length(beta_para))
for(b in 1:length(beta_para)){for(r in 1:3){for(i in 1:iteration){
  NewInf_beta=matrix(NA,nrow((RR[[r]][[i]])),sim)
  for(t in 1:sim){for(n in 1:nrow((RR[[r]][[i]]))){NewInf_beta[n,t]=rnbinom(1,size=beta_para[b],mu=RR[[r]][[i]]$TP[n])}
    Infestation_beta[t,i,r,b]=sum(NewInf_beta[,t])}};print(r)};print(b)}

Abun_beta=matrix(NA,3,length(beta_para)); for(b in 1:length(beta_para)){for(r in 1:3){Abun_beta[r,b]=mean(Infestation_beta[,,r,b])}}
Minf_beta=melt(Infestation_beta)
colnames(Minf_beta)=c('Simulation','Iteration','Route','Beta','Abundance')
Minf_beta$Route=factor(Minf_beta$Route,levels=c('1','2','3'))
Minf_beta$Route=ifelse(Minf_beta$Route=='1', 'Route 1', ifelse(Minf_beta$Route=='2','Route 2','Route 3'))
Minf_beta$Iteration=as.factor(Minf_beta$Iteration); Minf_beta$Route=as.factor(Minf_beta$Route)
Minf_beta$Beta=as.factor(beta_para[Minf_beta$Beta])

ggplot(data=subset(Minf_beta, Abundance<10),mapping=aes(x=Abundance))+
  geom_histogram(binwidth=0.5)+facet_grid(Beta~Route)+theme(legend.position='none')+
  xlab('Count of sea lice per fish')+ylab('Proportion')+theme_bw()
