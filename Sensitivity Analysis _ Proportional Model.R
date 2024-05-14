rm(list=ls())

source(file='02_r_scripts/Raster map.R')

d=DistanceMatrix/1000 # to change m to km
Bandwidth=c(5,10,20,30,40,60,100)
sigma=Bandwidth/4

dSum=ip=matrix(NA,nrow(d),nrow(SL_unique)*length(Bandwidth)); dim(dSum)=dim(ip)=c(nrow(d),nrow(SL_unique),length(Bandwidth))
for(b in 1:length(Bandwidth)){W=function(d){(1/(sqrt(2*pi)*(Bandwidth[b]/4)))*exp((-d^2)/(2*(Bandwidth[b]/4)^2))} # Gaussian kernel density
for(i in 1:nrow(d)){for(j in 1:nrow(SL_unique)){dSum[i,j,b]=W(d[i,j])/1000}};print(b)}
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(j in 1:nrow(SL_unique)){ip[i,j,b]=ifelse(sum(dSum[,j,b])==0,0,dSum[i,j,b]/sum(dSum[,j,b]))}};print(i)}

LiceCountsPerFarm=62500000
LiceCounts=LiceCountsPerFarm*8*c(.18, .4, .7, 1, 2, 3, 4.12) # Lice count were updated 

# To calculate IP for each lice count values
iipp=matrix(NA,nrow(d),nrow(SL_unique)*length(Bandwidth)*length(LiceCounts)); dim(iipp)=c(nrow(d),nrow(SL_unique),length(Bandwidth),length(LiceCounts))
IP=matrix(NA,nrow(d),length(Bandwidth)*length(LiceCounts)); dim(IP)=c(nrow(d),length(Bandwidth),length(LiceCounts))
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(j in 1:nrow(SL_unique)){for(s in 1:length(LiceCounts)){iipp[i,j,b,s]=ip[i,j,b]*LiceCounts[s]}}};print(i)}
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(s in 1:length(LiceCounts)){IP[i,b,s]=sum(iipp[i,,b,s])}};print(i)}
PI=list(); for(s in 1:length(LiceCounts)){ppii=list(); for(b in 1:length(Bandwidth)){ppii[[b]]=cbind(PixelInfo,IP=IP[,b,s])}; PI[[s]]=ppii}
PI_SA=PI
#save(dSum,ip,iipp,IP,Bandwidth,LiceCounts,PI_SA, file='01_tidy_data/PI_SA.RData')

PI_SA
#####################################################

load(file='01_tidy_data/PI_SA.RData') # PI: IP, TP
source(file='02_r_scripts/MigrationFunction.R')
load(file="01_tidy_data/beta.RData")
load(file="01_tidy_data/SA_PR.RData")
length(RR); length(RR[[2]]); length(RR[[2]][[1]]); length(RR[[2]][[1]][[1]]); head(RR[[2]][[1]][[1]])
gamma=rev(c(0.20, 0.47, 0.73, 1.00, 1.17, 1.33, 1.5))

SA_infes=function(r,c,t,g,b){
  Trawl=function(IP){exp(-5.49+log(0.0032141*gamma[g])+0.425*log(IP))}
  if(r==1){b=4}
  NewInf=numeric()
  i=sample(1:iteration,1)
  for(n in 1:nrow(RR[[r]][[i]][[b]])){TP_SA=Trawl(PI_SA[[c]][[t]][which(RR[[r]][[i]][[b]][n,3]==PI_SA[[c]][[t]][,3] & RR[[r]][[i]][[b]][n,4]==PI_SA[[c]][[t]][,4]),]$IP)
  NewInf[n]=rnbinom(1,size=0.921,mu=TP_SA)}
  return(sum(NewInf))
}

ProportionalModel=function(r,c,t,g,b){
  Trawl=function(IP){exp(-11.6+log(0.0032141*gamma[g])+log(IP))}
  if(r==1){b=4}
  NewInf=numeric()
  i=sample(1:iteration,1)
  for(n in 1:nrow(RR[[r]][[i]][[b]])){TP_SA=Trawl(PI_SA[[c]][[t]][which(RR[[r]][[i]][[b]][n,3]==PI_SA[[c]][[t]][,3] & RR[[r]][[i]][[b]][n,4]==PI_SA[[c]][[t]][,4]),]$IP)
  NewInf[n]=rnbinom(1,size=0.921,mu=TP_SA)}
  return(sum(NewInf))
}
############################################################
# pearson SA excluding Route
N=3000 # repetition
c=sample(1:length(LiceCounts),N, replace=TRUE)
t=sample(1:length(Bandwidth),N, replace=TRUE)
g=sample(1:length(gamma),N, replace=TRUE)
b=sample(1:nrow(beta),N, replace=TRUE)

r=1;
Output=numeric(); for(i in 1:N){Output[i]=ProportionalModel(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Number of sea lice","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df1 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
save(N, r, c, t, g, b, df1, Output, file='01_tidy_data/SA_result_1_ProportionalModel.RData')

r=2
Output=numeric(); for(i in 1:N){Output[i]=ProportionalModel(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Number of sea lice","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df2 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
save(N, r, c, t, g, b, df2, Output, file='01_tidy_data/SA_result_2_ProportionalModel.RData')

r=3
Output=numeric(); for(i in 1:N){Output[i]=ProportionalModel(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Number of sea lice","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df3 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
save(N, r, c, t, g, b, df3, Output, file='01_tidy_data/SA_result_3_ProportionalModel.RData')

# pearson SA including Route
r=sample(1:3,N, replace=TRUE)

Output=numeric(); for(i in 1:N){Output[i]=ProportionalModel(r[i], c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_r1_2 <- cor(r[which(r==1 | r==2)], Output[which(r==1 | r==2)], method="pearson")
correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson"); correlation_coefficient_r1_3 <- cor(r[which(r==1 | r==3)], Output[which(r==1 | r==3)], method="pearson")
input_names <- c("Number of sea lice","Bandwidth","Swimming speed", "Progression rate", "Route 2", "Route 3")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b, correlation_coefficient_r1_2, correlation_coefficient_r1_3)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
save(N, r, c, t, g, b, df, Output, file='01_tidy_data/SA_result_ProportionalModel.RData')

#############################################################################################
# results
load(file='01_tidy_data/SA_result_ProportionalModel.RData')
df$Variable <- reorder(df$Variable, abs(df$Coefficient))
load(file='01_tidy_data/SA_result_1_ProportionalModel.RData') 
df1=df1[-4,]; df1$Variable <- reorder(df1$Variable, abs(df1$Coefficient))
load(file='01_tidy_data/SA_result_2_ProportionalModel.RData') 
df2$Variable <- reorder(df2$Variable, abs(df2$Coefficient))
load(file='01_tidy_data/SA_result_3_ProportionalModel.RData')
df3$Variable <- reorder(df3$Variable, abs(df3$Coefficient))

library(ggplot2); library(gridExtra)
R0=ggplot(df,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.75,.75))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R1=ggplot(df1,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df1$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.75,.75))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R2=ggplot(df2,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df2$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.75,.75))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R3=ggplot(df3,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df3$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.75,.75))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))

library(patchwork)
tiff(filename = "03_plots/Fig 8.tiff", width = 19.05, height = 10, units = 'cm', res = 600)
R0+(R1/R2/R3)+plot_annotation(tag_levels='A')
dev.off()
