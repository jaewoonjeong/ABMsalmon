############ Sensitivity Analysis

source(file='Raster map.R')
source(file='MigrationFunction.R')

d=DistanceMatrix/1000 # to change m to km
Bandwidth=c(5,10,20,30,40,60,100)
sigma=Bandwidth/4

dSum=ip=matrix(NA,nrow(d),nrow(SL_unique)*length(Bandwidth)); dim(dSum)=dim(ip)=c(nrow(d),nrow(SL_unique),length(Bandwidth))
for(b in 1:length(Bandwidth)){W=function(d){(1/(sqrt(2*pi)*(Bandwidth[b]/4)))*exp((-d^2)/(2*(Bandwidth[b]/4)^2))} # Gaussian kernel density
for(i in 1:nrow(d)){for(j in 1:nrow(SL_unique)){dSum[i,j,b]=W(d[i,j])/1000}};print(b)}
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(j in 1:nrow(SL_unique)){ip[i,j,b]=ifelse(sum(dSum[,j,b])==0,0,dSum[i,j,b]/sum(dSum[,j,b]))}};print(i)}

LiceCountsPerFarm=62500000
LiceCounts=LiceCountsPerFarm*8*c(0.5,1.5,3,6,8,10,12)/6 # Lice count were updated 

# To calculate IP for each lice count values
iipp=matrix(NA,nrow(d),nrow(SL_unique)*length(Bandwidth)*length(LiceCounts)); dim(iipp)=c(nrow(d),nrow(SL_unique),length(Bandwidth),length(LiceCounts))
IP=matrix(NA,nrow(d),length(Bandwidth)*length(LiceCounts)); dim(IP)=c(nrow(d),length(Bandwidth),length(LiceCounts))
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(j in 1:nrow(SL_unique)){for(s in 1:length(LiceCounts)){iipp[i,j,b,s]=ip[i,j,b]*LiceCounts[s]}}};print(i)}
for(i in 1:nrow(d)){for(b in 1:length(Bandwidth)){for(s in 1:length(LiceCounts)){IP[i,b,s]=sum(iipp[i,,b,s])}};print(i)}
PI=list(); for(s in 1:length(LiceCounts)){ppii=list(); for(b in 1:length(Bandwidth)){ppii[[b]]=cbind(PixelInfo,IP=IP[,b,s])}; PI[[s]]=ppii}
PI_SA=PI

#####################################################
# PR1= -beta1*55+87
pr1=rev(c(31.1, 31.1, 31.1, seq(31.1, 31.1/2, length.out=4)))
beta1=(86.1-pr1)/55
RR1=list()
for(r in 1){II=list()
for(i in 1:iteration){BB=list()
for(b in 1:length(beta1)){BB[[b]]=SA_function(MaxTime, r, beta1[b])}
II[[i]]=BB
print(i)}
RR1[[r]]=II}

# PR2= -beta2*15+33
pr2=c(12.5, 14.4, 16.2, 18.1, 19.8, 21.4, 23.1)
beta2=(33.1-pr2)/15
RR2=list()
for(r in 2){II=list()
for(i in 1:iteration){BB=list()
for(b in 1:length(beta2)){BB[[b]]=SA_function(MaxTime, r, beta2[b])}
II[[i]]=BB; print(i)}
RR2[[r]]=II}

# PR3= -beta3*20+30
pr3=c(8.7, 9.7, 10.8, 11.8, 13.2, 14.5, 15.9)
beta3=(33.7-pr3)/21.6
RR3=list()
for(r in 3){II=list()
for(i in 1:iteration){BB=list()
for(b in 1:length(beta3)){BB[[b]]=SA_function(MaxTime, r, beta3[b])}
II[[i]]=BB; print(i)}
RR3[[r]]=II}

RR=list(RR1[[1]], RR2[[2]], RR3[[3]])

beta=data.frame(beta1, beta2, beta3)
#################################################################
gamma=rev(c(0.20, 0.47, 0.73, 1.00, 1.17, 1.33, 1.5))

SA_infes=function(r,c,t,g,b){
  Trawl=function(IP){exp(-5.49+log(0.0032141*gamma[g])+0.425*log(IP))}
  if(r==1){b=4}
  NewInf=numeric()
  i=sample(1:iteration,1)
  for(n in 1:nrow(RR[[r]][[i]][[b]])){TP_SA=Trawl(PI_SA[[c]][[t]][which(RR[[r]][[i]][[b]][n,3]==PI_SA[[c]][[t]][,3] & RR[[r]][[i]][[b]][n,4]==PI_SA[[c]][[t]][,4]),]$IP)
  NewInf[n]=rnbinom(1,size=2.733,mu=TP_SA)}
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
Output=numeric(); for(i in 1:N){Output[i]=SA_infes(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Sea lice count","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df1 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
df1=df1[-4,]; df1$Variable <- reorder(df1$Variable, abs(df1$Coefficient))

r=2
Output=numeric(); for(i in 1:N){Output[i]=SA_infes(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Sea lice count","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df2 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
df2$Variable <- reorder(df2$Variable, abs(df2$Coefficient))

r=3
Output=numeric(); for(i in 1:N){Output[i]=SA_infes(r, c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson")
input_names <- c("Sea lice count","Bandwidth","Swimming speed", "Progression rate")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df3 <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
df3$Variable <- reorder(df3$Variable, abs(df3$Coefficient))

# pearson SA including Route
r=sample(1:3,N, replace=TRUE)

Output=numeric(); for(i in 1:N){Output[i]=SA_infes(r[i], c[i], t[i], g[i], b[i]);print(i)}
correlation_coefficient_c <- cor(c, Output, method="pearson"); correlation_coefficient_t <- cor(t, Output, method="pearson"); correlation_coefficient_r1_2 <- cor(r[which(r==1 | r==2)], Output[which(r==1 | r==2)], method="pearson")
correlation_coefficient_g <- cor(g, Output, method="pearson"); correlation_coefficient_b <- cor(b, Output, method="pearson"); correlation_coefficient_r1_3 <- cor(r[which(r==1 | r==3)], Output[which(r==1 | r==3)], method="pearson")
input_names <- c("Sea lice count","Bandwidth","Swimming speed", "Progression rate", "Route 2", "Route 3")
coefficients <- c(correlation_coefficient_c, correlation_coefficient_t, correlation_coefficient_g, correlation_coefficient_b, correlation_coefficient_r1_2, correlation_coefficient_r1_3)
df <- data.frame(Variable = input_names, Coefficient = coefficients)
df <- df[order(abs(df$Coefficient), decreasing = TRUE), ]  # Order by absolute value for ranking
df$Variable <- reorder(df$Variable, abs(df$Coefficient))
#############################################################################################
# results
R0=ggplot(df,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.5,.5))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R1=ggplot(df1,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df1$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.5,.5))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R2=ggplot(df2,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df2$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.5,.5))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
R3=ggplot(df3,aes(x=Variable,y=Coefficient))+geom_col(fill=ifelse(df3$Coefficient>0,"gray30","gray60"),width=0.6)+coord_flip()+labs(title="",x="",y="Pearson Correlation Coefficient")+theme_bw()+ggtitle('')+scale_y_continuous(limits=c(-.5,.5))+theme(axis.title.x = element_text(size = 7), axis.text = element_text(size = 7),plot.margin = unit(c(0, .1, 0, .1), "in"))
# Figure 7
R0+(R1/R2/R3)+plot_annotation(tag_levels='A')



