########## Before Model Fitting #########
#### Load neccessary packages ####
library(manipulate)
library(tidyverse)
library(readr)
library(deSolve)
library(readxl)
library(FME)
library(gridExtra) 

#### Load in data and do some cleaning ####

Phos<-read.csv("~/MathBio/Data/PhosphateAssayFullSecond.csv") #Read in phosphate assay
Phos$Date<-Phos$Date %>% as.Date("%m_%d") #correct date variable
start_dateP=Phos$Date[1] #Baseline date
start_dateP
Phos<-Phos %>% mutate(day=as.numeric(difftime(Date, start_dateP, units="days"))+1) #Days from start variable

Phos_5dc=Phos %>% #Create dataset for just 5 diluted and combined data
  filter(Sample=="5")%>%
  filter(Rep=="DandC")%>%
  select(ConcP,"day") 
Phos_5dc

Phos_17dc=Phos %>% #Create dataset for just 1.7 diluted and combined data
  filter(Sample=="1.7")%>%
  filter(Rep=="DandC")%>%
  select(ConcP,"day")
Phos_17dc
Phos_17dc$ConcP[c(1,2,3)]<-Phos_17dc$ConcP[c(1,2,3)]-0.4*10^(-3)


Phos_340dc=Phos %>% #Create dataset for just .34 diluted and combined data
  filter(Sample=="340")%>%
  filter(Rep=="DandC")%>%
  select(ConcP,"day")
Phos_340dc



algae_dilluted_data<- read_excel("~/MathBio/Data/Algae Data Diluted.xlsx") #Loading in count data
start_date=algae_dilluted_data$Date[1] #Count start date

algae_dilluted_data=algae_dilluted_data%>% #Create average count variable and days variable
  mutate(Scene_Average=(Scene_Count_1+Scene_Count_2+Scene_Count_3+Scene_Count_4
                        +Scene_Count_5+Scene_Count_6+Scene_Count_7+Scene_Count_8)/8,
         Clamy_Average=(Clamy_Count_1+Clamy_Count_2+Clamy_Count_3+Clamy_Count_4
                        +Clamy_Count_5+Clamy_Count_6+Clamy_Count_7+Clamy_Count_8)/8)%>%
  mutate(day=as.numeric(difftime(Date, start_date, units="days"))+1)

algae_dilluted_data_1=algae_dilluted_data%>% #reduce dataset dimension
  filter(Flask=="5DC")%>%
  select(Scene_Average, Clamy_Average, day)


algae_dilluted_data_17DC<-algae_dilluted_data %>% #1.7 diluted and combined subset
  filter(Flask=="1.7DC") %>%
  select(Scene_Average,Clamy_Average, day)


algae_dilluted_data_340DC<-algae_dilluted_data %>% #.34 diluted and combined subset
  filter(Flask=="0.34DC") %>%
  select(Scene_Average,Clamy_Average, day)

#### manipulate function to aid model fitting process ####

 manipulate({
   
######Define the Diff Eqs
   compderivs=function(t,pop,parms) {
     with(as.list(c(pop,parms)), {
       #"Known" Quantities
       S=pop[1]
       C=pop[2]
       QS=pop[3]
       QC=pop[4]
       Pf=pop[5]
       
      #Parameters
       muS=3.2
       QminS=parms[1]
       alphaS= parms[2]
       VmaxS=parms[3]
       KmS=parms[4]
       QmaxS=parms[5]
       muC=3.2
       QminC=parms[6]
       alphaC= 0.071
       VmaxC=parms[7]
       KmC=parms[8]
       QmaxC=parms[9]
       
       
       VS=VmaxS*(Pf/(KmS+Pf))*(QmaxS-QS)/(QmaxS-QminS)
       VC=VmaxC*(Pf/(KmC+Pf))*(QmaxC-QC)/(QmaxC-QminC)
       
       #Differential Equations
       dS=muS*S*(1-(QminS/QS))-alphaS*S
       dC=muC*C*(1-(QminC/QC))-alphaC*C
       dQS= VS-muS*(QS-QminS)
       dQC= VC-muC*(QC-QminC)
       dPf=-VS*S-VC*C
 
       
       #Return List
       list(c(dS,dC,dQS,dQC,dPf))
     }
     )}
   
   #Initial Conditions
   #Units: 
   yinitcomp=c(0.1*10^4,0.1*10^4,9*10^(-15),5*10^(-15), 5.642936e-09)
   #Parameters determined by manipulate
   parms=c(QminSs,alphaSs,VmaxSs,KmSs,QmaxSs,
           QminCs,VmaxCs,KmCs,QmaxCs)
   t=seq(0,19, length=1000) #Time sequence
   
   #Solve the Diff Eq with this inputs!
   outcomp=ode(y=yinitcomp,times=t,func = compderivs, parms=parms)
   
   
   dayC=algae_dilluted_data_1$day
   averageC=algae_dilluted_data_1$Clamy_Average
   
   dayS=algae_dilluted_data_1$day
   averageS=algae_dilluted_data_1$Scene_Average
   
   par(mfrow=c(2,2))
   plot(outcomp[,1],outcomp[,2],type='l',ylim = c(0,15*10^4), main="Scene Counts")
   points(dayS,averageS*10^4)
   plot(outcomp[,1],outcomp[,3],type='l',ylim = c(0,5*10^4), main="Chlamy Counts")
   points(dayC,averageC*10^4)
   
   plot(outcomp[,1],outcomp[,6],type='l',ylim = c(0,6.5e-09), main="Free Phosphorus")
   points(Phos_5dc$day, Phos_5dc$ConcP*10^(-3))
   
   print(parms)
   
 },QminSs=slider(1*10^(-15),5*10^(-15),initial =1.6*10^(-15),step=0.01*10^(-15)), 
 alphaSs=slider(0,5,initial=1,step=0.01),
 VmaxSs=slider(1*10^(-15),5*10^(-15),initial =1.25*10^(-15),step=0.01*10^(-15)),
 KmSs=slider(0.1*10^(-9),10*10^(-9),initial =5*10^(-9),step=0.01*10^(-9)),
 QmaxSs=slider(5*10^(-15),50*10^(-15),initial =9.1*10^(-15),step=0.01*10^(-15)),
 QminCs=slider(1*10^(-15),10*10^(-15),initial =1.6*10^(-15),step=0.01*10^(-15)), 
 VmaxCs=slider(1*10^(-15),5*10^(-15),initial = 1.25*10^(-15),step=0.01*10^(-15)),
 KmCs=slider(0.1*10^(-9),10*10^(-9),initial =5*10^(-9),step=0.01*10^(-9)),
 QmaxCs=slider(3*10^(-15),50*10^(-15),initial =5*10^(-15),step=0.01*10^(-15))
 )

########## Model Fitting Code #########
#### Fitting the DE model ####


dayC=algae_dilluted_data_1$day #create day vector (clamy)
averageC=algae_dilluted_data_1$Clamy_Average*10^4 #create properly scaled average count vector (clamy)

dayS=algae_dilluted_data_1$day #create day vector (scene)
averageS=algae_dilluted_data_1$Scene_Average*10^4 #create properly scaled average count vector (scene)

#### Define the Diff Eqs ######
compderivs=function(t,pop,parms) {
  with(as.list(c(pop,parms)), {
    #"Known" Quantities
    S=pop[1] #scene count
    C=pop[2] #clamy count
    QS=pop[3] #scene uptake (imputed)
    QC=pop[4] #clamy uptake (imputed)
    Pf=pop[5] #free phosphourus 
    
    #Parameters (See DE model in paper)
    muS=parms[1] 
    QminS=parms[2]
    alphaS= parms[3]
    VmaxS=parms[4]
    KmS=parms[5]
    QmaxS=parms[6]
    muC=parms[7]
    QminC=parms[8]
    alphaC= parms[9]
    VmaxC=parms[10]
    KmC=parms[11]
    QmaxC=parms[12]
    
    
    #Differential Equations (in a little different order so that R doesn't bug out)
    
    VS=VmaxS*(Pf/(KmS+Pf))*(QmaxS-QS)/(QmaxS-QminS) 
    VC=VmaxC*(Pf/(KmC+Pf))*(QmaxC-QC)/(QmaxC-QminC)
    
    dS=muS*S*(1-(QminS/QS))-alphaS*S
    dC=muC*C*(1-(QminC/QC))-alphaC*C
    dQS= VS-muS*(QS-QminS)
    dQC= VC-muC*(QC-QminC)
    dPf=-VS*S-VC*C
    
    
    #Return List
    list(c(dS,dC,dQS,dQC,dPf))
  }
  )}

#initial starting value for S, C, Qs, Qc, Pf
yinitcomp=c(0.1*10^4,0.1*10^4,9*10^(-15),5*10^(-15), 5.642936e-09) 


tcomp=unique(algae_dilluted_data$day)
tcomp

#### function to calculate squared errors ####
# multiply by a normalized standard error between compartments. 
sse.comp=function(parms0){  
  muS=parms0[1] 
  QminS=parms0[2]
  alphaS= parms0[3]
  VmaxS=parms0[4]
  KmS=parms0[5]
  QmaxS=parms0[6]
  muC=parms0[7]
  QminC=parms0[8]
  alphaC= parms0[9]
  VmaxC=parms0[10]
  KmC=parms0[11]
  QmaxC=parms0[12]
  out1comp=ode(y = yinitcomp, times = tcomp, func = compderivs, parms=c(muS,QminS,alphaS,VmaxS,KmS,QmaxS,
                                                                        muC,QminC,alphaC,VmaxC,KmC,QmaxC),
               method="ode45")
  # This is where model fits from the line directly above are compared to the actual data
  return(c(((out1comp[,2]-averageS)^2),(out1comp[,3]-averageC)^2*2,((out1comp[,6]-Phos_5dc$ConcP*10^(-3))*10^13.5)^2))   #sum of squared errors
}

########## Actually Fitting Models ########
#### Fit a model for 5 DC ####

# initial parameter estimates for this flask (drawn from modFit exploration)
parms5=c(3,1.42*10^(-15),0.5,4.09*10^(-15),1.97*10^(-9),10*10^(-15),
         3,9.28*10^(-15), 0.071,4.63*10^(-15),6.09*10^(-9),10*10^(-15))

# fitting the model from these estimates
fitcomp_5=modFit(f=sse.comp, p=parms5,lower=c(0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),5*10^(-15),
                                              0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),3*10^(-15)),
                 upper=c(10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15),
                         10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15)))



fitcomp_5$par
t=seq(0,19, length=1000)
out2comp_5=ode(y=yinitcomp,times=t,func = compderivs,parms=fitcomp_5$par)


#### Plotting model fits from above section ####
par(mfrow=c(2,2))
plot(out2comp_5[,1],out2comp_5[,2],type='l',ylim = c(0,15*10^4), main="Scene Counts")
points(dayS,averageS)
plot(out2comp_5[,1],out2comp_5[,3],type='l',ylim = c(0,5*10^4), main="Chlamy Counts")
points(dayC,averageC)

plot(out2comp_5[,1],out2comp_5[,6],type='l',ylim = c(0,6.5e-09), main="Free Phosphorus")
points(Phos_5dc$day, Phos_5dc$ConcP*10^(-3))

plot(out2comp_5)

#### Fitting a model for 1.7 DC data ####



dayC=algae_dilluted_data_17DC$day
averageC=algae_dilluted_data_17DC$Clamy_Average*10^4

dayS=algae_dilluted_data_17DC$day
averageS=algae_dilluted_data_17DC$Scene_Average*10^4

###### Define the Diff Eqs
compderivs=function(t,pop,parms) {
  with(as.list(c(pop,parms)), {
    #"Known" Quantities
    S=pop[1]
    C=pop[2]
    QS=pop[3]
    QC=pop[4]
    Pf=pop[5]
    
    #Parameters
    muS=parms[1]
    QminS=parms[2]
    alphaS= parms[3]
    VmaxS=parms[4]
    KmS=parms[5]
    QmaxS=parms[6]
    muC=parms[7]
    QminC=parms[8]
    alphaC= parms[9]
    VmaxC=parms[10]
    KmC=parms[11]
    QmaxC=parms[12]
    
    
    VS=VmaxS*(Pf/(KmS+Pf))*(QmaxS-QS)/(QmaxS-QminS)
    VC=VmaxC*(Pf/(KmC+Pf))*(QmaxC-QC)/(QmaxC-QminC)
    
    #Differential Equations
    dS=muS*S*(1-(QminS/QS))-alphaS*S
    dC=muC*C*(1-(QminC/QC))-alphaC*C
    dQS= VS-muS*(QS-QminS)
    dQC= VC-muC*(QC-QminC)
    dPf=-VS*S-VC*C
    
    
    #Return List
    list(c(dS,dC,dQS,dQC,dPf))
  }
  )}

yinitcomp=c(0.3755*10^4,0.25*10^4,9*10^(-15),5*10^(-15), 2.2e-06)

tcomp=unique(algae_dilluted_data$day)
tcomp
sse.comp=function(parms0){  #function to calculate squared errors
  muS=parms0[1]
  QminS=parms0[2]
  alphaS= parms0[3]
  VmaxS=parms0[4]
  KmS=parms0[5]
  QmaxS=parms0[6]
  muC=parms0[7]
  QminC=parms0[8]
  alphaC= parms0[9]
  VmaxC=parms0[10]
  KmC=parms0[11]
  QmaxC=parms0[12]
  out1comp=ode(y = yinitcomp, times = tcomp, func = compderivs, parms=c(muS,QminS,alphaS,VmaxS,KmS,QmaxS,
                                                                        muC,QminC,alphaC,VmaxC,KmC,QmaxC),
               method="ode45")
  return(c(((out1comp[,2]-averageS)^2*averageC/averageS),(out1comp[,3]-averageC)^2,((out1comp[,6]-Phos_17dc$ConcP*10^(-3)))^2*mean(averageC)/(Phos_17dc$ConcP*10^(-3))))   #sum of squared errors
}
parms17=c(3,1.42*10^(-15),0.5,1.09*10^(-15),1.97*10^(-9),10*10^(-15),
          3,9.28*10^(-15), 0.071,4.63*10^(-15),6.09*10^(-9),10*10^(-15))
fitcomp_17=modFit(f=sse.comp, p=parms17,lower=c(0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),5*10^(-15),
                                                0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),3*10^(-15)),
                  upper=c(10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15),
                          10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15)))




fitcomp_5$par
t=seq(0,19, length=1000)
out2comp_17=ode(y=yinitcomp,times=t,func = compderivs,parms=fitcomp_17$par)
out2comp_17<- out2comp_17 %>% as.matrix() %>% as.data.frame()
out2comp_17<- out2comp_17 %>% rename(Scene='1',Chlamy='2',QS='3',QC='4',Pf='5')

#### Plotting model for above section ####


ggP<-ggplot() +
  geom_point(data=Phos_17dc,aes(x=day,y=ConcP*10^(-3))) +
  geom_line(data=out2comp_17,aes(x=time,y=Pf)) 

ggCell<-ggplot() +
  geom_point(data=algae_dilluted_data_17DC,aes(x=day,y=Scene_Average*10^4),color="red") +
  geom_line(data=out2comp_17,aes(x=time,y=Scene),color="red") +
  geom_point(data=algae_dilluted_data_17DC,aes(x=day,y=Clamy_Average*10^4),color="blue") +
  geom_line(data=out2comp_17,aes(x=time,y=Chlamy),color="blue")

ggQ<- ggplot() +
  geom_line(data=out2comp_17,aes(x=time,y=QS),color="yellow") +
  geom_line(data=out2comp_17,aes(x=time,y=QC),color="green")


grid.arrange(ggCell,ggQ,ggP)


par(mfrow=c(2,2))
plot(out2comp_5[,1],out2comp_5[,2],type='l',ylim = c(0,15*10^4), main="Scene Counts")
points(dayS,averageS)
plot(out2comp_5[,1],out2comp_5[,3],type='l',ylim = c(0,20*10^4), main="Chlamy Counts")
points(dayC,averageC)

plot(out2comp_5[,1],out2comp_5[,6],type='l',ylim = c(0,3e-6), main="Free Phosphorus")
points(Phos_17dc$day, Phos_17dc$ConcP*10^(-3))


#### Fitting model for 0.34 data ####



dayC=algae_dilluted_data_340DC$day
averageC=algae_dilluted_data_340DC$Clamy_Average*10^4

dayS=algae_dilluted_data_340DC$day
averageS=algae_dilluted_data_340DC$Scene_Average*10^4

yinitcomp340=c(0.3755*10^4,0.25*10^4,9*10^(-15),5*10^(-15), 400e-09)



sse.comp=function(parms0){  #function to calculate squared errors
  muS=parms0[1]
  QminS=parms0[2]
  alphaS= parms0[3]
  VmaxS=parms0[4]
  KmS=parms0[5]
  QmaxS=parms0[6]
  muC=parms0[7]
  QminC=parms0[8]
  alphaC= parms0[9]
  VmaxC=parms0[10]
  KmC=parms0[11]
  QmaxC=parms0[12]
  out1comp=ode(y = yinitcomp340, times = tcomp, func = compderivs, parms=c(muS,QminS,alphaS,VmaxS,KmS,QmaxS,
                                                                           muC,QminC,alphaC,VmaxC,KmC,QmaxC),
               method="ode45")
  return(c(((out1comp[,2]-averageS)^2*averageC/averageS),(out1comp[,3]-averageC)^2,((out1comp[,6]-Phos_340dc$ConcP*10^(-3)))^2*mean(averageC)/(Phos_340dc$ConcP*10^(-3))))   #sum of squared errors
}
parms340=c(3,1.42*10^(-15),0.5,1.09*10^(-15),1.97*10^(-9),10*10^(-15),
           3,9.28*10^(-15), 0.071,4.63*10^(-15),6.09*10^(-9),10*10^(-15))
fitcomp_340=modFit(f=sse.comp, p=parms340,lower=c(0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),5*10^(-15),
                                                  0,1*10^(-15),0,1*10^(-15),0.1*10^(-9),3*10^(-15)),
                   upper=c(10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15),
                           10,10*10^(-15),5,10*10^(-15),10*10^(-9),50*10^(-15)))




fitcomp_340$par
t=seq(0,19, length=1000)
out2comp_340=ode(y=yinitcomp340,times=t,func = compderivs,parms=fitcomp_340$par)
out2comp_340<- out2comp_340 %>% as.matrix() %>% as.data.frame()
out2comp_340<- out2comp_340 %>% rename(Scene='1',Chlamy='2',QS='3',QC='4',Pf='5')

#### Plotting models for 0.34 data ####

ggP<-ggplot() +
  geom_point(data=Phos_340dc,aes(x=day,y=ConcP*10^(-3))) +
  geom_line(data=out2comp_340,aes(x=time,y=Pf)) 

ggCell<-ggplot() +
  geom_point(data=algae_dilluted_data_340DC,aes(x=day,y=Scene_Average*10^4),color="red") +
  geom_line(data=out2comp_340,aes(x=time,y=Scene),color="red") +
  geom_point(data=algae_dilluted_data_340DC,aes(x=day,y=Clamy_Average*10^4),color="blue") +
  geom_line(data=out2comp_340,aes(x=time,y=Chlamy),color="blue")

ggQ<- ggplot() +
  geom_line(data=out2comp_340,aes(x=time,y=QS),color="yellow") +
  geom_line(data=out2comp_340,aes(x=time,y=QC),color="green")


grid.arrange(ggCell,ggQ,ggP)





