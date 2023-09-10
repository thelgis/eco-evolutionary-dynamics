# @ author: Loukas Theodosiou
# @ comments: I think you can fix it and start thinking of data. 

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

#                       A1    
#                     .
#                    / V\
#                  / `  /
#                 <<   |
#                 /    |
#               /      |      PREY and PREDARTOR  - NO EVOLUTION
#             /        |                          - NO NETWORKS
#            /    \  \ /                          - NO PLASTICITY
#           (      ) | |
#   ________|   _/_  | |
#   <__________\______)\__)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

rm(list = ls())      # clean the history

######################################################################################################
# libraries that I need 
######################################################################################################

# CRAN Installations
if (!require("deSolve", quietly = TRUE)) install.packages("deSolve")
if (!require("WaveletComp", quietly = TRUE)) install.packages("WaveletComp")
if (!require("pracma", quietly = TRUE)) install.packages("pracma")
if (!require("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!require("changepoint", quietly = TRUE)) install.packages("changepoint")
if (!require("igraph", quietly = TRUE)) install.packages("igraph")
if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")

# Uncomment if needed:
# Bioconductor Installations
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!require("", quietly = TRUE)) BiocManager::install("")

library(deSolve)
library(WaveletComp) #Fourier transformation
library(pracma)      #Find the period and amplitude
library(tidyverse)
library(changepoint) # package to find transient 
library(igraph)      # package for the network to find the clusters
library(cowplot)

######################################################################################################
# parameters, the orders of the parameters is important
######################################################################################################
#### c has to be bigger than 1/2 .... c>1/2 ######

theta=0.3
phi= 2.5
c=1
parameters=list(theta, phi, c)

######################################################################################################
# time
######################################################################################################

time2 <- seq(0,100, by=0.1)

######################################################################################################
# initial values
######################################################################################################

yini=c(H=1, P=0.1)

######################################################################################################
# A. function -  Predator and Prey with no Evolution 
######################################################################################################

library(deSolve)

pp1 <- function(t, yini, parameters) {
  with(as.list(c(yini, parameters)), {
    
    dH=H*(1-theta*H)-(P*H/(1+H))
    
    dP=((phi*P*H)/(1+H))-(c*P)
    
    return(list(c(dH,dP)))  
  })
}

######################################################################################################
# solution
######################################################################################################
#out <- ode(y = yini, times = time2, func = coexist, parms = parameters, method="lsoda", atol=10^-15,rtol=10^-15)
out <- ode(y = yini, times = time2, func = pp1, parms = parameters, method="ode45")

######################################################################################################
# analysis
######################################################################################################

out.df = as.data.frame(out)
head(out.df)
tail(out.df)

##### check if you have any negative values in your data frame #######

out.df[out.df$H<0,"H"]
out.df[out.df$P<0,"P"]

######################################################################################################
# plot the predator-prey population dynamics
######################################################################################################

#dat1 <- out.df[1:2501,] #I select only the time points from 1 to 25
dat1=out.df
colnames(dat1) <- c("time","prey","predator")
data <- reshape2::melt(dat1, id.vars="time")
colnames(data)=c("time", "species","pop.size")

ggplot(data, aes(x=time, y=log(pop.size), color=species), fill=species) +
  geom_line(size=1.5)+
  scale_color_manual(values= c("#588300","#1380A1"), labels=c("prey","predator"),
                     guide= guide_legend(title="Species"))+
  theme_bw()+  
  ylab("Population Size in log scale")+
  ggtitle("Predator-Prey Dynamics (no prey evolution)")+
  #ylim(-25,50)+
  #scale_x_continuous(name="Time(Day)",breaks=c(1,2,4,5,7),labels=c(1,2,4,5,7))+
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold", family="Avenir"), 
        legend.title.align=0.5,
        panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
        axis.ticks.x = element_line(colour = "#333333"),
        axis.ticks.y = element_line(colour = "#333333"),
        axis.ticks.length =  unit(0.26, "cm"))


#######################################################################################################
### WAVELET ANALYSIS ###
#          |   _   _        
#    . | . x .|.|-|.|   1.Smooth time series and transform the time-series to Fourier 
# |\ ./.\-/.\-|.|.|.|   2.Do wavelet transformation (Fourier transformation) and find max power
# |.|_|.|_|.|.|.|_|.|   3.For the max-power find the period and amplitude of prey and predator

#######################################################################################################

prey_smooth= smooth.spline(out.df$time,(out.df$H),spar=0.5)      #smooth the host population dynamics  
predator_smooth= smooth.spline(out.df$time,(out.df$P),spar=0.5)  #smooth the prey population dynamics

#######################################################################################################
#For the prey| Transform to Fourier and find the max power 
#######################################################################################################

plot(prey_smooth$x,prey_smooth$y,type ='l')                      #plot the smoothed data over time 
prey.sm <- data.frame(day=prey_smooth$x,popul=prey_smooth$y)     #make the smoothed data to data.frame

#Fourier
wave.prey =analyze.wavelet(prey.sm,"popul")          #Do the Fourier transformation with analyze.wavelet
wt.image(wave.prey)

#find the max power
max.power.prey=which.max(wave.prey[[9]])            # asks for the position that we see max power
period.prey=(wave.prey$Period[max.power.prey])      # the actual period
period.prey=period.prey/10
max.power.prey.amp=(wave.prey$Ampl[,max.power.prey]) # ask for the amplitudes at the dominat period 
ampl.prey=max(max.power.prey.amp) # gives you the amplitude value
ampl.prey = ampl.prey*(prey_smooth$y)
ampl.prey=max(ampl.prey)

#######################################################################################################
#For the predator| Transform to Fourier and find the max power 
#######################################################################################################

plot(predator_smooth$x,predator_smooth$y,type ='l')                      #plot the smoothed data over time 
predator.sm <- data.frame(day=predator_smooth$x,popul=predator_smooth$y)     #make the smoothed data to data.frame

#Fourier
wave.predator =analyze.wavelet(predator.sm,"popul")          #Do the Fourier transformation with analyze.wavelet
wt.image(wave.predator)

#find the max power
max.power.predator=which.max(wave.predator[[9]])            # asks for the position that we see max power
period.predator=(wave.predator$Period[max.power.predator])      # the actual period
period.predator=period.predator/10
max.power.predator.amp=(wave.predator$Ampl[,max.power.predator]) # ask for the amplitudes at the dominat period 
ampl.predator=max(max.power.predator.amp) # gives you the amplitude value
ampl.predator = ampl.predator*(predator_smooth$y)
ampl.predator=max(ampl.predator)

#######################################################################################################
# put the amplitude and period together 
#######################################################################################################

type=c("prey", "predator")
period=c(period.prey, period.predator)
amplitude=c(ampl.prey, ampl.predator)

eco.params=data.frame(type, period, amplitude)
eco.params

#################################################################################
# phase plane plot | plot(prey, predator)
#################################################################################

ggplot(out.df, aes(out.df$H, out.df$P)) + 
  geom_point(size=0.5, color="#A426A6") + #use a geom_point with a nice purple color
  theme_bw()+                             #theme_bw leaves some grids inside but minimla
  theme_cowplot()+                        #theme_cowplot removes all the grid lines           
  xlab("Prey density")+                   #to name x-axis
  ylab("Predator density")+               #to name y-axis
  ggtitle("Phase predator-prey (with no evolution)")+
  theme(text=element_text(family = "Avenir"), #choose the font for all elements of the graph
        panel.grid.major = element_blank(),   # no panel grids
        panel.grid.minor = element_blank(),   
        panel.border = element_rect(colour = "black", size=0.5), #control the border of the graph
        plot.title = element_text(hjust=0.5, size=14,face="bold")) #control the title

#########################################################################################################

#           'x|`
#         '|xx| `          '|x|                 USE the phaseR package:
# `   '    |xx|    `   '    |x|`                1.set up the equation to the packages needs
#          |xx|             |x|                 2.find velocity field
#   ============|===============|===--          3.nullclines
#     ~~~~~|xx|~~~~~~~~~~~~~|x|~~~ ~~  ~   ~    4.trajectories of multiple initial points

#########################################################################################################

library(phaseR)   #this library requires in a way to resolve the equation but in a different way
y=c(H=1, P=0.1)   #y are the dependent values that are changing over time
parameters=list(  #the parameters
  theta, phi, c)

pred.prey <- function(t,y,parameters){
  x <- y[1]       #rewrite the function of the differential equations
  y <- y[2]       #denote as x the H and as y the P  
  theta <- parameters[1] #parameters
  phi<- parameters[2]
  c <- parameters[3]
  dy <- numeric(2) #I do not know what is this yet but I got it from this page 
  # (https://journal.r-project.org/archive/2014/RJ-2014-023/RJ-2014-023.pdf) pp.48
  dy[1]<-x*(1-theta*x)-(y*x/(1+x)) #here I write the equations
  dy[2]<-((phi*y*x)/(1+x))-(c*y)
  list(dy)
  
}

#########################################################################################################
# plot the velocity field, the nullclines, and trajectory from different starting points
#########################################################################################################

#solve and find the velocity field
pred.prey.field <- flowField(pred.prey, xlim = c(0, 3), ylim = c(0, 5),parameters=c(0.3, 2.5, 1), points=19,
                             add=FALSE, xlab="prey", ylab="predator", family="Avenir")
title("velocity field")

#identify and plot the nullclines
pred.prey.nullclines <-nullclines(pred.prey, xlim = c(-1, 3), ylim = c(-1, 5),
                                  parameters = c(0.3, 2.5, 1), points = 500, 
                                  col=c("#588300","#1380A1"), family="Avenir",lwd=c(2,2))

#plot the trajectory from different starting points
pred.prey.trajectory <-trajectory(pred.prey, y0 = c(1,1), tlim =c(0, 5),
                                  parameters = c(0.3, 2.5, 1), col="#A426A6", family="Avenir", lwd=2)


#########################################################################################################
# find the stability points
#########################################################################################################

pred.prey.stability.1 <-stability(pred.prey, ystar = c(1, 0.2), parameters=c(0.3, 2.5, 1))

#########################################################################################################

#  _   |~  _
# [_]--'--[_]   BIFURCATION diagram
# |'|""`""|'|   
# | | /^\ | |   How different θ values alter the amplitude and period of predator and prey
# |_|_|I|_|_|

#########################################################################################################

bifur=data.frame()                                     # opening a data frame to store my results
for(i in 1:100){
  theta.seq <- seq(0, 1, length=100)                   #here I create a sequence  
  theta <- c(sample(theta.seq, size=1, replace=TRUE))  #here I choose one sample from the sequence
  out <- ode(y = yini, times = time2, func = pp1, parms = parameters, method="ode45")
  out.df = as.data.frame(out)                                      #run the diff.equation and make it a data.frame
  ## prey 
  prey_smooth= smooth.spline(out.df$time,(out.df$H),spar=0.5)      #smooth the host population dynamics  
  prey.sm <- data.frame(day=prey_smooth$x,popul=prey_smooth$y)     #make the smoothed data to data.frame
  wave.prey =analyze.wavelet(prey.sm,"popul")          #Do the Fourier transformation with analyze.wavelet
  max.power.prey=which.max(wave.prey[[9]])            # asks for the position that we see max power
  period.prey=(wave.prey$Period[max.power.prey])      # the actual period
  period.prey=period.prey/10
  max.power.prey.amp=(wave.prey$Ampl[,max.power.prey]) # ask for the amplitudes at the dominat period 
  ampl.prey=max(max.power.prey.amp) # gives you the amplitude value
  ampl.prey = ampl.prey*(prey_smooth$y)
  ampl.prey=max(ampl.prey)
  ## predator
  predator_smooth= smooth.spline(out.df$time,(out.df$P),spar=0.5)  #smooth the prey population dynamics
  predator.sm <- data.frame(day=predator_smooth$x,popul=predator_smooth$y)     #make the smoothed data to data.frame
  wave.predator =analyze.wavelet(predator.sm,"popul")          #Do the Fourier transformation with analyze.wavelet
  max.power.predator=which.max(wave.predator[[9]])            # asks for the position that we see max power
  period.predator=(wave.predator$Period[max.power.predator])      # the actual period
  period.predator=period.predator/10
  max.power.predator.amp=(wave.predator$Ampl[,max.power.predator]) # ask for the amplitudes at the dominat period 
  ampl.predator=max(max.power.predator.amp) # gives you the amplitude value
  ampl.predator = ampl.predator*(predator_smooth$y)
  ampl.predator=max(ampl.predator)
  
  bifur[i,1] <- theta
  bifur[i,2] <- ampl.prey
  bifur[i,3] <- period.prey
  bifur[i,4] <- ampl.predator
  bifur[i,5] <- period.predator
  
}

write.csv(bifur, file = "bifurcation.csv")

bifur=read.csv("bifurcation.csv",header=TRUE)
bifur=bifur[,-1];colnames(bifur) <- c("theta", "ampl.prey", "period.prey", "ampl.predator", "period.predator")
bifur=subset(bifur, theta!=0)






ggplot() + 
  geom_point(data=bifur, aes(theta, ampl.prey, color='ampl.prey')) +
  geom_point(data=bifur, aes(theta, ampl.predator, colour='ampl.predator'))+
  geom_point(data=bifur, aes(theta, period.prey, colour='period.prey'))+
  geom_point(data=bifur, aes(theta, period.predator, colour='period.predator'))+
  scale_color_manual(name="Trigonometry", values=c("ampl.prey"="#588300", "ampl.predator"="#1380A1","period.prey"="#8bd000", "period.predator"="#1cb6e4"))+
  theme_bw()+
  theme_cowplot()+
  xlab("theta (self regulation)")+                   #to name x-axis
  ylab("abudance")+               #to name y-axis
  ggtitle("The effect of theta in the Amplitude and Period")+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold", family="Avenir"), 
        legend.title.align=0.5,
        panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
        panel.border = element_rect(colour = "#333333", size=0.5), #control the border of the graph
        axis.ticks.x = element_line(colour = "#333333"),
        axis.ticks.y = element_line(colour = "#333333"),
        axis.ticks.length =  unit(0.26, "cm"),
        legend.position = "bottom")

