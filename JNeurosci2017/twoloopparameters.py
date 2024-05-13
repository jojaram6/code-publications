# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:14:15 2015

@author: User
"""
###### parameters
import numpy as np

exp = 'wm'#hanks'#'wm_correlation'#'locallongrange'#'locallongrange'#'hanks'#locallongrange'


if exp=='dm': 
   maxtime=2500
elif exp=='wm' or exp=='locallongrange':
     maxtime= 23500
elif exp=='wm_correlation':
     maxtime = 60000  

elif exp=='hanks':
    maxtime=900
elif exp=='luminance':
    maxtime=700
elif exp=='huk':
    maxtime = 1900
pfcdisconnect='no'
    


stepsize=1
#print np.shape(maxtime)
time=np.arange(0,maxtime,stepsize)
totaltrials = 1#1000#80#1000#500#1000#00#1050
trialnumber = np.arange(0,totaltrials)

alltrials=range(0,totaltrials)
allcoh1= [0,3.2,6.4,12.8,25.6, 51.2]#[-51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6, 51.2]
allinh = [0,1]#[0.4, 0.6, 0.8, 1]

totalcoh=len(allcoh1)

tau_s=60 ###niyogi value tau_s=100, johnthesis value tau_s=60
gE=0.3725
gI=0.1137
gext=0.0011
Ib=0.33#4#7#7#47#0.3347# johnthesis Ib=0.3347
I_e=0.012#0.012#0.0165
gam=0.641
mu=30.0
sigma= 0#0.003#0.003#0.009#0.009##0.009#0.009#0.003 for autocorrelation
tau0= 2.0


### parameters for PPC rate
stepJ=0.02
J_m11_min = 0.15
J_m11_max = 0.45
J_m21_min = 0
J_m21_max = 0.10

####parameters for newfig/parameter sweep
if exp == 'wm':
   if pfcdisconnect == 'yes':
        J1_set = np.arange(0.34,0.43,0.01)#np.arange(0.34,0.42,0.03)#[0.35]
        J2_set = np.arange(0.34,0.43,0.01)#[0]#np.arange(0.34,0.42,0.03)#[0.42]
        Imem_set = np.arange(0,0.1,0.002)#np.arange(0,0.86,0.02)
        r_all_vec = np.zeros((len(Imem_set),4,len(time)))

        
                     # for small/higher resolution
   else:
        J1_set = np.arange(0.34,0.43,0.01)#[0.35]
        J2_set = np.arange(0.34,0.43,0.01)#[0.42]
        Imem_set = np.arange(0,2,0.002)#np.arange(0,0.86,0.02)
        r_all_vec = np.zeros((len(Imem_set),4,len(time)))

if exp == 'wm_correlation':      
     J1_set = [0.35]
     J2_set = [0.4182]
     Imem_set = [0.0]
     r_all_vec = np.zeros((4,len(time)))

     
elif exp == 'dm':
    J1_set = np.arange(0.34,0.43, 0.01)#[0.34,0.36,0.38,0.40]#np.arange(0.35,0.42,0.02)#[0.35]
    J2_set = np.arange(0.34,0.43, 0.01)#[0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41]#np.arange(0.35,0.42,0.02)#[0.42]
    Imem_set = [3.2,6.4,12.8,25.6]#[0,3.2,6.4,12.8,25.6,51.2]#allcoh1
else:
    J1_set = np.arange(0.34,0.43, 0.01)#[0.34,0.36,0.38,0.40]#np.arange(0.35,0.42,0.02)#[0.35]
    J2_set = np.arange(0.34,0.43, 0.01)#[0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41]#np.arange(0.35,0.42,0.02)#[0.42]
    #Imem_set = [3.2,6.4,12.8,25.6]#
tdistractor_revfig = 1500

    
t_beg = 1000
t_end = 3000
mem_crit = 15.0
dis_crit = 10.0
dis_mem_vec = np.zeros((len(J1_set),len(J2_set),len(Imem_set)))
min_mem_vec = np.zeros((len(J1_set),len(J2_set),len(Imem_set)))
I_all_vec = np.zeros((len(Imem_set),4,totaltrials, len(time)))


####parameters for hanks simulation

c_inh=0# paramter for hanks experiment, strength of current to the other population
poissontotal=30
amp_clicks = 0.8
motion_offset = 0#-10 
t_target = 0
pulseduration=50
number_realizations = 60#60
norm=1
independent='no' #no means each acc and rate generate their own time vectors
points = [200,220,400]#[200,250,300,350]#np.array([200,250, 300,350])### time points to caclulate sigmoid encoding of accumulator
binnumber = 8 #number of bins for accumulator
binsize = 2
delay=[0,140]#120
amp_pulse = I_e
tolerance=1



s1_0=[0.1,0.1]
noise_0=Ib

currentcue=0.1
currentdistractor=0.1#0.1
timecue=200
timedistractor=1500#1500
cueduration = 100
distractorduration = cueduration
Ivector=np.arange(0,0.37,0.002)

t_luminance=100
t_reward=100 #190 for pulse
tau_luminance=300
tau_reward=80

#reactiontimes=np.zeros(totaltrials)
#discriminationtimes=np.zeros(totaltrials)
#onsettimes=np.zeros(totaltrials)
#inhpercent=1
xlimit=600#maxtimeylimit=3
allcurrent=['IL1', 'IR1', 'IL2', 'IR2']


##target current parameters
visual=90 #90 for 34
tau_visual=10.0

#tau_rise=10.0
#tau_decay=20.0
A_risedecay=0.115 #0.12 for 35 spk/s peak
tau_rise=9.0#10.0
tau_decay=10.0#11.0
#t_target=30.0
t_motion=30.0
##FIcurve parameters
a=270.0 # Hz/na
b=108.0 #Hz
d=0.154 #s
if exp=='huk':
    ahuk = 270.0
else:
    ahuk = 270.0
#J_m11vec=[0.35, 0.32, 0.3]
#    J_m22vec=[0.35, 0.32, 0.35]
#    J_m12vec=[0.085, 0.12, 0.085]
#    J_m21vec=[0.05, 0.085, 0.11]
######conductances
actualsetting= 'other'#'ppcrate'#discrimination'# 'bargraphs' for bargraphs: rev_fig for review figure
def weights(setting, localindex):
    J_m11vec= [0.375, 0.2, 0.32]  #[0.35, 0.31, 0.31]
    J_m22vec=  [0.375, 0.2, 0.38] #[0.35, 0.31, 0.35]
    J_m12vec= [0.05, 0.22, 0.065]  #[0.074, 0.12, 0.085]
    J_m21vec= [0.035, 0.2, 0.065]  #[0.055, 0.09, 0.08]
    
    
    if setting=='bargraphs':
        J_m11=J_m11vec[localindex]#0.36127#original value 0.36127
        J_m22=J_m22vec[localindex]#0.42127# original value 0.40127
        J_m12=J_m12vec[localindex]#0.085# original value 0.085
        J_m21=J_m21vec[localindex]# original value 0.05
    elif setting=='ppcrate':
         loc1,loc2=localindex
         
         J_m11vec = np.arange(J_m11_min, J_m11_max, stepJ)
         J_m21vec = np.arange(J_m21_min, J_m21_max, stepJ/2)
         
         J_m11 =J_m11vec[loc1]#0.36127#original value 0.36127
         J_m22= 0.4182#0.41#0.42127# original value 0.40127
         J_m12=0.15# original value 0.05      
         J_m21=J_m21vec[loc2]#0.085# original value 0.085
    elif setting=='discrimination':
        loc1 = localindex
        J_m11vec = np.arange(0.35, 0.42, 0.01)
        J_m11 =J_m11vec[loc1]
        J_m22 = 0
        J_m12 = 0
        J_m21 = 0
    elif setting=='rev_fig':
         loc1,loc2=localindex
         
         J_m11vec = J1_set
         J_m22vec = J2_set
         
         J_m11 =J_m11vec[loc1]#0.36127#original value 0.36127
         if pfcdisconnect=='yes':
             J_m21= 0.0#0.4182#0.42127# original value 0.40127
         else:
             J_m21= 0.04#0.4182#0.42127# original value 0.40127
         J_m12=0.15# original value 0.05      
         J_m22=J_m22vec[loc2]#0.085# original value 0.085
        
         
    else:
        J_m11=0.35# value for DM 0.38127#value fortime constants#0.34127#original value 0.36127
        J_m22=0.4#82#82#0.4182# 
        J_m12=0.15#0.15# original value 0.085
        if pfcdisconnect=='yes':
           J_m21=0
        else:
           
           J_m21=0.04# original value 0.04
    
       
        
    J_p11=0.28#3#87#0.28387#original value 0.28387
    J_p22=J_p11#0.28387# original value 0.28387
    J_p12=0# original value 0
    
    J_p21=0# original value 0
    J_L1L1= (J_p11+J_m11)/2
    J_R1R1=J_L1L1
    J_R1L1=(J_p11-J_m11)/2
    J_L1R1=J_R1L1
    J_L2L2=(J_p22+J_m22)/2
    J_R2R2=J_L2L2
    J_R2L2=(J_p22-J_m22)/2
    J_L2R2=J_R2L2
    J_L1L2=(J_p12+J_m12)/2
    J_R1R2=J_L1L2
    J_L1R2=(J_p12-J_m12)/2
    J_R1L2=J_L1R2
    J_L2L1=(J_p21+J_m21)/2
    J_R2R1=J_L2L1
    J_L2R1=(J_p21-J_m21)/2
    J_R2L1=J_L2R1
    
    #J_L1L1 = 0
    #J_L1L2 = 0
    #J_L1R1 = 0
    #J_L1R2 = 0
    
    
    return np.array([J_L1L1,J_L1R1, J_L2L2, J_R2L2, J_L1L2, J_L1R2, J_L2L1, J_L2R1, J_m11, J_m22, J_m12, J_m21])


#onset time, reaction times, discriminationtimes
RTthreshold=40
DTthreshold=12
OTthreshold=7

 
Ib1=np.zeros(len(time))
Ib1[0]=Ib
Ib2=np.zeros(len(time))
Ib2[0]=Ib
Ib3=np.zeros(len(time))
Ib3[0]=Ib
Ib4=np.zeros(len(time))
Ib4[0]=Ib

#Ib_p=np.ones(len(time))*Ib
#Ib2=Ib1
#Ib3=Ib1
#Ib4=Ib1
  
sL1=np.zeros(len(time))
sL1[0]=0.08#0.08
sR1=np.zeros(len(time))
sR1[0]=0.08#0.08
sL2=np.zeros(len(time))
sL2[0]=0.08#0.08
sR2=np.zeros(len(time))
sR2[0]=0.08#0.08


#
alpha = 0.6
colors_set1 = {'red':(228./255, 26./255, 28./255),
          'blue':(55./255, 126./255, 184./255),
          'cyan':(118./255, 202./255, 215./255),#          'cyan':(158./255, 202./255, 225./255),

          'green':(77./255, 175./255, 74./255),
          'lightgreen':(1-(1-77./255)*alpha, 1-(1-175./255)*alpha, 1-(1-74./255)*alpha),
          'lightorange':(1-(1-255./255)*alpha, 1-(1-127./255)*alpha, 1-(1-0./255)*alpha),


          'purple':(152./255, 78./255, 163./255),
          'orange':(255./255, 127./255, 0./255),
          'yellow':(255./255, 235./255, 151./255),
          'brown':(166./255, 86./255, 40./255),
          'pink':(247./255, 129./255, 191./255),
          'gray':(153./255,153./255,153./255)}
cmap=[colors_set1['purple'], colors_set1['blue'], colors_set1['green'], colors_set1['yellow'], colors_set1['orange'], colors_set1['red']]
inhcmap= [colors_set1['blue'], colors_set1['green'], colors_set1['orange'], colors_set1['red']]
cMap = [colors_set1['blue'], colors_set1['green'], 'cyan', colors_set1['orange'], colors_set1['red']]
 
#alpha = 0.2
#clight = [1-(1-ci)*alpha for ci in c]
#$print weights(actualsetting, [0,3]) 