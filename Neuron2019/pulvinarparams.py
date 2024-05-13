import numpy as np

'''DM or WM
Conflict or no conflict
lateral Lesion
fixed or variable pulvinar
newvalues?
for confidence: no lesion, no pulvinar, loop-like
'''
#################### Simulation parameters ##############
#########################################################

exp = 'dm'#'hanks'#locallongrange'


if exp=='dm':
   maxtime = 800#1500 
elif exp=='wm' or exp=='locallongrange':
   maxtime= 2000
elif exp=='hanks':
    maxtime=900
elif exp=='luminance':
    maxtime=700

pfcdisconnect='no'
cx_impl = 'no'#'yes'    
usecluster = 'no'
exptype = 'confidence'#'conflict'#'wmgate'#'wmgate'#'wmgate'#'confidence'#'confidence'
allcoh1 = [0, 3.2, 6.4, 12.8, 25.6, 51.2]
#allcoh1 = [2, 9, 23]

trials_transition = 10#000#1000#500#500#
if exptype == 'confidence':
   coh1vec = [2,5,8]#[2,5,8]#[2, 9, 23]
   Ibconf = 0.36 #
   sigma = 0.007#0.008#0.004#
   pfcdisconnect = 'yes'
   changef = -0.008  ##### neural activity after action/saccade
   tchange = 670 ###### time of action/saccade
else:
    coh1vec = [20]#[4]
    Ibconf = 0.32
    sigma = 0.006#0.006#0.008#008#0.008	
    changef = 0#-0.008  ##### neural activity after action/saccade
    tchange = 700 ###### time of action/saccade				
trials_half = int(trials_transition/2)
pulse_change  = 0#-0.08
t_pulsechange = 130
transitionplot = 0
totaltrials=1050
alltrials=range(0,totaltrials)
allinh = [0,1]#[0.4, 0.6, 0.8, 1]
totalcoh =  len(coh1vec)#len(allcoh1)--
inhpercenttest=1
sim = 0
endtrial = 800
decision= 'conflict'#'conflict'#'visual'
RTepsilon=2
tolerance = RTepsilon
J_ratio = np.arange(0,5,0.2)#[0.5,1,1.5,2,2.5,3]
J_p = np.arange(0,2,0.2)
lambdas =260+45*np.arange(0,1.1,0.2)#[260, 275,290,305]# np.arange(60,320,40)#np.arange(200,320,20)#np.arange(60,320,40)##[0.0,300.0]#np.arange(250.0,335.0,15)#[200,250,300,350,400]
actualsetting=['loop-like','confidence']# 'bargraphs' for bargraphs 'loop_like' for confidence
if transitionplot==1:
    lenp = len(lambdas)
else:
    lenp = 1#len(lambdas)
if exptype=='confidence':
   lenp = 1
##########################################################
########Cortical circuit parameters ######################
##########################################################
##FIcurve parameters

a = 270.0#270.0 # Hz/na
b = 108.0 #Hz
d = 0.154 #s
gam = 0.641
Ib = 0.334#0.3335
tau0 = 2.0 # noise current
tau_s = 60 ### NMDA time constant wongwang, niyogi value tau_s=100, johnthesis value tau_s=60
a_mod1 = 1.0*a
a_mod2 = 1.0*a #1.3*a for supplementary figures cortical modulation
################################################################
########  Cortical input currents and parameters ################
################################################################
I_e= 0.0156#0.007#0.012#0.007 for confidence
currentcue = 0.11#0.08#080#.1
currentdistractor = currentcue#0.08# 0.1#0.1
timecue = 30##30
timedistractor = 800#timecue
cueduration=200
distractorduration=cueduration
motionduration = 2700


Ivector=np.arange(0,0.37,0.002)

t_luminance=100
t_reward=100 #190 for pulse
tau_luminance=300
tau_reward=80
##target current parameters
visual=90 #90 for 34
tau_visual=10.0
A_risedecay=0.115 #0.12 for 35 spk/s peak
tau_rise=13.0#9.0
tau_decay=14.0#10.0
t_target=timecue#30.0
t_motion=timecue
#onset time, reaction es, discriminationtimes
RTthreshold=40#30
DTthreshold=12
OTthreshold=8

#####################################################
################Pulvinar circuit parameters##########
####################################################
tau_p = 2.0 



les = 1 ## les = 1 is with intact pulvinar (no lesion)
pulvinars = 0 ## pulvinars = 1 means that there are two segregated pulvinar populations
p_number  = 1-pulvinars
newvalues = 0

#####Params for pulvinar transient######
s = 0.3
amp_conf = 6000
decay_conf = 30
offset_conf = 0#0.015#0.05
off_pulvmotion = 2000

####################################################################
############Confidence experiment parameters #####
####################################################################
disctime = 660
time_ini = disctime-200
time_end = disctime
norm_maxblack = 1/15.0
absfunction = 'no'#'yes' ## yes if using absolute function to calculate pulvinar activity
#escape_epsilon = 5#8#4#6#5#5
stepsize_switch = 1.0#.001
stepsize_switch = 1.0#.001
time_p = np.arange(0,maxtime,stepsize_switch)
lentime = int(maxtime/stepsize_switch)
Ib_p=np.ones(len(time_p))*1*Ibconf
conf_coupling = 0.05#0.15
stepsize=1.0
time=np.arange(0,maxtime,stepsize)


tau_pp=tau_p/1000.0


############################################################################################
##########Cortico-thalamic parameters for short-term plasticity#############################
############################################################################################
F_in = 0.0#initial probability
p = 0.45#0.45#/1000
tau_D = 0.6*1000
a_F = 0.35 #0.35
tau_F = 0.5*1000
tau_E = 0.005*1000
tau_I = 0.020*1000
rate_offset = 5#10#10.0
Iinh = 0#9
Iexc = 0#0.02#0.2

J_exc = 2.85#2.5#use 2.5 for previous results
J_inh = -2.6#-3.0#-1.6#-J_exc/1.31#-1.8#use -1.8 for previous results
exc_plasticity = 'yes'
inh_plasticity = 'yes'

S_rates = np.zeros((2,len(time)))
S_exc = np.zeros((2,len(time)))
S_inh = np.zeros((2,len(time)))
D = np.zeros((2,len(time)))
F = np.zeros((2,len(time)))
D[:,0] = 0.7
F[:,0] = 0.1
S_exc[0] = 0.0
S_inh[0] = 0.0

###### Kagan parameters
Jfix = 0#-.5#-20#-0.1
I_h = 0.35
t_shift = 0
slopeR = 0.002
slopeL = 0.002
tau_fix = 20.0
tau_h_dep = 20
tau_h_hyp = 100
microduration = 200
t_microstart = t_target+80
Ifix_amp = 0.2
micro_amp =  0#0.1
ipsi = 0
ipsi_amp =0

############################################################
######### Parameters for structured connections ############
############################################################

rec = 1 ## rec = 1 is excitatory recurrency

frac_Ib_Lp = 1
frac_Ib_Rp = frac_Ib_Lp

frac_Ib_Lpp =1
frac_Ib_Rpp = frac_Ib_Lpp

frac_reward=18#8 for plusreward condition
inh_p = -.81#-1#-.3
k_p = 1

k_1p = 1.8#1.2 #1.6#1 
k_p1 = 0.2#0.4
k_2p = 0.1#0.4
k_p2 = 1.8#1.2#1.6#1

coh1vec_ = [40]# not important, just the size

xlimit=600#maxtimeylimit=3
allcurrent=['IL1', 'IR1', 'IL2', 'IR2']

def weights(setting, localindex):
    setting_LR, setting_p = setting
    local_LR, local_lambda = localindex
    J_p = np.arange(0,2,0.2)
    

    if setting_LR == 'symmetric_variable' or setting_LR=='symmetric_fixed':
       J_m11= 0.38#0.37#0.4862#for conflict
       J_m22= 0.37#0.37#0.4862#  for conflict
       
       
       J_m21 = 0.02#0.05
       if setting_LR == 'symmetric_variable':
           J_m12 = J_ratio[localindex]*J_m21
       else:
           J_m12 = 0.05#0.05#0.05 #0.01 for conflict
       J_p11 = 0.2588#original value 0.28387
       J_p22 = 0.2588# original value 0.28387
       J_p12=0# original value 0
       J_p21=0# original value 0
       
       
         
    elif setting_LR == 'loop-like':
        J_m11 = 0.34#35# value for DM 0.38127#value fortime constants#0.34127#original value 0.36127
        J_m22 = 0.40#4182#82#0.492127# original value 0.40127
        if pfcdisconnect=='yes':
           J_m12=0
           J_m21=0
        else:
            if cx_impl=='yes':
               J_m12 = 0
            else:
               J_m12 = 0.03#0.1#14#5#0.15 (original loop model)
            J_m21 = 0.04#0.02#0.01#0.05
        J_p11 = 0.2588#0.28387#0.28387#0.28387#original value 0.28387
        J_p22 = 0.2588#0.28387# original value 0.28387
        J_p12 = 0# original value 0
        J_p21 = 0# original value 0
       
    
    
    if setting_p == 'yes_p_variable':
        lambda_p = lambdas[local_lambda]        
        #base_p = J_p[local_p]*0.4#0.5
    #else: 
    base_p = 0.28#0.28 #wilke# for conflict DM#0.5 #0.48 for confidence?
    # pulvinar recurrent inhibition
    J_LpLp = k_p*base_p
    J_RpRp = J_LpLp
    J_LpRp = 0#inh_p*J_LpLp
    J_RpLp = J_LpRp
    ##module 1 and pulvinar
    J_L1Lp = k_1p*base_p
    J_L1Rp = inh_p*J_L1Lp
    J_R1Lp = J_L1Rp
    J_R1Rp = J_L1Lp
    
    J_LpL1 = k_p1*base_p
    J_RpL1 = inh_p*J_LpL1 # was inh_p*J_L1Lp
    J_LpR1 = J_RpL1
    J_RpR1 = J_LpL1
    
    ##module 2 and pulvinar
    J_L2Lp = k_2p*base_p
    J_L2Rp = inh_p*J_L2Lp# was inh_p*J_L1Lp
    J_R2Lp = J_L2Rp
    J_R2Rp = J_L2Lp
    
    J_LpL2 = k_p2*base_p
    J_RpL2 = inh_p*J_LpL2
    J_LpR2 = J_RpL2
    J_RpR2 = J_LpL2
    
    
    ##########Individual cortical weights
     
    J_L1L1 = (J_p11+J_m11)/2# + delW_11[0,0]
    J_R1R1 = (J_p11+J_m11)/2# + delW_11[1,1]#J_L1L1
    J_R1L1 = (J_p11-J_m11)/2# + delW_11[0,1]
    J_L1R1 = (J_p11-J_m11)/2# + delW_11[1,0]#J_R1L1
    
    J_L2L2 = (J_p22+J_m22)/2# + delW_22[0,0]
    J_R2R2 = (J_p22+J_m22)/2# + delW_22[1,1]#J_L2L2
    J_R2L2 = (J_p22-J_m22)/2# + delW_22[0,1]
    J_L2R2 = (J_p22-J_m22)/2# + delW_22[1,0]#J_R2L2
    
    J_L1L2 = (J_p12+J_m12)/2# + delW_12[0,0]
    J_R1R2 = (J_p12+J_m12)/2# + delW_12[1,1]#J_L1L2
    J_L1R2 = (J_p12-J_m12)/2# + delW_12[1,0]
    J_R1L2 = (J_p12-J_m12)/2# + delW_12[0,1]#J_L1R2
    
    J_L2L1 = (J_p21+J_m21)/2# + delW_21[0,0]
    J_R2R1 = (J_p21+J_m21)/2# + delW_21[1,1]#J_L2L1
    J_L2R1 = (J_p21-J_m21)/2# + delW_21[1,0]
    J_R2L1 = (J_p21-J_m21)/2# + delW_21[0,1]#J_L2R1
    
    J_LR =  [J_L1L1,J_R1R1,J_R2R2,J_R2R1,J_R1R2,J_L1R1, J_L2L2, J_R2L2, J_L2R2,J_L1L2, J_L1R2, J_L2L1, J_L2R1, J_R1L1,J_R1L2, J_R2L1]   
    if setting_p == 'yes_p_variable' or setting_p == 'yes_p_fixed':
        J_p = [J_LpLp, J_RpRp,J_LpL1,J_RpL1,J_L1Lp,J_L1Rp, J_LpL2, J_RpL2, J_L2Lp, J_L2Rp,J_LpR1,J_RpR1,J_R1Lp, J_R1Rp,J_LpR2,J_RpR2,J_R2Lp,J_R2Rp,J_LpRp,J_RpLp]
        
    if setting_p == 'yes_p_fixed':
       lambda_p = 230#300#300#210#600#300#300#300
    elif setting_p == 'confidence':
        J_p = np.zeros(20)
        if les==1:
           J_LpL1 = conf_coupling#0#0.18#0.15#
        elif les==0:
             J_LpL1 = 0			
        J_LpR1 =J_LpL1
        J_p[2] = J_LpL1
        J_p[10] = J_LpR1
        lambda_p = 300
    elif setting_p == 'pcompetition':
         J_p = np.zeros(20)
         J_L1Lp = 0.154
	   
         J_LpL1 = 0.18
	    
         J_LpRp = -0.841
	   
         J_RpLp = J_LpRp
	 
         J_L1Rp = 0.14
	   
         J_RpL2 = 1.5#0.154
	   
         J_p[2] = J_LpL1
	   
         J_p[4] = J_LpL1
	    
         J_p[5] = J_L1Rp
	   
         J_p[7] = J_RpL2
	   
         J_p[18] = J_LpRp
	   
         J_p[19] = J_p[18]			
	   
         lambda_p = 300

        
				
    else:
        J_p = np.zeros(20)
        lambda_p = 0#280.0
    
    return [J_LR, J_p,lambda_p]

#Ibfix = 0
#############################################
################# Vector initialization ######
#############################################
rL1all = np.zeros((trials_transition,lenp,len(time)))
rL2all = np.zeros((trials_transition,lenp,len(time)))
rR1all = np.zeros((trials_transition,lenp,len(time)))
rR2all = np.zeros((trials_transition,lenp,len(time)))
rLpall = np.zeros((trials_transition,lenp,len(time)))
rRpall = np.zeros((trials_transition,lenp,len(time)))
rPall = np.zeros((trials_transition,len(time)))
rfixall = np.zeros((trials_transition,len(time)))
sL1all = np.zeros((trials_transition,len(time)))
sL2all = np.zeros((trials_transition,len(time)))
sLpall = np.zeros((trials_transition,len(time)))

sR1all = np.zeros((trials_transition,len(time)))
sR2all = np.zeros((trials_transition,len(time)))
sfixall = np.zeros((trials_transition,len(time)))

vec_trial = np.zeros((trials_transition))
winner_vec = np.zeros((trials_transition))

all_latencies = np.zeros((trials_transition))
vec_ratio = np.zeros((1,lenp,len(coh1vec)))

a_I = 1.1
Ib1=np.zeros(len(time))
Ib1[0] = Ib#weights(actualsetting, [0,0])[3][0]#Ib
Ib2=np.zeros(len(time))
Ib2[0] = Ib#weights(actualsetting, [0,0])[3][1]#Ib
Ib3=np.zeros(len(time))
Ib3[0] = Ib#weights(actualsetting, [0,0])[4][0]#Ib
Ib4=np.zeros(len(time))
Ib4[0] = Ib#weights(actualsetting, [0,0])[4][1]#Ib

Ib5 = np.zeros(len(time))
Ib5[0] = Ib/1.#1.5
Ib6 = np.zeros(len(time))
Ib6[0] = Ib/1.#1.5
Ib7 = np.zeros(len(time))
Ib7[0] = Ib
Ib8=np.zeros(len(time))
Ib8[0] = Ib
Ibfix=np.zeros(len(time))
Ibfix[0]=Ib
##For supplementary figure
Icomp_L = 0#.1##0.15#0.1
Icomp_R = 0#0.1#0.1
L_comp = 1#0.6#1#0.3#1.3
R_comp = 1#1.3#0.6 
pcompetition_shift = 500
  
sL1=np.zeros(len(time))
sL1[0]=0.08
sR1=np.zeros(len(time))
sR1[0]=0.08
sL2=np.zeros(len(time))
sL2[0]=0.08
sR2=np.zeros(len(time))
sR2[0]=0.08
sLp=np.zeros(len(time))
sLp[0]=0.08
sRp=np.zeros(len(time))
sRp[0]=0.08
sLpp=np.zeros(len(time))
sLpp[0]=0.08
sRpp=np.zeros(len(time))
sRpp[0]=0.08
sfix = np.zeros(len(time))
sfix[0] = 0.1



s1_0=[0.1,0.1]
noise_0=Ib


###############################
####### COLORS#################
###############################



alpha = 0.3 
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
 
inhcmap= [colors_set1['blue'], colors_set1['green'], colors_set1['orange'], colors_set1['red']]
