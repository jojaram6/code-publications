# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 18:01:02 2016

@author: jorgejaramillo
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:57:26 2015

@author: User
"""
import numpy as np
from pulvinarparams import *####twoloopparametersLR
#import thalamswitch as ts
from scipy.interpolate import interp1d
from scipy.signal import gaussian
from scipy.ndimage import filters
from scipy.optimize import curve_fit

import pylab as py


def step(t, t_shift):
    return (np.sign(t-t_shift)+1)/2
    
def FIcurve (I,a):
    #if np.abs(a*I-b)<10**-8:
     #   return 1/d
    #else:
        return (a*I-b)/(1-np.e**(-d*(a*I-b)))
def FIcurve_2(I,a_p):
        #a_p=300.0# confidence values
        b_p=112#112 confidence values
        d_p=0.2#0.2 confidence values
        return (a_p*I-b_p)/(1-np.e**(-d_p*(a_p*I-b_p)))#alphaL*I-betaL)
def FIcurve_h (I,h):
    if I<I_h:
        a_h = a+100*h
    else: 
        a_h = a
    return (a_h*I-b)/(1-np.e**(-d*(a_h*I-b)))
        
def J_eq(rate,shift,amp):
    return amp*(4*(1-0.5**np.e**((rate/shift)**1.7))-3)


def Ipulse(t,t_start, t_stop):
    return step(t,t_start)-step(t,t_stop)
def Iramp(t,t_start,slope):
    return slope*(t-t_start)*step(t,t_start)        
def Itarget(t, exp, coh):
    if exp=='wm':
        return 0#currentdistractor/10*(t>=timecue)*(t<timecue+100)
    elif exp=='dm':
         return Irisedecay(t,coh, tau_decay)
         
    else: 
        return 0#Irisedecay(t,coh)#np.e*(gext*visual-Imotion(t,coh))*(t-t_target)/tau_visual*np.e**(-(t-t_target)/tau_visual)*step(t, t_target)*(t>t_target)#*(t<t_motion)#+gext*(6+44*np.e**(-(t-t_motion)/tau_ad))*step(time,t_motion)*(t>=t_motion)
def Imotion(t,coh):
    return I_e*(1+coh/100.0)*(t>t_motion)-0.12*(t>t_motion+motionduration)*(t<t_motion+motionduration+100)
def Imotion2(t):
    return 1#(t>t_motion)*(t<t_motion+800)

def Irisedecay(t,coh,tau_decay):
    t_peak=t_target+tau_decay*tau_rise/(tau_decay-tau_rise)*np.log(tau_decay/tau_rise)
    Irisedecay.tpeak=t_peak
    normfactor=1/(np.e**(-(t_peak-t_target)/tau_decay)-np.e**(-(t_peak-t_target)/tau_rise))
    Irisedecay.normfactor=normfactor

    if exptype=='confidence':
        return (A_risedecay-Imotion(t,coh))*normfactor*(np.e**(-(t-t_target)/tau_decay)-np.e**(-(t-t_target)/tau_rise))*step(t, t_target)*(t>t_target)
    else:
        return (A_risedecay-Imotion(t,coh))*normfactor*(np.e**(-(t-t_target)/tau_decay)-np.e**(-(t-t_target)/tau_rise))*step(t, t_target)*(t>t_target)
def Irisedecay_pul(t,amp,tau_decay):
    t_peak=t_target+tau_decay*tau_rise/(tau_decay-tau_rise)*np.log(tau_decay/tau_rise)
    normfactor=1/(np.e**(-(t_peak-t_target)/tau_decay)-np.e**(-(t_peak-t_target)/tau_rise))
    return amp*normfactor*(np.e**(-(t-t_target)/tau_decay)-np.e**(-(t-t_target)/tau_rise))*step(t, t_target)*(t>t_target)
    

def Icue(t,currentcue):
    return currentcue*(t>=timecue)*(t<timecue+200)
def Idistractor(t,currentdistractor):
    return currentdistractor*(t>=timedistractor)*(t<timedistractor+200)
def Iwhite(t):
    return np.random.normal(0,1)#gaussian
def Inoise1(noise,t,Ib_ini):
    #Ib1=noise[0]
    #Ib2=noise[1]
    Iwhitevector=np.vectorize(Iwhite)
    return -(noise-Ib_ini)/tau0 + Iwhitevector(t)*np.sqrt(tau0)/tau0*sigma/np.sqrt(stepsize)

def Iext1(t, exp, coh,tau_decay):
    if exp=='dm':
        return Imotion(t,coh)*(t>t_motion)#Irisedecay(t,coh,tau_decay)+Imotion(t,coh)*(t>t_motion)
    else: 
        return Icue(t,coh)

def Iext2(t, exp,coh,tau_decay):
    if exp=='dm':
        return Imotion(t,coh)*(t>t_motion)#(Irisedecay(t,coh,tau_decay)+Imotion(t,coh)*(t>t_motion))#*(t<t_motion+500)
    else: 
        return Idistractor(t,coh)

def Iext3(t, exp, coh,tau_decay):
    if exp=='dm':
        if exptype=='confidence':
            return 0
        elif exptype=='conflict':
            return Imotion(t,coh)
    else: 
        return Icue(t,coh)

def Iext4(t, exp,coh,tau_decay):
    if exp=='dm':
        if exptype=='confidence':
            return 0
        elif exptype=='conflict':
            return Imotion(t,coh)
    else: 
        return Idistractor(t,coh)

def Iext5(t, exp, coh,tau_decay):
    
    if exp=='dm':
        return 0#-1*Ipulse(t, 200,250)#Irisedecay_pul(t,coh,tau_decay)#1.3*(Irisedecay(t,coh,tau_decay)+Imotion(t,coh))#+Imotion(t,coh))
    else:
        return 0

def synL1(sL1pre, sR1pre, sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent,localindex):
    allweightsLR, allweightsP = weights(actualsetting, localindex)[0:2]
    [J_L1L1,J_R1L1, J_L2L1, J_R2L1] = [allweightsLR[0],allweightsLR[5], allweightsLR[11], allweightsLR[15]]
    [J_LpL1,J_RpL1] = [allweightsP[2],allweightsP[3]]
    [J_LppL1,J_RppL1] = [J_LpL1,J_RpL1]
#    
    IL1 = J_L1L1*sL1pre+J_R1L1*sR1pre+J_L2L1*sL2pre+J_R2L1*sR2pre+J_LpL1*sLp+J_RpL1*sRp+ J_LppL1*sLpp+J_RppL1*sRpp +  Iext1(t,exp,coh,10)+Itarget(t,exp,coh)+Ib1[np.where(time==t)]
    #IL1 = Iext1(t,exp,coh,10)+Itarget(t,exp,coh)+Ib1[np.where(time==t)]
    return (FIcurve(IL1,a_mod1)*gam/1000*(1-sL1pre)-sL1pre/tau_s+changef*step(t,tchange)), IL1
def synR1(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh, inhpercent,localindex):
    allweightsLR, allweightsP = weights(actualsetting, localindex)[0:2]
    [J_R1R1,J_L1R1, J_R2R1, J_L2R1]= [allweightsLR[0],allweightsLR[5], allweightsLR[11], allweightsLR[15]]
    [J_LpR1,J_RpR1] = [allweightsP[10],allweightsP[11]]
    [J_LppR1,J_RppR1] = [J_LpR1,J_RpR1]
    IR1 = J_R1R1*sR1pre+J_L1R1*sL1pre+J_R2R1*sR2pre+J_L2R1*sL2pre+J_LpR1*sLp+J_RpR1*sRp+J_LppR1*sLpp+J_RppR1*sRpp+Iext2(t,exp,coh,10)+Itarget(t,exp,coh)+Ib2[np.where(time==t)]# use Idistractor(t)/Imotion(t,coh2) for working memory/decision making
    #IR1 = 0##print IR1
    return FIcurve(IR1,a_mod1)*gam/1000*(1-sR1pre)-sR1pre/tau_s+changef*step(t,tchange), IR1
    
def synL2(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent, localindex):
    allweightsLR, allweightsP = weights(actualsetting, localindex)[0:2]
    [J_L2L2,J_R2L2, J_L1L2, J_R1L2]= [allweightsLR[2],allweightsLR[7], allweightsLR[9], allweightsLR[14]]
    [J_LpL2,J_RpL2] = [allweightsP[6],allweightsP[7]]
    [J_LppL2,J_RppL2] = [J_LpL2,J_RpL2]
    
    IL2 = J_L2L2*sL2pre+J_R2L2*sR2pre+J_L1L2*sL1pre+inhpercent*J_R1L2*sR1pre+J_LpL2*sLp+J_RpL2*sRp + J_LppL2*sLpp+J_RppL2*sRpp+ Ib3[np.where(time==t)]+Iext3(t,exp,coh,10)
    
    return FIcurve(IL2,a_mod2)*gam/1000*(1-sL2pre)-sL2pre/tau_s, IL2
    
def synR2(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent,localindex):
    allweightsLR, allweightsP = weights(actualsetting, localindex)[0:2]
    [J_R2R2,J_L2R2, J_R1R2, J_L1R2]= [allweightsLR[2],allweightsLR[7], allweightsLR[9], allweightsLR[14]]
    [J_LpR2,J_RpR2] = [allweightsP[14],allweightsP[15]]
    [J_LppR2,J_RppR2] = [J_LpR2,J_RpR2]
    
    IR2 = J_R2R2*sR2pre+J_L2R2*sL2pre+J_R1R2*sR1pre+inhpercent*J_L1R2*sL1pre+J_LpR2*sLp+J_RpR2*sRp+J_LppR2*sLpp+J_RppR2*sRpp+Ib4[np.where(time==t)]+Iext4(t,exp,coh,10)
    
    return FIcurve(IR2,a_mod2)*gam/1000*(1-sR2pre)-sR2pre/tau_s, IR2

def synLp(S,sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,s_fix, t,coh,inhpercent, localindex):
    allweightsLR, allweightsP,lambda_p = weights(actualsetting, localindex)[0:3]   
    [J_LpLp,J_RpLp,J_L2Lp,J_R2Lp,J_L1Lp,J_R1Lp] = [allweightsP[0], allweightsP[19],allweightsP[8], allweightsP[16], allweightsP[4], allweightsP[12]]
    [J_L2Lpp, J_R2Lpp] = [J_L2Lp, J_R2Lp]
        
    
    ILp = p_number*(rec*J_LpLp*sLp + J_RpLp*sRp + J_L2Lp*sL2pre + J_R2Lp*sR2pre) + J_L1Lp*sL1pre + J_R1Lp*sR1pre+frac_Ib_Lp*Ib5[np.where(time==t)]+Icomp_L*step(t,pcompetition_shift)+Jfix*s_fix + (1-ipsi)*ipsi_amp*Ipulse(t,t_microstart,t_microstart+microduration)+Iext5(t,exp,coh,30)
    ILpp = J_L2Lpp*sL2pre + J_R2Lpp*sR2pre +frac_Ib_Lpp*Ib7[np.where(time==t)]#+Iext3(t,exp,coh)
    if exptype == 'confidence':
        
        ILp = get_pulvinarcurrent(S)+(Iext5(t, 'dm', amp_conf,decay_conf)+offset_conf)*(t<t_motion+off_pulvmotion)#Iext5(t,'dm',800,30)-0.05
        
    if absfunction=='yes':
        return -sLp/tau_p+np.abs(FIcurve(IL1,a)-FIcurve(IR1))/1000,ILp,-sLpp/tau_p+FIcurve(ILpp,a)/1000,ILpp
    else:
        if cx_impl == 'yes':
            return -sLp/tau_s + FIcurve_2(ILp,L_comp*lambda_p)*gam/1000*(1-sLp)+pulse_change*step(t,t_pulsechange),ILp,-sLpp/tau_p+FIcurve(ILpp,a)/1000,ILpp
        else:
            
            return -sLp/tau_p + 2*FIcurve_2(ILp,L_comp*lambda_p)/1000+pulse_change*step(t,t_pulsechange),ILp,-sLpp/tau_p+FIcurve(ILpp,a)/1000,ILpp



def synRp(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,s_fix, t,coh,inhpercent,localindex):
    allweightsLR, allweightsP, lambda_p = weights(actualsetting, localindex)[0:3]        
    [J_RpRp,J_LpRp,J_R2Rp,J_L2Rp,J_R1Rp,J_L1Rp] = [allweightsP[1], allweightsP[18],allweightsP[17], allweightsP[9], allweightsP[13], allweightsP[5]]
    [J_L2Rpp, J_R2Rpp] = [J_L2Rp, J_R2Rp] 
    
    IRp = p_number*(rec*J_RpRp*sRp+J_LpRp*sLp+J_R2Rp*sR2pre+J_L2Rp*sL2pre)+J_R1Rp*sR1pre+J_L1Rp*sL1pre+frac_Ib_Rp*Ib6[np.where(time==t)]+Icomp_R*step(t,pcompetition_shift) + Jfix*s_fix + (ipsi)*micro_amp*Ipulse(t,t_microstart,t_microstart+microduration)+Iext5(t,exp,coh,30)
    IRpp = J_R2Rpp*sR2pre+J_L2Rpp*sL2pre+frac_Ib_Rpp*Ib8[np.where(time==t)]#+Iext4(t,exp,coh)
   
    if cx_impl == 'yes':
            return -sRp/tau_s + FIcurve_2(IRp,R_comp*lambda_p)*gam/1000*(1-sRp)+pulse_change*step(t,t_pulsechange),IRp,-sRpp/tau_p+FIcurve(IRpp,a)/1000,IRpp
    else:
            
            return -sRp/tau_p + FIcurve_2(IRp,R_comp*lambda_p)/1000+pulse_change*step(t,t_pulsechange),IRp,-sRpp/tau_p+FIcurve(IRpp,a)/1000,IRpp

    #return -sRp/tau_p+FIcurve_2(IRp,R_comp*lambda_p)/1000,IRp,-sRpp/tau_p+FIcurve(IRpp,a)/1000,IRpp



def synfix(s_fix,sRp, sLp,t): 
#    
    Ifix = Ibfix[np.where(time==t)]+Ifix_amp*(step(t,0))-Iramp(t,t_target,(Ifix_amp-fix_level)/fix_decay)+ Iramp(t,t_target+fix_decay,(Ifix_amp-fix_level)/fix_decay) + fix_level*step(t,t_target+fix_decay) + micro_amp*Ipulse(t,t_microstart,t_microstart+microduration)
    
    return -s_fix/tau_fix+FIcurve(Ifix)/1000, Ifix
#
#

def get_confrates(sL1pre, sR1pre, sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent,localindex):
    IL1_conf = synL1(sL1pre, sR1pre, sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent,localindex)[1]
    IR1_conf = synR1(sL1pre, sR1pre, sL2pre, sR2pre,sLp, sRp,sLpp,sRpp,t,coh,inhpercent,localindex)[1]
    r1_conf = FIcurve(IL1_conf,a)
    r2_conf = FIcurve(IR1_conf,a)
    #print r1_conf, r2_conf			
    return r1_conf,r2_conf
def get_pulvinarcurrent(S):
    SL,SR = S
    #py.plot(SL+SR)
    Ip = SL+SR+Ibconf
    #print FIcurve_2(ILp)
    return Ip




#def synLpp(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,t,coh,inhpercent, localindex):
#    allweightsLR, allweightsP = weights(actualsetting, localindex)   
#    [J_LpLp,J_RpLp,J_L2Lp,J_R2Lp,J_L1Lp,J_R1Lp] = [allweightsP[0], allweightsP[19],allweightsP[8], allweightsP[16], allweightsP[4], allweightsP[12]]
#    ILp = rec*J_LpLp*sLp + J_RpLp*sRp + J_L2Lp*sL2pre + J_R2Lp*sR2pre + J_L1Lp*sL1pre + J_R1Lp*sR1pre+frac_Ib_L*Ib5[np.where(time==t)]#+Iext3(t,exp,coh)
#    return -sLp/tau_p+FIcurve(ILp)/1000,ILp #FIcurve(ILp)*gam/1000*(1-sLp)-sLp/tau_p, ILp
#    
#def synRpp(sL1pre,sR1pre,sL2pre, sR2pre,sLp, sRp,t,coh,inhpercent,localindex):
#    allweightsLR, allweightsP = weights(actualsetting, localindex)        
#    [J_RpRp,J_LpRp,J_R2Rp,J_L2Rp,J_R1Rp,J_L1Rp] = [allweightsP[1], allweightsP[18],allweightsP[17], allweightsP[9], allweightsP[13], allweightsP[5]]
#    IRp = rec*J_RpRp*sRp+J_LpRp*sLp+J_R2Rp*sR2pre+J_L2Rp*sL2pre+J_R1Rp*sR1pre+J_L1Rp*sL1pre+frac_Ib_R*Ib6[np.where(time==t)]#+Iext4(t,exp,coh)
#    return -sRp/tau_p+FIcurve(IRp)/1000,IRp#FIcurve(IRp)*gam/1000*(1-sRp)-sRp/tau_p, IRp
#



#J_LR =  [J_L1L1,J_R1R1,J_R2R2,J_R2R1,J_R1R2,J_L1R1, J_L2L2, J_R2L2, J_L2R2,J_L1L2, J_L1R2, J_L2L1, J_L2R1, J_R1L1,J_R1L2, J_R2L1]   
   
#J_p = [J_LpLp, J_RpRp,J_LpL1,J_RpL1,J_L1Lp,J_L1Rp, J_LpL2, J_RpL2, J_L2Lp, J_L2Rp,J_LpR1,J_RpR1,J_R1Lp, J_R1Rp,J_LpR2,J_RpR2,J_R2Lp,J_R2Rp,J_LpRp,J_RpLp]
  

def allplots(nrows,ncolumns,width,height):
    fig, axes = py.subplots(nrows,ncolumns,figsize=(width,height)) #(6,5) for most plots   # Hide the right and top spines
    for ax in axes:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')
    return fig,axes

def plotticks(ax):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')
def plotticks2(ax):
        #ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('right')
        ax.xaxis.set_ticks_position('bottom')
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')


def winnercounter(timevector1, timevector2):
    if np.size(timevector2)==0 and np.size(timevector1)!=0:
       winnercount=2
       winner=timevector1
    elif np.size(timevector2)!=0 and np.size(timevector1)!=0:
         if np.min(timevector1)<np.min(timevector2):
            winnercount=2
            winner=timevector1
         else:
             winnercount=1
             winner=timevector2
    elif np.size(timevector2)!=0 and np.size(timevector1)==0:
         winnercount=1
         winner=timevector2
    elif np.size(timevector2)==0 and np.size(timevector1)==0:
         #print "george"
         winnercount=0
         winner=0
    return [winnercount,winner]
def quintile(index, timevector,stepr):
    return np.mean(timevector[(timevector>np.min(timevector)+index*stepr)*(timevector>np.min(timevector)+(index+1)*stepr)]) 
    #else:
     #    continue
def filter_t(timevector,sigma):
     
    b = gaussian(100, sigma)
    filtered_timevector = filters.convolve1d(timevector, b/b.sum())
    return filtered_timevector         
    
