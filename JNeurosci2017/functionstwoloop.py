# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:57:26 2015

@author: User
"""

#from twoloopparameters import *
import numpy as np
import pylab as py
from scipy.signal import gaussian
from scipy.ndimage import filters
from scipy.optimize import curve_fit


def step(t, t_shift):
    return (np.sign(t-t_shift)+1)/2
def pulse(t, t_start, t_duration):
    return step(t,t_start)-step(t,t_start+t_duration)
    
def FIcurve (I,a_param):
#    if np.abs(a*I-b)<10**-8:
#        return 1/d
#    else:
        
        return (a_param*I-b)/(1-np.e**(-d*(a_param*I-b)))

def Itarget(t, exp, coh):
    if exp=='dm':
        return Irisedecay(t,coh)#currentdistractor/10*(t>=timecue)*(t<timecue+100)
    else:
        return 0#0.04*(t>t_target)*(0.6+1.0*np.exp(-(t-t_target)/10.0))#np.e*(gext*visual-Imotion(t,coh))*(t-t_target)/tau_visual*np.e**(-(t-t_target)/tau_visual)*step(t, t_target)*(t>t_target)#*(t<t_motion)#+gext*(6+44*np.e**(-(t-t_motion)/tau_ad))*step(time,t_motion)*(t>=t_motion)
def Imotion(t,coh):
    return I_e*(1+coh/100.0)*(t>t_motion)#*(t<t_motion+500)
def Irisedecay(t,coh):
    t_peak=t_target+tau_decay*tau_rise/(tau_decay-tau_rise)*np.log(tau_decay/tau_rise)
    Irisedecay.tpeak=t_peak
    normfactor=1/(np.e**(-(t_peak-t_target)/tau_decay)-np.e**(-(t_peak-t_target)/tau_rise))
    Irisedecay.normfactor=normfactor

    return (A_risedecay-Imotion(t,coh))*normfactor*(np.e**(-(t-t_target)/tau_decay)-np.e**(-(t-t_target)/tau_rise))*step(t, t_target)*(t>t_target)


def I_luminance(t,signal):
    return signal*(np.e**-((t-t_luminance)/tau_luminance))*step(t,t_luminance)#
#signal*(step(t,t_luminance)-step(t,t_luminance+luminanceduration))#
def Ireward(t,ramplitude):
    return ramplitude*(1-np.e**-((t-t_reward)/tau_reward))*step(t,t_reward)



######### for hanks simulations#########
def gen_clicks(rate):
    t_poisson=np.empty(poissontotal+1)
    t_poisson[0]=np.random.uniform(0,2)
    for index in range(0,poissontotal):
        t_poisson[index+1]=t_poisson[index]-np.log(np.random.uniform(0,1))/rate
    return t_poisson[1:len(t_poisson)]

def clicks_L(t, Lvec):
    clicks_L.Lvec=Lvec

    if t in np.round(Lvec):
       return np.array([amp_clicks/stepsize,1])
    else:
       return np.array([0,0])
def clicks_R(t,Rvec):
    clicks_R.Rvec=Rvec
    if t in np.round(Rvec):
       return np.array([amp_clicks/stepsize,1])
    else:
       return np.array([0,0])

def I_clicks_L(t, Lvec):
    current_old=0
    current_new=0
    for a in np.round(Lvec):
        #print a
        current_new = amp_pulse*pulse(t,a,pulseduration)+current_old
        current_old=current_new
    return current_new
def I_clicks_R(t, Rvec):
    current_old=0
    current_new=0
    for a in np.round(Rvec):
        current_new = amp_pulse*pulse(t,a,pulseduration)+current_old
        current_old=current_new
    return current_new
    
     

def Icue(t,currentcue):
    return currentcue*(t>=timecue)*(t<timecue+cueduration)
def Idistractor(t,currentdistractor,tdistractor):
    return currentdistractor*(t>=tdistractor)*(t<tdistractor+distractorduration)
def Iwhite(t):
    return np.random.normal(0,1)#gaussian
def Inoise1(noise,t):
    #Ib1=noise[0]
    #Ib2=noise[1]
    Iwhitevector=np.vectorize(Iwhite)
    return -(noise-Ib)/tau0 + Iwhitevector(t)*np.sqrt(tau0)/tau0*sigma/np.sqrt(stepsize)

def Iext1(t, exp, coh):
    if exp=='dm':
        return Imotion(t,coh)
    elif exp=='wm' or exp=='locallongrange': 
        return Icue(t,coh)
    elif exp=='hanks':
        return I_clicks_L(t+motion_offset,coh[0])-c_inh*I_clicks_R(t+motion_offset,coh[1])#+0.6*I_e
    elif exp=='luminance':
        return I_luminance(t,coh)
    else:
        return 0


def Iext2(t, exp,coh,tdistractor):
    if exp=='dm':
        return Imotion(t,coh)
    elif exp=='wm': 
        if actualsetting == 'rev_fig':
            return Idistractor(t,coh,tdistractor_revfig)
        else:
            return Idistractor(t,coh,tdistractor)
    elif exp=='hanks':
        return I_clicks_R(t+motion_offset,coh[1])-c_inh*I_clicks_L(t+motion_offset,coh[0])#+0.6*I_e
    elif exp=='luminance':
        return I_luminance(t,coh)
    elif exp=='locallongrange':
        return 0#I_luminance(t,coh)
    else:
        return 0
    
def Ihuk(t,exp,coh):
    if exp=='huk':
        return Imotion(t,coh)
    else:
        return 0 


def synL1(sL1pre, sR1pre, sL2pre, sR2pre,t,coh,inhpercent,localindex):
    #1[np.where(time==t)]
    allweights=weights(actualsetting, localindex)
    [J_L1L1,J_R1L1, J_L2L1, J_R2L1]= [allweights[0],allweights[1], allweights[6], allweights[7]]
    IL1=J_L1L1*sL1pre+J_R1L1*sR1pre+J_L2L1*sL2pre+inhpercent*J_R2L1*sR2pre+Iext1(t,exp,coh)+Itarget(t,exp,coh)+Ib1[np.where(time==t)]
    # use Icue(t)/Imotion(t,coh1) for working memory/decision making
    return (FIcurve(IL1,ahuk)*gam/1000*(1-sL1pre)-sL1pre/tau_s), IL1

def synR1(sL1pre,sR1pre,sL2pre, sR2pre,t,coh, inhpercent,localindex):
    allweights=weights(actualsetting, localindex)
    [J_R1R1,J_L1R1, J_R2R1, J_L2R1]= [allweights[0],allweights[1], allweights[6], allweights[7]]
    IR1=J_R1R1*sR1pre+J_L1R1*sL1pre+J_R2R1*sR2pre+inhpercent*J_L2R1*sL2pre+Iext2(t,exp,coh,localindex)+Itarget(t,exp,coh)+Ib2[np.where(time==t)]# use Idistractor(t)/Imotion(t,coh2) for working memory/decision making
    return FIcurve(IR1,a)*gam/1000*(1-sR1pre)-sR1pre/tau_s, IR1
    
def synL2(sL1pre,sR1pre,sL2pre, sR2pre,t,coh,inhpercent, localindex):
    allweights=weights(actualsetting, localindex)
    [J_L2L2,J_R2L2, J_L1L2, J_R1L2]= [allweights[2],allweights[3], allweights[4], allweights[5]]
    IL2=J_L2L2*sL2pre+J_R2L2*sR2pre+J_L1L2*sL1pre+inhpercent*J_R1L2*sR1pre+Ib3[np.where(time==t)]+ Itarget(t,exp,coh) + Ihuk(t,exp,coh)# + Iext1(t,exp,coh)
    #return FIcurve(IR2,a)*gam/1000*(1-sR2pre)-sR2pre/tau_s, IR2
    #print J_R1L2
    #print J_R1L2*inhpercent
    #print J_R1L2+J_L1L2
    #print J_L1L2
    return FIcurve(IL2,a)*gam/1000*(1-sL2pre)-sL2pre/tau_s, IL2
    
def synR2(sL1pre,sR1pre,sL2pre, sR2pre,t,coh,inhpercent,localindex):
    allweights=weights(actualsetting, localindex)
    [J_R2R2,J_L2R2, J_R1R2, J_L1R2]= [allweights[2],allweights[3], allweights[4], allweights[5]]
    IR2=J_R2R2*sR2pre+J_L2R2*sL2pre+J_R1R2*sR1pre+inhpercent*J_L1R2*sL1pre+Ib4[np.where(time==t)]+ Itarget(t,exp,coh) + Ihuk(t,exp,coh)# + Iext2(t,exp,coh,localindex)
    return FIcurve(IR2,a)*gam/1000*(1-sR2pre)-sR2pre/tau_s, IR2

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


def winnercounter(timevector1, timevector2):
    #print np.size(timevector1)
    #print np.size(timevector2)
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
    return winnercount,winner
def quintile(index, timevector,stepr):
    #return np.mean(timevector[(timevector>np.min(timevector)+index*stepr)*(timevector<np.min(timevector)+(index+1)*stepr)]) 
    timevectorsort=np.sort(timevector)
    return np.mean(timevectorsort[np.floor(stepr*index):np.floor((index+1)*stepr)-1])
         
def errormemory(vector1, vector2,timedecision):
    if vector1[timedecision]>vector2[timedecision]:
       errmem=0
    else:
        errmem=1
    return errmem 


def autocorr(timevector):
    timevector=timevector/np.max(timevector)
    #timevector = timevector[1000:34000]
    #[result.size/2:]
    result = np.correlate(timevector-np.mean(timevector), timevector-np.mean(timevector), mode='full')
    return result[result.size/2:]/np.max(result)

def autocorr_2(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result






def filter_autocorr(timevector,sigma):
     
    b = gaussian(100, sigma)
    filtered_timevector = filters.convolve1d(timevector, b/b.sum())
    return filtered_timevector
    
def linearfit(vec, vec_x):
    def func(x, a, b):
        return a*x+b
    
    a1=vec
    
    
# 
    Jfit= np.arange(0.35,0.42,0.01)

    
    
    poptR1, pcovR1 = curve_fit(func, vec_x, a1)
    #poptR2, pcovR2 = curve_fit(func, , a2)
    
    
    return Jfit,func(Jfit, *poptR1)   