# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:14:05 2015

@author: jorgejaramillo
"""
from __future__ import division
import numpy as np
import pylab as py
import matplotlib as mp
import random
from twoloopparameters import *
from functionstwoloop import *
from maintwoloop import simulation

font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 7}

mp.rc('font', **font)

newdata = 'yes'
def plotticks(ax):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')

def J_eq(rate,shift,amp):
    return amp*(4*(1-0.5**np.e**((rate/shift)**1.7))-3)

def dmfig():
    #from maintwoloop import simulation
    markersize = 0.7
    fig, axes=allplots(2,1, 2.2,4.6)
    cohlowhigh=[6.4, 51.2]
    colorcoh=[colors_set1['orange'], colors_set1['green'], colors_set1['lightorange'],colors_set1['lightgreen']]
    if newdata=='yes':
        I=np.empty((4,1,len(time),len(cohlowhigh)))
        #S=np.empty((4,len(time),len(cohlowhigh)))

    
    fig.text(-0.008, 0.55, 'Firing rate (spikes/s)', ha='center', va='center', rotation='vertical')
    for cohindex in range(0,len(cohlowhigh)):
        coh1=cohlowhigh[cohindex]
        coh2=-1*coh1
        #I=np.load('Idmfig.npy')
        #I=Idm
        if newdata=='yes':
            I[...,cohindex]=simulation(1,coh1,coh2,1,1)[0]
        else:
            I=I2
        #S[...,cohindex]=simulation(1,coh1,coh2,1,1)[1]

        
        Irecent = I
        #print np.shape(Irecent)
        rL1=FIcurve(Irecent[0,0,...,cohindex],a)
        rR1=FIcurve(Irecent[1,0,...,cohindex],a)
        rL2=FIcurve(Irecent[2,0,...,cohindex],a)
        rR2=FIcurve(Irecent[3,0,...,cohindex],a)
        ax=axes[0]
        
        #sL1 = (S[0,...,cohindex])
        #sR1 = (S[1,...,cohindex])
        #sL2 = (S[2,...,cohindex])
        #sR2 = (S[3,...,cohindex])
        
        #J_LP = 0.3
        #J_RP = -0.3
        
        amp =0.25
        rateshift = 30
        #J_LP_eq = -amp*step(rL1,0)+2*amp*step(rL1,rateshift)
        #J_RP_eq = -amp*step(rR1,0)+2*amp*step(rR1,rateshift)

        #J_LP_eq = J_eq(rL1, rateshift,amp)
        #J_RP_eq = J_eq(rR1, rateshift,amp)

        
        #rP = FIcurve(J_LP_eq*sL1+J_RP_eq*sR1+Ib)
        
        
        
    
        ax.plot(time,rR1, color=colorcoh[cohindex+2], linewidth=markersize)
        ax.plot(time,rL1, color=colorcoh[cohindex],linewidth=markersize)
       
        ax.set_xlim(0,800)
        ax.set_xticks([0,400,800])

        ax.set_ylim(0,60)
        #fig, ax2=allplots()
        ax2=axes[1]
        time_th1 = time[0:340]
        time_th2 = time[0:640]
        if cohindex==1:
            timeth=time_th1
        else:
            timeth=time_th2
        ax2.plot(time,rR2, color=colorcoh[cohindex+2], linewidth=markersize)
        ax2.plot(timeth,rL2[timeth], color=colorcoh[cohindex],linewidth=markersize)

        #ax2.plot(time,rL1, linewidth=markersize)
        #print colorcoh[cohindex]
        ax2.plot(time, 40*np.ones(len(time)), '--', color='black')
        ax2.set_xlabel('Time (ms)')

        #ax2.set_ylabel('Firing rate (spikes/s)')
        ax2.set_xlim(0,800)
        ax2.set_xticks([0,400,800])
        ax2.set_ylim(0,60)
        #ax3 =axes[2]
        #ax3.plot(time,rP, color=colorcoh[cohindex], ms=5)
        #ax3.set_xlim(0,900)
        
        
    ax.set_xticklabels([])
    
    
    for ax in axes:
        ax.text(620, 18, "easy\n(contrast= 51.2\%)", fontsize=7, color= colors_set1['green'])
        ax.text(620, 9, "hard\n(contrast= 6.4\%)", fontsize=7, color=colors_set1['orange'])
    #fig.savefig("DMwithtargetnov2016.pdf")
    return rL1,rL2
dmfig()
#I2=np.load('easyhardapril.npy')