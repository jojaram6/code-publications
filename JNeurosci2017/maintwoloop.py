# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:56:10 2015

Jorge J
"""
from __future__ import division
import numpy as np
import pylab as py
import random
from twoloopparameters import *

from functionstwoloop import *


def simulation(totaltrials,coh1,coh2,inhpercent, localindex):  
     
    IL1=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IR1=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IL2=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IR2=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]

    for trialnumber in np.arange(0, totaltrials):
        #print trialnumber
    
        #for index in np.arange(0, np.ceil(maxtime/stepsize)-1):
            
            
        
    
    
        for indexs in np.arange(0, np.ceil(maxtime/stepsize)-1):
            indexs=int(indexs)
            
            sL1[indexs+1]=sL1[indexs]+stepsize*synL1(sL1[indexs], sR1[indexs],sL2[indexs],sR2[indexs], indexs*stepsize, coh1, inhpercent, localindex)[0] 
            sR1[indexs+1]=sR1[indexs]+stepsize*synR1(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], indexs*stepsize, coh2,inhpercent, localindex)[0]
            sL2[indexs+1]=sL2[indexs]+stepsize*synL2(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], indexs*stepsize,coh1,inhpercent,localindex)[0]
            sR2[indexs+1]=sR2[indexs]+stepsize*synR2(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], indexs*stepsize,coh2, inhpercent,localindex)[0]
            
            Ib1[indexs+1] = Ib1[indexs]+stepsize*Inoise1(Ib1[indexs],indexs*stepsize)  
            Ib2[indexs+1] = Ib2[indexs]+stepsize*Inoise1(Ib2[indexs],indexs*stepsize)
            Ib3[indexs+1] = Ib3[indexs]+stepsize*Inoise1(Ib3[indexs],indexs*stepsize)
            Ib4[indexs+1] = Ib4[indexs]+stepsize*Inoise1(Ib4[indexs],indexs*stepsize)
            
    
        
        IL1[trialnumber]=synL1(sL1, sR1, sL2, sR2,time, coh1, inhpercent, localindex)[1] # J_L1L1*sL1+J_R1L1*sR1+J_L2L1*sL2+J_R2L1*sR2+Imotion(time,coh1)+Itarget(time)+Ib1
        IR1[trialnumber]=synR1(sL1, sR1, sL2, sR2,time, coh2,inhpercent, localindex)[1]#J_R1R1*sR1+J_L1R1*sL1+J_R2R1*sR2+J_L2R1*sL2+Imotion(time,coh2)+Itarget(time)+Ib2
        IL2[trialnumber]=synL2(sL1, sR1, sL2, sR2,time,coh1,inhpercent,localindex)[1]#J_L2L2*sL2+J_R2L2*sR2+J_L1L2*sL1+J_R1L2*sR1+Ib3
        IR2[trialnumber]=synR2(sL1, sR1, sL2, sR2,time,coh2,inhpercent,localindex)[1]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
    #print sL1[0]
    #print sL1[1]
#    print coh1
#    print coh2
    I=np.array([IL1, IR1, IL2, IR2])
    S=np.array([sL1, sR1, sL2, sR2])
    simulation.testnoise1=Ib1
    simulation.testnoise2=Ib2
    return I#,S

#  
def allcohtrials(inhpercent):
    start_time = clocktime.time()
    resultscoh=np.empty((4,totaltrials,len(time), totalcoh))
    #print resultscoh
    for cohindex, coh1 in enumerate(allcoh1):
        #global coh1
        #global coh2
        coh2=-1*coh1
        print coh1
        resultssim=np.asarray(simulation(totaltrials,coh1,coh2,inhpercent,0))
        resultscoh[:,:,:,cohindex]=resultssim#stores all simulation trials
    print start_time-clocktime.time()
    return resultscoh
#resultssimcoh=allcohtrials()
#np.save('allcohtarget', resultssimcoh )
def inhtrials():
    allinhtrials = np.empty((4,totaltrials,len(time), totalcoh, len(allinh)))
    for inhindex, inhpercent in enumerate (allinh):
        #J_L1R2=(J_p12-J_m12)/2*inhpercent
        #J_R1L2=J_L1R2
        #J_L2R1=(J_p21-J_m21)/2*inhpercent
        #J_R2L1=J_L2R1
        print inhindex
        allinhtrials[...,inhindex]=allcohtrials(inhpercent)
#        print J_L1R2
#        print J_L1R2+J_L1L2
    return allinhtrials

#allcohtarget=inhtrials()
#testtrials = allcohtrials(1)
