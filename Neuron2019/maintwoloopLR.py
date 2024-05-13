# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:56:10 2015

Jorge J
"""
from __future__ import division
import numpy as np
import pylab as py
import random
from pulvinarparams import *     ###twoloopparametersLR
import thalamswitch as ts
from functionstwoloopLR_pnumber import *



def simulation(totaltrials,coh1,coh2,coh3,coh4,inhpercent, localindex):  
     
    #print coh1
    #print coh3
    IL1=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IR1=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IL2=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    IR2=np.empty((totaltrials, len(time))) #[0 for tindex in range(totaltrials)]
    ILp=np.empty((totaltrials, len(time))) 
    IRp=np.empty((totaltrials, len(time))) 
    ILpp=np.empty((totaltrials, len(time))) 
    IRpp=np.empty((totaltrials, len(time))) 
    Ifix=np.empty((totaltrials, len(time))) 
    
    for trialnumber in np.arange(0, totaltrials):
       
        for indexs in np.arange(0, np.ceil(maxtime/stepsize)-1):
            indexs=int(indexs)
        
            
            sL1[indexs+1]= sL1[indexs]+stepsize*synL1(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs], indexs*stepsize, coh1, inhpercent, localindex)[0] 
            sR1[indexs+1]= sR1[indexs]+stepsize*synR1(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs], indexs*stepsize, coh2,inhpercent, localindex)[0]
            sL2[indexs+1]=sL2[indexs]+stepsize*synL2(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],indexs*stepsize,coh3,inhpercent,localindex)[0]
            sR2[indexs+1]=sR2[indexs]+stepsize*synR2(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],indexs*stepsize,coh4,inhpercent,localindex)[0]
            #sfix[indexs+1]=sfix[indexs]+stepsize*synfix(sfix[indexs], sRp[indexs], sLp[indexs],indexs*stepsize)[0]

                        
            if exptype=='confidence':
               
               bothrate_ext = get_confrates(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],indexs*stepsize,coh1,inhpercent,localindex)
               
               for index_rate,rate_ext in enumerate(bothrate_ext):
                   
                   S_exc[index_rate,indexs + 1]= S_exc[index_rate, indexs] + stepsize_switch*ts.s_eq('exc',(rate_ext-rate_offset)*(rate_ext>rate_offset),indexs*stepsize,S_exc[index_rate, indexs],F[index_rate,indexs],D[index_rate,indexs])
                   S_inh[index_rate,indexs + 1]= S_inh[index_rate, indexs] + stepsize_switch*ts.s_eq('inh',rate_ext,indexs*stepsize,S_inh[index_rate, indexs],F[index_rate,indexs],D[index_rate,indexs])
                   F[index_rate,indexs+1]= F[index_rate,indexs]+ stepsize_switch*ts.fac_arg(F[index_rate,indexs],indexs*stepsize,(rate_ext-rate_offset)*(rate_ext>rate_offset)/1000.0)
                   D[index_rate,indexs + 1]= D[index_rate,indexs] + stepsize_switch*ts.dep_arg(D[index_rate,indexs],indexs*stepsize,rate_ext/1000.0)
                   S_rates[index_rate,indexs] = J_exc*S_exc[index_rate, indexs]+J_inh*S_inh[index_rate, indexs]
                   
            
            
                   sRp[indexs+1]=0
            else:
               #sRp[indexs+1]=0
                   sRp[indexs+1] = sRp[indexs]+stepsize*synRp(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs], sfix[indexs], indexs*stepsize,coh1,inhpercent,localindex)[0]

            
            
            if les==1:
                sLp[indexs+1]=sLp[indexs]+stepsize*synLp(S_rates[:,indexs],sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],sfix[indexs],indexs*stepsize,coh1,inhpercent,localindex)[0]
            else:
                sLp[indexs+1]=0
            
            
            
            if pulvinars==1:
                sLpp[indexs+1]=sLpp[indexs]+stepsize*synLp(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],indexs*stepsize,coh4,inhpercent,localindex)[2]            
                sRpp[indexs+1]=sRpp[indexs]+stepsize*synRp(sL1[indexs], sR1[indexs], sL2[indexs],sR2[indexs], sLp[indexs], sRp[indexs], sLpp[indexs], sRpp[indexs],indexs*stepsize,coh4,inhpercent,localindex)[2]
            else:
                sLpp[indexs+1]=0
                sRpp[indexs+1]=0
            
            
            Ib1[indexs+1]=Ib1[indexs]+stepsize*Inoise1(Ib1[indexs],indexs*stepsize,Ib1[0])  
            Ib2[indexs+1]=Ib2[indexs]+stepsize*Inoise1(Ib2[indexs],indexs*stepsize,Ib2[0])
            Ib3[indexs+1]=Ib3[indexs]+stepsize*Inoise1(Ib3[indexs],indexs*stepsize,Ib3[0])
            Ib4[indexs+1]=Ib4[indexs]+stepsize*Inoise1(Ib4[indexs],indexs*stepsize,Ib4[0])
            Ib5[indexs+1]=Ib5[indexs]+stepsize*Inoise1(Ib5[indexs],indexs*stepsize,Ib5[0])
            Ib6[indexs+1]=Ib6[indexs]+stepsize*Inoise1(Ib6[indexs],indexs*stepsize,Ib6[0])
            Ib7[indexs+1]=Ib7[indexs]+stepsize*Inoise1(Ib7[indexs],indexs*stepsize,Ib7[0])
            Ib8[indexs+1]=Ib8[indexs]+stepsize*Inoise1(Ib8[indexs],indexs*stepsize,Ib8[0])
            #Ibfix[indexs+1]=Ibfix[indexs]+stepsize*Inoise1(Ibfix[indexs],indexs*stepsize,Ibfix[0])
 
            
            
        IL1[trialnumber]=synL1(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, time, coh1, inhpercent, localindex)[1] # J_L1L1*sL1+J_R1L1*sR1+J_L2L1*sL2+J_R2L1*sR2+Imotion(time,coh1)+Itarget(time)+Ib1
        IR1[trialnumber]=synR1(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, time, coh2,inhpercent, localindex)[1]#J_R1R1*sR1+J_L1R1*sL1+J_R2R1*sR2+J_L2R1*sL2+Imotion(time,coh2)+Itarget(time)+Ib2
        IL2[trialnumber]=synL2(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, time, coh3,inhpercent,localindex)[1]#J_L2L2*sL2+J_R2L2*sR2+J_L1L2*sL1+J_R1L2*sR1+Ib3
        IR2[trialnumber]=synR2(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, time, coh4,inhpercent,localindex)[1]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
        ILp[trialnumber]=synLp(S_rates,sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, sfix,time, coh1,inhpercent,localindex)[1]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
        #print ILp
        IRp[trialnumber]=synRp(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, sfix,time, coh1,inhpercent,localindex)[1]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
        #Ifix[trialnumber]=synfix(sfix, sRp, sLp,time)[1]             
        if pulvinars==1:
            ILpp[trialnumber]=synLp(S_rates, sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, sfix, time, coh4,inhpercent,localindex)[3]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
            IRpp[trialnumber]=synRp(sL1, sR1, sL2, sR2,sLp,sRp, sLpp, sRpp, time, coh4,inhpercent,localindex)[3]#J_R2R2*sR2+J_L2R2*sL2+J_R1R2*sR1+J_L1R2*sL1+Ib4
        else:
            ILpp[trialnumber]=0
            IRpp[trialnumber]=0

    
    #print sL1[0]
    #print sL1[1]
#    print coh1
#    print coh2
    I=np.array([IL1, IR1, IL2, IR2, ILp, IRp, ILpp, IRpp, Ifix])
    #print Ifix
    S=np.array([sL1, sR1, sL2, sR2, sLp, sRp, sfix])
    simulation.testnoise1=Ib1
    simulation.testnoise2=Ib2
    return I,S

#  
def allcohtrials(inhpercent):
    resultscoh=np.empty((4,totaltrials,len(time), totalcoh))
    #print resultscoh
    for cohindex, coh1 in enumerate(allcoh1):
        #global coh1
        #global coh2
        coh2=-1*coh1
        print (coh1)
        resultssim=np.asarray(simulation(totaltrials,coh1,coh2,inhpercent,0))
        resultscoh[:,:,:,cohindex]=resultssim#stores all simulation trials
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
        print (inhindex)
        allinhtrials[...,inhindex]=allcohtrials(inhpercent)
#        print J_L1R2
#        print J_L1R2+J_L1L2
    return allinhtrials

#allcohtarget=inhtrials()
#
