# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:11:42 2015

@author: jorgejaramillo
"""
from __future__ import division
import numpy as np
import pylab as py
import matplotlib as mp
import random
import pulvinarparamsWM2 as pp
#from twoloopparameters import *
import functionstwoloop as ftwo
from scipy.optimize import curve_fit

#def wmdistractorfig():
from maintwoloop import simulation
font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 7}

mp.rc('font', **font)
mp.rc('legend',**{'fontsize':7})
mp.rc('text', usetex=True)
def memoryplusdistractors(trial):
    figwidth = 1.6
    figheight = 1.2
    ylimit=80
    maxtime=3500
    if trial=='errortrial':
        tdistractor=timecue+cueduration+500
        currentcue = 0.04#0.04#0.04#use 0.04 for errortrial
        currentdistractor = 0.067#use 0.067 for errortrial
    else:
        tdistractor = 1200#1000#1500
        currentcue = 0.09#use 0.09 for correct trial
        currentdistractor = currentcue

    
    if pfcdisconnect=='yes':
       fig_dis=py.figure(frameon=False, figsize=(figwidth,figheight))
       ax1 = fig_dis.add_subplot(111)
       ax1.set_xlabel('Time (ms)')
       ax1.set_xlim(0,maxtime)
       ax1.set_xticks([200,1200,2200])
       ax1.set_xticklabels([0,1000,2000])

       
       plotticks(ax1)
       
    
    else:
       fig_top = py.figure(frameon=False, figsize=(figwidth,figheight))#,axes_top=allplots(2,1,1.6,2.6)
       ax1 = fig_top.add_subplot(111)
       
       fig_bottom = py.figure(frameon=False, figsize=(figwidth,figheight))
       ax2 = fig_bottom.add_subplot(111)
    
   

    
    
    
    #ax1=axes[0]
    I=simulation(1,currentcue,currentdistractor,1,tdistractor)
    rL1=FIcurve(I[0,0,...],a)
    rR1=FIcurve(I[1,0,...],a)
    rL2=FIcurve(I[2,0,...],a)
    rR2=FIcurve(I[3,0,...],a)
    r_all = np.array([rL1, rR1,rL2, rR2])               
    ax1.plot(time,rL1,color=colors_set1['blue'])
    ax1.plot(time,rR1, color=colors_set1['red'])
    ax1.set_ylim(0,ylimit)
    #ax1.set_yticks([0,40,80])

    if pfcdisconnect=='no':
        #ax2=axes[1]
        ax2.plot(time,rL2,color=colors_set1['cyan'])
        ax2.plot(time,rR2, color=colors_set1['purple'])
        ax2.set_xlabel('Time (ms)')
        #ax2.set_ylabel('Firing rate (spikes/s)')
        ax2.set_xlim(0,maxtime)
        ax2.set_xticks([200,1200,2200])
        ax2.set_ylim(0,ylimit)
        #ax2.set_yticks([0,40,80])
        
    if pfcdisconnect=='yes':
       fig_dis.savefig("WMdisconnectsep.pdf")
    else:
       ax1.set_xlim(0,maxtime)
       ax2.set_xlim(0,maxtime)
       ax1.set_xticks([200,1200,2200])
       ax1.set_xticklabels([])
       ax2.set_xticklabels([0,1000,2000])
       plotticks(ax1)
       plotticks(ax2)
#       if trial == 'errortrial':
#           fig_top.savefig("WMdistractortop_errorsep.pdf")
#           fig_bottom.savefig("WMdistractorbottom_errorsep.pdf")  
#       else:   
#           fig_top.savefig("WMdistractortop_correctsep.pdf")
#           fig_bottom.savefig("WMdistractorbottom_correctsep.pdf")  

    return rL1, rR1, rL2, rR2
        
   
    
#    
    
def distraction_time_data():
    firsttime = timecue + cueduration
    tdistractortimes = np.array([firsttime,firsttime+50, firsttime+100,firsttime+200,2500])
    trialsmemory = 1#600 for errorrate figure#100#3
    dis_sd = 10**-8# 0.04 #for error trials 0.04 #10**-8#0.04#0.02
    cue_sd = dis_sd
    cue_mean = 0.09#0.08#0.06
    dis_mean = cue_mean
    errorall = np.zeros((trialsmemory))
    error = np.zeros(len(tdistractortimes))
    rL1_alldistractor = np.zeros((len(tdistractortimes),len(time)))
    rL2_alldistractor = np.zeros((len(tdistractortimes),len(time)))
    rR1_alldistractor = np.zeros((len(tdistractortimes),len(time)))
    rR2_alldistractor = np.zeros((len(tdistractortimes),len(time)))

    r_all_old = 0 
    r_error_old = 0
    r_correct_old = 0
    count_errortrials = 0
    count_correcttrials = 0
    #r_correct_new = []
    
    for tdindex, tdistractor in enumerate(tdistractortimes):
        r_all_old = 0 
        r_error_old = 0
        r_correct_old = 0
        count_errortrials = 0
        count_correcttrials = 0
        for trials in np.arange(0,trialsmemory):
        
            print trials
            currentcue = np.random.normal(cue_mean, cue_sd)
            currentdistractor = np.random.normal(dis_mean,dis_sd)#currentcue#0.1#currentcue+0.03
            
            I=simulation(1,currentcue,currentdistractor,1,tdistractor)
            rL1=FIcurve(I[0,0,...],a)
            rR1=FIcurve(I[1,0,...],a)
            rL2=FIcurve(I[2,0,...],a)
            rR2=FIcurve(I[3,0,...],a)
            r_all = np.array([rL1, rR1,rL2, rR2])               
            
            if pfcdisconnect=='no':
               errorall[trials] = errormemory(rL2, rR2,2600)
            elif pfcdisconnect=='yes':
               errorall[trials] = errormemory(rL1, rR1,2600) 
            if errorall[trials] ==1:#errormemory(rL2,rR2,2600)==1:
               r_error = np.array([rL1, rR1,rL2, rR2])
               r_error_new = r_error + r_error_old
               r_error_old = r_error_new
               
               count_errortrials = count_errortrials + 1
            if errorall[trials] ==0:#errormemory(rL2,rR2,2600)==0:
               r_correct = np.array([rL1, rR1,rL2, rR2])
               r_correct_new = r_correct + r_correct_old
               r_correct_old = r_correct_new
               count_correcttrials = count_correcttrials + 1
            
            r_all_new = r_all + r_all_old
            r_all_old = r_all_new
                
        if r_correct_new !=[]:
            norm_target = 1#np.max(r_correct_new[0])
            #print norm_target
            rL1_alldistractor[tdindex, :] = r_correct_new[0]/norm_target#(count_correcttrials*norm_target)
            rL2_alldistractor[tdindex, :] = r_correct_new[2]/count_correcttrials
            rR1_alldistractor[tdindex, :] = r_correct_new[1]/norm_target#(count_correcttrials*norm_target)
            rR2_alldistractor[tdindex, :] = r_correct_new[3]/count_correcttrials


        error[tdindex] = np.mean(errorall, axis=0)
    return rL1_alldistractor, rR1_alldistractor, rL2_alldistractor, rR2_alldistractor, error #r_correct_new[0]

def distraction_time_fig(run):
#    if run=='new':
#        alldistractor = distraction_time_data()
    #alldistractor = distraction_time_data()
    
    firsttime = timecue+cueduration    
    errorfig=0
    figwidth = 1.6
    figheight = 1.2   
    xlim = 1000
    
    
    colors = [colors_set1['red'], colors_set1['orange'],colors_set1['lightorange'], colors_set1['pink'],colors_set1['red'],colors_set1['red']]
    colors_bottom = colors#[colors_set1['purple'], 'lightgreen', 'yellow','gray', 'gray']
    tdistractortimes = np.array([firsttime,firsttime+50, firsttime+100,firsttime+200])#),2500])



    if errorfig==1:
        fig2=py.figure(frameon=False, figsize=(1.6,1.2))
        ax= fig2.add_subplot(111) 
    
        r_distractor = alldistractor_nopfc[0:4]
        error = alldistractor[4]
        error_nopfc = alldistractor_nopfc[4]
    else:
        r_distractor = alldistractor[0:4]
        fig_corr_top = py.figure(frameon=False, figsize=(figwidth,figheight))
        fig_corr_bottom = py.figure(frameon=False, figsize=(figwidth,figheight))
        axes_corr_top = fig_corr_top.add_subplot(111)
        axes_corr_bottom = fig_corr_bottom.add_subplot(111)
    
        

    
    
    if errorfig==0:
        left, bottom, width, height = [0.62, 0.34, 0.4, 0.3]#[0.75, 0.65, 0.15, 0.2]
        ax_inset = fig_corr_bottom.add_axes([left, bottom, width, height])
        plotticks(ax_inset)
        ax_inset.get_xaxis().set_visible(False)
        ax_inset.get_yaxis().set_visible(False)
        ax_inset.axis('off')
        for tdindex, tdistractor in enumerate(tdistractortimes):
        
            
            axes_corr_top.plot(time,r_distractor[1][tdindex,:], color = colors[tdindex])
            axes_corr_bottom.plot(time,r_distractor[3][3-tdindex,:], color = colors_bottom[3-tdindex])  
            ####insetplot####
            ax_inset.plot(time,r_distractor[3][3-tdindex,:], color = colors_bottom[3-tdindex]) 
            ax_inset.set_ylim(0,0.6)
            ax_inset.set_xlim(250,900)
            ax_inset.set_yticks([0,0.6])
            ax_inset.set_xticks([250,1000])

            


            if tdindex == 0:
                axes_corr_top.plot(time,r_distractor[0][4-tdindex,:],color=colors_set1['blue'])
                axes_corr_bottom.plot(time,r_distractor[2][4-tdindex,:],color=colors_set1['cyan'])
   
        
        
        axes_corr_bottom.set_xlabel('Time (ms)')
        axes_corr_top.set_xlim(100,xlim)
        axes_corr_top.set_xticks([200,400,600,800,1000])

        axes_corr_top.set_xticklabels([])
        axes_corr_top.set_yticks([0,40,80])
        axes_corr_bottom.set_yticks([0,40,80])
        axes_corr_bottom.set_xlim(100,xlim)
        axes_corr_bottom.set_xticks([200,400,600,800,1000])
        axes_corr_bottom.set_xticklabels([0,200,400,600,800])

        axes_corr_top.set_ylim(0,80)
        axes_corr_bottom.set_ylim(0,80)
        plotticks(axes_corr_top)
        plotticks(axes_corr_bottom)
        #fig_corr_top.savefig("gottfig2_top_sep.pdf")
        #fig_corr_bottom.savefig("gottfig2_bottom_sep.pdf")
    
   
    if  errorfig==1:
        ax.plot(tdistractortimes-200, error*100, '.',color = colors_set1['blue'])
        ax.plot(tdistractortimes-200, error_nopfc*100, '.',color = colors_set1['red'])
        ax.plot(np.arange(100,400,10), fit(alldistractor[4],'exponential',[36,0.03,0])[1], color = colors_set1['blue'],label = 'control')
        ax.plot(np.arange(100,400,10), fit(alldistractor_nopfc[4],'exponential',[36,0.03,80])[1],color = colors_set1['red'] , label = 'no PFC')
        ax.set_ylim(-0.1,105)
        ax.set_xlim(50,330)
        ax.set_xticks([100,200,300])
        plotticks(ax)
        ax.set_ylabel(r'Error rate (\%)')
        ax.set_xlabel('TDOA (ms)')
        py.legend(bbox_to_anchor=(1.2, 1.12), frameon=False,numpoints=1)
        #fig2.savefig("errorrate_sep.pdf")
    
    #py.show()

    
    
def fit(vectorfit,type_fit,p0fit):
    firsttime = timecue+cueduration
    initialtimefit = 100
    def func(x, a, b, c):
        if type_fit == 'exponential':
            return c-a* np.exp(-b * (x-initialtimefit)) 
            
        elif type_fit == 'linear':
            return a * (x-initialtimefit) + b
    
    tdistractortimes = np.array([firsttime,firsttime+50, firsttime+100,firsttime+200,2500])

    timeforfit = (tdistractortimes[0:4]-200)
    yforfit = vectorfit[0:4]*100
    popt, pcov = curve_fit(func, timeforfit, yforfit,p0 = p0fit )
    
    return popt,func(np.arange(initialtimefit,400,10), *popt)

def timeconstantfit():
    initialtimefit = 200
    endtimefit = 470
    totaltrialsfit = 1
    I=simulation(totaltrialsfit,0.1,0,1,0)
    meanI=np.mean(I, axis=1)
    rL1=FIcurve(meanI[0,...],a)
    rR1=FIcurve(meanI[1,...],a)
    rL2=FIcurve(meanI[2,...],a)
    rR2=FIcurve(meanI[3,...],a)
   
    def func(x, a, b, c):
        return a * np.exp(-b * (x-initialtimefit)) + c
    
   
    timeforfit = time[initialtimefit:endtimefit]#time[np.where(time==initialtimefit)[0]: np.where(time==endtimefit)[0]]
    yR1forfit=rR1[initialtimefit:endtimefit]#[np.where(time==initialtimefit)[0]: np.where(time==endtimefit)[0]]
    yR2forfit=rR2[initialtimefit:endtimefit]#[np.where(time==initialtimefit)[0]: np.where(time==endtimefit)[0]]
    
    x=timeforfit
    
    
    poptR1, pcovR1 = curve_fit(func, timeforfit, yR1forfit)
    poptR2, pcovR2 = curve_fit(func, timeforfit, yR2forfit)
    
    
    fig=py.figure(frameon=False, figsize=(1.6,1.2))
    axes = fig.add_subplot(111)

    axes.plot(timeforfit, func(timeforfit, *poptR1), 'k--', linewidth=1.5)#, label="PPC exponential fit")
    axes.plot(timeforfit, yR1forfit, color=colors_set1['red'], label=r"PPC" )

    axes.plot(timeforfit, func(timeforfit, *poptR2), 'k--', linewidth=1.5)#), label="PFC exponenential fit")
    axes.plot(timeforfit, yR2forfit, color=colors_set1['purple'], label=r"PFC")
    axes.set_xlabel('Time (ms)')
    axes.set_ylabel('Firing rate (spikes/s)')
    
    axes.set_xticks([200,350,500])
    axes.set_xticklabels([0,150,300])
    axes.set_yticks([0,1.5,3])
    
    print round(1/poptR1[1],2)
    print round(1/poptR2[1],2)
    #axes.legend(loc='upper right', frameon=False)
    #axes.show()

    left, bottom, width, height = [0.62, 0.54, 0.28, 0.3]#[0.75, 0.65, 0.15, 0.2]
    ax_inset = fig.add_axes([left, bottom, width, height])
    plot_inset(ax_inset, 'suppression',poptR1, poptR2)
    ax_inset.set_xticks([])
    ax_inset.set_ylim([0,80])

    ax_inset.set_yticks([0,40,80])
    ax_inset.set_title(r"$\tau_{\rm sup}$ (ms)" )
    plotticks(axes)
    
    
    
    
    
    #fig.savefig("WMtimeconstantssep.pdf")
    
    return poptR1, poptR2, rL1, rL2    

def plot_inset(ax_inset, typeplot,popt1,popt2):
    bar_width = 0.03
    index = 0
    
    if typeplot=='suppression':
        slope1 = round(1/popt1[1])
        slope2 = round(1/popt2[1])
        col1=colors_set1['red']
        col2=colors_set1['purple']
    elif typeplot=='intrinsic':
        slope1 = round(1/popt1[1])
        slope2 = round(1/popt2[1])
        col1=pp.colors_set1['blue']
        col2=pp.colors_set1['cyan']
        #print slope1
        #print slope2
        
        
    
    rects1 = ax_inset.bar(index, slope1, bar_width,
                 
                 color = col1,
                 
                 label = 'PPC')

    rects2 = ax_inset.bar(index + bar_width, slope2, bar_width,
                 
                 color=col2,
                 
                 label='PFC')
    
    ftwo.plotticks(ax_inset)
    #py.show()
def correlation(timeseries1,timeseries2):
    totaltrialsfit = 1
    #I= simulation(totaltrialsfit,0,0,1,0)
    
    #meanI =np.mean(I, axis=1)

    rL1 = timeseries1#FIcurve(meanI[0,...],a)
    #rR1=FIcurve(meanI[1,...],a)
    rL2 = timeseries2#FIcurve(meanI[2,...],a)
    #rR2=FIcurve(meanI[3,...],a)
    
    
    end_corr = 7900#8000#pp.time[-1]#34000
    ini_corr = 1000#2000#1000
    corr_range = (pp.time<end_corr)*(pp.time>ini_corr)
    print np.shape(rL1)
    rL1_filt = ftwo.filter_autocorr(rL1,20)
    rL2_filt = ftwo.filter_autocorr(rL2,20)
    #print np.size(rL1_filt)
    #print np.size(corr_range)
    rL1_filt[corr_range]
    rL1_autocorr = ftwo.autocorr(rL1_filt[corr_range])
    rL2_autocorr = ftwo.autocorr(rL2_filt[corr_range])
    
    
    return corr_range,rL1_autocorr, rL2_autocorr,rL1_filt, rL2_filt

def correlationfig(autocorr1,autocorr2,corr_range):
    rL1_autocorr = autocorr1
    rL2_autocorr = autocorr2

    
    fig2=py.figure(frameon=False, figsize=(1.6,1.2))# (3,1.65)
    axes = fig2.add_subplot(111)
    
   
    
   
    axes.set_xlabel('Time (ms)')
    axes.set_ylabel('Autocorrelation')
    axes.set_xticks([0,500,1000])
    axes.set_yticks([-0.5,0,0.5,1])
    py.ylim(-0.5,1.05)
    #
    ftwo.plotticks(axes)
    
    
    initialtimefit = 0
    endtimefit  = 1000
    fit_range = (pp.time<endtimefit)*(pp.time>initialtimefit)
    timeforfit = pp.time[fit_range]
    def func(x, a, b, c):
            return (a * np.exp(-b * (x-initialtimefit)) + c)#*step(x-initialtimefit,0)
        
       
    yL1forfit = rL1_autocorr[fit_range]#rR1[np.where(time==initialtimefit)[0]: np.where(time==endtimefit)[0]]
    yL2forfit = rL2_autocorr[fit_range]#rR2[np.where(time==initialtimefit)[0]: np.where(time==endtimefit)[0]]

   
    
    poptL1, pcovL1 = curve_fit(func, timeforfit, yL1forfit)
    poptL2, pcovL2 = curve_fit(func, timeforfit, yL2forfit)
        
    axes.plot(timeforfit, func(timeforfit, *poptL1), 'k--', linewidth = 1) 
    axes.plot(timeforfit, func(timeforfit, *poptL2), 'k--', linewidth = 1)  
    #axes.plot(timeforfit, guess_R1[0]*np.exp(-guess_R1[1] * (timeforfit-initialtimefit)) + guess_R1[2],'red')
    axes.set_xlim(0,1000)
    
    axes.plot(rL1_autocorr, color = pp.colors_set1['blue'], label='PPC')
    axes.plot(rL2_autocorr, color = pp.colors_set1['cyan'], label='PFC')
    
    #print round(1/poptL1[1],2)
    #print round(1/poptL2[1],2)
    left, bottom, width, height = [0.62, 0.54, 0.28, 0.3]#[0.75,0.58, 0.12, 0.3]#[0.75, 0.65, 0.15, 0.2]
    ax_inset = fig2.add_axes([left, bottom, width, height])
    
    plot_inset(ax_inset, 'intrinsic',poptL1,poptL2)
    ax_inset.set_xticks([])
    #ax_inset.yaxis.set_label_coords(1.65, 0.5) 
    #ax_inset.yaxis.tick_right()

    ax_inset.set_ylim([0,500])

    ax_inset.set_yticks([0,250,500])
    ax_inset.set_title(r"$\tau_{ fluct}$ (ms)")
    return round(1/poptL1[1],2), round(1/poptL2[1],2)
    #fig2.savefig("autocorrelationjune.pdf")
    

def load_data():
    import os
    os.chdir('/Users/jorgejaramillo/Documents/Work/loopmodel/loopdata')

    alldistractor=np.load('alldistractornov20.npy')
    alldistractor_nopfc = np.load('alldistractor_nopfc.npy')
    corr_range,rL1_autocorr, rL2_autocorr= np.load('intrinsiccorrelation.npy')#correlation() # intrinsic time sacles


    
    return corr_range,rL1_autocorr, rL2_autocorr#alldistractor
#alldistractor = load_data()
#corr_range, autocorr1, autocorr2 = correlation()#load_data()    
#alldistractor= distraction_time_data()
distraction_time_fig('new')    
    
    #
#correlationfig() #intrinsic timescales
#popt,popt2,rL1, rL2=timeconstantfit()   



memoryplusdistractors('correcttrial')

#alldistractor_nopfc = distraction_time_data()
#

#alldistractor = np.load('alldistractornov20.npy')
