from __future__ import division
import numpy as np
import pylab as py
import matplotlib as mp
import random
from pulvinarparams import *  ### Iext5(time,'dm', 0, 30)
from functionstwoloopLR_pnumber import *
from maintwoloopLR import simulation
from thalamswitch import *
from reaction_inhibition import inhanalysis
from reaction_inhibition import analysisonset
from reaction_inhibition import plot_reactionfig
import thalamusloop as tl


from scipy.integrate import odeint
font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 7}

mp.rc('font', **font)


def winner(coh,trial_type,vecs,epsilon):
    
    """ returns a vector with the trials where a correct choice was made. The trials are divided into two
   'first' (first half of the trials)
   'second' (second half of the trials)"""
    
    
    rL1all,rL2all,rR1all,rR2all,rLpall,rRpall = vecs
    #winnerall = np.zeros((len(coh1vec),trials_half))
    winnerall = np.zeros(trials_half)
    
    if trial_type == 'first':
        rL= rL1all[:,0:trials_half,:]
        rR = rR1all[:,0:trials_half,:]

    elif trial_type == 'last':
        rL = rL1all[:,trials_half:trials_transition,:]
        rR = rR1all[:,trials_half:trials_transition,:]
    for trial_index in np.arange(0,trials_half):
		rL2test = filter_t(rL[coh,trial_index, :], 20)[disctime]
		rR2test = filter_t(rR[coh,trial_index,:],20)[disctime]
		

            
		if rL2test>rR2test and get_disc(rL2test,rR2test) > epsilon:
		   winnerall[trial_index] = 2
                
		elif rL2test<rR2test and get_disc(rL2test,rR2test) > epsilon:
			winnerall[trial_index] = 1


    
    
    '''for  coh_index in np.arange(0,len(coh1vec)):
        for trial_index in np.arange(0,trials_half):
            rL2test = rL[coh_index,trial_index, dectime]
            rR2test = rR[coh_index,trial_index,dectime]

            
            if rL2test>rR2test and get_disc(rL2test,rR2test) > epsilon:
                winnerall[coh_index, trial_index] = 2
                
            elif rL2test<rR2test and get_disc(rL2test,rR2test) > epsilon:
                winnerall[coh_index, trial_index] = 1'''
    return winnerall#[coh,:]
    
    
def escape(coherence,trial_type,vecs,epsilon):
    rL1all,rL2all,rR1all,rR2all,rLpall,rRpall = vecs
    trials_no = np.zeros(trials_half)
    if trial_type == 'first':
        rL = rL1all[:,0:trials_half,:]
        rR = rR1all[:,0:trials_half,:]

    elif trial_type == 'last':
        rL = rL1all[:,trials_half:trials_transition,:]
        rR = rR1all[:,trials_half:trials_transition,:]

	rL = filter_t(rL,10)
	rR = filter_t(rR,10)
    for index_trial in np.arange(0,trials_half):
        
        if get_disc(rL[coherence,index_trial,disctime],rR[coherence,index_trial,disctime]) < epsilon:
           trials_no[index_trial]=index_trial 
         
    trials__no=trials_no[trials_no>0]      #index_no = index_no+1
    return (trials__no), trials_no>0
 
def RT_definition(vecs):
    
    rL= vecs[0]
    rR = vecs[2]				
    RTvec = np.zeros((len(coh1vec), trials_transition))

    choicevec = np.zeros((len(coh1vec),trials_transition))
    for coh_index, coh in enumerate(coh1vec):
        for trial_index in range(0,trials_transition):
            reactiontimes, choice = tl.get_TandW(rL[coh_index, trial_index,:],rR[coh_index, trial_index,:])[0:2]
            RTvec[coh_index,trial_index] = reactiontimes
            choicevec[coh_index,trial_index] = choice  
    return RTvec,choicevec    
def RT_plot(vecs1,vecs2): #vecs1, vecs2 are control and lesion arrays of shape (ratesineach area, totalcoh, trials, time)
    plotvec = np.zeros((2,totalcoh,2))## first dim = reactiontime or choice, second dim= coherence level, third dim= control or lesion
    vectuple = (vecs1,vecs2)
    for vec_index, vec in enumerate(vectuple):
        RTvec, choicevec = RT_definition(vec)			
        for coh_index in range(0,totalcoh):
		  RTvec_ = RTvec[coh_index][RTvec[coh_index]>0]
		  choicevec_ = choicevec[coh_index][choicevec[coh_index]>0]	
		  plotvec[0,coh_index,vec_index] = np.mean(RTvec_)
		  plotvec[1,coh_index,vec_index]	= np.mean(choicevec_-1)
    return plotvec
                 
"""def RT_plot(RT, choice):
    fig=py.figure(frameon=False, figsize=(1.5,1.5))
    axes = allplots(2,1, 5,5)
    axes[0].plot(coh, choice)
    axes[1].plot(coh, reaction_time)
#    
"""
def confidence_behavior_plot(main,epsilon,plot_type):

    """ This is a plot of behavioral choice / performance as a function of contrast and trial type (correct, error, escape)
    as in Komura Fig: 3B """ 


    fig=py.figure(frameon=False, figsize=(1,1.52))
    axes = fig.add_subplot(111)#allplots(3,1, 5,5)
    wins_vec = np.zeros((3,2))
    escapes_vec = np.zeros((3,2))
    errors_vec = np.zeros((3,2))


    pointsize = 3
    cohuniform = [55,25,0, 45, 75, 100]
    cohuniform2 = [0,25, 45,55, 75, 100]


    for coh_index in np.arange(0, len(coh1vec)):
        w_first = winner(coh_index,'first',main,epsilon)
        wins =len(w_first[w_first>1])/(trials_half) #np.count_nonzero(winner(coh_index,'first',main))/(trials_half)
        escapes = np.size(escape(coh_index,'first',main,epsilon)[0])/(trials_half)
        errors = 1-wins-escapes
        
        w_last = winner(coh_index,'last',main,epsilon)
        wins_last = len(w_last[w_last>1])/(trials_half)#np.count_nonzero(winner(coh_index,'last',main))/(trials_half)
        escapes_last = np.size(escape(coh_index,'last',main,epsilon)[0])/(trials_half)
        errors_last = 1-wins_last-escapes_last

        wins_vec[coh_index,:] = wins, wins_last
        escapes_vec[coh_index,:] = escapes, escapes_last
        errors_vec[coh_index,:] = errors, errors_last
        

        if plot_type == 'muscimol':
		  axes.plot(cohuniform[coh_index], escapes, 'co', ms = pointsize)
		  axes.plot(cohuniform[coh_index+3], escapes_last, 'co', ms = pointsize)


        else:
		  			
		  axes.plot(cohuniform[coh_index], wins, 'ko', ms = pointsize)
		  axes.plot(cohuniform[coh_index], escapes, 'co', ms = pointsize)
		  axes.plot(cohuniform[coh_index], errors, color = colors_set1['pink'],marker = 'o', ms = pointsize)
		  axes.plot(cohuniform[coh_index+3], wins_last, 'ko', ms = pointsize)
		  axes.plot(cohuniform[coh_index+3], escapes_last, 'co', ms = pointsize)
		  axes.plot(cohuniform[coh_index+3], errors_last, color = colors_set1['pink'], marker = 'o',ms = pointsize)
        
    newwins_vec=[wins_vec[2,0],wins_vec[1,0], wins_vec[0,1], wins_vec[0,0],wins_vec[1,1],wins_vec[2,1]]
    newescapes_vec=[escapes_vec[2,0],escapes_vec[1,0], escapes_vec[0,1], escapes_vec[0,0],escapes_vec[1,1],escapes_vec[2,1]]
    newerrors_vec=[errors_vec[2,0],errors_vec[1,0], errors_vec[0,1], errors_vec[0,0],errors_vec[1,1],errors_vec[2,1]]

    if plot_type == 'muscimol':
        axes.plot(cohuniform2, newescapes_vec, 'c', ms = pointsize)
        axes.plot(cohuniform2, escape_lesion, 'b', ms = pointsize)
        axes.plot(cohuniform2, escape_lesion, 'bo', ms = pointsize)
        axes.set_ylim(-0.05,0.6)
        plotticks(axes)
        axes.set_xlim(-1.5,103)
        axes.set_xlabel('up-down ratio%')
        axes.set_ylabel('Escape choice')
        #fig.savefig('confidencefigs/confidence_muscimol3.pdf')


					     
    else:
	   axes.plot(cohuniform2, newwins_vec, 'k', ms = pointsize)
	   axes.plot(cohuniform2, newescapes_vec, 'c', ms = pointsize)
	   axes.set_ylim(-0.05,1.05)
	   axes.plot(cohuniform2, newerrors_vec, color = colors_set1['pink'], ms = pointsize)
	   axes.set_ylabel('Behavioral choice')
	   plotticks(axes)
	   axes.set_xlim(-1.5,103)
	   axes.set_xlabel('up-down ratio%')
	   fig.savefig('confidencefigs/confidence_behaviordecay.pdf')

    
             
      
      
      
    
    #axes.set_xticks(coh1u)
    #axes.set_xticklabels([0,25,50])
    #axes.set_ylim(-0.05,0.65)
    
    
    return newescapes_vec
	
    
def confidence_rates_plot(rtype,vecs,epsilon,pos):
    
    '''This is a 3 panel plot firing_rate as a function of 
    time for correct,error, and escape trials and three different contrast conditions
    as in Komura Fig: 3C '''    
    
    
    fig, axes=allplots(1,3,2.8,0.8)
    LW = 1.6 #linewidth
    rL1all,rL2all,rR1all,rR2all,rLpall,rRpall = vecs
    #rLpall = rPP#calc_rp(rtype)
    if rtype=='absolute':
        rLpall = np.abs (rL1all- rR1all) 
    
    #rLpall = rLpall[:,int(trials_transition/2):int(trials_transition),:]
    if pos == 'first':
       rLpall = rLpall[:,0:int(trials_transition/2.0),:]
    else:
	  rLpall = rLpall[:,int(trials_transition/2):int(trials_transition),:]
								
				

    
    
    for coh_index in np.arange(0,3):
        axes[2-coh_index].plot(time, np.mean(rLpall[coh_index, winner(coh_index,pos,vecs,epsilon)==2,:], axis=0), color = 'black', linewidth =LW)
        if coh_index <7:
            axes[2-coh_index].plot(time, np.mean(rLpall[coh_index, winner(coh_index,pos,vecs,epsilon)==1,:], axis=0), color = colors_set1['pink'], linewidth = LW)
            axes[2-coh_index].plot(time, np.mean(rLpall[coh_index, escape(coh_index,pos,vecs,epsilon)[1],:], axis=0), color = colors_set1['cyan'], linewidth = LW)
        
        axes[coh_index].set_xticks([0,250,500,750])
        
        axes[coh_index].set_ylim(0,22)
        axes[coh_index].set_yticks([])

        axes[coh_index].set_xlim(0,750)
    axes[0].set_ylabel('Firing rate (Hz)') 
    axes[0].set_yticks([0,10,20])

    fig.savefig('confidencefigs/confidence_ratesdecay.pdf')
def mean_rate(vecs,pos,coh_index,epsilon,rate_index, trialtype):
    rp = vecs[rate_index]
    if pos == 'first':
	  rp = rp[:,0:int(trials_half),:]
    else:
       rp = rp[:,trials_half:trials_transition,:]
				
    if trialtype=='correct':
	  rp_w = rp[coh_index, winner(coh_index,pos,vecs,epsilon)==2,:]
    elif trialtype=='error':
	    rp_w = rp[coh_index, winner(coh_index,pos,vecs,epsilon)==1,:]			
	
    elif trialtype =='escape':	
	    rp_w = rp[coh_index, escape(coh_index,pos,vecs,epsilon)[1],:]
					
    mean = np.mean(rp_w,axis =0) 				
    return mean#,rp_w			

def plot_av(epsilon,vecs, coh_index,trialtype, pos):
    fig2,axes = allplots(2,1, 4,6) #width= 1.4, height = kk
    ax1,ax2 = axes   
    #for ax_index, ax in enumerate(axes):
    rL1 = vecs[0][coh_index,winner(coh_index,pos,vecs,epsilon)==1,:]
    rR1 = vecs[2][coh_index,winner(coh_index,pos,vecs,epsilon)==1,:]
    lenrp = np.shape(rL1)[0]			
    rp = np.zeros((lenrp,len(time)))			
    """for i_trials in range(0,lenrp):
	   print i_trials
	   #py.plot(time, rL1[i_trials,:],time, rR1[i_trials,:])
	   rp[i_trials,:] = rate_abs(rL1[i_trials,:],rR1[i_trials,:])[1]				
	"""			
    meanL1 = 	mean_rate(vecs,pos,coh_index,epsilon, 0, trialtype)
    meanR1 = mean_rate(vecs,pos,coh_index,epsilon, 2, trialtype)
    ax1.plot(time,meanL1,'b', time, np.mean(rL1, axis = 0) ,'c')
    ax1.plot(time, meanR1,'r', time, np.mean(rR1, axis = 0) ,'m')
    ax1.set_ylim(0,30)
    ax2.plot(time, mean_rate(vecs,pos,coh_index,epsilon, 4, trialtype),'g')
    ax2.plot(time, np.abs(meanL1-meanR1), 'k')
    ax2.plot(time, rate_abs(meanL1,meanR1)[1], 'r')	
    #ax2.plot(time, np.mean(rp, axis =0),'b')					

    ax2.set_ylim(0,20)
    for ax in axes:
	   ax.set_xlim(0,disctime)

	
	

			
    






def rates_behavior_plot(rtype,vecs,epsilon): 
    '''This is a plot of (firing_rate) normalized activity as a function of 
    contrast and trial type (correct, error, escape)
    as in Komura Fig: 3D'''




    
    wins_vec = np.zeros((len(coh1vec),2))
    escapes_vec = np.zeros((len(coh1vec),2))
    errors_vec = np.zeros((len(coh1vec),2))
    fig=py.figure(frameon=False, figsize=(1.,1.52))
    axes = fig.add_subplot(111)#allplots(3,1, 5,5)
    rL1all,rL2all,rR1all,rR2all,rLp_all,rRpall = vecs
    #rLp_all = rPP#calc_rp(rtype)
    if rtype=='absolute':
        rLp_all = np.abs (rL1all- rR1all) 
    
    elif rtype=='abstest':
		      for coh_index in np.arange(0,len(coh1vec)):
			     for trial_index in np.arange(0,trials_transition):
		              print trial_index
		              rL1coh = rL1all[coh_index,trial_index,:]
		              rR1coh = rR1all[coh_index, trial_index, :]
		              rp = rate_abs(rL1coh,rR1coh)[1]
		              rLp_all[coh_index,trial_index,:] = rp	
				
    rp_first = rLp_all[:,0:int(trials_transition/2.0),:]
    rp_last = rLp_all[:,int(trials_transition/2.0):int(trials_transition),:]
    cohuniform=[55,25,0, 45, 75, 100]
    cohuniform2 = [0,25, 45,55, 75, 100]

    if rtype=='absolute':
        norm = 0.07
    normfirst = 1.0/np.mean(np.mean(rp_first[2, winner(2,'first',vecs,epsilon)==2,time_ini:time_end], axis=0))
    normlast = 1.0/np.mean(np.mean(rp_last[2, winner(2,'last',vecs,epsilon)==2,time_ini:time_end], axis=0))

    pointsize = 3
    for coh_index in np.arange(0,len(coh1vec)):
        						
									
							
        axes.plot(cohuniform[coh_index], normfirst*np.mean(np.mean(rp_first[coh_index, winner(coh_index,'first',vecs,epsilon)==2,time_ini:time_end], axis=0)), color = 'black', marker = 'o', ms = pointsize)
        axes.plot(cohuniform[coh_index], normfirst*np.mean(np.mean(rp_first[coh_index, escape(coh_index,'first',vecs,epsilon)[1],time_ini:time_end], axis=0)), color = colors_set1['cyan'], marker = 'o', ms = pointsize)
        axes.plot(cohuniform[coh_index], normfirst*np.mean(np.mean(rp_first[coh_index, winner(coh_index,'first',vecs,epsilon)==1,time_ini:time_end], axis=0)), color = colors_set1['pink'], marker = 'o', ms = pointsize)        
        axes.plot(cohuniform[coh_index+3], normlast*np.mean(np.mean(rp_last[coh_index, winner(coh_index,'last',vecs,epsilon)==2,time_ini:time_end], axis=0)), color = 'black', marker = 'o', ms = pointsize)
        axes.plot(cohuniform[coh_index+3], normlast*np.mean(np.mean(rp_last[coh_index, escape(coh_index,'last',vecs,epsilon)[1],time_ini:time_end], axis=0)), color = colors_set1['cyan'], marker = 'o', ms = pointsize)
        axes.plot(cohuniform[coh_index+3], normlast*np.mean(np.mean(rp_last[coh_index, winner(coh_index,'last',vecs,epsilon)==1,time_ini:time_end], axis=0)), color = colors_set1['pink'], marker = 'o', ms = pointsize)
    
        wins_vec[coh_index,:] = normfirst*np.mean(np.mean(rp_first[coh_index, winner(coh_index,'first',vecs,epsilon)==2,time_ini:time_end], axis=0)), normlast*np.mean(np.mean(rp_last[coh_index, winner(coh_index,'last',vecs,epsilon)==2,time_ini:time_end], axis=0))
        escapes_vec[coh_index,:] = normfirst*np.mean(np.mean(rp_first[coh_index, escape(coh_index,'first',vecs,epsilon)[1],time_ini:time_end], axis=0)), normlast*np.mean(np.mean(rp_last[coh_index, escape(coh_index,'last',vecs,epsilon)[1],time_ini:time_end], axis=0))
        errors_vec[coh_index,:] = normfirst*np.mean(np.mean(rp_first[coh_index, winner(coh_index,'first',vecs,epsilon)==1,time_ini:time_end], axis=0)), normlast*np.mean(np.mean(rp_last[coh_index, winner(coh_index,'last',vecs,epsilon)==1,time_ini:time_end], axis=0))
    
    
    newwins_vec=[wins_vec[2,0],wins_vec[1,0], wins_vec[0,1], wins_vec[0,0],wins_vec[1,1],wins_vec[2,1]]
    newescapes_vec=[escapes_vec[2,0],escapes_vec[1,0], escapes_vec[0,1], escapes_vec[0,0],escapes_vec[1,1],escapes_vec[2,1]]
    newerrors_vec=[errors_vec[2,0],errors_vec[1,0], errors_vec[0,1], errors_vec[0,0],errors_vec[1,1],errors_vec[2,1]]

    axes.plot(cohuniform2, newwins_vec, 'k', ms = pointsize)
    axes.plot(cohuniform2, newescapes_vec, color = colors_set1['cyan'] , ms = pointsize)
    axes.plot(cohuniform2, newerrors_vec, color = colors_set1['pink'], ms = pointsize)
    
    
    axes.set_xlim(-2,103)
    axes.set_ylim(-0.01,1.1)
    axes.set_xlabel('up-down ratio%')
    axes.set_ylabel('Normalized activity')
    #axes.set_ylim(0,1.5)
    plotticks(axes)
    fig.savefig('confidencefigs/rates_behaviordecay2.pdf')    
    
def Jeqplot():
    rtest = np.arange(0,50,0.01)
    fig=py.figure(frameon=False, figsize=(2,2))
    axes = fig.add_subplot(111)
    axes.plot(rtest,J_eq(rtest,rateshift,amp))
    axes.set_ylim(-0.5,0.5)
    axes.set_xlabel('Rate (spikes/s)')
    axes.set_ylabel('Effective weight (nA)')
    plotticks(axes)
    #fig.savefig('Jeqplot.pdf')
    

    
def calc_rp(rtype):
    rL1all,rL2all,rR1all,rR2all,rLpall,rRpall,sL1all,sL2all,sR1all,sR2all = vecs_p_withs
    if rtype=='plasticity':
       for cohindex in np.arange(0,len(coh1vec)):
           print cohindex
           for trialindex in np.arange(0,trials_transition):
               rL1 = rL1all[cohindex,trialindex,:] 
               rR1 = rR1all[cohindex,trialindex,:] 
               rPall[cohindex,trialindex,:] = rate_total(rL1,rR1)[0]## with excitatory and inhibitory plasticity
       return rPall
    else:
        
        J_LP_eq = J_eq(rL1all, rateshift,amp)
        J_RP_eq = J_eq(rR1all, rateshift,amp)
        rP = FIcurve(J_LP_eq*sL1all+J_RP_eq*sR1all+Ib*Imotion2(time)+Iext5(time,'dm',400)) 
        return rP
    ## with simplified Jeq transfer function

def get_disc(rL,rR):
    return abs(rL-rR)   
    
def load_data():
    import os
    os.chdir('/Users/jorgejaramillo/Documents/Work/Thalamus-Pulvinar/pulvinardata')

    vecs= np.load('confidence_Dataoct31night.npy')
    
    return vecs
#vecs_p_withs = load_data()
#vecs_p = vecs_p_withs[0:6]
    
