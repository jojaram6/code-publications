from __future__ import division
import numpy as np
import pylab as py
import matplotlib as mp
import random
#from pulvinarparamsWM import *  ### Iext5(time,'dm', 0, 30)
from functionstwoloopLR_pnumber import *
from maintwoloopLR import simulation
#from thalamswitch import *
#import workingmemoryfigure as wm


from scipy.integrate import odeint
font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 7}

mp.rc('font', **font)



def main():
    
    ''' calculates all firing rates as a function of contrast, without using the cluster'''
    allmain = np.zeros((6,len(coh1vec),trials_transition, lenp,len(time)))
    for coh_index, coh in enumerate(coh1vec):  
        allmain[:, coh_index, :] = main_for_cluster(coh_index)
    return allmain


def  main_for_cluster(coh_index):
     
     if exp == 'wm':
         
         coh1 = currentcue#0.08
         coh2 = currentdistractor# ############### 0 if there is no distractor
         coh4 = 0#coh1
         coh3 = 0#-coh4
     elif exp == 'dm':
         
         if decision=='conflict':
             coh1 = coh1vec[coh_index]
             
             coh2 = -coh1
             coh4 = coh1
             coh3 = -coh4
         elif decision=='visual':
             coh1 = coh1vec[coh_index]
             coh2 = -coh1
             coh4 = -coh1
             coh3 = -coh4
         elif decision=='visualandreward':
              coh1 = coh1vec[coh_index]
              coh2 = -coh1
              coh4 = -frac_reward
              coh3 = -coh4
     for ratio_LR in np.arange(0,1):
         for ratio_p in np.arange(0,lenp):
              ratio = [ratio_LR, ratio_p]
              
                
              for trial_index in np.arange(0,trials_transition):
                    print(trial_index)
                    print(coh1)
                
                    Itest, Stest=simulation(1,coh1, coh2,coh3,coh4, inhpercenttest, ratio)
                    
                    
                    
                    #def FIcurve(I):
                        #return (I)
                    
                    rL1test=FIcurve(Itest[0,sim,...],a_mod1)#,coherenceindex, 0])
                    rR1test=FIcurve(Itest[1,sim,...],a_mod1)#,coherenceindex,0])
                    rL2test=FIcurve(Itest[2,sim,...],a_mod2)#,coherenceindex,0])
                    rR2test=FIcurve(Itest[3,sim,...],a_mod2)#,coherenceindex,0])
                    rfixtest=FIcurve(Itest[8,sim,...],a)
                    #print np.shape(rL1test)
                    
                    if transitionplot==1:
                        lambda_p = weights(actualsetting,ratio)[2]
                    else:
                        lambda_p = weights(actualsetting, ratio)[2]
                    #print lambda_p    
                    if cx_impl == ' yes':
                        rLptest = FIcurve_2(Itest[4,sim,...],lambda_p)#,coherenceindex,0])
                        rRptest = FIcurve_2(Itest[5,sim,...],lambda_p)
                    else:
                        rLptest=FIcurve_2(Itest[4,sim,...],lambda_p)#,coherenceindex,0])
                        rRptest=FIcurve_2(Itest[5,sim,...],lambda_p)#,coherenceindex,0])
                    
                    rLpptest=FIcurve(Itest[6,sim,...],a)#,coherenceindex,0])
                    rRpptest=FIcurve(Itest[7,sim,...],a)#,coherenceindex,0])
                    
                    
                    sL1test = Stest[0,...]
                    sR1test = Stest[1,...]
                    sL2test = Stest[2,...]
                    sR2test = Stest[3,...]
                    sLptest = Stest[4,...]
                    sRptest = Stest[5,...]
                    sfixtest = Stest[6,...]  
                    
                    rL1all[trial_index,ratio_p,:] = rL1test
                    rL2all[trial_index,ratio_p,:] = rL2test
                    rR1all[trial_index,ratio_p,:] = rR1test
                    rR2all[trial_index,ratio_p,:] = rR2test
                    sL1all[trial_index,:] = sL1test
                    sL2all[trial_index,:] = sL2test
                    sR1all[trial_index,:] = sR1test
                    sR2all[trial_index,:] = sR2test
                    sLpall[trial_index,:] = sLptest																
                    sfixall[trial_index,:] = rfixtest

                    
                    
                    rLpall[trial_index,ratio_p,:] = rLptest
                    rRpall[trial_index,ratio_p,:] = rRptest
                    rfixall[trial_index,:] = rfixtest
                
        
    
     r_all_vector = [rL1all,rL2all,rR1all,rR2all,rLpall,rRpall]
     if usecluster=='yes':
         np.save('pulvinardata/confidence1/rvector'+str(coh_index), r_all_vector)           


     if transitionplot==1:
        return rL1all,rL2all,rR1all,rR2all,rLpall,rRpall, rfixall, sfixall, vec_ratio
     else:
        return rL1all,rL2all,rR1all,rR2all,rLpall,rRpall#, sLpall#, vec_ratio ##rL1all,rL2all,rR1all,rR2all,rLpall,rRpall,sL1all,sL2all,sR1all,sR2all use for confidence

def get_clusterfiles():
    allcluster = np.zeros((6,len(coh1vec),trials_transition, lenp, len(time)))
    for coh_index, coh_value in enumerate(coh1vec):
        clustfile = np.load('pulvinarclusterdata/desimone/control/rvector'+str(coh_index)+'.npy')
        allcluster[:, coh_index,:,0, :] = clustfile
    return allcluster
    
def plot_desimone(vecscont, vecsles,vecsles_distractor):
    fig=py.figure(frameon=False, figsize=(1.2,1))

    axes = fig.add_subplot(111)#allplots(3,1, 5,5)
    bar_width = 0.06
    index = 0
    winner_vec = np.zeros((3, trials_transition))
				
    colors = [colors_set1['blue'], colors_set1['purple'],colors_set1['green']]
    for vecs_index, vecs in enumerate([vecscont, vecsles, vecsles_distractor]):
			
        for trial_index in range(0,np.shape(vecscont)[2]):
            rL2 = vecs[1,0,trial_index,0,:]
            rR2 = vecs[3,0,trial_index,0,:]
            winner_vec[vecs_index, trial_index] = get_TandW(rL2,rR2)[1]
            #print winner_vec#[trial_index]									
        if vecs_index==2:
               correct = ( np.sum(winner_vec[vecs_index,:])/np.count_nonzero(winner_vec[vecs_index,:])-1)*100
               #print correct
        else:
               correct = (2-np.sum(winner_vec[vecs_index,:])/np.count_nonzero(winner_vec[vecs_index,:]))*100   
               #print correct

        axes.bar(index+0.08*vecs_index+bar_width, correct, bar_width,color = colors[vecs_index])
    axes.set_ylim(0,50)#0,60
    axes.set_ylabel(r'Errors %')
    axes.set_xticks([])
    plotticks(axes)
    fig.savefig('lesionfigs/desimonebarplotapril30.pdf')

    return winner_vec
 

def counter_transition(tvector1, tvector2, tvector3, tvector4,transition_type):
    if transition_type =='L':
    
        if tvector1[endtrial]>tvector3[endtrial]:
           if tvector2[endtrial]>tvector4[endtrial]:
               counter = 2
           else: 
               counter = 1
        else:
            counter=0
    elif transition_type=='congruent':
        if tvector1[endtrial]>tvector3[endtrial] and tvector2[endtrial]>tvector4[endtrial]:
           counter =2
        elif tvector1[endtrial]<tvector3[endtrial] and tvector2[endtrial]<tvector4[endtrial]:
            counter = 1
        else:
            counter = 0
    return counter

def transition():
    vec_transition = np.zeros((lenp,len(coh1vec),trials_transition))
    	#rL1all,rL2all,rR1all,rR2all,rLpall,rRpall	 = alltrialstransition(:,:)	
    for i_p in np.arange(0,len(lambdas)):
        for i_coh in np.arange(0,len(coh1vec)):
            for i_trial in np.arange(0,trials_transition):
                tvector1,tvector2,tvector3,tvector4 = alltrialstest[0:4,i_coh,i_trial,i_p,:]	
                vec_transition[i_p,i_coh,i_trial] = counter_transition(tvector1, tvector2, tvector3, tvector4,'congruent')					
    vec_transition_mean = np.sum(vec_transition, axis=2)/np.count_nonzero(vec_transition,axis=2)-1
    return vec_transition, vec_transition_mean
def plot_transition():
    vec_transition, vec_transition_mean = transition()
    #vecs_ptransition = vecs_p[6]
    fig=py.figure(frameon=False, figsize=(2,1.6))
    axes = fig.add_subplot(111)
    axes.plot(lambdas, vec_transition_mean[:,0], 'gray')
    axes.plot(lambdas, vec_transition_mean[:,1], 'k' )
    #axes.plot(lambdas, lowcoh[0,:,0], 'gray')
    #axes.plot(lambdas, highcoh[0,:,0], 'k' )
    axes.set_xlabel('Pulvinar gain $\lambda$ (Hz/nA)')#('Pulvinar structure')
    axes.set_ylabel('Prob. Cx1 switches Cx2')
    axes.set_ylim(0,1)
    axes.set_xlim(220,300)
    plotticks(axes)
    fig.savefig("switchprob_may1st.pdf")  

    
    
    
    
   
def plot_allrates():
    fig2,axes2 = allplots(3,1, 3,3)
    disc_trials = np.arange(0,trials_transition)#slow_fast(0)
    for i_coh in np.arange(0,len(coh1vec)):
        for i_trials in disc_trials:#np.arange(0,trials_transition):
            plot_rates(i_coh,i_trials)
            
def plot_rates(cohindex,trial,i_lambda,vecs, det):
    
    """" Returns a plot of cortex1, cortex2 and pulvinar, for confidence, conflict, and working memory"""
    #vec_ratio, rL1test, rR1test, rL2test, rR2test, rLptest, rRptest, rLpptest, rRpptest, winner_vec,all_latencies,rp_analytical = main()
    if transitionplot==1:
        rL1all,rL2all,rR1all,rR2all,rLpall,rRpall, rfixall, sfixall, vec_ratio= vecs #main()
    else:
        rL1all,rL2all,rR1all,rR2all,rLpall,rRpall = vecs 
    
    
    #############################################################    
    ############ Definition of firing-rate variables#############  
    #############################################################
    rL1test = rL1all[cohindex,trial,i_lambda,:]
    rL2test = rL2all[cohindex,trial,i_lambda,:]
    rR1test = rR1all[cohindex,trial,i_lambda,:]
    rR2test = rR2all[cohindex,trial,i_lambda,:]
    if det=='absolute':
       rLptest = np.abs(rL1test-rR1test)
    else:
       rLptest = rLpall[cohindex,trial,i_lambda,:]
       rRptest = rRpall[cohindex,trial,i_lambda,:]
    #rfixtest = rfixall[cohindex,trial,:]
    #sfixtest = sfixall[cohindex,trial,:]

    #rPPtest = rPP[cohindex,trial,:]
    
    
    if exptype=='confidence':
        fig2,axes = allplots(2,1, 1,2.8) #width= 1, height = 2.8
        ax1,ax2 = axes    
    				
        
        ax1.plot(time, rL1test, color = colors_set1['blue'])
        ax1.plot(time, rR1test, color = colors_set1['red'])
        #ax1.set_ylim(8,15)
        for ax in axes:
            ax.set_xlim(0,730)
            ax.set_xticks([0,650])
        ax1.set_ylim(0,40)
        ax1.set_xticklabels([])
        ax1.set_yticks([0,20,40])
        #ax2.plot(time, rL2test,color = colors_set1['blue'])
        #ax2.plot(time, rR2test,color = colors_set1['red'])
        ax2.set_ylim(0,22)
        ax2.set_yticks([0,10,20])
        #ax2.set_xticklabels([])
    
    
    #ax2.plot(time,RTthreshold*np.ones(len(time)),'k--')
        ax2.plot(time, rLptest,color = colors_set1['green'])# colors_set1['green'])
        #ax2.plot(time, rRptest,color = colors_set1['orange'])# colors_set1['green'])
        #ax3.plot(time, rate_abs(rL1test,rR1test)[1], 'k')							

        
    #ratax3.set_ylim(0,5)
    #ax3.set_yticks([0,50,100])
    
        #ax3.set_ylim(-0.03,0.03)
        #ax3.set_yticks([0,10,20])
        print (get_TandW(rL1test,rR1test)[0])
    
    
        ax2.set_xlabel('Time (ms)')
        fig2.savefig('confidencefigs/error_decay.pdf')
    elif exptype=='conflict':#'gainmod':
        fig2,axes = allplots(3,1,1.4,2.4) #width= 1.4, height = 2.4
        ax1,ax2,ax3 = axes    
        for ax in axes:
	       
            ax.set_xlim(0,800)							
        
        ax1.plot(time, rL1test, color = colors_set1['blue'])
        ax1.plot(time, rR1test, color = colors_set1['red'])
        #ax1.set_ylim(0,50)
        ax1.set_xticks([0,250,500,750])
        ax1.set_xticklabels([])
        ax1.set_yticks([0,15,30])
        ax2.plot(time, rL2test,color = colors_set1['blue'])
        ax2.plot(time, rR2test,color = colors_set1['red'])
        ax2.set_yticks([0,15,30])
        #ax2.set_ylim(0,70)
        ax2.set_xticks([0,250,500,750])
        ax2.set_xticklabels([])
    
    
        ax3.plot(time, rLptest,color = colors_set1['green'])# colors_set1['green'])
        #ax3.set_xlim(100,260)
        ax2.set_xticklabels([0,250,500,750])

    
        ax3.plot(time, rRptest, color = colors_set1['orange'])
    
        #ax3.set_ylim(0,12)
        ax3.set_yticks([0,10])#30,60])
        reactiontimes, winner = get_TandW(rL2test,rR2test)
        print (reactiontimes)
        fig2.savefig('gainmodfigs/conflictWMnomod.pdf')#fig2.savefig('gainmodfigs/conflictlargelambdahighconflict.pdf')		
	#elif exptype=='c			
				
				
				
    elif exptype=='wmgate':
        fig2,axes = allplots(3,1, 1.4,2.4) #width= 1.4, height = 2.4
        ax1,ax2,ax3 = axes    
        xticks = [0,1000,2000]
        labels =[]
        ylim = [0,134]#0,104

        for ax in axes:
            ax.set_xticks(xticks)
            ax.set_xticklabels(labels)
            ax.set_ylim(ylim)
        yticks = [0,40,80]#[0,2,4]#
        ax1.plot(time, rL1test, color = colors_set1['blue'])
        ax1.plot(time, rR1test, color = colors_set1['red'])
        
        ax1.set_yticks(yticks)
        #ax1.set_xticks(xticks)
        ax1.set_xticklabels(labels)
        ax2.plot(time, rL2test,color = colors_set1['blue'])
        ax2.plot(time, rR2test,color = colors_set1['red'])
        ax2.set_yticks(yticks)
        ax3.plot(time, rLptest,color = colors_set1['green'])# colors_set1['green'])
        ax3.set_xticklabels([0,1000,2000])
        ax3.plot(time, rRptest, color = colors_set1['orange'])
        ax3.set_yticks([0,40,80])#yticks)

    
        

        
    
        ax3.set_xlabel('Time (ms)')
        fig2.savefig('gainmodfigs/conflictWMmod2.pdf')        
    

    return 0#np.transpose(rp_analytical[0]),np.transpose(rp_analytical[1])
def get_TandW(rL,rR): ## rL and rR are time vectors
    hit_threshold1 = np.where(abs(rL-RTthreshold)<RTepsilon)[0]
    hit_threshold2 = np.where(abs(rR-RTthreshold)<RTepsilon)[0]
    #print hit_threshold1				
				
    if len(hit_threshold1)!=0:
       hit_threshold1 = hit_threshold1[hit_threshold1>100]
    if len(hit_threshold2)!=0:
       hit_threshold2 = hit_threshold2[hit_threshold2>100]

				
    winner, hit_winner = winnercounter(hit_threshold1, hit_threshold2)
    reactiontimes= time[np.min(hit_winner)]-timecue
    return reactiontimes, winner#, hit_threshold1

    
def print_latency(rL,rR):
    T = get_TandW(rL,rR)[0]
    #print T
    return T
def get_latency(vecscont, vecsles):
    latency_vec = np.zeros((2,trials_transition))
    for vecs_index, vecs in enumerate([vecscont, vecsles]):
			
        for trial_index in range(0,np.shape(vecscont)[2]):
            rL2 = vecs[1,0,trial_index,0,:]
            rR2 = vecs[3,0,trial_index,0,:]
            latency_vec[vecs_index, trial_index] = get_TandW(rL2,rR2)[0]				
    latencymean_vec = np.mean(latency_vec, axis =1)		
    latencystd_vec = np.std(latency_vec, axis =1)		

    return latencymean_vec, latencystd_vec	
def latency_plot(meanlatency,stdlatency):
    fig=py.figure(frameon=False, figsize=(0.6,1))
    ax = fig.add_subplot(111)#allplots(3,1, 5,5)
    N = 1
    control, lesion = meanlatency
    controlStd, lesionStd = stdlatency			
    #control = np.mean(vec_control)#[18, 35, 30, 35, 27]
    #controlStd =   np.std (vec_control)#[2, 3, 4, 1, 2]
    #lesion = np.mean(vec_lesion)
    #lesionStd =   np.std(vec_lesion)

## necessary variables
    ind = np.arange(N)                # the x locations for the groups
    width = 0.0081#0.006                   # the width of the bars
    rects1 = ax.bar(ind, control, width,
                color=colors_set1['blue'],
                yerr=controlStd,
                error_kw=dict(elinewidth=2,ecolor='k'))

    rects2 = ax.bar(ind+1.3*width, lesion, width, color= colors_set1['purple'], yerr=lesionStd, error_kw=dict(elinewidth=2,ecolor='k'))
    ax.set_xticks([])
    ax.set_ylim(100,370)
    plotticks(ax)
    fig.savefig("latencyplotapril302018.pdf")
# axes and labels




def plot_choice():
    
    """ plot that reproduces Wilke Experiment: contralesional effects and how reward
    ameliorates them"""
    fig=py.figure(frameon=False, figsize=(3,3))
    axes = fig.add_subplot(111)#allplots(3,1, 5,5)
    ones_control = np.mean(vecs_pcontrol[3])*100#np.count_nonzero(vecs_pcontrol[3])/250*100
    zeros_control= 100-ones_control
    ones_les = np.mean(vecs_ples[3])*100#np.count_nonzero(vecs_ples[3])/250*100
    zeros_les= 100-ones_les
    ones_plusreward = np.mean(vecs_pplusreward[3])*100#np.count_nonzero(vecs_pplusreward[3])/250*100
    zeros_plusreward = 100-ones_plusreward
    bar_width = 0.06
    index = 0
    
   
    axes.bar(index, ones_control, bar_width,color = colors_set1['blue'],label='PPC')

    axes.bar(index + bar_width, ones_les, bar_width, color = colors_set1['purple'], label='PFC')
                 
    axes.bar(index+2*bar_width, ones_plusreward, bar_width, color = colors_set1['green'])
             
    axes.bar(index+4*bar_width, zeros_control, bar_width, color = colors_set1['blue'])

    axes.bar(index + 5*bar_width, zeros_les, bar_width,color = colors_set1['purple'])
     
    axes.bar(index + 6*bar_width, zeros_plusreward, bar_width,color = colors_set1['green'])           
              
          
    
#    axes.hist(vecs_p[9], color ='red')
#    axes.hist(vecs_pcontrol[9]+0.3, color ='blue')
#    #axes.hist(vecs_pcontrol[7]+0.3, color ='blue')


    axes.set_xticks([])
    axes.set_xlim([0,.5])
    axes.set_ylim(0,100)
    plotticks(axes)
    #fig.savefig("choice_rewardAug.pdf")

    return ones_control, zeros_control, ones_les, zeros_les,ones_plusreward, zeros_plusreward

def calc_corr (allfile):
    rL1, rL2 = allfile[0:2]
    t_1 = np.zeros(lenp)
    t_2 = np.zeros(lenp)			
    #print 'RL1 isd%', np.shape(rL1)

    
    for i_lambda in range(0,lenp):
        rL1_old = 0
        rL2_old = 0
        rL1_mean = 0
        rL2_mean = 0					
        for i in range(0,trials_transition):
            corr_range, rL1_autocorr, rL2_autocorr, rL1_filt, rL2_filt = wm.correlation(rL1[:,i,i_lambda].flatten(),rL2[:,i,i_lambda].flatten())
            rL1_mean = rL1_autocorr+rL1_old
            rL1_old = rL1_mean
            rL2_mean = rL2_autocorr+rL2_old
            rL2_old = rL2_mean
        rL1_mean = rL1_mean/trials_transition
        rL2_mean = rL2_mean/trials_transition
    #print np.shape(rL2_mean)
        t_1[i_lambda], t_2[i_lambda] = wm.correlationfig(rL1_mean, rL2_mean,corr_range)
    return t_1, t_2							
def calc_corr_fig(times1,times2):    
    lw = 2.4
    fromjorgem = [4.8212,  5.0189, 5.2530, 5.5353, 5.7011, 5.6723]
    fig=py.figure(frameon=False, figsize=(3,3))
    axes = fig.add_subplot(111)#allplots(3,1, 5,5)	
    font = {'family' : 'sans-serif',
    'weight' : 'normal',
    'size'   : 7}
    mp.rc('font', **font)
			
				
    axes.plot(np.arange(0,len(lambdas)), times2-times1, '.', color = colors_set1['green'] )
    axes.plot(np.arange(0,len(lambdas)), times2-times1, color = colors_set1['green'], linewidth = lw)

    #axes.plot(lambdas,0.7*t_2, 'g.')
    #axes.plot(lambdas,1.3*t_1, 'r.')
			
    ax2 = axes.twinx()
    ax2.plot(np.arange(0,len(lambdas)), fromjorgem, '.', color = colors_set1['blue'])	
    ax2.plot(np.arange(0,len(lambdas)), fromjorgem, color = colors_set1['blue'], linewidth = lw)		
    #axes.set_xlim(250,310)
    axes.set_ylim(0,350)
    ax2.set_ylabel('V4-IT hierarchical distance ', color = colors_set1['blue'], fontsize = 7)		
    axes.set_ylabel('V4-IT intrinsic-timescale difference (ms)',color = colors_set1['green'] )		
    axes.tick_params(axis='y', colors=colors_set1['green'])
    ax2.tick_params(axis='y', colors=colors_set1['blue'] )

    axes.set_xlim(-0.3,5.3)
    axes.set_xlabel('Pulvinar gain (normalized)')				
    plotticks(axes)
    plotticks2(ax2)
    fig.savefig('hierarchy2.pdf')			
				
#vecs_pall = np.load('vecs_pall.npy')
#timescales = np.load('timescales_jorgemapril8.npy')
#vecs_pplusreward = main()
#vecs_p = main(0)

