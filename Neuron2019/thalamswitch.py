import numpy as np
import pylab as py
from functionstwoloopLR_pnumber import *
import matplotlib as mp
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 7}

mp.rc('font', **font)




def step(t, t_shift):
    return (np.sign(t-t_shift)+1)/2
    
    
def s_eq(plas_type,rate_ext,t_index,S,F,D):
    if plas_type == 'exc':
       tau = tau_E
       Seq = -S/tau+F*rate_ext/1000.0+Iexc
    if plas_type == 'inh':
       tau = tau_I
       Seq = -S/tau+p*D*rate_ext/1000.0+Iinh
    return Seq    
def dep_eq(t_index,rate_ext):
    D = np.zeros(lentime)
    D[0]=0.7
    #for t_index in np.arange(0, np.ceil(maxtime/stepsize_switch)-1):
    D[t_index + 1]= D[t_index] + stepsize_switch*dep_arg(D[t_index],t_index,rate_ext/1000.0)
    if inh_plasticity == 'yes':
        return D[t_index]
    else: 
        return D[0]*np.ones(len(time_p))
def dep_arg(D_, t_index,rate_ext):
    return -p*D_*rate_ext+(1-D_)/tau_D
    
def fac_eq(t_index,rate_ext):
    F = np.zeros(maxtime/stepsize_switch)
    F[0]=0.1
    #for t_index in np.arange(0, np.ceil(maxtime/stepsize_switch)-1):
    F[t_index+1]= F[t_index]+ stepsize_switch*fac_arg(F[t_index],t_index,rate_ext/1000.0)
    if exc_plasticity == 'yes':
        return F[t_index]
    else:
        return 0.2*np.ones(len(time_p))
def fac_arg(F,t_index,rate_ext):
    return a_F*(1-F)*(rate_ext)+(F_in-F)/tau_F #### change to U-F so that F is not zero when r is zero(what is U?)


def all_rates_abs():
    """ returns Itot, a function of rate and time"""
    lenrate = 45
    Itot = np.zeros((lenrate,len(time_p)))
    for rate_i, rate in enumerate(np.arange(0,lenrate)):
        rate = rate*np.ones(len(time))
        #print rate

        Iexc = J_exc*thalam_abs('exc', rate)[0]#[0]
        Iinh = J_inh*thalam_abs('inh', rate)[0]#[0]
        Itot[rate_i,:] = Iexc+Iinh
    fig, ax = py.subplots(figsize = [4,10])
    #ax.plot(ratetest, Itot[:,500])
    ax.set_xlim(0,30)
    plotticks(ax)
    return Itot
        

def thalam_abs(plas_type,rate_ext):
    S = np.zeros(lentime)
    S[0]=0
    D = np.zeros(lentime)
    F = np.zeros(lentime)
    F_ = np.zeros((2,lentime))
    D_ = np.zeros((2,lentime)	)			

    D[0]=0.7
    F[0]=0.1
    D_[0,0] = 0.7	
    F_[0] = 0.1				
    #F = fac_eq(rate_ext)
    #D = dep_eq(rate_ext)
				
    for t_index in np.arange(0, lentime-1):
        indexs = t_index
        index_rate = 0
        rate_ext_ = rate_ext[indexs] 							
        S_exc[index_rate,indexs + 1]= S_exc[index_rate, indexs] + stepsize_switch*s_eq('exc',(rate_ext_-rate_offset)*(rate_ext_>rate_offset),indexs*stepsize,S_exc[index_rate, indexs],F_[index_rate,indexs],D_[index_rate,indexs])
        S_inh[index_rate,indexs + 1]= S_inh[index_rate, indexs] + stepsize_switch*s_eq('inh',rate_ext_,indexs*stepsize,S_inh[index_rate, indexs],F_[index_rate,indexs],D_[index_rate,indexs])
        F_[index_rate,indexs+1]= F_[index_rate,indexs]+ stepsize_switch*fac_arg(F_[index_rate,indexs],indexs*stepsize,(rate_ext_-rate_offset)*(rate_ext_>rate_offset)/1000.0)
        D_[index_rate,indexs + 1]= D_[index_rate,indexs] + stepsize_switch*dep_arg(D_[index_rate,indexs],indexs*stepsize,rate_ext_/1000.0)
        S_rates[index_rate,indexs] = J_inh*S_inh[index_rate, indexs]+J_exc*S_exc[index_rate, indexs]#+J_inh*S_inh[index_rate, indexs]
                   				
					
					
					
        if plas_type == 'exc':
           tau = tau_E
           Seq = (-S[t_index]/tau+F[t_index]*(rate_ext[t_index]-rate_offset)*(rate_ext[t_index]>rate_offset)/1000.0)
           Feq = fac_arg(F[t_index],t_index,(rate_ext[t_index]-rate_offset)*(rate_ext[t_index]>rate_offset)/1000.0)
           
           S[t_index + 1]= S[t_index] + stepsize_switch*Seq 
           F[t_index+1]= F[t_index]+ stepsize_switch*Feq
        if plas_type == 'inh':
           tau = tau_I
           S[t_index + 1]= S[t_index] + stepsize_switch*(-S[t_index]/tau+p*D[t_index]*rate_ext[t_index]/1000.0+Iinh) 
           D[t_index + 1]= D[t_index] + stepsize_switch*dep_arg(D[t_index],t_index,rate_ext[t_index]/1000.0)
    
    return S, S_rates#,F,D, F_,D_





def rate_abs(rate1,rate2):
    #rate1 = rate1/1000.0
    #rate2 = rate2/1000.0				
    Iexc_L = J_exc*thalam_abs('exc', rate1)[0]
    Iinh_L = J_inh*thalam_abs('inh', rate1)[0]
    Iexc_R = J_exc*thalam_abs('exc',rate2)[0]
    Iinh_R = J_inh*thalam_abs('inh',rate2)[0]
    I_L = Iexc_L+Iinh_L
    I_R = Iexc_R+Iinh_R
    I_syn = I_L + I_R
    
    I_total = I_syn+Ibconf + Iext5(time, 'dm', amp_conf,decay_conf) + offset_conf
    return I_syn, FIcurve_2(I_total,300)#r,
    
   
def rate_analytical(rate1,rate2):
    rate1 = rate1#/1000.0
    rate2 = rate2#/1000.0
    tau_FF = tau_F/1000.0
    tau_DD = tau_D/1000.0
    tau_EE = tau_E/1000.0
    tau_II = tau_I/1000.0
    I_L =  J_exc*tau_EE*(tau_FF*a_F*(rate1**2)/(1+tau_FF*a_F*rate1)+Iexc)+J_inh*tau_II*(p*rate1/(1+tau_DD*p*rate1)+Iinh) 
    I_R =  J_exc*tau_EE*(tau_FF*a_F*(rate2**2)/(1+tau_FF*a_F*rate2)+Iexc)+J_inh*tau_II*(p*rate2/(1+tau_DD*p*rate2)+Iinh)  
    return I_L+I_R, FIcurve_2(I_L+I_R+Ib_p,300)

def rate_totalplot(rate1,rate2):  
    '''
    calculates a plot of the absolute value as a function of time given two
    input rates, 
    
    rate1: input rate 1
    rate2: input rate 2
    
    '''
    #currents = Iexc_L+Iinh_L
    rate1 = 2.5*(step(time_p,0)-step(time_p,t_motion))+rate1*step(time_p,t_motion)+180*Iext5(time_p, 'dm', 80,30)
    rate2 = 2.5*(step(time_p,0)-step(time_p,t_motion))+rate2*step(time_p,t_motion)+180*Iext5(time_p, 'dm', 80,30)
    r = rate_abs(rate1,rate2)[1]
    
    r_an = rate_analytical(rate1,rate2)[1]
    #r_high = rate_hightotal(rate1,rate2, 0.0012,8,0.03)
    
    fig = py.figure(figsize=[3,3])    
    ax = fig.add_subplot(1,1,1)    
    ax.plot(time_p, r,'k')
    ax.plot(time_p,r_an, 'g')
    #ax.plot(time_p,r_high, 'c')
    ax.plot(time_p, rate1,'b')
    ax.plot(time_p,rate2, 'r')
    #ax.set_ylim(-5,40)
    plotticks(ax)
    ax.set_ylim(0,10)
    #fig.savefig('blueredlarge.pdf')







    
def get_J(rate):
    tau_FF = tau_F/1000.0
    tau_DD = tau_D/1000.0
    tau_EE = tau_E/1000.0
    tau_II = tau_I/1000.0
    a1 = tau_FF*a_F
    a2 = tau_DD*p
    tau_ast = tau_II/tau_EE
    J_ast = p*tau_ast*(1+a1*rate)/(a1*rate+a1*a2*rate**2)
    return J_ast



def calc_ratep(rate1,rate2):
    rPall = np.zeros((len(coh1vec), int(trials_transition), lentime))
    for cohindex in np.arange(0,len(coh1vec)):
           print cohindex
           for trialindex in np.arange(0,trials_transition):
               rL1 = rate1[cohindex,trialindex,:] 
               rR1 = rate2[cohindex,trialindex,:] 
               rPall[cohindex,trialindex,:] = rate_abs(rL1,rR1)[1]## with excitatory and inhibitory plasticity
    return rPall






def rate_total(rate1,rate2,t_index):
    
    #rate1 = rate1*step(time_p,0)
    #rate2 = rate2*step(time_p,0)

    Iexc_L = J_exc*s_eq('exc', rate1,t_index)#[0]
    Iinh_L = J_inh*s_eq('inh', rate1,t_index)#[0]
    Iexc_R = J_exc*s_eq('exc',rate2,t_index)#[0]
    Iinh_R = J_inh*s_eq('inh',rate2,t_index)#[0]
    I_L = Iexc_L+Iinh_L
    I_R = Iexc_R+Iinh_R

    
    I_total = Iexc_L+Iinh_L + Iexc_R+Iinh_R+Ib_p#+1*Iext5(time,'dm',400)
    #r = 0#FIcurve_2(I_total)
    return (I_total)#r, I_L,Iexc_L,Iinh_L,I_R, I_total
    
    #rP = FIcurve(J_LP_eq*sL1+J_RP_eq*sR1+Ib*Imotion2(time)+Iext5(time,'dm',400),a)    



def steady_stateplot(plas, J_exc, J_inh,p):
    ''' plot of current as a function of input rate, that mimicks an absolute-value function
         plas: include plasticity, 'yes' or 'no' 
         J_exc: value of excitatory synaptic weight
         J_inh: value of inhibitory synaptic weight
         p: 
    
    '''
    ratetest = np.arange(0,31,0.01)
    ratetest = ratetest/1000.0
    offset1 = rate_offset/1000.0
    offset2 = rate_offset/1000.0
    if plas=='yes':
        s_inh = p*ratetest*tau_I/(1+tau_D*p*ratetest)
        s_exc = a_F*((ratetest-offset1)**2)*tau_F*tau_E/(1+a_F*(ratetest-offset2)*tau_F)*(ratetest>offset2)

    elif plas=='stf': 
        s_inh = p*ratetest*tau_I*D[0,0]
        s_exc = a_F*((ratetest-offset1)**2)*tau_F*tau_E/(1+a_F*(ratetest-offset2)*tau_F)*(ratetest>offset2)
    elif plas=='std':
        s_inh = p*ratetest*tau_I/(1+tau_D*p*ratetest)
        s_exc = 	tau_E*ratetest*F[0,0]
					
    fig = py.figure(figsize=[1,1.4])    
    ratetest = ratetest*1000
    ax = fig.add_subplot(1,1,1) 
    #ax.plot(ratetest,s_exc)
    #ax.plot(ratetest,-s_inh,'r')
    I_steady = J_exc*s_exc+J_inh*s_inh+Iexc				
    ax.plot(ratetest,I_steady,'k')
    ax.plot(ratetest,J_exc*s_exc, color = colors_set1['lightgreen'] )	
    ax.plot(ratetest,J_inh*s_inh, color = colors_set1['lightgreen'] )				
			
    ratelow = ratetest[0:1400]
    ratezero = ratetest[1400:1402]
    ratehigh = ratetest[1402:-1]
    #py.xlim(0,45)
    #ax.plot(ratezero, [-0.105,0], 'g')
    #ax.plot(ratelow, -0.008*ratelow, 'g')
    #ax.plot(ratehigh, 0.01*(ratehigh-13), 'g')
    ax.set_xlabel('cortical firing rate (Hz)')
    ax.set_ylabel('cortico-thalamic current (nA)')
    ax.set_xticks([0,15,30])#,45])
    ax.plot()
    ax.set_ylim(-0.08,0.12)
    ax.set_yticks([-0.05,0,0.05,0.1])
    zeroindex = (np.where(np.abs(I_steady)<0.0001))
    mean = np.max(ratetest[zeroindex])
    print mean
    plotticks(ax)
    fig.savefig('confidencefigs/steadystateabsolute_decay.pdf')
    return I_steady, ratetest 
        
    
def abs_comparison():  
    '''produces a plot of real absolute value vs approximated absolute value'''
    fig, ax = py.subplots(figsize = [1,1.2])
    offset1 = rate_offset/1000.0
    offset2 = rate_offset/1000.0
    rangenumber = 0.7			
    eps_mean = .1
    I_steady,ratetest = steady_stateplot('yes', J_exc, J_inh,p)    
    zeroindex = (np.where(np.abs(I_steady)<0.0001))
    mean = np.max(ratetest[zeroindex])
    allvalues = np.arange(mean-np.ceil(rangenumber*mean),mean+np.ceil(rangenumber*mean),0.1)/1000.0#np.random.normal(mean,4,80)/1000.0
    mean = mean/1000.0
    off = 0/1000.0
    realrate = np.zeros((len(allvalues), len(allvalues)))
    fakerate = np.zeros((len(allvalues), len(allvalues)))

    valuesg20 = allvalues#[allvalues>mean+off]
    valuesl20 = allvalues#[allvalues<mean-off]
    for indexa, valuesa in enumerate(valuesg20):
        for indexb, valuesb in enumerate(valuesl20):
            valuecondition = valuesa+valuesb-2*mean
            if valuecondition<eps_mean/1000.0:
                realabs = (valuesa-valuesb)
                #print realabs
                s_inh_a = p*valuesa*tau_I/(1+tau_D*p*valuesa)
                s_exc_a = a_F*((valuesa-offset1)**2)*tau_F*tau_E/(1+a_F*(valuesa-offset2)*tau_F)*(valuesa>offset2)
                s_inh_b = p*valuesb*tau_I/(1+tau_D*p*valuesb)
                s_exc_b = a_F*((valuesb-offset1)**2)*tau_F*tau_E/(1+a_F*(valuesb-offset2)*tau_F)*(valuesb>offset2)
           
                fakeabsrate = FIcurve_2(J_exc*(s_exc_a+s_exc_b)+J_inh*(s_inh_a+s_inh_b)+Ibconf,300)
                #fakeabscurrent = J_exc*(s_exc_a+s_exc_b)+J_inh*(s_inh_a+s_inh_b)
                ax.plot(realabs*1000,fakeabsrate, '.', markersize= 0.3,color = colors_set1['lightgreen'])
                realrate[indexa, indexb] = realabs
                fakerate[indexa,indexb] = fakeabsrate
    real1d = realrate.flatten()
    fake1d = fakerate.flatten()
    f = interp1d(1000*real1d, fake1d)  
    print f        
    xnew = np.arange(-2*rangenumber*(mean*1000),2*rangenumber*(mean*1000),0.1)
    faked1d_filt = filter_t(f(xnew),20)
    ax.plot(xnew, faked1d_filt, color = colors_set1['green'])
    
    plotticks(ax)
    ax.set_xticks([-20,-10,0,10,20])
    ax.set_ylim(-2,16)
    ax.set_yticks([0,5,10,15])
    ax.set_xlabel('real diff.(Hz)')
    ax.set_ylabel('estimated diff.(Hz)')

    fig.savefig('confidencefigs/abscomparison_decay.pdf')
    #ax.set_ylim(-0.2,0.2)
    return mean*1000#realrate, fakerate
    
    
    

'''def plot_plasticitycurrents(rate1,rate2):
    rate1 = rate1*step(time_p,0)
    Iexc = rate_total(rate1,rate2)[2]
    Iinh = rate_total(rate1,rate2)[3]
    I_R = rate_total(rate1,rate2)[4]
    fig = py.figure(figsize=[2,2])    
    ax = fig.add_subplot(1,1,1) 
    ax.plot(time_p, Iexc,'b')
    ax.plot(time_p, Iinh, 'r')
    ax.plot(time_p, Iexc + Iinh+I_R, 'k')
    return Iexc+Iinh+I_R
'''    



    
"""   
def rate_high(alpharate,offrate1,offrate2): 
    rate_vec = np.arange(0,40,1)# approximation to absolute value with high firing rate
    py.plot(rate_vec, alpharate*(rate_vec-offrate1)**1.5-offrate2)

    #py.plot(rate_vec, alpharate*(rate_vec-offrate))
def rate_hightotal(rate1,rate2,alpharate,offrate1,offrate2):
    I1 = alpharate*(rate1-offrate1)**1.5-offrate2
    I2 = alpharate*(rate2-offrate1)**1.5-offrate2
    return FIcurve_2(I1+I2+Ib_p,300)
"""
"""def plot_abs():
    tau_FF = tau_F/1000.0
    tau_DD = tau_D/1000.0
    tau_EE = tau_E/1000.0
    tau_II = tau_I/1000.0
    factor = 1000
    rate1_vec = np.arange(0,40,1)
    rate2_vec = np.arange(0,40,1)
    Inh = 0
    for rate1 in rate1_vec:
        for rate2 in rate2_vec:
            rate_diff = rate1-rate2
            I_L =  J_exc*tau_EE*tau_FF*a_F*(rate1**2)/(1+tau_FF*a_F*rate1)+J_inh*tau_II*(p*rate1/(1+tau_DD*p*rate1)+Iinh) 
            I_R =  J_exc*tau_EE*tau_FF*a_F*(rate2**2)/(1+tau_FF*a_F*rate2)+J_inh*tau_II*(p*rate2/(1+tau_DD*p*rate2)+Inh)  
            py.plot(rate_diff, np.abs(rate_diff), 'k.')
            py.plot(rate_diff, factor*(I_L+I_R), 'g.')
    
  """