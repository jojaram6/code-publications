clear VS
clear VD
clear IS
clear ID
clear m
clear mslow
clear mpersistent
clear n
clear h
clear IappS
clear IappD
clear Iapp
clear t
clear t

% PR.m
% This model is the two-compartment model of Pinsky and Rinzel (1994)
dt = 0.00005;
tmax=1;

iclamp_flag = 1; % if this is 1, run under current clamp conditions
vclamp_flag = 0; % otherwise this should be 1, for voltage clamp conditions

idendrite_flag =1; %if this is one, apply current in dendrite, otherwise in soma
IappS = 0.0;
IappD = 0.0;

isine_flag = 0;  % if this is one send in an oscillatory current
freq = 9;        % frequency of current oscillations

%%%%%%%%%%%%%%%%%%%%%%%
%%%SOMA input THETA%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

jitter=0;
initialphaseinsoma=0;
phasesoma=jitter+initialphaseinsoma;
%phasesoma=0;
dcsoma=0;
amplitudesoma=20e-9;
soma=amplitudesoma*cos(2*pi*8*(tsmall-phasesoma))+dcsoma;
soma2=zeros(size(tsmall));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



istart = 0.0; % time applied current starts
ilength=1;   % length of applied current pulse
Ie=0e-9;     % magnitude of applied current pulse
%Ie = 0.78e-7; % threshold for constant spiking with no A-current

vstart = 0.25;  % time to step voltage
vlength = 0.5;  % length of voltage step
V0 = -0.080;    % initial voltage before step
Ve = -0.000;    % value of votage stepped to

E_L = -0.065;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.090;   % reversal for potassium channels

g_LD = 1.8e-6;     % specific leak conductance in Siemens per mm-square
g_LS = 1.8e-6;     % specific leak conductance in Siemens per mm-square
g_Na = 0.55e-3;   % specific sodium conductance
g_K = 0.28e-3;    % specific potassium conductance
g_Ks=5.4e-5; %kamondi values: 1.4e-5;      % slow potassium
g_Nap=5e-7; %kamondi value 5e-7      % persistent sodium



g_Link = 1e-5; % kamondi value 1 conductance linking dendrite and soma divided by total membrane area
%g_Link = 0e-6; % conductance linking dendrite and soma divided by total membrane area


S_frac = 0.15;  % fraction of total membrane area that is soma
D_frac = 1-S_frac; % rest of area is dendritic

g_S_Link = g_Link/S_frac; %link conductance divided by somatic area
g_D_Link = g_Link/D_frac; % link conductance divided by dendritic area

cm = 10e-9;     % specific membrane capacitance in Farads per mm-square

t=0:dt:tmax;        % time vector
VS=zeros(size(t));  % somatic voltage vector
VD=zeros(size(t));  % dendritic voltage vector
VD(1) = E_L;
VS(1) + E_L;



I_LD= zeros(size(t));   % leak current in dendrite
I_LS= zeros(size(t));   % leak current in soma
I_Na = zeros(size(t));
I_K = zeros(size(t));
I_Ks = zeros(size(t));
I_Nap = zeros(size(t));

if ( iclamp_flag ) % i.e. if in current-clamp mode 
    VS(1) = E_L;    % set the inititial value of somatic voltage     
    VD(1) = E_L;    % set the inititial value of dendritic voltage     
end

n=zeros(size(t));   % n: potassium activation gating variable
m=zeros(size(t));   % m: sodium activation gating variable
h=zeros(size(t));   % h: sodim inactivation gating variable
mpersistent=(size(t));
mslow=zeros(size(t)); % CaT current activation gating variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************************************
%****EPSCS DUE TO INHOMOGENOUS POISSON PROCESS***************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I0=15e-9;
normal=10;
testpulse=2.2e-8*(sign(t-0.1)+1)/2;

if noiseinputflag==1
    Iapp=noisetocell;
else
    Iapp=conshift; %%%%%%%% it should be conshift;
end

if normal==1
Iapp = zeros(size(t));
if ( iclamp_flag )   % i.e. if in current-clamp mode 
    for i = 1:round(istart/dt)
        Iapp(i) = I0;
    end
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        if ( isine_flag ) % if the pulse is oscillatory
            Iapp(i) = Ie - (Ie-I0)*cos((t(i)-istart)*freq*2*pi);
            
           
        else
            Iapp(i) = Ie;
        end
    end
    for i = round((istart+ilength)/dt):length(Iapp)
        Iapp(i) = I0;
    end
end
end




Vapp=zeros(size(t)); % Applied voltage, relevant in voltage-clamp mode
if ( vclamp_flag )
    for i = 1:round(vstart/dt)          % % make V0 before pulse
        Vapp(i) = V0;
    end
    for i=round(vstart/dt)+1:round((vstart+vlength)/dt) % make Ve for duration of voltage pulse
        Vapp(i) = Ve;
    end
    for i=round((vstart+vlength)/dt):length(Vapp) % make V0 following pulse
        Vapp(i) = V0;
    end

    Vapp = V0+(Ve-V0)*t/tmax;
end

Itot=zeros(size(t)); % in case we want to plot and look at the total current

for i = 2:length(t); % now see how things change through time
    I_LS(i) = g_LS*(E_L-VS(i-1));
    I_LD(i) = g_LD*(E_L-VD(i-1));
   
    Vm = VS(i-1)*1000; % converts voltages to mV
    VmD = VD(i-1)*1000; % converts voltages to mV
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. 
    
    %%%%%%%%m for sodium currnet I_Na
    if ( Vm == -31 ) 
        alpha_m = 1;
    else 
        alpha_m = 0.1*(Vm+31)/(1-exp(-0.1*(Vm+31)));
    end

    beta_m = 4*exp(-(Vm+56)/18);
    
    alpha_h = 0.07*exp(-(Vm+47)/20);
    beta_h = 1/(1+exp(-0.1*(Vm+17)));
    

    if ( Vm == -34 ) 
       alpha_n = 0.01/0.1;
    else
        alpha_n = 0.01*(Vm+34)/(1-exp(-0.1*(Vm+34)));
    end
    beta_n = 0.125*exp(-(Vm+44)/80);
     
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1e-3/(alpha_m+beta_m);      % time constant converted from ms to sec
    m_inf = alpha_m/(alpha_m+beta_m);
 
    if( m_inf < 0 ) 
        m_inf = 0;
    end
    if ( m_inf > 1 ) 
        m_inf = 1;
    end
        
    tau_h = 1e-3/(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h+beta_h);
    if ( h_inf < 0 ) 
        h_inf = 0;
    end
    if ( h_inf > 1 )
        h_inf = 1;
    end
    
    tau_n = 1e-3/(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n+beta_n);   
    if ( n_inf < 0 ) 
        n_inf = 0;
    end
    if ( n_inf > 1 ) 
        n_inf = 1;
    end
    
    
    m_infp = 1/(1+exp(-(VmD+57.7)/7.7));
    m_infs = 1/(1+exp(-(VmD+35)/6.5));
    tau_ms = 200/(exp(-(VmD+55)/30)+exp((Vm+55)/30));
    testtau(i)=tau_ms;
    m(i) = m_inf - (m_inf-m(i-1))*exp(-3.33*dt/tau_m);    
    
    mpersistent(i) = m_infp;  
    
    mslow(i) = m_infs - (m_infs-mslow(i-1))*exp(-3.33*dt/tau_ms);  

    
    h(i) = h_inf - (h_inf-h(i-1))*exp(-3.33*dt/tau_h);    % Update h
    
    n(i) = n_inf - (n_inf-n(i-1))*exp(-3.33*dt/tau_n);    % Update n
    
     
    g_Na_now = g_Na*m(i)*m(i)*m(i)*h(i);
    I_Na(i) = g_Na_now*(E_Na-VS(i-1)); % sodium current in soma
    
    g_K_now = g_K*n(i)*n(i)*n(i)*n(i);
    I_K(i) = g_K_now*(E_K-VS(i-1)); % potassium delayed rectifier current, soma
    
    g_Nap_now= g_Nap*m_infp*m_infp*m_infp;
    I_Nap(i) = g_Nap_now*(E_Na-VS(i-1));
    
    g_Ks_now = g_Ks*mslow(i);
    I_Ks(i) = g_Ks_now*(E_K-VS(i-1));
    
   
    
    I_Link(i) = g_Link*(VD(i-1)-VS(i-1));
    
    if ( idendrite_flag==1 ) 
        IappD = Iapp(i);
    else
        IappS = Iapp(i);
    end
    
    IappS=soma(i);
    %IappS=0;
    IS(i) = I_LS(i)+I_Na(i)+I_K(i)+I_Link(i)/S_frac+IappS; % total current in soma
    testcurrent(i)=IappS;
    testIS(i)=I_LS(i)+I_Na(i)+I_K(i)+I_Link(i)/S_frac+testcurrent(i);
    ID(i) = I_LD(i)+I_Nap(i)+I_Ks(i)-I_Link(i)/D_frac +IappD; % total current in dendrite
    
    gS_Tot = g_LS+g_Na_now+g_K_now+g_S_Link;
    VS_inf = (g_LS*E_L + g_Na_now*E_Na + g_K_now*E_K ...
            + VD(i-1)*g_S_Link + IappS)/gS_Tot;
                   
    gD_Tot = g_LD+ g_Nap_now+ g_Ks_now+ g_D_Link;
    VD_inf = (g_LD*E_L + g_Nap_now*E_Na + g_Ks_now*E_K + ...
            + VS(i-1)*g_D_Link + IappD )/gD_Tot;
                   
    VS(i) = VS_inf - (VS_inf-VS(i-1))*exp(-dt*gS_Tot/cm);  % Update the membrane potential, V.
    VD(i) = VD_inf - (VD_inf-VD(i-1))*exp(-dt*gD_Tot/cm);  % Update the membrane potential, V.
   
    
    if ( vclamp_flag )      % if we are using voltage clamp
        VS(i) = Vapp(i);     % ignore the voltage integration and set V to be the applied voltage
        VD(i) = Vapp(i);     % ignore the voltage integration and set V to be the applied voltage
    end
   
    Idendrite(i)=IappD;
    Isoma(i)=IappS;
end
%figure;
%plot(t,1000*(VS),'r')
%clear Iapp

%clear convolutionreduced
%figure
% plot(t,testpulse)
% hold all
% plot(t,conshift)

 