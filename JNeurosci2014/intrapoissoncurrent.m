%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%poisson process, spike generation
%%% GENERATION OF EPSCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%theta wave components%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%clear
clear t
clear n
clear m
clear i
clear tsmall

t=-.8:0.005:1.4;
tsmall=0:0.00005:1;
the=1*1e-3; %theta amplitude
ftheta=8;%theta frequency
phasetheta=0; %phase of theta 
theta=the*sin(2*pi*ftheta*t+phasetheta);% theta
A=10;%%DC amplitude of oscillation
B=5;%%% AC amplitude of oscillation
Pmax=2*1e-11;
tscurrent=0.006;
foscillation=9;
period=1/foscillation;
Nneuron=200;
rmax=1.1*(A+B);%maximum firing rate
xrand1=rand(1);
tspike(1)=0;
bin=t;%0:0.001:1.1;
binperiod=0:0.007:period;
oscillation=A+B*cos(2*pi*foscillation*t);
m=0;
p=1;
tspike(1)=0;
amplambda=1; %amplitude of firing rate
placewidth=0.47;
delaylambda=0.5;
%spiketime=spikeposition*0.68;%changing normalized position to time in phaseplot
for numberneurons=1:Nneuron
    numberneurons;
    for n=1:150
        tspike(n+1)=tspike(n)-log(rand(1))/rmax;
        %line([tspike(n),tspike(n)],[0.45,0.85]);
        %hold on
        tspikepool(numberneurons,n)=tspike(n);
         %if tspike(n)>0.7
          % break
         %end
    end
    for n=2:length(tspike)-1
        %fin=find((bin-tspike(n)).^2==min((bin-tspike(n)).^2));
         %fin=find((spiketime-tspike(n)).^2==min((spiketime-tspike(n)).^2));
        %finn(n)=fin;
         %
        if  amplambda*exp(-power((tspike(n)-delaylambda)/placewidth,2))*(A+B*cos(2*pi*foscillation*tspike(n)))/rmax>rand(1);% firing rate, gaussian
             
            %((1+sign(tspike(n)))/2-(1+sign(tspike(n)-1))/2)
            %
            %
            %((spikecount(fin))/4)/rmax>rand(1)  
            
            
             %%%%%%expression for firing rate
             
            
            total(p)=tspike(n);
            %line([total(p),total(p)],[1.6,2]);
            totalpool(numberneurons,p)=total(p);
             
         %if total(p)>0.69
          
         %break
        %end
        else
            m=m+1;
        end

        
        p=n+1-m;
       
    end
    
    m=0;
    p=1;


end

      %      
randomcell=ceil(Nneuron.*rand(1,1));      %%% generation of a random number from the Ninputs
%randomCA3cellspikes=totalpool(randomcell, :); 
pooled=totalpool;
      
%figure;
%hold on
%axis([-0.85 1.4 0 12]);
%plot(t,0.1*fire,'r');
%plot(spikeposition*0.68, spikecount/15, 'r');
%hold off
%axis([-0.1 1 0 12]);
%%%%%%%%%%%%%%%%%%
%%%histogram%%%%%%
%%%%%%%%%%%%%%%%%%
phasebin=mod(2*pi*foscillation*bin,period);
phasetotal=mod(totalpool,period);
mm=0;
xx=size(totalpool)
% for nn=1:length(binperiod)-1
%     for ii=1:Nneuron
%         for mo=1:xx(2)
%         if phasetotal(ii,mo)>binperiod(nn) & phasetotal(ii,mo)<binperiod(nn+1)
%             mm=mm+1;
% 
%         end
%     end
%     
%     end
%     spikecount(nn)=mm;
%         mm=0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%convolution of spike train with EPSP%%%%%%%%%%%%%%%%%%%
%epsp=Pmax/ts*(bin).*exp(1-bin)/ts;
tspikevectorr=totalpool(:);
o=1;
for i=1:length(tspikevectorr)
    if tspikevectorr(i)~=0
        tspikevector(o)=tspikevectorr(i);
        o=o+1;
    else
        o=o;
    end

end
o=1;
pp=1;
for ff=1:length(bin)
    for g=1:length(tspikevector)
        %lines of CA3 spikes
%         axis([0 1 0 3])
%         line([tspikevector(g),tspikevector(g)], [1.6,2])
%         hold on
        
        if tspikevector(g)~=0
           argument=bin(ff)-tspikevector(g);
           arg(g,ff)=argument;
          epsc(g,ff)=Pmax/tscurrent*(argument).*exp(1-argument/tscurrent)*(1+sign(argument))/2;
        end
           end

end
convolutions=sum(epsc);
convolution=convolutions;
range=180;%%range of jitter is +- this value in degrees
jitter=range/(8*360)*(2*rand(1)-1);%jitter is in time
jitter2=range/(8*360)*(2*rand(1)-1);
conshift2=interp1(t, convolution, tsmall-jitter);
convolutionreduced1=interp1(t,convolution, tsmall);
area= 0.0415; %%%%% in mm^2
convolutionreduced=convolutionreduced1/area;
conshift=conshift2/area;
mean(convolutionreduced);

 figure
 plot(tsmall,convolutionreduced*10e8,'g')
 hold on
% plot(t, amplambda*exp(-power((t-delaylambda)/placewidth,2)).*(A+B*cos(2*pi*foscillation*t)))

clear arg
clear g
clear tspikevectorr
clear tspikevector
clear tspikepool
clear total
clear totalpool
clear convolutions
