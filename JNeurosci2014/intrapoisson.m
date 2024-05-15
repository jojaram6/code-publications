%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%poisson process, spike generation%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%theta wave components%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%clear
baseline=-68;
threshold=-55;

t=-0.8:0.005:1.4;
the=3*1e-3; %theta amplitude
ftheta=8;%theta frequency
phasetheta=0; %phase of theta 
theta=the*sin(2*pi*ftheta*t+phasetheta);% theta
A=20;%%DC amplitude of oscillation
B=18;%%% AC amplitude of oscillation
c=(A-B)/(A+B);
Pmax=1.5*1e-4;
ts=0.01;
foscillation=9;
period=1/foscillation;
Nneuron=200;

rmax=1.1*(A+B);%maximum firing rate
%rmax=20;%only for phaseplot.m


xrand1=rand(1);
tspike(1)=0;
bin=t;%0:0.001:1.1;
binperiod=0:0.007:period;
oscillation=A+B*cos(2*pi*foscillation*t);
m=0;
p=1;
tspike(1)=-0.81;
amplambda=1; %amplitude of firing rate
placewidth=0.4;
delaylambda=0.5;
%spiketime=spikeposition;%changing normalized position to time in phaseplot
for numberneurons=1:Nneuron
    %numberneurons
    for n=1:150
        tspike(n+1)=tspike(n)-log(rand(1))/rmax;
        %line([tspike(n),tspike(n)],[5,7]);
        %hold on
%         if tspike(n+1)>1
%             break
%          end
        tspikepool(numberneurons,n)=tspike(n);
         
        
         
         
    end
    for n=2:length(tspike)-1
        %fin=find((bin-tspike(n)).^2==min((bin-tspike(n)).^2));
         %fin=find((spiketime-tspike(n)).^2==min((spiketime-tspike(n)).^2));
        %finn(n)=fin;
         %
%         if  ((realsmooth(fin)))/rmax>rand(1)
            if amplambda*exp(-power((tspike(n)-delaylambda)/placewidth,2))*(A+B*cos(2*pi*foscillation*tspike(n)))/rmax>rand(1);% firing rate, gaussian
            
            % 
            %((1+sign(tspike(n)))/2-(1+sign(tspike(n)-1))/2)
            %
            %
              
            
            
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
pp=1;
for f=1:length(bin)
    for g=1:length(tspikevector)
        if tspikevector(g)~=0
           argument=bin(f)-tspikevector(g);
           arg(g,f)=argument;
          epsp(g,f)=Pmax/ts*(argument).*exp(1-argument/ts)*(1+sign(argument))/2;
        end
        
end
 end
convolution=sum(epsp);
convolutioninplacefield= convolution(161:361);
convolutioninplacefieldpluszeros=[zeros(size([-0.8:0.005:0])) convolutioninplacefield zeros(size([1:0.005:1.4-0.01]))];
cpz=convolutioninplacefieldpluszeros;
poissonmembrane=(cpz+theta)*1000+baseline;
figure;
plot(bin,poissonmembrane);% multiply times 1000 to get millivolts!
%title (['normalized Membrane potential for N=', num2str(Nneuron), ' A=', num2str(A), ' B=', num2str(B)])
P=findpeaks(bin,poissonmembrane,10, threshold,10,10,1);
spiketimes=P(:,2);

xlabel('time (s)')
ylabel('potential (mV)')
box off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%amplitude of convolution at t=0 is larger. try either converting
%%the sparse matrix and elminating 0s, or using a break inthe for loop to
%%eliminate the 0s
deltavmax=20*1e-3;
rho=B*sqrt((ts*Nneuron)/A);
SNRmin=SNR*0.7;
low=(SNRmin^2*A)/(8*ts*B^2);
high=deltavmax/(ts*Pmax*(A+B)*exp(1));

            


