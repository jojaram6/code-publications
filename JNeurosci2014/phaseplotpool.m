%%%%%%%
noiseinputflag=1;
%%%%%%%%
speed=40;
numberofrealizations=1;
spikeincrement=1;
tspikeallpool=[];
phasepool=[];
oldcorr1=0;
for number=1:numberofrealizations
    if noiseinputflag==1
        noiseinput
    else
    intrapoissoncurrent
    end
%     jitter=range/(8*360)*(2*rand(1)-1);
         number
%     conshift2=interp1(bin, convolution, tsmall-jitter);
%     conshift=conshift2/area;
    %mean(conshift)
    prca1newparameters
    %%%%%%%FIGURES OF MEMBRANE POTENTIAL AND INPUT CURRENTS%%%%%%%%%
    
      %figure
      hold all
      %plot(t,1000*VS+60)
     %hold on
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Plots for noisetocell.m (Sponantenous/Poisson random input%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     plot(t,1000*(12e-9/amplitudesoma*7e4*Isoma-0.064),'r')
%     plot(t,1e9*noisetocell-50, 'g')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    spikecount;
   
    spikesperrealization(number)=length(tspikeall);
    
    phase=360*ftheta*(tspikeall-phasesoma);
    if length(tspikeall)~=0
       
    for xspike=1:length(tspikeall)
        %test=56
        tspikeallpool(spikeincrement)=tspikeall(xspike);
        phasepool(spikeincrement)=phase(xspike);
        spikeincrement=spikeincrement+1;
    end
    end
%    distance=speed*tspikeallpool;
   
if tspikeall==37
    allspikes=37
    number
else
  allspikes=tspikeallpool;  
end 
    
    
    
    
    %tpeakperrealization(number)=tpeak;
   %clear tspikeall
   clear tspikenoise
   clear noisetocell
   clear epspmany
%    allspikes=tspikeallpool; 
%    allphases=mod(phasepool,360);
jittertotal(number)=jitter;
  
clear tspikeallpool % clear tspikeallpool for crosshistogram.mnoiseinput.m/intrapoisson.m
  %  clear phasepool


%clear soma
clear VS_inf
clear gS_Tot
clear VD_inf
clear gD_Tot
%VS=0
clear Iapp
clear tspikevectorr
clear tspikevector
clear tspikepool
%clear VS
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
jitter


  crosshistogram %Crosshistogram computes the correlation function between CA3 and CA1 spikes
  spikeincrement=1; %mset for  crosshistogram.m/intrapoisson.m
  newcorr1=oldcorr1+corr1;% summation of correlation 
  oldcorr1=newcorr1;% summation of correlation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correlation plot%%%%%%%%%%
  figure
  plot(lags*deltabin, newcorr1/Nneuron,'r')% plot of sum/pooled correlations over CA1 neurons

  
  %  allphases=mod(phasepool,360);
 %figure
 %plot(distance/max(distance),allphases, '.')
 
 %hold on
 
  %hist(mod(phasepool,360),100)
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 %histogram of spike phases
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  figure 
%  hist(allphases,40)
%  hold on
%  hist(allphases+360,40)
%  figure
%  hist(allspikes, 40)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure
%  axis([0 1 0 3])
%   for indexrandominput=2:length(randomCA3cellspikes)
%      line([randomCA3cellspikes(indexrandominput),randomCA3cellspikes(indexrandominput)],[1.6,2])
%   end
%   hold on
%   for indexallspikes=1:length(allspikes)
%      line([allspikes(indexallspikes),allspikes(indexallspikes)],[0.6,1])
%  end
 

 
 
 
 
 %plot(distance/max(distance), realphase+360, '.')
% figure
 %plot(t,VS)
 %plot(distance, realphase-360, '.')
%%%%%%%%%%%
%%circular phase precession correlation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%extract maxslope



% realphase=allphases;
% phi=realphase/360*2*pi;
% slope=-10:0.01:0;
% for slopenumber=1:length(slope)
% firstterm(slopenumber)=1/length(realphase).*sum(cos(phi-2*pi*slope(slopenumber)*distance));
% secondterm(slopenumber)=1/length(realphase).*sum(sin(phi-2*pi*slope(slopenumber)*distance));
% R(slopenumber)=sqrt(firstterm(slopenumber)^2+secondterm(slopenumber)^2);
% end
% maxslope=slope(find(R==max(R)));
% 
% 
% 
% theta=mod(2*pi*abs(maxslope)*distance, 2*pi);
% phibar=angle(sum(exp(j*phi)));
% thetabar=angle(sum(exp(j*theta)));
% rhoc=sum(sin(phi-phibar).*sin(theta-thetabar))./sqrt(sum(sin(phi-phibar).^2)*sum(sin(theta-thetabar).^2));
% S=sum(sin(phi-2*pi*maxslope*distance));
% C=sum(cos(phi-2*pi*maxslope*distance));
% offset=atan(S/C);
% if C<0
%     realoffset=offset+pi;
% end
% if C>0&&S<0
%     realoffset=offset+2*pi;
% end




 
 