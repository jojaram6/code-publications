%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% count number of spikes%%%%%%
Vth=-0.01;
deltatmin=1.5e-3;
tlist=t(find(VS>Vth));%%%%%%%%%% tlist is the list of times which contain the spike times somewhere in the middle
if tlist~=0
    tmin=tlist(1);
    i=0;
    for listindex=1:length(tlist)-1
        if tlist(listindex+1) - tlist(listindex)>deltatmin
            tmax=tlist(listindex);
            timeofspike=(tmin+tmax)/2;
            tmin=tlist(listindex+1);
            i=i+1;
            tspikeall(i)=timeofspike;
        end
        tspikeall(i+1)=(tlist(length(tlist))+tmin)/2;
    end
else
     
      %tspikeallpool=[];
%      phasepool=[];
end
