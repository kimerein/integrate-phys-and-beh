function pelletTouchAlignedReaching(alltbt,out,metadata)

% note that pellet is only present after cue
% withinRange=[3 6];  % take reaches in this range, time is in seconds from onset of cue
withinRange=[0.25 6];  % take reaches in this range, time is in seconds from onset of cue
nIndsToTake=200;

timeStep=mode(diff(nanmean(alltbt.times,1)));
cueInd=find(nanmean(alltbt.cueZone_onVoff,1)>0,1,'first');

withinRange_inds=[cueInd+ceil(withinRange(1)/timeStep) cueInd+ceil(withinRange(2)/timeStep)];

% temp=alltbt.reachBatch_success_reachStarts+alltbt.reachBatch_drop_reachStarts;
temp=alltbt.reachBatch_drop_reachStarts;
touchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);

% temp=alltbt.reachBatch_miss_reachStarts+alltbt.pelletmissingreach_reachStarts;
temp=alltbt.reachBatch_miss_reachStarts;
reachNoTouchAligned=alignToEvent(temp,withinRange_inds,alltbt.all_reachBatch,nIndsToTake);

figure(); 
plot(0:timeStep:timeStep*(nIndsToTake-1),nanmean(touchAligned,1),'Color','k');
hold on;
plot(0:timeStep:timeStep*(nIndsToTake-1),nanmean(reachNoTouchAligned,1),'Color','m');
xlabel('Time (seconds)');
ylabel('Reaching');
legend({'Aligned to touch of pellet','Aligned to reach but no touch of pellet'});

figure();
plot(0:timeStep:timeStep*(nIndsToTake-1),nanmean(touchAligned,1)./nanmean(reachNoTouchAligned,1),'Color','b');
xlabel('Time (seconds)');
ylabel('Fractional suppression of reaching after pellet touch');

figure();
plot(0:timeStep:timeStep*(nIndsToTake-1),nanmean(touchAligned,1)-nanmean(reachNoTouchAligned,1),'Color','b');
xlabel('Time (seconds)');
ylabel('Difference in reaching after pellet touch');
[mi,nm]=nanmin(nanmean(touchAligned,1)./nanmean(reachNoTouchAligned,1));

% figure();
% plot(0:timeStep:timeStep*(nIndsToTake-1),nanmean(touchAligned,1)./nanmean(reachNoTouchAligned,1),'Color','b');
% hold on;
% x=0:0.001:10;
% mu=0.5;
% sigma=1;
% y=1-lognpdf(x,mu,sigma);
% y=y-nanmin(y);
% y=y./nanmax(y);
% plot(x,y,'Color','r');
% xlabel('Time (seconds)');
% ylabel('Fractional suppression of reaching after pellet touch -- lognormal fit');

[mi,nm]=nanmin(nanmean(touchAligned,1)-nanmean(reachNoTouchAligned,1));
figure();
plot(0:timeStep:timeStep*(nIndsToTake-nm),nanmean(touchAligned(:,nm:end),1)-nanmean(reachNoTouchAligned(:,nm:end),1),'Color','b');
hold on;
x=0:0.001:7;
tau=0.5;
tau2=5;
a=0.7;
b=0.3;
y=-nanmax(abs(nanmean(touchAligned,1)-nanmean(reachNoTouchAligned,1)))*(a*exp(-x/tau)+b*exp(-x/tau2));
plot(x,y,'Color','r');
xlabel('Time (seconds)');
ylabel('Difference in reaching after pellet touch -- exponential fit');

% [mi,nm]=nanmin(nanmean(touchAligned,1)./nanmean(reachNoTouchAligned,1));
% figure();
% plot(0:timeStep:timeStep*(nIndsToTake-nm),nanmean(touchAligned(:,nm:end),1)./nanmean(reachNoTouchAligned(:,nm:end),1)-mi,'Color','b');
% hold on;
% x=0:0.001:10;
% tau=6.5;
% y=1-exp(-x/tau);
% plot(x,y,'Color','r');
% xlabel('Time (seconds)');
% ylabel('Fractional suppression of reaching after pellet touch -- plot from minimum, exponential fit');

end

function alignedOut=alignToEvent(event,withinRange_inds,whatToAlign,howMuchToTake)

% for each trial, get first pellet touch in this range
% align to this pellet touch
% either reachBatch_success_reachStarts or reachBatch_drop_reachStarts

event=event(:,withinRange_inds(1):withinRange_inds(2));

eventInds=nan(1,size(event,1));
for i=1:size(event,1)
    if any(event(i,:)>0.5,2)
        f=find(event(i,:)>0.5,1,'first');
        eventInds(i)=f;
    end
end
eventInds=eventInds+withinRange_inds(1)-1;

j=1;
alignedOut=nan(sum(~isnan(eventInds)),howMuchToTake);
for i=1:size(event,1)
    if ~isnan(eventInds(i))
        alignedOut(j,:)=whatToAlign(i,eventInds(i):eventInds(i)+howMuchToTake-1);
        j=j+1;
    end
end

end


