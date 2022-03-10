function alltbt=checkForOptoEnhancedReach(alltbt,metadata,trialTypes,reachType,ledOnCond,cueFieldName,aroundCueWindow,percentGreaterCutoff)

u=unique(metadata.sessid);
ledIsOn=eval(ledOnCond);
[~,fi]=nanmax(nanmean(alltbt.(cueFieldName),1));
times=nanmean(alltbt.times,1);
startInd=fi+floor(aroundCueWindow(1)/mode(diff(nanmean(alltbt.times,1))));
endInd=fi+floor(aroundCueWindow(2)/mode(diff(nanmean(alltbt.times,1))));
optoenhance=nan(1,length(u));
alltbt.opto_enhanced_reach=nan(size(alltbt.sessid));
for i=1:length(u)
    currsess=u(i);
    temp=alltbt.(reachType);
    sub_ledIsOn=ledIsOn(metadata.sessid==currsess);
    sub_temp=temp(metadata.sessid==currsess,:);
%     figure();
%     plot(times,nanmean(sub_temp(sub_ledIsOn==1,:),1),'Color','r'); hold on;
%     plot(times,nanmean(sub_temp(sub_ledIsOn==0,:),1),'Color','k');
    optoOnReaching=nanmean(nanmean(sub_temp(sub_ledIsOn==1,startInd:endInd),1),2);
    optoOffReaching=nanmean(nanmean(sub_temp(sub_ledIsOn==0,startInd:endInd),1),2);
    if optoOnReaching>optoOffReaching+(percentGreaterCutoff/100)*optoOffReaching
        % opto-enhanced reach
        optoenhance(i)=true;
    else
        optoenhance(i)=false;
    end
    alltbt.opto_enhanced_reach(alltbt.sessid==currsess)=optoenhance(i);
end