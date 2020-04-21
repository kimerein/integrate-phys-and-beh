function [dprimes,hit_rates,FA_rates,out]=get_dprime_per_session(tbt,out,metadata,whichReach,nameOfCue,settings)

% signal in this case is cue
% response if mouse performed a "cued reach" -- definition in
% RTanalysis_settings.m
% so hit is when mouse reached to cue
% miss is when mouse failed to reach to cue
% false alarm is when mouse performed uncued reach (take a time window
% before cue, for example)
% and correct rejection is when mouse did not reach without cue (for
% example, in pre-cue window)

if isempty(settings)
    settings=reachExpt_analysis_settings();
end
hitWindow_start=settings.reachAfterCueWindow_start; % wrt cue onset
hitWindow_end=settings.reachAfterCueWindow_end; % wrt cue onset
FAWindow_start=settings.preCueWindow_start; % wrt trial onset
FAWindow_end=settings.preCueWindow_end; % wrt trial onset

% if settings.excludePawOnWheelDuringCue
%     f=fieldnames(tbt);
%     paw_during_wheel=out.paw_during_wheel;
%     for i=1:length(f)
%         temp=tbt.(f{i});
%         tbt.(f{i})=temp(paw_during_wheel==0,:);
%     end
%     f=fieldnames(metadata);
%     for i=1:length(f)
%         temp=metadata.(f{i});
%         if isnumeric(temp)
%             metadata.(f{i})=temp(paw_during_wheel==0,:);
%         else
%             % is cell
%             metadata.(f{i})=temp(paw_during_wheel==0);
%         end
%     end
%     f=fieldnames(out);
%     for i=1:length(f)
%         temp=out.(f{i});
%         out.(f{i})=temp(paw_during_wheel==0);
%     end
% end

% calculate hit rates per session

% Convert time window wrt cue onset into indices into data
cueInd=find(nanmean(tbt.(nameOfCue),1)>settings.lowThresh,1,'first');
startInds=floor(abs(hitWindow_start)/mode(diff(nanmean(tbt.times,1))));
if hitWindow_start<0
    startInds=-startInds;
end
endInds=floor(abs(hitWindow_end)/mode(diff(nanmean(tbt.times,1))));
if hitWindow_end<0
    endInds=-endInds;
end
startInds=cueInd+startInds;
endInds=cueInd+endInds;
if startInds<1
    startInds=1;
end
if endInds>size(tbt.times,2)
    endInds=size(tbt.times,2);
end

temp=tbt.(whichReach);
hits=any(temp(:,startInds:endInds),2);
% throw out all nan trials
allnantrials=all(isnan(tbt.times),2);
hit_rates=trials_per_session(metadata,hits==1,allnantrials==0);

% calculate false alarm rate

% Convert time window wrt trial onset into indices into data
if ~iscell(FAWindow_start)
    useInds=zeros(1,size(temp,2));
    startInds=floor(FAWindow_start/mode(diff(nanmean(tbt.times,1))));
    endInds=floor(FAWindow_end/mode(diff(nanmean(tbt.times,1))));
    if startInds<1
        startInds=1;
    end
    if endInds>size(tbt.times,2)
        endInds=size(tbt.times,2);
    end
    useInds(startInds:endInds)=1;
else
    useInds=zeros(1,size(temp,2));
    for i=1:length(FAWindow_start)
        currStretch=FAWindow_start{i};
        startInds=floor(currStretch(1)/mode(diff(nanmean(tbt.times,1))));
        endInds=floor(currStretch(2)/mode(diff(nanmean(tbt.times,1))));
        if startInds<1
            startInds=1;
        end
        if endInds>size(tbt.times,2)
            endInds=size(tbt.times,2);
        end
        useInds(startInds:endInds)=1;
    end
end

FAs=any(temp(:,useInds==1),2);
FA_rates=trials_per_session(metadata,FAs==1,allnantrials==0);

dprimes=dprime(hit_rates,FA_rates);

% dprime for session, organize per trial
out.dprime=nan(size(out.led,1),1);
for i=1:size(out.led,1)
    out.dprime(i)=dprimes(metadata.sessid(i));
end

end

function fractionTrialsInGroup=trials_per_session(metadata,in_group_trials,total_trials)

[sessid,sessStartInd]=unique(metadata.sessid);
fractionTrialsInGroup=nan(1,length(sessid));
sessStartInd=[sessStartInd; length(metadata.sessid)+1];
for i=1:length(sessid)
    curr_sessStartInd=sessStartInd(i);
    sessinds=curr_sessStartInd:sessStartInd(i+1)-1;
    fractionTrialsInGroup(i)=sum(in_group_trials(sessinds)==1)/sum(total_trials(sessinds)==1);
end

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end