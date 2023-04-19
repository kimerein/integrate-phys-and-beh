function [isreaching_out,dprimes]=getCuedResponseVsSuppression(alltbt,metadata,out,nameOfCue,reachName,useTrials,firstSess,reachAfterCueWindow_start,reachAfterCueWindow_end,doPlot,doRawReachRates)

% this code assumes that all sessions are from the same mouse
if isempty(doPlot)
    doPlot=1;
end
doQuiver=1;
settingsForDp=settingsForDprimes(alltbt,nameOfCue,false);
preCueWindow_start1=settingsForDp.preCueWindow_start1; % define start of time window from trial onset, in seconds -- for first window, assuming that trial onset is 0 sec
preCueWindow_end1=settingsForDp.preCueWindow_end1; % define end of time window from trial onset, in seconds -- for first window
% preCueWindow_start2=3.81; % define start of time window from trial onset, in seconds -- for second window
% preCueWindow_end2=5.31; % define end of time window from trial onset, in seconds -- for second window
% [~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
% cuetimeat=mode(diff(nanmean(alltbt.times,1)))*ma;
preCueWindow_start2=settingsForDp.preCueWindow_start2; 
preCueWindow_end2=settingsForDp.preCueWindow_end2;

% doRawReachRates=0;
mean_cued_reach_rate=[];
max_non_cued=[]; 
if doRawReachRates==false
    
    % useTrials only
    if ~isempty(useTrials)
        f=fieldnames(alltbt);
        for i=1:length(f)
            temp=alltbt.(f{i});
            if size(temp,1)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1,:);
            alltbt.(f{i})=temp;
        end
        f=fieldnames(out);
        for i=1:length(f)
            temp=out.(f{i});
            if length(temp)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1);
            out.(f{i})=temp;
        end
        f=fieldnames(metadata);
        for i=1:length(f)
            temp=metadata.(f{i});
            if length(temp)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1);
            metadata.(f{i})=temp;
        end
    end
    
    % fix sessids in metadata
    u=unique(metadata.sessid);
    for i=1:length(u)
        metadata.sessid(metadata.sessid==u(i))=i;
    end
    
    isreaching_out=countMouseSuccessDropMiss(alltbt,metadata,out);
    
    % Various methods for calculating dprimes
    settings=RTanalysis_settings();
    settings.preCueWindow_start=preCueWindow_start1; % define start of time window from trial onset, in seconds
    settings.preCueWindow_end=preCueWindow_end1; % define end of time window from trial onset, in seconds
    settings.reachAfterCueWindow_start=reachAfterCueWindow_start; % in sec, wrt cue onset
    settings.reachAfterCueWindow_end=reachAfterCueWindow_end; % in sec, wrt cue onset
    disp(['preCueWindow 1 is ' num2str(settings.preCueWindow_start) ' to ' num2str(settings.preCueWindow_end) ' secs from beginning of trial']);
    [dprimes_preCue,hit_rates,FA_rates_preCue]=get_dprime_per_session(alltbt,out,metadata,reachName,nameOfCue,settings);
    
    settings=RTanalysis_settings();
    settings.preCueWindow_start=preCueWindow_start2; % define start of time window from trial onset, in seconds
    settings.preCueWindow_end=preCueWindow_end2; % define end of time window from trial onset, in seconds
    [~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
    cuetimeat=mode(diff(nanmean(alltbt.times,1)))*ma;
    disp(['preCueWindow 2 is ' num2str(settings.preCueWindow_start-cuetimeat) ' to ' num2str(settings.preCueWindow_end-cuetimeat) ' secs from cue onset']);
    settings.reachAfterCueWindow_start=reachAfterCueWindow_start; % in sec, wrt cue onset
    settings.reachAfterCueWindow_end=reachAfterCueWindow_end; % in sec, wrt cue onset

    if settingsForDp.postCue_onlyWithDistractorTrials==true
        [dprimes_postCue,~,FA_rates_postCue]=distractVNoDistract_dprime(alltbt,metadata,trialTypes);
    else
        [dprimes_postCue,~,FA_rates_postCue]=get_dprime_per_session(alltbt,out,metadata,reachName,nameOfCue,settings);
    end
    
    max_FA=max([FA_rates_preCue; FA_rates_postCue],[],1);
    
    isreaching_out.dprimes_preCue=dprimes_preCue;
    isreaching_out.dprimes_postCue=dprimes_postCue;
    isreaching_out.max_FA=max_FA;
    dprimes=min([dprimes_preCue; dprimes_postCue],[],1);
    
    if ~isempty(firstSess)
        if length(firstSess)==1
            useSess=isreaching_out.nth_session>=firstSess;
            isreaching_out.nth_session=isreaching_out.nth_session(useSess);
            max_FA=max_FA(useSess);
            hit_rates=hit_rates(useSess);
        elseif length(firstSess)>1
            % use these sessions
            useSess=firstSess;
            isreaching_out.nth_session=isreaching_out.nth_session(useSess);
            max_FA=max_FA(useSess);
            hit_rates=hit_rates(useSess);
        end
    end
    
    % Plot hit rate vs false alarm
    if doPlot==1
        figure();
        [~,i]=sort(isreaching_out.nth_session);
        scatter(max_FA(i(1)),hit_rates(i(1)),[],'k');
        hold on;
        plot(max_FA(i),hit_rates(i),'Color','k');
        xlabel('Fraction trials with false alarm');
        ylabel('Fraction trials with hit');
    end
else
    
    dprimes=[];
    
    % useTrials only
    if ~isempty(useTrials)
        f=fieldnames(alltbt);
        for i=1:length(f)
            temp=alltbt.(f{i});
            if size(temp,1)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1,:);
            alltbt.(f{i})=temp;
        end
        f=fieldnames(out);
        for i=1:length(f)
            temp=out.(f{i});
            if length(temp)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1);
            out.(f{i})=temp;
        end
        f=fieldnames(metadata);
        for i=1:length(f)
            temp=metadata.(f{i});
            if length(temp)~=length(useTrials)
                continue
            end
            temp=temp(useTrials==1);
            metadata.(f{i})=temp;
        end
    end
    
    % fix sessids in metadata, must be 1,2,3,4,etc.
    u=unique(metadata.sessid);
    metasessid=metadata.sessid;
    for i=1:length(u)
        metadata.sessid(metasessid==u(i))=i;
    end
    
    isreaching_out=countMouseSuccessDropMiss(alltbt,metadata,out);
    % Various methods for calculating cued and noncued reach rates
    settings=RTanalysis_settings();
    settings.preCueWindow_start=preCueWindow_start1; % define start of time window from trial onset, in seconds
    settings.preCueWindow_end=preCueWindow_end1; % define end of time window from trial onset, in seconds
    settings.reachAfterCueWindow_start=reachAfterCueWindow_start; % in sec, wrt cue onset
    settings.reachAfterCueWindow_end=reachAfterCueWindow_end; % in sec, wrt cue onset
    disp(['preCueWindow 1 is ' num2str(settings.preCueWindow_start) ' to ' num2str(settings.preCueWindow_end) ' secs from beginning of trial']);
    [mean_cued_reach_rate,mean_noncued_reach_rate]=reachRateCuedVsNonCued(alltbt,out,metadata,reachName,nameOfCue,settings);
    
    settings=RTanalysis_settings();
    settings.preCueWindow_start=preCueWindow_start2; % define start of time window from trial onset, in seconds
    settings.preCueWindow_end=preCueWindow_end2; % define end of time window from trial onset, in seconds
    [~,ma]=nanmax(nanmean(alltbt.(nameOfCue),1));
    cuetimeat=mode(diff(nanmean(alltbt.times,1)))*ma;
    disp(['preCueWindow 2 is ' num2str(settings.preCueWindow_start-cuetimeat) ' to ' num2str(settings.preCueWindow_end-cuetimeat) ' secs from cue onset']);
    settings.reachAfterCueWindow_start=reachAfterCueWindow_start; % in sec, wrt cue onset
    settings.reachAfterCueWindow_end=reachAfterCueWindow_end; % in sec, wrt cue onset
    [mean_cued_reach_rate_postCue,mean_noncued_reach_rate_postCue]=reachRateCuedVsNonCued(alltbt,out,metadata,reachName,nameOfCue,settings);
    
    max_non_cued=max([mean_noncued_reach_rate; mean_noncued_reach_rate_postCue],[],1);

    isreaching_out.cued_reach_rate=mean_cued_reach_rate;
    isreaching_out.noncued_reach_rate=max_non_cued;
    
    if ~isempty(firstSess)
        if length(firstSess)==1
            useSess=isreaching_out.nth_session>=firstSess;
            isreaching_out.nth_session=isreaching_out.nth_session(useSess);
            max_non_cued=max_non_cued(useSess);
            mean_cued_reach_rate=mean_cued_reach_rate(useSess);
        elseif length(firstSess)>1
            % use these sessions
            useSess=firstSess;
            isreaching_out.nth_session=isreaching_out.nth_session(useSess);
            max_non_cued=max_non_cued(useSess);
            mean_cued_reach_rate=mean_cued_reach_rate(useSess);
        end
    end

    % Plot noncued reaching vs cued reaching
    if doPlot==1
        if doQuiver==1
            figure(); 
            cmap=colormap('cool');
            k=1;
            kstep=ceil(size(cmap,1)/length(isreaching_out.nth_session));
            [~,i]=sort(isreaching_out.nth_session);
            scatter(max_non_cued(i(1)),mean_cued_reach_rate(i(1)),[],'k');
            hold on;
            for i=1:length(max_non_cued)-1
                quiver(max_non_cued(i),mean_cued_reach_rate(i),max_non_cued(i+1)-max_non_cued(i),mean_cued_reach_rate(i+1)-mean_cued_reach_rate(i),'Color',cmap(k,:));
                k=k+kstep;
                if k>size(cmap,1)
                    k=size(cmap,1);
                end
            end
            xlabel('Non-cued reaching rate (Hz)');
            ylabel('Cued reaching rate (Hz)');
        else
            figure();
            [~,i]=sort(isreaching_out.nth_session);
            scatter(max_non_cued(i(1)),mean_cued_reach_rate(i(1)),[],'k');
            hold on;
            plot(max_non_cued(i),mean_cued_reach_rate(i),'Color','k');
            xlabel('Non-cued reaching rate (Hz)');
            ylabel('Cued reaching rate (Hz)');
        end
    end
end

end

function [out_dprime,out_hit,out_FA]=distractVNoDistract_dprime(alltbt,metadata,trialTypes)

backups.alltbt=alltbt; backups.metadata=metadata; backups.trialTypes=trialTypes;
afterCueBins=floor((settings.preCueWindow_start-cuetimeat)./mode(diff(nanmean(alltbt.times,1))));
useTrials=any(alltbt.movie_distractor(:,ma:ma+afterCueBins)>0.5,2);
if ~isempty(useTrials)
    f=fieldnames(alltbt);
    for i=1:length(f)
        temp=alltbt.(f{i});
        if size(temp,1)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1,:);
        alltbt.(f{i})=temp;
    end
    f=fieldnames(out);
    for i=1:length(f)
        temp=out.(f{i});
        if length(temp)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1);
        out.(f{i})=temp;
    end
    f=fieldnames(metadata);
    for i=1:length(f)
        temp=metadata.(f{i});
        if length(temp)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1);
        metadata.(f{i})=temp;
    end
end
[dprimes_wDistract,hit_rates_wDistract,FA_rates_wDistract]=get_dprime_per_session(alltbt,out,metadata,reachName,nameOfCue,settings);

alltbt=backups.alltbt; metadata=backups.metadata; trialTypes=backups.trialTypes;
useTrials=~any(alltbt.movie_distractor(:,ma:ma+afterCueBins)>0.5,2);
if ~isempty(useTrials)
    f=fieldnames(alltbt);
    for i=1:length(f)
        temp=alltbt.(f{i});
        if size(temp,1)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1,:);
        alltbt.(f{i})=temp;
    end
    f=fieldnames(out);
    for i=1:length(f)
        temp=out.(f{i});
        if length(temp)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1);
        out.(f{i})=temp;
    end
    f=fieldnames(metadata);
    for i=1:length(f)
        temp=metadata.(f{i});
        if length(temp)~=length(useTrials)
            continue
        end
        temp=temp(useTrials==1);
        metadata.(f{i})=temp;
    end
end
[dprimes_noDistract,hit_rates_noDistract,FA_rates_noDistract]=get_dprime_per_session(alltbt,out,metadata,reachName,nameOfCue,settings);

out_hit=hit_rates_noDistract-hit_rates_wDistract;
out_FA=FA_rates_noDistract-FA_rates_wDistract;
out_dprime=dprime(out_hit,out_FA);

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end

function [mean_cued_reach_rate,mean_noncued_reach_rate,out]=reachRateCuedVsNonCued(tbt,out,metadata,whichReach,nameOfCue,settings)

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
    settings=RTanalysis_settings();
end
hitWindow_start=settings.reachAfterCueWindow_start; % wrt cue onset
hitWindow_end=settings.reachAfterCueWindow_end; % wrt cue onset
FAWindow_start=settings.preCueWindow_start; % wrt trial onset
FAWindow_end=settings.preCueWindow_end; % wrt trial onset

% calculate reach rate after cue
if ~isfield(settings,'lowThresh')
    tempset=reachExpt_analysis_settings;
    settings.lowThresh=tempset.lowThresh;
end
% Convert time window wrt cue onset into indices into data
cueInd=find(nanmean(tbt.(nameOfCue),1)>settings.lowThresh,1,'first');
startInds=floor(abs(hitWindow_start)/mode(diff(nanmean(tbt.times,1))))-1;
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

timeBin=mode(diff(nanmean(tbt.times,1)));
temp=tbt.(whichReach);
cuedWindowTimeBin=double(length(startInds:endInds).*timeBin);
temp(temp<settings.lowThresh)=0; temp(temp>=settings.lowThresh)=1;
cued_reach_rate=sum(temp(:,startInds:endInds),2,'omitnan'); %./cuedWindowTimeBin; % divide through later
% throw out all nan trials
allnantrials=all(isnan(tbt.times),2);
cued_reach_rate(allnantrials==1)=nan;
mean_cued_reach_rate=meanValue_per_session(metadata,cued_reach_rate);
mean_cued_reach_rate=mean_cued_reach_rate./cuedWindowTimeBin;

% calculate reach rate non-cued

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

FAWindowTimeBin=nansum(useInds==1).*timeBin;
FA_reach_rate=sum(temp(:,useInds==1),2,'omitnan')./FAWindowTimeBin;
% throw out all nan trials
FA_reach_rate(allnantrials==1)=nan;
mean_noncued_reach_rate=meanValue_per_session(metadata,FA_reach_rate);

% cued and non-cued reach rates for session, organize per trial
out.cued_reach_rate=nan(size(out.led,1),1);
for i=1:size(out.led,1)
    out.cued_reach_rate(i)=mean_cued_reach_rate(metadata.sessid(i));
end

out.noncued_reach_rate=nan(size(out.led,1),1);
for i=1:size(out.led,1)
    out.noncued_reach_rate(i)=mean_noncued_reach_rate(metadata.sessid(i));
end

end

function meanValForTrialsInGroup=meanValue_per_session(metadata,values)

[sessid,sessStartInd]=unique(metadata.sessid);
meanValForTrialsInGroup=nan(1,length(sessid));
sessStartInd=[sessStartInd; length(metadata.sessid)+1];
for i=1:length(sessid)
    curr_sessStartInd=sessStartInd(i);
    sessinds=curr_sessStartInd:sessStartInd(i+1)-1;
    meanValForTrialsInGroup(i)=nanmean(values(sessinds));
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

function out=countMouseSuccessDropMiss(tbt,metadata,ttclass)

% gets the # of times mouse grabbed pellet successfully, dropped pellet or
% missed pellet and as a function of pellet availability

thresh=0.5;
% Get trials where pellet arrived before a cue
%[~,indForPresented]=max(nanmean(tbt.pelletPresented,1)); % when was pellet presented, in indices
[~,indForCue]=max(nanmean(tbt.cueZone_onVoff,1)); % when cue on, indices
pelletPresentOnThisTrial=any(tbt.pelletPresent(:,indForCue-25:indForCue)>thresh,2);

sesstypes=unique(metadata.sessid);
out.hasSuccess=nan(1,length(sesstypes));
out.hasDrop=nan(1,length(sesstypes));
out.hasMiss=nan(1,length(sesstypes));
out.hasReach=nan(1,length(sesstypes));
out.pelletPresent=nan(1,length(sesstypes));
out.touched_pellet=nan(1,length(sesstypes));
if ~isfield(tbt,'reachBatch_success_reachStarts_pawOnWheel')
    disp('discarding paw on wheel successes, etc.');
    disp('calculating successes, drops, etc. AFTER CUE');
    for i=1:length(sesstypes)
        currsessid=sesstypes(i);
        out.hasSuccess(i)=sum(any(tbt.reachBatch_success_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasDrop(i)=sum(any(tbt.reachBatch_drop_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasMiss(i)=sum(any(tbt.reachBatch_miss_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasReach(i)=sum(any(tbt.all_reachBatch(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.pelletPresent(i)=sum(pelletPresentOnThisTrial(metadata.sessid==currsessid));
        out.touched_pellet(i)=sum(ttclass.touched_pellet(metadata.sessid==currsessid));
        out.totalTrials(i)=sum(metadata.sessid==currsessid);
        out.consumedPelletAndNoPawOnWheel(i)=sum(ttclass.consumed_pellet(metadata.sessid==currsessid) & ~ttclass.paw_during_wheel(metadata.sessid==currsessid));
        if isfield(metadata,'optoOnHere')
            out.optoOnHere(i)=nanmean(metadata.optoOnHere(metadata.sessid==currsessid));
        end
        if isfield(metadata,'nth_session')
            out.nth_session(i)=nanmean(metadata.nth_session(metadata.sessid==currsessid));
        end
    end
else
    disp('using paw on wheel successes, etc.');
    disp('calculating successes, drops, etc. AFTER CUE');
    for i=1:length(sesstypes)
        currsessid=sesstypes(i);
        out.hasSuccess(i)=sum(any(tbt.reachBatch_success_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2))+sum(any(tbt.reachBatch_success_reachStarts_pawOnWheel(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasDrop(i)=sum(any(tbt.reachBatch_drop_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2))+sum(any(tbt.reachBatch_drop_reachStarts_pawOnWheel(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasMiss(i)=sum(any(tbt.reachBatch_miss_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2))+sum(any(tbt.reachBatch_miss_reachStarts_pawOnWheel(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.hasReach(i)=sum(any(tbt.all_reachBatch(metadata.sessid==currsessid,indForCue:end)>thresh,2));
        out.pelletPresent(i)=sum(pelletPresentOnThisTrial(metadata.sessid==currsessid));
        out.touched_pellet(i)=sum(ttclass.touched_pellet(metadata.sessid==currsessid));
        out.totalTrials(i)=sum(metadata.sessid==currsessid);
        out.consumedPelletAndNoPawOnWheel(i)=sum(ttclass.consumed_pellet(metadata.sessid==currsessid) & ~ttclass.paw_during_wheel(metadata.sessid==currsessid));
        if isfield(metadata,'optoOnHere')
            out.optoOnHere(i)=nanmean(metadata.optoOnHere(metadata.sessid==currsessid));
        end
        if isfield(metadata,'nth_session')
            out.nth_session(i)=nanmean(metadata.nth_session(metadata.sessid==currsessid));
        end
    end
end

end
