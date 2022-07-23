function [T,signals,signal_names,signal_phys_times]=bernardoTrialTable(tbt, bin1, bin2, bin3, signalSpikeFields, spikeFieldsNames, tbt_photo, signalPhotoFields, photoFieldsNames)

% bin1, bin2, bin3 are with respect to cue onset
% Time delay Arduino to cue onset

% Re-expand alignment into single row for signals
% Using alignment of photometry and physiology
% and times in tbt's
backup_tbt=[];
if ~isempty(signalPhotoFields) && isempty(signalSpikeFields)
    % this experiment had photometry but not physiology
    backup_tbt=tbt;
    tbt=tbt_photo; % use behavior tbt associated with photometry
elseif isempty(signalPhotoFields) && ~isempty(signalSpikeFields)
    % this experiment had physiology but not photometry
    % tbt is behavior tbt associated with physiology
elseif ~isempty(signalPhotoFields) && ~isempty(signalSpikeFields)
    % this experiment had both
    % only populate bernardoTrialTable for trials with both phys and photo
    % throw out trials that don't exist in both
    [tbt, tbt_photo, tookForPhys, tookForPhoto]=trimTbtsToMatch(tbt, tbt_photo);
    signalSpikeFields=takeOnlyTheseTrials(signalSpikeFields,tookForPhys);
    signalPhotoFields=takeOnlyTheseTrials(signalPhotoFields,tookForPhoto);
else
    % no neural data, only behavior data
end
[tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimes(tbt, tbt_photo, signalSpikeFields, signalPhotoFields);

% Add signals and names of signals
[signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signal_times_wrt_thisTrialStart,cuechunkoffset,signal_phys_times]=getLinedUpCueSignal(tbt, signalSpikeFields);
signal_cue_times_wrt_thisTrialCueStart(signal_cue_times_wrt_thisTrialCueStart<-100)=nan;
signal_cue_times(signal_cue_times<-100)=nan;
signal_times_wrt_thisTrialStart(signal_cue_times<-100)=nan;
% Double check
nextSignalTbt=tbt.cueZone_onVoff;
nextSignalTimes=tbt.times_wrt_trial_start;
nextSignalCue=tbt.cueZone_onVoff;
%sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
signal_ind=1;
signals{signal_ind}=sig;
signal_names{signal_ind}='cueZone_onVoff';
figure();
plot(signals{signal_ind});
hold on;
scatter(find(signal_cue>0.5),ones(size(find(signal_cue>0.5))));
% Now fill in other signals
nextSignalTbt=tbt.all_reachBatch;
nextSignalTimes=tbt.times_wrt_trial_start;
nextSignalCue=tbt.cueZone_onVoff;
%sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
signal_ind=signal_ind+1;
signals{signal_ind}=sig;
signal_names{signal_ind}='all_reachBatch';
% then go through phys and photo signals, adding them
% single units
if ~isempty(signalSpikeFields) 
    % get opto from phys acquisition
    nextSignalTbt=signalSpikeFields.opto;
    nextSignalTimes=signalSpikeFields.times_wrt_trial_start;
    nextSignalCue=signalSpikeFields.cue;
    %sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
    sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}='opto';
else
    % get opto from behavior acquisition
    nextSignalTbt=tbt.optoOn;
    nextSignalTimes=tbt.times_wrt_trial_start;
    nextSignalCue=tbt.cueZone_onVoff;
    %sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
    sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}='opto';
end
for i=1:length(spikeFieldsNames)
    nextSignalTbt=signalSpikeFields.(spikeFieldsNames{i});
    nextSignalTimes=signalSpikeFields.times_wrt_trial_start;
    nextSignalCue=signalSpikeFields.cue;
    resampledSignalCue=nan(size(nextSignalTbt));
    resampledPhysTimepoints=nan(size(nextSignalTbt));
    for j=1:size(nextSignalCue,1)
        if ~isempty(regexp(spikeFieldsNames{i},'unit'))
            cti=signalSpikeFields.unitTimes(j,:);
        else 
            cti=nextSignalTimes(j,:);
        end
        resampledSignalCue(j,:)=fillInFromMatchingTimes(nextSignalCue(j,:), nextSignalTimes(j,:), cti);
        resampledPhysTimepoints(j,:)=fillInFromMatchingTimes(signalSpikeFields.phys_timepoints(j,:), nextSignalTimes(j,:), cti);
        temp=resampledPhysTimepoints(j,:);
        endofgood=find(temp==nanmax(signalSpikeFields.phys_timepoints(j,:)),1,'first');
        temp(endofgood+1:end)=nan;
        temp(1:endofgood)=linspace(temp(1),temp(endofgood),endofgood);
        resampledPhysTimepoints(j,:)=temp;
    end
    %sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
    if ~isempty(regexp(spikeFieldsNames{i},'unit'))
        signalSpikeFields.resampPhysTimepoints=resampledPhysTimepoints;
        sig=fillInOtherSignalsPhysiology(nextSignalTbt,nextSignalTimes,resampledSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signalSpikeFields,spikeFieldsNames{i},signal_cue,signal_cue_times);
    else
        sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    end
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}=spikeFieldsNames{i};
end
% photometry channels
for i=1:length(photoFieldsNames)
    nextSignalTbt=signalPhotoFields.(photoFieldsNames{i});
    nextSignalTimes=signalPhotoFields.times_wrt_trial_start;
    % for photometry, sampling rate of cue might not match sampling rate of
    % processed photometry signal, so resample cue to match processed
    % photometry sampling rate
    nextSignalCue=signalPhotoFields.cue;
    resampledSignalCue=nan(size(nextSignalTbt));
    for j=1:size(nextSignalCue,1)
        resampledSignalCue(j,:)=fillInFromMatchingTimes(nextSignalCue(j,:), signalPhotoFields.cue_times(j,:), signalPhotoFields.green_time(j,:));
    end
    %sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial);
    %sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,resampledSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial);
    sig=fillInOtherSignalsPhotometry(nextSignalTbt,nextSignalTimes,resampledSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signalPhotoFields,photoFieldsNames{i},signal_cue,signal_cue_times);
    signal_ind=signal_ind+1;
    signals{signal_ind}=sig;
    signal_names{signal_ind}=photoFieldsNames{i};
end


% Trials
Trials=(1:size(tbt.cue,1))';
% TrialStart
% Pellet presentation wheel begins to move
TrialStart=[1 find(diff(signal_which_trial)==1)+1];
% Trial start seconds
TrialStart_s=signal_cue_times(TrialStart);
% Cue onset
Cue_onset=[0+1+cuechunkoffset find(diff(signal_which_trial)==1)+1+cuechunkoffset];
% Opto manipulation onset
whichIsCurrent=find(ismember(signal_names,'opto'));
temp=signals{whichIsCurrent};
optoStarts=[];
optoEnds=[];
for i=1:length(temp)-1
    if temp(i)<0.5 && temp(i+1)>0.5
        optoStarts=[optoStarts i+1];
    elseif temp(i)>0.5 && temp(i+1)<0.5
        optoEnds=[optoEnds i];
    end
end
% assign optos to trials
Opto_onset=nan(1,length(Cue_onset));
% Opto manipulation end
Opto_offset=nan(1,length(Cue_onset));
u=unique(signal_which_trial);
for i=1:length(u)
   f=find(signal_which_trial==u(i));
   thisopto=find(ismember(optoStarts,f(1:floor(end/2))));
   % opto should only ever be in first half of trial
   if ~isempty(thisopto)
       Opto_onset(i)=optoStarts(thisopto);
       Opto_offset(i)=optoEnds(thisopto);
   end
end
% First reach batch timing
% First reach batch type
[ReactionTimes,FirstCuedReachType]=getFirstReachesAndReachTypes(tbt);
% Bin 1 start
% find indices matching bin 1 start
[Bin1starts, Bin1ends]=getWindowIndices(bin1, signal_which_trial, signal_cue_times_wrt_thisTrialCueStart);
% Bin 1 end
% Number of reaches bin 1 (before the cue)
whichIsCurrent=find(ismember(signal_names,'all_reachBatch'));
ReachesInBin1=findReachesBetween(signals{whichIsCurrent}, Bin1starts, Bin1ends);
% Bin 2 start
[Bin2starts, Bin2ends]=getWindowIndices(bin2, signal_which_trial, signal_cue_times_wrt_thisTrialCueStart);
% Bin 2 end
% Number of reaches bin 2 (after the cue)
ReachesInBin2=findReachesBetween(signals{whichIsCurrent}, Bin2starts, Bin2ends);
% Bin 3 start
[Bin3starts, Bin3ends]=getWindowIndices(bin3, signal_which_trial, signal_cue_times_wrt_thisTrialCueStart);
% Bin 3 end
% Number of reaches bin 3 (long time after the cue)
ReachesInBin3=findReachesBetween(signals{whichIsCurrent}, Bin3starts, Bin3ends);
TrialStart=TrialStart';
TrialStart_s=TrialStart_s';
Cue_onset=Cue_onset';
Opto_onset=Opto_onset';
Opto_offset=Opto_offset';
ReactionTimes=ReactionTimes';
FirstCuedReachType=FirstCuedReachType';
Bin1starts=Bin1starts';
Bin1ends=Bin1ends';
ReachesInBin1=ReachesInBin1';
Bin2starts=Bin2starts';
Bin2ends=Bin2ends';
ReachesInBin2=ReachesInBin2';
Bin3starts=Bin3starts';
Bin3ends=Bin3ends';
ReachesInBin3=ReachesInBin3';
T=table(Trials,TrialStart,TrialStart_s,Cue_onset,Opto_onset,Opto_offset,ReactionTimes,...
        FirstCuedReachType,Bin1starts,Bin1ends,ReachesInBin1,Bin2starts,Bin2ends,ReachesInBin2,...
        Bin3starts,Bin3ends,ReachesInBin3,...
        'VariableNames',...
       {'Trials','TrialStart','TrialStart_s','Cue_onset','Opto_onset','Opto_offset','ReactionTimes',...
       'FirstCuedReachType','Bin1starts','Bin1ends','ReachesInBin1','Bin2starts','Bin2ends','ReachesInBin2',...
       'Bin3starts','Bin3ends','ReachesInBin3'});

end

function countReaches=findReachesBetween(reaches, binstarts, binends)

countReaches=nan(1,length(binstarts));
for i=1:length(binstarts)
    countReaches(i)=sum(reaches(binstarts(i):binends(i))>0.5,'all','omitnan');
end

end

function [binStarts, binEnds]=getWindowIndices(binWindow, signal_which_trial, signal_times_wrt_thisTrialStart)

t=unique(signal_which_trial(~isnan(signal_which_trial)));
binStarts=nan(1,length(t));
binEnds=nan(1,length(t));
for i=1:length(t)
    temp=signal_times_wrt_thisTrialStart;
    temp(signal_which_trial~=t(i))=nan;
    [~,binStarts(i)]=min(abs(temp-binWindow(1)),[],2,'omitnan');
    [~,binEnds(i)]=min(abs(temp-binWindow(2)),[],2,'omitnan');
end

end

function [reactionTimes,reachTypes]=getFirstReachesAndReachTypes(tbt)

temp=tbt.times_wrt_trial_start; 
temp=temp'; 
temp2=diff(temp(1:end)); 
timestep=mode(temp2(temp2>0));
[~,fcue]=max(nanmean(tbt.cueZone_onVoff,1),[],'all','linear','omitnan');
reactionTimes=nan(1,size(tbt.times_wrt_trial_start,1));
reachTypes=cell(1,size(tbt.times_wrt_trial_start,1));
for i=1:size(tbt.all_reachBatch,1)
    temp=tbt.all_reachBatch(i,:);
    firstreachind=find(temp(fcue:end)>0.5,1,'first')+fcue-1;
    reachAfterCueInd=find(temp(fcue:end)>0.5,1,'first');
    if isempty(reachAfterCueInd)
        % no reach, leave reactionTimes as nan
        reachTypes{i}='no reach';
        continue
    end
    reactionTimes(i)=reachAfterCueInd*timestep;
    if tbt.reachBatch_success_reachStarts(i,firstreachind)>0.5
        % success
        thisreachtype='success';
    elseif tbt.reachBatch_drop_reachStarts(i,firstreachind)>0.5
        % drop
        thisreachtype='drop';
    elseif tbt.reachBatch_miss_reachStarts(i,firstreachind)>0.5
        % miss
        thisreachtype='miss';
    elseif tbt.reachBatch_success_reachStarts_pawOnWheel(i,firstreachind)>0.5
        thisreachtype='success_pawOnWheel';
    elseif tbt.reachBatch_drop_reachStarts_pawOnWheel(i,firstreachind)>0.5
        thisreachtype='drop_pawOnWheel';
    elseif tbt.reachBatch_miss_reachStarts_pawOnWheel(i,firstreachind)>0.5
        thisreachtype='miss_pawOnWheel';
    else
        thisreachtype='unknown';
    end
    reachTypes{i}=thisreachtype;
end

end

function [tbt, tbt_photo, signalSpikeFields, signalPhotoFields]=addTimes(tbt, tbt_photo, signalSpikeFields, signalPhotoFields)

% have already trimmed all tbt's to have same number of rows in each field
if ~isempty(tbt)
    if ~isfield(tbt,'times_wrt_trial_start') 
        tbt=addTimesWrtTrialStarts(tbt,'times','cueZone_onVoff');
    end
end
if ~isempty(tbt_photo)
    if ~isfield(tbt_photo,'times_wrt_trial_start')
        tbt_photo=addTimesWrtTrialStarts(tbt_photo,'times','cueZone_onVoff');
    end
end
if ~isempty(signalSpikeFields)
    if ~isfield(signalSpikeFields,'times_wrt_trial_start')
        signalSpikeFields.phys_timepoints(signalSpikeFields.phys_timepoints==0)=nan;
        signalSpikeFields=addTimesWrtTrialStarts(signalSpikeFields,'phys_timepoints','cue');
        signalSpikeFields.times_wrt_trial_start(signalSpikeFields.times_wrt_trial_start<0)=nan;
    end
end
if ~isempty(signalPhotoFields)
    if ~isfield(signalPhotoFields,'times_wrt_trial_start')  
        signalPhotoFields=addTimesWrtTrialStarts(signalPhotoFields,'green_time','cue');
    end
end

end

function [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signal_times_wrt_thisTrialStart,cuechunkoffset,signal_phys_times]=getLinedUpCueSignal(tbt, signalSpikeFields)

% could use timesfromarduino to line trials back up 
% or movieframeinds*timestep sec
temp=tbt.times_wrt_trial_start; 
temp=temp'; 
temp2=diff(temp(1:end)); 
timestep=mode(temp2(temp2>0));
% each signal has its own time
% use cues to align across signals
% Plus need to resample each signal
% I think the best way to do this is going to be to put the cues at fixed
% times, then just fill in the other time points using the value from the
% closest time in the existing data structure -- yes

% If is physiology, use cue times from physiology
% Else use cue times from behavior tbt
if ~isempty(signalSpikeFields)
    temp=diff(signalSpikeFields.cue_times(1,:)); 
    ts=mode(temp(temp~=0));
    cuechunkoffset=floor(1.015/ts);
    [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signal_times_wrt_thisTrialStart,signal_phys_times]=lineUpCueChunks(signalSpikeFields, 'cue', 'cue_times', 'cuetimes_wrt_trial_start', cuechunkoffset);
else
    temp=diff(tbt.times_wrt_trial_start(1,:)); 
    ts=mode(temp(temp~=0));
    cuechunkoffset=floor(1.015/ts);
    [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signal_times_wrt_thisTrialStart,signal_phys_times]=lineUpCueChunks(tbt, 'cueZone_onVoff', 'times_wrt_trial_start', 'times_wrt_trial_start', cuechunkoffset);
    signal_cue_times=0:timestep:(length(signal_cue_times)-1)*timestep;
end


end

function sig=fillInOtherSignalsSkipCue(nextSignalTbt,nextSignalTimes,signal_times_wrt_thisTrialStart,signal_which_trial)

sig=nan(size(signal_times_wrt_thisTrialStart)); % must be the same size as other signals
% just use nearest neighbor to interpolate or downsample
for i=1:size(nextSignalTbt,1)
    temp=nextSignalTbt(i,:);
    times=nextSignalTimes(i,:);
    % relative to the start of this trial
    % resample temp
    f=find(signal_which_trial==i);
    timeswrtcuestart=signal_times_wrt_thisTrialStart(f);
    timeswrtcuestart_nextsignal=times;
    if max(timeswrtcuestart,[],'all','omitnan')>max(timeswrtcuestart_nextsignal,[],'all','omitnan')
        % fill in times
        tempdiff=diff(timeswrtcuestart_nextsignal);
        ts=mode(tempdiff(~isnan(tempdiff)));
        fin=find(~isnan(timeswrtcuestart_nextsignal),1,'last')+1;
        timeswrtcuestart_nextsignal(fin:end)=timeswrtcuestart_nextsignal(fin-1)+ts:ts:timeswrtcuestart_nextsignal(fin-1)+length(timeswrtcuestart_nextsignal(fin:end))*ts;
    end
    for j=1:length(timeswrtcuestart)
        if isnan(timeswrtcuestart(j))
            sig(f(j))=nan;
            continue
        end
        % find nearest neighbor
        [~,mi]=nanmin(abs(timeswrtcuestart_nextsignal-timeswrtcuestart(j)));
        sig(f(j))=temp(mi);
    end
    %if nansum(diff(sig)==0)>50
    %    a=1;
    %end
end

end

function sig=fillInOtherSignals(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial)

sig=nan(size(signal_cue_times_wrt_thisTrialCueStart)); % must be the same size as other signals
% just use nearest neighbor to interpolate or downsample
for i=1:size(nextSignalTbt,1)
    temp=nextSignalTbt(i,:);
    times=nextSignalTimes(i,:);
    tempcue=nextSignalCue(i,:);
    % relative to the timing of the cue in this trial
    % resample temp
    f=find(signal_which_trial==i);
    timeswrtcuestart=signal_cue_times_wrt_thisTrialCueStart(f);
    fcue=find(tempcue>0.5,1,'first');
    timeswrtcuestart_nextsignal=times-times(fcue);
    if max(timeswrtcuestart,[],'all','omitnan')>max(timeswrtcuestart_nextsignal,[],'all','omitnan')
        % fill in times
        tempdiff=diff(timeswrtcuestart_nextsignal);
        ts=mode(tempdiff(~isnan(tempdiff)));
        fin=find(~isnan(timeswrtcuestart_nextsignal),1,'last')+1;
        timeswrtcuestart_nextsignal(fin:end)=timeswrtcuestart_nextsignal(fin-1)+ts:ts:timeswrtcuestart_nextsignal(fin-1)+length(timeswrtcuestart_nextsignal(fin:end))*ts;
    end
    for j=1:length(timeswrtcuestart)
        if isnan(timeswrtcuestart(j))
            sig(f(j))=nan;
            continue
        end
        % find nearest neighbor
        [~,mi]=nanmin(abs(timeswrtcuestart_nextsignal-timeswrtcuestart(j)));
        sig(f(j))=temp(mi);
    end
    %if nansum(diff(sig)==0)>50
    %    a=1;
    %end
end

end

function sig=fillInOtherSignalsPhotometry(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,photo_tbt,whichCh,signal_cue,signal_cue_times)

% make single continuous photometry signal out of unique times
green_times=[];
green_ch_data=[];
cue_data=[];
photodata=photo_tbt.(whichCh);
for i=1:size(photo_tbt.green_time,1)
    f=find(~isnan(photo_tbt.green_time(i,:)) & photo_tbt.green_time(i,:)>0);
    if i~=1
        [~,mi]=nanmin(abs(photo_tbt.green_time(i,:)-lastTimeTaken));
        fend=f(end);
        f=mi:fend;
    end
    green_times=[green_times photo_tbt.green_time(i,f)];
    green_ch_data=[green_ch_data photodata(i,f)];
    cue_data=[cue_data nextSignalCue(i,f)];
    lastTimeTaken=photo_tbt.green_time(i,f(end));
end

firstcuephoto=green_times(find(cue_data>0.5,1,'first'));
firstcuephys=signal_cue_times(find(signal_cue>0.5,1,'first'));
figure(); plot(green_times-green_times(find(cue_data>0.5,1,'first')),cue_data,'Color','k'); 
hold on; plot(signal_cue_times-signal_cue_times(find(signal_cue>0.5,1,'first')),signal_cue,'Color','b');
answer=questdlg('Does alignment look good?');
switch answer
    case 'Yes'
        useAlignment=true;
    case 'No'
        useAlignment=false;
    case 'Cancel'
        error('Stop script');
    case ''
        error('Stop script');
end

sig=nan(size(signal_cue_times_wrt_thisTrialCueStart)); % must be the same size as other signals
if useAlignment==false
    % just use nearest neighbor to interpolate or downsample
    for i=1:size(nextSignalTbt,1)
        temp=nextSignalTbt(i,:);
        times=nextSignalTimes(i,:);
        tempcue=nextSignalCue(i,:);
        greentimesreference=photo_tbt.green_time(i,:);
        % relative to the timing of the cue in this trial
        % resample temp
        f=find(signal_which_trial==i);
        timeswrtcuestart=signal_cue_times_wrt_thisTrialCueStart(f);
        fcue=find(tempcue>0.5,1,'first');
        timeswrtcuestart_nextsignal=times-times(fcue);
        if max(timeswrtcuestart,[],'all','omitnan')>max(timeswrtcuestart_nextsignal,[],'all','omitnan')
            % fill in times
            tempdiff=diff(timeswrtcuestart_nextsignal);
            ts=mode(tempdiff(~isnan(tempdiff)));
            fin=find(~isnan(timeswrtcuestart_nextsignal),1,'last')+1;
            timeswrtcuestart_nextsignal(fin:end)=timeswrtcuestart_nextsignal(fin-1)+ts:ts:timeswrtcuestart_nextsignal(fin-1)+length(timeswrtcuestart_nextsignal(fin:end))*ts;
        end
        for j=1:length(timeswrtcuestart)
            if isnan(timeswrtcuestart(j))
                sig(f(j))=nan;
                continue
            end
            % find nearest neighbor
            [~,mi]=nanmin(abs(timeswrtcuestart_nextsignal-timeswrtcuestart(j)));
            % mi is index into this row
            % use green times as reference
            [~,mi2]=nanmin(abs(green_times-greentimesreference(mi)));
            %sig(f(j))=temp(mi);
            sig(f(j))=green_ch_data(mi2);
        end
        %if nansum(diff(sig)==0)>50
        %    a=1;
        %end
    end
else
    % use alignment
    green_times=green_times-green_times(find(cue_data>0.5,1,'first'))+signal_cue_times(find(signal_cue>0.5,1,'first'));
    figure(); plot(green_times,cue_data,'Color','k'); 
    hold on;
    plot(signal_cue_times,signal_cue,'Color','b');
    % times in green_times now match times in signal_cue_times
    % so just get closest value from green_ch for each time in
    % signal_cue_times
    for i=1:length(signal_cue_times)
        [~,mi]=nanmin(abs(green_times-signal_cue_times(i)));
        sig(i)=green_ch_data(mi);
    end
end

end

function sig=fillInOtherSignalsPhysiology(nextSignalTbt,nextSignalTimes,nextSignalCue,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,photo_tbt,whichCh,signal_cue,signal_cue_times)

% make single continuous photometry signal out of unique times
green_times=[];
green_ch_data=[];
cue_data=[];
photodata=photo_tbt.(whichCh);
for i=1:size(photo_tbt.resampPhysTimepoints,1)
    f=find(~isnan(photo_tbt.resampPhysTimepoints(i,:)) & photo_tbt.resampPhysTimepoints(i,:)>0);
    if i~=1
        [~,mi]=nanmin(abs(photo_tbt.resampPhysTimepoints(i,:)-lastTimeTaken));
        fend=f(end);
        f=mi:fend;
    end
    green_times=[green_times photo_tbt.resampPhysTimepoints(i,f)];
    if i<size(photo_tbt.resampPhysTimepoints,1)
        % find next trial's cue
        nextcue=find(nextSignalCue(i+1,:)>0.5,1,'first');
        takebegininds=1:nextcue-1;
        photodata(i,f(end)-length(takebegininds)+1:f(end))=photodata(i+1,takebegininds);
    end
    green_ch_data=[green_ch_data photodata(i,f)];
    cue_data=[cue_data nextSignalCue(i,f)];
    lastTimeTaken=photo_tbt.resampPhysTimepoints(i,f(end));
end

firstcuephoto=green_times(find(cue_data>0.5,1,'first'));
firstcuephys=signal_cue_times(find(signal_cue>0.5,1,'first'));
figure(); plot(green_times-green_times(find(cue_data>0.5,1,'first')),cue_data,'Color','k'); 
hold on; plot(signal_cue_times-signal_cue_times(find(signal_cue>0.5,1,'first')),signal_cue,'Color','b');
answer=questdlg('Does alignment look good?');
switch answer
    case 'Yes'
        useAlignment=true;
    case 'No'
        useAlignment=false;
    case 'Cancel'
        error('Stop script');
    case ''
        error('Stop script');
end
useAlignment=true;

sig=nan(size(signal_cue_times_wrt_thisTrialCueStart)); % must be the same size as other signals
if useAlignment==false
    % just use nearest neighbor to interpolate or downsample
    for i=1:size(nextSignalTbt,1)
        temp=nextSignalTbt(i,:);
        times=nextSignalTimes(i,:);
        tempcue=nextSignalCue(i,:);
        greentimesreference=photo_tbt.green_time(i,:);
        % relative to the timing of the cue in this trial
        % resample temp
        f=find(signal_which_trial==i);
        timeswrtcuestart=signal_cue_times_wrt_thisTrialCueStart(f);
        fcue=find(tempcue>0.5,1,'first');
        timeswrtcuestart_nextsignal=times-times(fcue);
        if max(timeswrtcuestart,[],'all','omitnan')>max(timeswrtcuestart_nextsignal,[],'all','omitnan')
            % fill in times
            tempdiff=diff(timeswrtcuestart_nextsignal);
            ts=mode(tempdiff(~isnan(tempdiff)));
            fin=find(~isnan(timeswrtcuestart_nextsignal),1,'last')+1;
            timeswrtcuestart_nextsignal(fin:end)=timeswrtcuestart_nextsignal(fin-1)+ts:ts:timeswrtcuestart_nextsignal(fin-1)+length(timeswrtcuestart_nextsignal(fin:end))*ts;
        end
        for j=1:length(timeswrtcuestart)
            if isnan(timeswrtcuestart(j))
                sig(f(j))=nan;
                continue
            end
            % find nearest neighbor
            [~,mi]=nanmin(abs(timeswrtcuestart_nextsignal-timeswrtcuestart(j)));
            % mi is index into this row
            % use green times as reference
            [~,mi2]=nanmin(abs(green_times-greentimesreference(mi)));
            %sig(f(j))=temp(mi);
            sig(f(j))=green_ch_data(mi2);
        end
        %if nansum(diff(sig)==0)>50
        %    a=1;
        %end
    end
else
    % use alignment
    green_times=green_times-green_times(find(cue_data>0.5,1,'first'))+signal_cue_times(find(signal_cue>0.5,1,'first'));
    figure(); plot(green_times,cue_data,'Color','k'); 
    hold on;
    plot(signal_cue_times,signal_cue,'Color','b');
    % times in green_times now match times in signal_cue_times
    % so just get closest value from green_ch for each time in
    % signal_cue_times
    % note that phys is up-sampled wrt behavior, so need to average to
    % prevent aliasing issues
    % resample green_times and green_ch_data
    green_times=downSampAv(green_times,ceil(length(green_times)/length(signal_cue_times)));
    green_ch_data=downSampAv(green_ch_data,ceil(length(green_ch_data)/length(signal_cue_times)));
    for i=1:length(signal_cue_times)
        [~,mi]=nanmin(abs(green_times-signal_cue_times(i)));
        sig(i)=green_ch_data(mi);
    end
end

end

function out=fillInFromMatchingTimes(signal, signal_times, match_to_times)

out=nan(1,length(match_to_times));
for i=1:length(match_to_times)
    % find nearest neighbor
    [~,mi]=nanmin(abs(signal_times-match_to_times(i)));
    out(i)=signal(mi);
end

end

function [signal_cue,signal_cue_times,signal_cue_times_wrt_thisTrialCueStart,signal_which_trial,signal_times_wrt_thisTrialStart,signal_phys_times]=lineUpCueChunks(tbt, cuefield, cuetimes, timeswrttrialstart, cuechunkoffset)

cue=tbt.(cuefield);
cue_times=tbt.(cuetimes);
timeswrt=tbt.(timeswrttrialstart);
if isfield(tbt, 'phys_timepoints')
    phystimes=tbt.phys_timepoints;
end
% Each chunk begins at the start of each cue onset, chunk continues until
% the start of the next chunk onset
signal_cue=[];
signal_cue_times=[];
signal_cue_times_wrt_thisTrialCueStart=[];
signal_times_wrt_thisTrialStart=[];
signal_which_trial=[];
signal_phys_times=[];
% if this is physiology_tbt, can also use phys_timepoints as a reference
for i=1:size(cue,1)
    % find first cue onset 
    temp=cue(i,:);
    f=find(temp>0.5,1,'first');
    fcue=f;
    fnan=find(isnan(temp(fcue:end)),1,'first')+fcue-1;
    temp(fnan:end)=nan;
    % find index before next cue onset (or last index of trial)
    f2=find(temp<0.5,1,'last');
    f=f-cuechunkoffset;
    f2=f2-cuechunkoffset;
    if i==1
        f=1; % include ITI before first trial
    else
        if isfield(tbt, 'phys_timepoints')
            % adjust according to last time taken
            [~,mi]=nanmin(abs(tbt.phys_timepoints(i,:)-lastTimeTaken));
            f=mi+1;
        end
    end
    if isfield(tbt, 'phys_timepoints')
        lastTimeTaken=tbt.phys_timepoints(i,f2);
    end
    chunk=temp(f:f2);
    signal_cue=[signal_cue chunk];
    if isfield(tbt, 'phys_timepoints')
        signal_phys_times=[signal_phys_times phystimes(i,f:f2)];
    end
    signal_cue_times=[signal_cue_times cue_times(i,f:f2)];
    signal_cue_times_wrt_thisTrialCueStart=[signal_cue_times_wrt_thisTrialCueStart cue_times(i,f:f2)-cue_times(i,fcue)];
    signal_times_wrt_thisTrialStart=[signal_times_wrt_thisTrialStart timeswrt(i,f:f2)];
    signal_which_trial=[signal_which_trial i*ones(size(chunk))];
end

end

function tbt=addTimesWrtTrialStarts(tbt,whichTimesField,cuename)

if ~isempty(whichTimesField)
    timesfield=tbt.(whichTimesField);
else
    % look for a times field
    % find it automatically by searching for token 'times'
    f=fieldnames(tbt);
    for i=1:length(f)
        if ~isempty(regexp(f{i},'times'))
            timesfield=tbt.(f{i});
            break
        end
    end
end

[~,f]=max(mean(tbt.(cuename),1,'omitnan'),[],2);
tbt.times_wrt_trial_start=timesfield-repmat(min(timesfield(:,1:f),[],2,'omitnan'),1,size(timesfield,2));

end

function [behphys_tbt, behphoto_tbt, tookForPhys, tookForPhoto]=trimTbtsToMatch(behphys_tbt, behphoto_tbt)

foundReference=false;
if isfield(behphys_tbt, 'this_is_which_beh')
    if behphys_tbt.this_is_which_beh==1
        if isfield(behphys_tbt, 'reference_into_beh2trialinds')
            if any(any(~isnan(behphys_tbt.reference_into_beh2trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_phys_to_photo=behphys_tbt.reference_into_beh2trialinds(:,1);
                % throw out trials that are missing
                tookForPhys=~isnan(ref_from_phys_to_photo);
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,~isnan(ref_from_phys_to_photo));
                % throw out trials that are missing 
                tookForPhoto=ref_from_phys_to_photo;
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,ref_from_phys_to_photo);
            else
                % this field is missing reference
            end
        end
    elseif behphys_tbt.this_is_which_beh==2
        if isfield(behphys_tbt, 'reference_into_beh1trialinds')
            if any(any(~isnan(behphys_tbt.reference_into_beh1trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_phys_to_photo=behphys_tbt.reference_into_beh1trialinds(:,1);
                % throw out trials that are missing 
                tookForPhys=~isnan(ref_from_phys_to_photo);
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,~isnan(ref_from_phys_to_photo));
                % throw out trials that are missing 
                tookForPhoto=ref_from_phys_to_photo;
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,ref_from_phys_to_photo);
            else
                % this field is missing reference
            end
        end
    end
elseif isfield(behphoto_tbt, 'this_is_which_beh')
    if behphoto_tbt.this_is_which_beh==1
        if isfield(behphoto_tbt, 'reference_into_beh2trialinds')
            if any(any(~isnan(behphoto_tbt.reference_into_beh2trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_photo_to_phys=behphoto_tbt.reference_into_beh2trialinds(:,1);
                % throw out trials that are missing
                tookForPhoto=~isnan(ref_from_photo_to_phys);
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,~isnan(ref_from_photo_to_phys));
                % throw out trials that are missing 
                tookForPhys=ref_from_photo_to_phys;
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,ref_from_photo_to_phys);
            else
                % this field is missing reference
            end
        end
    elseif behphoto_tbt.this_is_which_beh==2
        if isfield(behphoto_tbt, 'reference_into_beh1trialinds')
            if any(any(~isnan(behphoto_tbt.reference_into_beh1trialinds),2))
                foundReference=true;
                % use this set of trials and reference
                ref_from_photo_to_phys=behphoto_tbt.reference_into_beh1trialinds(:,1);
                % throw out trials that are missing 
                tookForPhoto=~isnan(ref_from_photo_to_phys);
                behphoto_tbt=takeOnlyTheseTrials(behphoto_tbt,~isnan(ref_from_photo_to_phys));
                % throw out trials that are missing 
                tookForPhys=ref_from_photo_to_phys;
                behphys_tbt=takeOnlyTheseTrials(behphys_tbt,ref_from_photo_to_phys);
            else
                % this field is missing reference
            end
        end
    end
end
if foundReference==false
    error('Could not find reference to align behavior data for photometry and physiology');
end
if size(behphys_tbt.cue,1)~=size(behphoto_tbt.cue,1)
    error('Tbts must be same size');
end

end

function tbt=takeOnlyTheseTrials(tbt,takeThese)

f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    if isvector(temp)
        continue
    end
    temp=temp(takeThese(~isnan(takeThese)),:);
    tbt.(f{i})=temp;
end

end