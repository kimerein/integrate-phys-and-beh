function [spikes,auxData]=convertSpikesToTrials_and_saveAuxChs(spikes,filename)

% Note that spikes before first trial will be discarded!

% seconds before cue to begin trial
cueRelativeToTrialStart=1; % in seconds

ledch='auxData193'; % WHISPER system distractor LED
cuech='auxData195invert'; % WHISPER system inverted cue
optoch='auxData194';

[cuetimes,cueData]=getEventsFromAnalogCh(filename,cuech);
[distractortimes,distractorData]=getEventsFromAnalogCh(filename,ledch);
[optotimes,optoData]=getEventsFromAnalogCh(filename,optoch);
auxData.cueData=cueData;
auxData.distractorData=distractorData;
auxData.distractorData.Values=smooth(auxData.distractorData.Values,3);
auxData.distractorData.Values=auxData.distractorData.Values>0;
auxData.optoData=optoData;

spikes=convertContinuousToTrialStructure(spikes,cuetimes-cueRelativeToTrialStart,cueRelativeToTrialStart);

% truncate spikes before first trial
throwOut=isnan(spikes.trials);
f=fieldnames(spikes);
numSpikesAtBeginning=length(spikes.trials);
for i=1:length(f)
    temp=spikes.(f{i});
    if isstruct(temp)
        % ignore
    elseif strcmp(f{i},'waveforms')
        temp=temp(throwOut==false,:,:);
    elseif length(temp)==numSpikesAtBeginning
        temp=temp(throwOut==false);
    end
    spikes.(f{i})=temp;
end

% plot trial durations histogram
[n,x]=histcounts(cuetimes(2:end)-cuetimes(1:end-1),10);
figure(); plot(nanmean([x(1:end-1); x(2:end)],1),n);
xlabel('trial duration (seconds)');
ylabel('count');
title('distribution of trial durations in this experiment');

end

function spikes=convertContinuousToTrialStructure(spikes,trialStartTimes,cueRelativeToTrialStart)

whichTrialStart=nan(1,length(spikes.spiketimes));
for i=1:length(trialStartTimes)
    temp=spikes.spiketimes-trialStartTimes(i);
    whichTrialStart(temp>=0)=i;
end

spikes.trials=whichTrialStart;
temp=1:length(trialStartTimes);
temp2=spikes.trials;
temp2(isnan(temp2))=1;
spikes.spiketimes=spikes.spiketimes-trialStartTimes(temp2);
spikes.info.cue_thisManyS_beforeTrialOnset=cueRelativeToTrialStart;

spikes.sweeps.trials=temp;
% note that we don't actually know last trial duration
spikes.sweeps.trialDuration=[trialStartTimes(2:end)-trialStartTimes(1:end-1) nan];

end

