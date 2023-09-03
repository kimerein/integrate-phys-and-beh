function [spikes,auxData,answersout]=convertSpikesToTrials_and_saveAuxChs(varargin)

% Note that spikes before first trial will be discarded!

if length(varargin)==2
    spikes=varargin{1};
    filename=varargin{2};
    answers=nan(1,6);
elseif length(varargin)==3
    spikes=varargin{1};
    filename=varargin{2};
    answers=varargin{3};
else
    error('Wrong number of arguments to convertSpikesToTrials_and_saveAuxChs.m');
end

% seconds before cue to begin trial
cueRelativeToTrialStart=1; % in seconds

% THERE IS AN OFFSET SUCH THAT SPIKEGLX READS "ch 192" 
% BUT THIS IS SAVED AND READ INTO MATLAB AS "ch 193"
ledch='auxData193'; % WHISPER system distractor LED
cuech='auxData195invert'; % WHISPER system inverted cue
optoch='auxData196invert';

answersout=nan(1,6);
[cuetimes,cueData,~,answerout]=getEventsFromAnalogCh(filename,cuech,answers(1:2));
answersout(1:2)=answerout;
[distractortimes,distractorData,~,answerout]=getEventsFromAnalogCh(filename,ledch,answers(3:4));
answersout(3:4)=answerout;
[optotimes,optoData,~,answerout]=getEventsFromAnalogCh(filename,optoch,answers(5:6));
answersout(5:6)=answerout;
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

