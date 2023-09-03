function [spiketimes,bincenters,N,edges]=SUresponseToCue(post_spikes,whichAssign,auxData,auxDataFieldForAlignment,timeWindowBeforeCueOnset,timeWindowAfterCueOnset,histogramBin)

% all times passed in must be in seconds

withOpto=nan;
if strcmp(auxDataFieldForAlignment,'cueWithoutOptoData') 
    auxDataFieldForAlignment='cueData';
    withOpto=false;
end
if strcmp(auxDataFieldForAlignment,'cuePlusOptoData')
    auxDataFieldForAlignment='cueData';
    withOpto=true;
end

optos=auxData.optoData.Values;

alignTo=auxData.(auxDataFieldForAlignment).Values;
if size(alignTo,1)>1 & size(alignTo,2)==1
    alignTo=alignTo';
end
unwrappedTimesVec=0:1/auxData.(auxDataFieldForAlignment).ADFreq:(1/auxData.(auxDataFieldForAlignment).ADFreq)*(length(alignTo)-1);
alignToOnsets=[0 diff(alignTo)==1];

spikes=filtspikes(post_spikes,0,'assigns',whichAssign);

alignOnsetsAt=find(alignToOnsets==1);
spiketimes=[];
for i=1:length(alignOnsetsAt)
    curronset=alignOnsetsAt(i);
    if isnan(withOpto)
        % take with or without opto
    elseif withOpto==false
        % take only cues without opto
        if optos(curronset)==1
            continue
        end
    elseif withOpto==true
        % take only cues with opto
        if optos(curronset)==0
            continue
        end
    end
    theseinds=spikes.unwrapped_times>=unwrappedTimesVec(curronset)-timeWindowBeforeCueOnset & spikes.unwrapped_times<=unwrappedTimesVec(curronset)+timeWindowBeforeCueOnset;
    spiketimes=[spiketimes spikes.unwrapped_times(theseinds)-unwrappedTimesVec(curronset)];
end

% histogram of spiketimes
bins=-timeWindowBeforeCueOnset:histogramBin:timeWindowAfterCueOnset;
[N,edges]=histcounts(spiketimes,bins);
bincenters=mean([edges(1:end-1); edges(2:end)],1,'omitnan');
figure(); 
plot(bincenters,N);

end