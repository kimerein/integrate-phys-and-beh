function loopThroughAndGetEachUnitCueResponse(dirLocation,saveDir,modulationTimeWindowBefore,modulationTimeWindowAfter)

% This was the channel ordering for Shi's spike sorting: [16 14 15 11 13 9 12 7 10 3 5 1 8 6 4 2 31 29 27 25 32 17 30 18 28 19 26 20 23 22 24 21]

files=findFiles(dirLocation);
% Loop through each spikes
answersout=nan(1,6);
for i=1:length(files)
    currname=files(i).name;
    number=str2double(regexp(currname, '(?<=spikes)\d+', 'match', 'once'));
    if isnan(number)
        number=1;
    end
    a=load(fullfile(dirLocation,currname));
    [post_spikes,auxData,answersout]=convertSpikesToTrials_and_saveAuxChs(a.spikes,fullfile(dirLocation,'auxData193.mat'),answersout);
    % find good units
    unittypes=post_spikes.labels(:,2);
    unitassignstouse=post_spikes.labels(unittypes==2,1);
    for j=1:length(unitassignstouse)
        chnumber=channelMap(number,post_spikes,unitassignstouse(j)); % 32 is most dorsal, 1 is most ventral
        saveName=fullfile(saveDir,['unit' num2str(unitassignstouse(j)) 'onCh' num2str(chnumber)]);
        [spiketimes,bincenters,N,edges,modout]=SUresponseToCue(post_spikes,unitassignstouse(j),auxData,'cueData',2,2,0.01,modulationTimeWindowBefore,modulationTimeWindowAfter);
        saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,'cueAligned',modout);
        [spiketimes,bincenters,N,edges,modout]=SUresponseToCue(post_spikes,unitassignstouse(j),auxData,'distractorData',2,2,0.01,modulationTimeWindowBefore,modulationTimeWindowAfter);
        saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,'distractorAligned',modout);
        [spiketimes,bincenters,N,edges,modout]=SUresponseToCue(post_spikes,unitassignstouse(j),auxData,'optoData',2,2,0.01,modulationTimeWindowBefore,modulationTimeWindowAfter);
        saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,'optoAligned',modout);
        [spiketimes,bincenters,N,edges,modout]=SUresponseToCue(post_spikes,unitassignstouse(j),auxData,'cueWithoutOptoData',2,2,0.01,modulationTimeWindowBefore,modulationTimeWindowAfter);
        saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,'cueWithoutOptoData',modout);
        [spiketimes,bincenters,N,edges,modout]=SUresponseToCue(post_spikes,unitassignstouse(j),auxData,'cuePlusOptoData',2,2,0.01,modulationTimeWindowBefore,modulationTimeWindowAfter);
        saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,'cuePlusOptoData',modout);
        close all;
    end
end

end

function saveUnitAlignmentData(spiketimes,bincenters,N,edges,saveName,addName,modulation)

alignData.spiketimes=spiketimes;
alignData.bincenters=bincenters;
alignData.N=N;
alignData.edges=edges;
alignData.modulation=modulation;
save([saveName addName '.mat'],'alignData');

end

function files=findFiles(directoryPath)

% Get list of all "_sorted.mat" files in the directory
files = dir(fullfile(directoryPath, '*_sorted.mat'));

end

function chnumber=channelMap(spikesnumber,spikes,currAssign)

% second row is Matlab index, first row is depth on probe, where 32 is most
% dorsall, 1 is most ventral
% for A1x32Edge
% was using 20 micron spacing between chs
         chmap=[1   21; ...
                2   24; ...
                3   22; ...
                4   23; ...
                5   20; ...
                6   26; ...
                7   19; ...
                8   28; ...
                9   18; ...
                10  30; ...
                11  17; ...
                12  32; ...
                13  25; ...
                14  27; ...
                15  29; ...
                16  31; ...
                17  2; ...
                18  4; ...
                19  6; ...
                20  8; ...
                21  1; ...
                22  5; ...
                23  3; ...
                24  10; ...
                25  7; ...
                26  12; ...
                27  9; ...
                28  13; ...
                29  11; ...
                30  15; ...
                31  14; ...
                32  16];

switch spikesnumber
    case 1
        trodeChsForSpikes=[16 14 15 11];
    case 2
        trodeChsForSpikes=[13 9 12 7];
    case 3
        trodeChsForSpikes=[10 3 5 1];
    case 4
        trodeChsForSpikes=[8 6 4 2];
    case 5
        trodeChsForSpikes=[31 29 27 25];
    case 6
        trodeChsForSpikes=[32 17 30 18];
    case 7
        trodeChsForSpikes=[28 19 26 20];
    case 8
        trodeChsForSpikes=[];
    case 9
        trodeChsForSpikes=[23 22 24 21];
    otherwise
        error('do not recognize name of spikes mat file');
end

% find ch where waveform biggest
amp=nan(1,size(spikes.waveforms,3));
for l=1:size(spikes.waveforms,3)
    amp(l)=abs(min(reshape(mean(spikes.waveforms(spikes.assigns==currAssign,:,l),1,'omitnan'),1,size(spikes.waveforms,2)),[],2,'omitnan'));
end
[~,si]=sort(amp,'ascend');
% sort trode chs for units
if any(si>length(trodeChsForSpikes))
    error(['problem indexing spike channel for unit']);
end
trodeChsForSpikes=trodeChsForSpikes(si(end));

chnumber=find(chmap(:,2)==trodeChsForSpikes);

end