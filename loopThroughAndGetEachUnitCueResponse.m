function loopThroughAndGetEachUnitCueResponse(dirLocation)

files=findFiles(dirLocation);
% Loop through each spikes
for i=1:length(files)
    currname=files(i).name;
    a=load(fullfile(dirLocation,currname));
    [post_spikes,auxData]=convertSpikesToTrials_and_saveAuxChs(a.spikes,fullfile(dirLocation,'auxData193.mat'));
    % find good units
end

end

function files=findFiles(directoryPath)

% Get list of all "_sorted.mat" files in the directory
files = dir(fullfile(directoryPath, '*_sorted.mat'));

end