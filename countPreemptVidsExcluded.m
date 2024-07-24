function countPreemptVidsExcluded(mainDir)

% Get a list of all subfolders
allSubFolders = genpath(mainDir);
remain = allSubFolders;
listOfFolderNames = {};

while true
    [singleSubFolder, remain] = strtok(remain, pathsep);
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder]; 
end

% Loop through each subfolder
numberPreemptVids=0;
totalNumberVids=0;
preemptCueValue=0;

% Initialize cell arrays to keep track of date format file names and preemptCue values
dateFileNames = {};
preemptCueValues = [];

for k = 1 : length(listOfFolderNames)
    thisFolder = listOfFolderNames{k};
    
    % Check for the preemptCue.mat file
    preemptCueFile = fullfile(thisFolder, 'preemptCue.mat');
    if exist(preemptCueFile, 'file')
        fprintf('Found preemptCue.mat in folder: %s\n', thisFolder);
        % Read the file here (example using textread)
        data = load(preemptCueFile); % Adjust according to file format
        preemptCueValue = data.preemptCue;
        totalNumberVids=totalNumberVids+1;
    else
        warning('preemptCue.mat missing!');
        pause;
    end
    numberPreemptVids=numberPreemptVids+preemptCueValue;
    % Get list of all files in the current folder
    filePattern = fullfile(thisFolder, '*.txt');
    txtFiles = dir(filePattern);
    
    % Loop through text files to find the one with a date format name
    for f = 1 : length(txtFiles)
        [~, fileName, ext] = fileparts(txtFiles(f).name);
        if length(fileName) == 8 && all(isstrprop(fileName, 'digit')) && strcmp(ext, '.txt')
            fprintf('Found date formatted file: %s in folder: %s\n', txtFiles(f).name, thisFolder);
            % Add the file name and preemptCue value to the lists
            dateFileNames{end+1} = fileName;
            preemptCueValues(end+1) = preemptCueValue; 
        end
    end
end

% Count the number of unique date file names
uniqueDateFileNames = unique(dateFileNames);
numUniqueDateFileNames = length(uniqueDateFileNames);

% Initialize a counter for days with preemptCue true
daysWithPreemptCueTrue = 0;

% Loop through unique dates to check preemptCue values
for i = 1:length(uniqueDateFileNames)
    date = uniqueDateFileNames{i};
    % Find indices of the current date in the list
    indices = find(strcmp(dateFileNames, date));
    % Check if all preemptCue values for this date are 1
    if all(preemptCueValues(indices) == 1)
        daysWithPreemptCueTrue = daysWithPreemptCueTrue + 1;
    end
end

% Display the results
fprintf('Total number of unique days: %d\n', numUniqueDateFileNames);
fprintf('Number of days with preemptCue true: %d\n', daysWithPreemptCueTrue);
fprintf('Total number of unique VIDEOS: %d\n', length(preemptCueValues));
fprintf('Number of VIDEOS with preemptCue true: %d\n', sum(preemptCueValues,'all','omitnan'));

end