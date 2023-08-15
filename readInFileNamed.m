function all_r2scores=readInFileNamed(directories,addtodirectories,filename)

% Example cell array of directories
% directories = {'/path/to/dir1', '/path/to/dir2', ...}; % replace with your actual paths

% Initialize an empty row vector to store the concatenated results
all_r2scores = [];

% Iterate through each directory
for i = 1:length(directories)
    dirPath = directories{i};
    
    % Check if r2scores.mat exists in the current directory
    matFilePath = fullfile(dirPath, addtodirectories, [filename '.mat']);
    
    if exist(matFilePath, 'file')
        % Load the content of r2scores.mat
        loadedData = load(matFilePath);
        
        % Assuming the loaded content contains a variable named r2scores
        % If the variable inside r2scores.mat has a different name, replace 'r2scores' with that name
        if isfield(loadedData, filename)
            temp=loadedData.(filename);
            all_r2scores = [all_r2scores, temp(:)']; % Ensuring it's a row vector and concatenating
        else
            warning(['Variable ' filename ' not found in %s'], matFilePath);
        end
    else
        warning('r2scores.mat not found in %s', dirPath);
    end
end

% At this point, all_r2scores contains the concatenated values from all r2scores.mat files
