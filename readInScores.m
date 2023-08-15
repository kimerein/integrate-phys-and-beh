function values=readInScores(directoryPath)

% Specify the directory
% directoryPath = '/path/to/your/directory'; % replace with your directory path

% Get a list of all .csv files in the directory
csvFiles = dir(fullfile(directoryPath, '*.csv'));

% Extract numbers from each filename
numbers = zeros(1, length(csvFiles));
for i = 1:length(csvFiles)
    name = csvFiles(i).name;
    numStr = regexp(name, 'neuron(\d+)_', 'tokens');
    
    % Check if we found a match
    if ~isempty(numStr)
        numbers(i) = str2double(numStr{1}{1});
    else
        warning('Filename format unexpected for: %s', name);
    end
end

% Sort files based on the extracted numbers
[~, idx] = sort(numbers);
csvFiles = csvFiles(idx);

% Initialize an array to store the values
values = zeros(1, length(csvFiles));

% Loop through each .csv file and extract the desired value
for i = 1:length(csvFiles)
    % Construct the full file path
    filePath = fullfile(directoryPath, csvFiles(i).name);
    
    % Read the .csv file
    data = readmatrix(filePath);
    
    % Extract the value from the fourth column and second row
    values(i) = data(2, 4);
end

% At this point, 'values' contains the values from the fourth column and second row of each .csv file

end