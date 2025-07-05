function processDAandTbtFiles(rootDir, dataBaseDir, branchNames, n, pre)
% processDAandTbtFiles  Recursively finds "tbt.mat" files and processes them
%   This version allows multiple branch names (e.g., 'cued_success',
%   'cued_failure', etc.) to be specified in a cell array.  Steps 5 & 6
%   from the original procedure are handled by the helper subfunction
%   processBranch.
%
%   USAGE:
%     branchNames = {'cued_success', 'cued_failure', 'other_branch', ...};
%     processDAandTbtFiles('./BehaviorData', './ProcessedData', branchNames, 100, 10);
%
%   INPUTS:
%     rootDir      – Root directory in which to search recursively for
%                    "tbt.mat".
%     dataBaseDir  – Higher‐level directory where the corresponding folders
%                    [datestring mousestring] live.
%     branchNames  – Cell array of strings, e.g. {'cued_success', 'cued_failure'}.
%                    Each entry corresponds to a subfolder under
%                    "photometry aligned to behNAcc".
%     n            – Number of indices to extract (total window length).
%     pre          – Number of indices before the event to begin the chunk.
%
%   For each "tbt.mat" found:
%     1) Load existing tbt struct.
%     2) In the same folder, find a .txt file whose name is all digits → datestring.
%     3) The parent folder name → mousestring.
%     4) Find corresponding folder = fullfile(dataBaseDir, [datestring mousestring]).
%     5) For each branchName in branchNames, call processBranch to:
%        • Navigate to photometry aligned to behNAcc/branchName/
%        • Load ch_<camelCase(branchName)>.mat, extract dataout & alignComp
%        • Run extractEventChunks(dataout,alignComp,origTbt,n,pre)
%        • Return a matrix of chunks
%     6) Store each branch's chunks in updatedTbt.(sprintf('%sChunks', branchName))
%     7) Save updatedTbt back into "tbt.mat".
%
%   NOTE: Requires extractEventChunks to be on the MATLAB path.

    % 1) Find all "tbt.mat" files under rootDir (recursive)
    fileList = dir(fullfile(rootDir, '**', 'tbt.mat'));
    
    for k = 1:numel(fileList)
        tbtPath = fullfile(fileList(k).folder, fileList(k).name);
        
        % Load existing tbt struct
        data = load(tbtPath, 'tbt');
        if ~isfield(data, 'tbt')
            warning('No variable "tbt" found in %s. Skipping.', tbtPath);
            continue;
        end
        origTbt = data.tbt; 
        
        % 2) Find .txt file with all digits in same folder as tbt.mat
        tbtFolder = fileList(k).folder;
        txtFiles = dir(fullfile(tbtFolder, '*.txt'));
        datestring = '';
        for i = 1:numel(txtFiles)
            [~, nameNoExt] = fileparts(txtFiles(i).name);
            if ~isempty(regexp(nameNoExt, '^\d+$', 'once'))
                datestring = nameNoExt;
                break;
            end
        end
        if isempty(datestring)
            warning('No all-numeric .txt file found in %s. Skipping.', tbtFolder);
            continue;
        end
        
        % 3) Determine mousestring as the parent-of-parent folder of tbtFolder
        parentDir = fileparts(tbtFolder);                 % parent of tbtFolder
        [~, mousestring] = fileparts(parentDir);           % parent of parentDir
        
        if isempty(mousestring)
            warning('Cannot determine mouse folder for %s. Skipping.', tbtFolder);
            continue;
        end
        
        % 4) Construct correspondingfolder path under dataBaseDir
        correspondingFolder = fullfile(dataBaseDir, datestring, mousestring); 
        if ~isfolder(correspondingFolder)
            warning('Corresponding folder "%s" does not exist. Skipping.', correspondingFolder);
            continue;
        end
        
        % Prepare to store updated fields in the new tbt struct
        updatedTbt = origTbt;
        
        % Base photometry directory
        photometryDir = fullfile(correspondingFolder, 'photometry aligned to behNAcc');
        
        % 5) Loop over each branchName and process
        for b = 1:numel(branchNames)
            branch = branchNames{b};
            try
                chunks = processBranch(branch, photometryDir, origTbt, n, pre);
                fieldName = sprintf('%sChunks', branch);
                updatedTbt.(fieldName) = chunks;
            catch ME
                warning('Error processing branch "%s" for "%s": %s', ...
                        branch, tbtPath, ME.message);
            end
        end

        % Base photometry directory
        photometryDir = fullfile(correspondingFolder, 'photometry aligned to beh');
        
        % 6) Loop over each branchName and process
        for b = 1:numel(branchNames)
            branch = branchNames{b};
            try
                chunks = processBranch(branch, photometryDir, origTbt, n, pre);
                fieldName = sprintf('%sChunks_pDMSt', branch);
                updatedTbt.(fieldName) = chunks;
            catch ME
                warning('Error processing branch "%s" for "%s": %s', ...
                        branch, tbtPath, ME.message);
            end
        end
        
        % 7) Save the updated tbt struct back to tbt.mat (overwrites)
        tbt = updatedTbt; 
        save(tbtPath, 'tbt');
        
        fprintf('Processed and updated "%s".\n', tbtPath);
    end
end

%% Subfunction: processBranch
function chunks = processBranch(branchName, photometryDir, origTbt, n, pre)
% processBranch  Handles steps 5 & 6 for a single branchName
%
%   INPUTS:
%     branchName    – e.g. 'cued_success' or 'cued_failure' or any folder under
%                     photometryDir.
%     photometryDir – full path to "photometry aligned to behNAcc" directory.
%     origTbt       – the original tbt struct loaded from file.
%     n             – number of samples to extract.
%     pre           – number of samples before event index.
%
%   OUTPUT:
%     chunks        – an M×n matrix of extracted data chunks, where M is
%                     the number of trials in origTbt (or origTbt.y length).

    % Construct branch-specific folder
    branchDir = fullfile(photometryDir, branchName);
    if ~isfolder(branchDir)
        warning('Branch folder "%s" does not exist. Returning empty.', branchDir);
        chunks = [];
        return;
    end
    
    % Convert branchName (e.g., 'cued_success') to camelCase ('cuedSuccess')
    camelName = toCamel(branchName);
    
    % Construct .mat filename: ch_<camelName>.mat
    matFileName = ['ch__' camelName '.mat'];
    matPath = fullfile(branchDir, matFileName);
    
    if ~isfile(matPath)
        warning('File "%s" not found under "%s". Returning empty.', matFileName, branchDir);
        chunks = [];
        return;
    end
    
    % Load dataout & alignComp
    tmp = load(matPath, 'dataout', 'alignComp');
    if ~isfield(tmp, 'dataout') || ~isfield(tmp, 'alignComp')
        warning('Variables "dataout" or "alignComp" missing in %s. Returning empty.', matPath);
        chunks = [];
        return;
    end
    
    % Call extractEventChunks using origTbt as the base
    resultTbt = extractEventChunks(tmp.dataout, tmp.alignComp, origTbt, n, pre);
    chunks = resultTbt.chunks;
end

%% Subfunction: toCamel
function camel = toCamel(snakeStr)
% toCamel  Converts 'snake_case' to 'camelCase'
    parts = split(snakeStr, '_');
    camel = parts{1};
    for i = 2:numel(parts)
        token = parts{i};
        if ~isempty(token)
            camel = [camel, upper(token(1)), token(2:end)]; %#ok<AGROW>
        end
    end
end
