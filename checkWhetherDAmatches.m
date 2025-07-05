function checkWhetherDAmatches(rootDir)
% checkWhetherDAmatches  Recursively checks tbt.mat files for mismatched field lengths
%
%   checkWhetherDAmatches(rootDir) searches all subfolders of rootDir for
%   files named "tbt.mat". For each tbt.mat, it loads the struct variable
%   tbt and compares the number of rows in tbt.cue to the number of rows in
%   any of the following fields (if they exist):
%       • cued_successChunks
%       • cued_failureChunks
%       • cued_successChunks_pDMSt
%       • cued_failureChunks_pDMSt
%
%   If any of those fields exists and its first dimension length does not
%   match size(tbt.cue,1), then the function prints the path to that tbt.mat
%   and looks in the same folder for a file named "humanchecked.txt" and
%   renames it to "fieldlengthsdontmatch.txt".
%
%   INPUT:
%     rootDir  – Root directory to begin the recursive search. If omitted,
%                uses the current working directory.
%
%   USAGE EXAMPLES:
%     % 1. Use the current folder:
%     checkWhetherDAmatches();
%
%     % 2. Specify a folder to search:
%     checkWhetherDAmatches('C:\MyDataFolder');
%

    if nargin < 1 || isempty(rootDir)
        rootDir = pwd;
    end

    % Find all "tbt.mat" files under rootDir (recursive)
    tbtFiles = dir(fullfile(rootDir, '**', 'tbt.mat'));

    % Define the list of DA‐related fields to check
    daFields = {...
        'cued_successChunks', ...
        'cued_failureChunks', ...
        'cued_successChunks_pDMSt', ...
        'cued_failureChunks_pDMSt' ...
    };

    % Loop through each found tbt.mat
    for k = 1:numel(tbtFiles)
        tbtPath = fullfile(tbtFiles(k).folder, tbtFiles(k).name);
        try
            % Load the variable 'tbt' from the .mat file
            data = load(tbtPath, 'tbt');
            if ~isfield(data, 'tbt')
                continue;
            end
            tbt = data.tbt;
        catch
            continue;
        end

        % Ensure the 'cue' field exists
        if ~isfield(tbt, 'cue')
            continue;
        end
        nCue = size(tbt.cue, 1);

        % Check each DA field for a row‐count mismatch
        mismatchFound = false;
        for iField = 1:numel(daFields)
            fieldName = daFields{iField};
            if isfield(tbt, fieldName)
                nField = size(tbt.(fieldName), 1);
                if nField ~= nCue
                    mismatchFound = true;
                    break;
                end
            end
        end

        % If any mismatch was found, print the tbt.mat location and rename file
        if mismatchFound
            fprintf('Mismatch found in: %s\n', tbtPath);
            folderPath = tbtFiles(k).folder;
            oldTxt = fullfile(folderPath, 'humanchecked.txt');
            newTxt = fullfile(folderPath, 'fieldlengthsdontmatch.txt');
            if isfile(oldTxt)
                try
                    movefile(oldTxt, newTxt);
                catch
                    warning('Could not rename "%s". Check file permissions.', oldTxt);
                end
            end
        end
    end
end
