function applyFidgetToAll(dataRootDir, behaviorByDateDir)
% applyFidgetToAll  Recursively process mouse/day folders to add fidgets
%   dataRootDir: root directory containing mouse-named folders
%   behaviorByDateDir: base path to behavioral data by date (e.g., 'Z:/MICROSCOPE/Kim/KER Behavior/By date/Low speed')

% List all mouse folders
mouseDirs = dir(dataRootDir);
mouseDirs = mouseDirs([mouseDirs.isdir] & ~startsWith({mouseDirs.name}, '.'));

for m = 1:numel(mouseDirs)
    mouseName = mouseDirs(m).name;
    mousePath = fullfile(dataRootDir, mouseName);
    
    % List day folders inside each mouse folder
    dayDirs = dir(mousePath);
    dayDirs = dayDirs([dayDirs.isdir] & ~startsWith({dayDirs.name}, '.'));
    
    for d = 1:numel(dayDirs)
        dayPath = fullfile(mousePath, dayDirs(d).name);
        
        % Load tbt.mat
        tbtFile = fullfile(dayPath, 'tbt.mat');
        if ~isfile(tbtFile)
            warning('Missing tbt.mat in %s', dayPath);
            continue;
        end
        S = load(tbtFile, 'tbt');
        tbt = S.tbt; clear S;
        
        % Find date .txt file (format YYYMMDD.txt)
        txtFiles = dir(fullfile(dayPath, '*.txt'));
        if isempty(txtFiles)
            warning('No .txt date file in %s', dayPath);
            continue;
        end
        dateStr = txtFiles(1).name(1:8);
        
        % Construct corresponding behavior folder
        behavFolder = fullfile(behaviorByDateDir, dateStr, mouseName);
        if ~isfolder(behavFolder)
            warning('Behavior folder not found: %s', behavFolder);
            continue;
        end
        o2Folder = fullfile(behavFolder, 'O2 output');
        if ~isfolder(o2Folder)
            warning('Missing O2 output in %s', behavFolder);
            continue;
        end
        
        % Load fidget
        fidgetFile = dir(fullfile(o2Folder, '*_fidget.mat'));
        if isempty(fidgetFile)
            warning('No *_fidget.mat in %s', o2Folder);
            continue;
        end
        S = load(fullfile(o2Folder, fidgetFile(1).name), 'fidget');
        fidget = S.fidget; clear S;
        
        % Load settings
        settingsFile = dir(fullfile(o2Folder, '*_autoReachSettings.mat'));
        if isempty(settingsFile)
            warning('No *_autoReachSettings.mat in %s', o2Folder);
            continue;
        end
        S = load(fullfile(o2Folder, settingsFile(1).name), 'settings');
        settings = S.settings; clear S;
        
        % Backup original tbt
        backupFile = fullfile(dayPath, 'backup_before_fidget_tbt.mat');
        save(backupFile, 'tbt');
        
        % Add back fidgets and overwrite tbt.mat
        tbt = addBackFidgets(tbt, settings, fidget, 1);
        save(tbtFile, 'tbt');
        
        fprintf('Processed %s on %s\n', mouseName, dateStr);
    end
end
end
