function results = concatenateFilesByName(dirList, fileList, concatDim)
% CONCATFILESBYNAME  Load identically-named arrays from multiple folders
%                    and concatenate along the requested dimension.
%
%   results = CONCATFILESBYNAME(dirList, fileList)   
%     loads each fileList{k}+'.mat' from every dir in dirList and
%     concatenates vertically (concatDim=1).
%
%   results = CONCATFILESBYNAME(dirList, fileList, concatDim)
%     where concatDim is 1 (vertical) or 2 (horizontal).
%
% Example:
%   dirs     = {'run1','run2','run3'};
%   files    = {'sessionA','sessionB'};
%   R        = concatFilesByName(dirs, files, 2);
%   % R.sessionA is [N x (M1+M2+M3)] if each run folder had a sessionA.mat
%

    if nargin<3, concatDim = 1; end
    assert(ismember(concatDim,[1 2]), 'concatDim must be 1 or 2');

    results = struct();
    for k = 1:numel(fileList)
        name = fileList{k};
        accum = [];
        for i = 1:numel(dirList)
            fn = fullfile(dirList{i}, [name '.mat']);
            if ~isfile(fn)
                warning('File not found: %s', fn);
                continue;
            end
            S = load(fn, name);
            A = S.(name);
            accum = cat(concatDim, accum, A);
        end
        results.(matlab.lang.makeValidName(name)) = accum;
    end
end
