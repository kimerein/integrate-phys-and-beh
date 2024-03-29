function [spikes sortvector] = filtspikes_without_sweeps(spikes,bOR,varargin)

if nargin < 3 || ( (nargin < 4) && ~iscell(varargin{1}) )
    error('Not enough arguments supplied')
end

if(~iscell(varargin{1}))
    numSortFields = length(varargin)/2;
    emptyspots = zeros(numSortFields);
    for i = 1:numSortFields
        sortfield(i) = {varargin{(i-1)*2+1}};
        sortvalue{i} = varargin{2*i};
        emptyspots(i) = ~isempty(sortfield{i});
    end
else
    numSortFields = length(varargin{1})/2;
    emptyspots = zeros(numSortFields);
    thecell = varargin{1};
    for i = 1:numSortFields
        sortfield{i} = thecell{(i-1)*2+1};
        sortvalue{i} = thecell{2*i};
        emptyspots(i) = ~isempty(sortfield{i});
    end
end
sortfield = sortfield(logical(emptyspots));
sortvalue = sortvalue(logical(emptyspots));
numSortFields = length(sortfield);
% Make logical vector to sort spikes
for i = 1:numSortFields
    if(isa(sortvalue{i}, 'function_handle'))
        fcn_hand = sortvalue{i};
        tempvector = arrayfun(fcn_hand, spikes.(sortfield{i}));
    else    
        tempvector = ismember(spikes.(sortfield{i}),sortvalue{i});     % Take trials = value
    end
        
    % Combine multiple conditions
    if (numSortFields > 1) && (i > 1)
        if bOR  % OR
            sortvector = sortvector | tempvector;
        else    % AND
            sortvector = sortvector & tempvector;
        end
    else
        sortvector = tempvector;
    end
end

% Find fields that have same length as spiketimes
reqLength = length(spikes.spiketimes);  
fieldList = fieldnames(spikes);

% Fix case where reqLength == 1
if reqLength == 1
    rmFields = {'params','info','sweeps','labels'};
    for i = 1:length(rmFields)
        fieldList(strcmp(rmFields{i},fieldList)) = '';
    end
end

% Use logical vector to extract spiketimes, trials, etc.
for i = 1:length(fieldList)
    if ismember(reqLength,size(spikes.(fieldList{i})));
        switch fieldList{i}
            case 'waveforms'
                spikes.(fieldList{i}) = spikes.(fieldList{i})(sortvector,:,:);
            otherwise
                spikes.(fieldList{i}) = spikes.(fieldList{i})(sortvector);      % Using dynamic field names
        end
    end
end