function tbt = extractEventChunks(dataout, alignComp, tbt, n, pre)
% extractEventChunks  Aligns and extracts fixed-length data chunks around events
%
%   tbt = extractEventChunks(dataout, alignComp, tbt, n, pre) searches each
%   trial (row) of alignComp.y for the first time the event turns “ON” (>0.05).
%   It then finds the matching timepoint in dataout.x and extracts exactly n
%   samples from dataout.y, starting pre indices before the event index.
%   The extracted chunks are stored in a matrix tbt.chunks of size [M × n],
%   where M is the number of trials. Trials with no event receive a row of
%   all NaNs.
%
%   INPUTS:
%     dataout   – struct with fields:
%                  • x : 1×T vector of timepoints
%                  • y : M×T matrix (M trials, T timepoints)
%     alignComp – struct with fields:
%                  • x : 1×U vector of timepoints (length U)
%                  • y : M×U matrix indicating event signal (NaN, 0, or >0.05)
%     tbt       – (possibly empty) struct.  This function will add/overwrite
%                  a field called “chunks” to tbt.
%     n         – scalar integer, total number of indices to extract per trial
%     pre       – scalar integer, number of indices BEFORE the event index to
%                  begin the chunk
%
%   OUTPUT:
%     tbt.chunks – an M×n matrix. For trial i:
%                    • If alignComp.y(i,:) is all NaN → row i is all NaNs
%                    • Otherwise, chunk begins pre samples before the event
%                      and spans n samples.  Out-of-bounds segments are
%                      padded with NaNs.
%
%   EXAMPLE USAGE:
%     % To extract 50 samples total, starting 10 indices before each event:
%     nSamples = 50;
%     preSamples = 10;
%     tbt = extractEventChunks(dataout, alignComp, struct(), nSamples, preSamples);
%

    % Number of trials (rows in dataout.y)
    [nTrials, T] = size(dataout.y);

    % Preallocate the chunks matrix with NaNs
    tbt.chunks = nan(nTrials, n);

    % Loop over each trial
    for i = 1:nTrials
        rowAlign = alignComp.y(i, :);

        % If all NaN → no event, leave row i as all NaNs
        if all(isnan(rowAlign))
            continue
        end

        % 1) Find the first index where the event is ON (>0.05)
        idxEventAlign = find(rowAlign > 0.05, 1, 'first');
        if isempty(idxEventAlign)
            % No ON event detected → leave row i as NaNs
            continue
        end

        % 2) Get the event time from alignComp.x
        timeEvent = alignComp.x(idxEventAlign);

        % 3) Find the closest index in dataout.x to this event time
        [~, idxEventData] = min(abs(dataout.x - timeEvent));

        % 4) Determine the start and end indices in dataout for the chunk
        startIdx = idxEventData - pre;
        endIdx = startIdx + n - 1;

        % Build an NaN-filled row of length n
        chunk = nan(1, n);

        % Determine which data indices fall within [1, T]
        dataInds = startIdx : endIdx;
        validMask = (dataInds >= 1) & (dataInds <= T);
        validDataInds = dataInds(validMask);

        % Corresponding positions within chunk
        chunkPositions = find(validMask);

        % Insert valid data into the chunk
        if ~isempty(validDataInds)
            chunk(chunkPositions) = dataout.y(i, validDataInds);
        end

        % 5 & 6) Store the chunk in row i of tbt.chunks
        tbt.chunks(i, :) = chunk;

        % Also return max and min of each chunk
    end
end
