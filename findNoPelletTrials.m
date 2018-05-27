function noPelletTrials=findNoPelletTrials(tbt)

if isfield(tbt,'pelletPresent')
    noPelletTrials=all(~(tbt.pelletPresent>0),2);
else
    % find trials where all reaches in trial have no pellet
    noPelletTrials=~any(tbt.reachStarts_pelletPresent,2);
end