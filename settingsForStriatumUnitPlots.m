function settings=settingsForStriatumUnitPlots()

% settings for alignToCompanion
settings.excludeHigherFR=false;
% will only use if row above is true
settings.excludeAboveFR=4; % in spikes / sec, exclude units with average firing rate above this
settings.cutAtTime=3; % stop plotting this many seconds after max of alignment companion
settings.ds=1; % spiking data was initially binned in 10 ms bins, further downsample by this integer
settings.normalizePSTHs=false;
settings.suppressPlots=true;
settings.padsize=1000;
settings.testForAlignment=false;
settings.unitbaseline=300; %150;
settings.maxUnitsPerSess=300; % can't get more than this many units per session
settings.keepAllSingleTrials=true;