function settings=settingsForStriatumPhotometryPlots()

% settings for alignToCompanion
settings.excludeHigherFR=false;
% will only use if row above is true
settings.excludeAboveFR=[]; % in spikes / sec, exclude units with average firing rate above this
settings.cutAtTime=3; % stop plotting this many seconds after max of alignment companion
settings.ds=1; % photo data was initially binned in 10 ms bins, further downsample by this integer
settings.normalizePSTHs=false;
settings.suppressPlots=true;
settings.padsize=1000;
settings.testForAlignment=false;
settings.unitbaseline=300; %150;
settings.maxUnitsPerSess=30; % can't get more than this many units per session
settings.keepAllSingleTrials=true;
settings.isPhotometry=true;

% generally only use the following for units, not photometry
settings.discardTrialsWhereOptoDuringCue=false; % never did photometry and silencing of pDMSt at the same time
settings.onlyTrialsWhereOptoDuringCue=false;
settings.discardTrialsIfAnyOpto=false;
settings.discardDrops=false;
settings.useSameTrainingSetForAllNeurons=false;
settings.useTestSet=false;
settings.makeTrainingSet=false;
settings.useTrainingSet=false;