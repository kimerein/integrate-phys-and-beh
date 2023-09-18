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
settings.unitbaseline=100; %600; %300; %150;
settings.maxUnitsPerSess=300; % can't get more than this many units per session
settings.keepAllSingleTrials=false;
settings.fracForTrainingSet=0.5; % fraction of trials to use for training set
settings.makeTrainingSet=false;
settings.useTestSet=true;
settings.useTheseTestSets={}; %{'cued_failureIntersect'}; %{'cued_failure','uncued_failure'};
settings.useTheseTestFilenames={}; %{'cuedFailure'}; %{'cuedFailure','uncuedFailure'};
settings.useTrainingSet=false;
settings.useSameTrainingSetForAllNeurons=true;
settings.discardDrops=false;
settings.dropFolderName='uncued_drop';
settings.dropFileName='uncuedDrop';
settings.discardTrialsWhereOptoDuringCue=true;
settings.onlyTrialsWhereOptoDuringCue=false;
settings.discardTrialsIfAnyOpto=false;
settings.minTrialLength=9.5; % in seconds

% basic checks
if settings.discardTrialsWhereOptoDuringCue==true && settings.onlyTrialsWhereOptoDuringCue==true
    error('Cannot both be true in settingsForStriatumUnitPlots.m: settings.discardTrialsWhereOptoDuringCue==true && settings.onlyTrialsWhereOptoDuringCue==true');
end