function scriptToOrganizeD1vD2unitResponses_wrapper(dd)

% settings for scriptToOrganizeD1vD2unitResponses
settings.pvalcutoff=[-0.001 0.2];
settings.responseType='uncued_success';
settings.timeWindow=[-0.22 3.32]; % relative to alignment companion onset, in seconds
settings.responseBaseline=[]; %[-1.05 -0.25]; % empty means don't subtract off baseline
settings.cueWindow=[0 0.5];
settings.beforeCueBaseline=[-1.05 -0.2];
settings.skipCueAlignment=false;
% settings for alignToCompanion
settings.excludeHigherFR=false;
settings.excludeAboveFR=4; % in spikes / sec, exclude units with average firing rate above this 
settings.cutAtTime=3; % stop plotting this many seconds after max of alignment companion
settings.ds=6; % spiking data was initially binned in 10 ms bins, further downsample by this integer
settings.normalizePSTHs=false;
settings.suppressPlots=true;
settings.padsize=1000;
settings.testForAlignment=false;
settings.unitbaseline=300; %150;
settings.beforeOptoBaseline=140;
settings.optoTagDuration=0.2; % in seconds

[D1tagged_cueResponse,D2tagged_cueResponse,cued_success_D1tagged,cued_success_D2tagged]=scriptToOrganizeD1vD2unitResponses(dd,settings);

end