function [alltbt,trialTypes,metadata]=findPelletMissingCues(alltbt,trialTypes,metadata)

[~,f]=nanmax(nanmean(alltbt.cueZone_onVoff,1));

alltbt.pelletMissingAtCue=zeros(size(alltbt.cue,1),1);
alltbt.pelletMissingAtCue=alltbt.pelletPresent(:,f)<0.5;
trialTypes.pelletMissingAtCue=alltbt.pelletPresent(:,f)<0.5;
metadata.pelletMissingAtCue=alltbt.pelletPresent(:,f)<0.5;

end