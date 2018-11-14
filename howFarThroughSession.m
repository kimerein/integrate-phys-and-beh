function [metadata,fractionThroughSess]=howFarThroughSession(metadata)

theseSess=unique(metadata.sessid);

metadata.fractionThroughSess=nan(size(metadata.sessid));
for i=1:length(theseSess)
    currSess=theseSess(i);
    nTrialsInSess=sum(metadata.sessid==currSess);
    metadata.fractionThroughSess(metadata.sessid==currSess)=1:nTrialsInSess;
    metadata.fractionThroughSess(metadata.sessid==currSess)=metadata.fractionThroughSess(metadata.sessid==currSess)./nTrialsInSess;
end

fractionThroughSess=metadata.fractionThroughSess;