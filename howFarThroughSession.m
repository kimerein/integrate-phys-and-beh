function [metadata,fractionThroughSess]=howFarThroughSession(metadata,excludeNoTouchTrialsFlag,trialTypes)

theseSess=unique(metadata.sessid);

metadata.fractionThroughSess=nan(size(metadata.sessid));
metadata.fractionThroughSess_adjusted=nan(size(metadata.sessid));
for i=1:length(theseSess)
    currSess=theseSess(i);
    if excludeNoTouchTrialsFlag==false
        nTrialsInSess=sum(metadata.sessid==currSess);
        metadata.fractionThroughSess(metadata.sessid==currSess)=1:nTrialsInSess;
        metadata.fractionThroughSess(metadata.sessid==currSess)=metadata.fractionThroughSess(metadata.sessid==currSess)./nTrialsInSess;
    else
        % session begins on the first trial when mouse touches pellet
        % session ends on the last trial when mouse touches pellet
        touches=trialTypes.touched_pellet(metadata.sessid==currSess);
        nTrialsInSess=sum(metadata.sessid==currSess);
        temp=1:nTrialsInSess;
        temp=temp./nTrialsInSess;
        % subtract off first trial with touch
        ftouch=temp(find(touches==true,1,'first'));
        temp=temp-ftouch;
        % rescale using last trial with touch
        ltouch=temp(find(touches==true,1,'last'));
        temp=temp./ltouch;
        metadata.fractionThroughSess_adjusted(metadata.sessid==currSess)=temp;
    end
end

if excludeNoTouchTrialsFlag==false
    fractionThroughSess=metadata.fractionThroughSess;
    metadata=rmfield(metadata,'fractionThroughSess_adjusted');
else
    fractionThroughSess=metadata.fractionThroughSess_adjusted;
    metadata=rmfield(metadata,'fractionThroughSess');
end