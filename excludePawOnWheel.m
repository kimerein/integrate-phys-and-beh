function [alltbt,metadata,out]=excludePawOnWheel(alltbt,metadata,out,cueZone)

% make a field in metadata where nans indicate which trials to exclude
% if paw on wheel before cue
[~,ma]=nanmax(nanmean(alltbt.(cueZone),1));
exclu=any(alltbt.reachBatch_all_pawOnWheel(:,1:ma-1),2) | any(alltbt.pawOnWheel(:,1:ma-1),2); % | any(alltbt.reachStarts(:,ma-29:ma-1),2);
metadata.metadata_nanField=zeros(size(metadata.sessid));
metadata.metadata_nanField(exclu)=nan;
[alltbt,metadata,out]=excludeTrials(alltbt,metadata,out,'metadata_nanField');
% remove filter field from metadata
metadata=rmfield(metadata,'metadata_nanField');

end