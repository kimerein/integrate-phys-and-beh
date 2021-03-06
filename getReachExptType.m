function [reachExptType,dbs_procData]=getReachExptType(continuingAnalysisDir,dbs_procData,behLog,writeToDir)

reachExptType=false(1,length(dbs_procData.vids_to_match_mouseIDs));
for i=1:length(dbs_procData.vids_to_match_mouseIDs)
    currVid=dbs_procData.vids_to_match_mouseIDs{i};
    indx=dbs_procData.indsInto_behLog(i);
    % Search for "absent" or "control"
    isControl=false;
    for j=1:length(behLog(indx,:))
        if ~isempty(regexpi(behLog(indx,j),'absent')) || ~isempty(regexpi(behLog(indx,j),'control'))
            isControl=true;
            break
        end
    end
    reachExptType(i)=isControl;
    if writeToDir==true && length(continuingAnalysisDir)==length(dbs_procData.vids_to_match_mouseIDs)
        if isControl==true
            save([continuingAnalysisDir{i} '\CONTROL_SAYS_BEH_LOG.mat']);
        else
            save([continuingAnalysisDir{i} '\REGULAR_SAYS_BEH_LOG.mat']);
        end
    end
end

dbs_procData.areControls=dbs_procData.areControls==1 && reachExptType==true;