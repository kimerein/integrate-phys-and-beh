function checkThatUnitExistsInAllFolders(data_loc_array)

checkTheseFolders={'cue','cue_followedby_success','cue_noReach','cued_failure','cued_failure_then_noReach','cued_success','uncued_failure','uncued_failure_then_noReach','uncued_reach','uncued_success'};
alignmentNames={'cueAligned','cueFollowedBySuccess','cueNoReach','cuedFailure','cuedFailureThenNoReach','cuedSuccess','uncuedFailure','uncuedFailureThenNoReach','uncuedReach','uncuedSuccess'};

if length(checkTheseFolders)~=length(alignmentNames)
    error('checkTheseFolders must have same length as alignmentNames in checkThatUnitExistsInAllFolders');
end

% make empty unit
dataout.x=[];
dataout.y=[];
alignComp.x=[];
alignComp.y=[];
phys_timepointsComp.x=[];
phys_timepointsComp.y=[];
for indIntoDataLoc=1:size(data_loc_array,1)
    dd=dir(data_loc_array{indIntoDataLoc,8});
    disp(['Processing ' data_loc_array{indIntoDataLoc,8}]);
    namesOfUnits={};
    countUnits=1;
    firstFolder=true;
    for i=1:length(dd)
        if ismember(dd(i).name,checkTheseFolders)
            subdd=dir([data_loc_array{indIntoDataLoc,8} sep dd(i).name]);
            for j=1:length(subdd)
                r=regexp(subdd(j).name,'unit');
                if ~isempty(r)
                    runder=regexp(subdd(j).name,'_');
                    getUnitName=subdd(j).name(1:runder-1);
                    if ~ismember(getUnitName,namesOfUnits)
                        if firstFolder==true
                            % add
                        else
                            % problem
                            % fill in empty unit
                            disp(['Uh oh ' getUnitName ' is missing from ' data_loc_array{indIntoDataLoc,8} sep dd(i).name ' -- filling in empty']);
                            pause;
                            unitnameWithTag=subdd(j).name(1:runder(end));
                            fillInAllFolders(data_loc_array,dd,unitnameWithTag,checkTheseFolders,alignmentNames,dataout,alignComp,phys_timepointsComp);
                        end
                        namesOfUnits{countUnits}=getUnitName;
                        countUnits=countUnits+1;
                    end
                end
            end
            firstFolder=false;
        end
    end
end

end

function fillInAllFolders(data_loc_array,dd,unitnameWithTag,checkTheseFolders,alignmentNames,dataout,alignComp,phys_timepointsComp)

for i=1:length(dd)
    if ismember(dd(i).name,checkTheseFolders)
        nameOfFile=[data_loc_array{indIntoDataLoc,8} sep dd(i).name sep unitnameWithTag alignmentNames{i} '.mat'];
        if ~exist(nameOfFile,'file')
            save(nameOfFile,'dataout','alignComp','phys_timepointsComp');
        end
    end
end

end
