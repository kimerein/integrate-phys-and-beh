function [out,dontUseTrials]=doUnitTest(foldername, unitname)

temp=getCriteriaForUnitsToPlot();
if ischar(temp)
    if strcmp(temp,'alltrue')
        out=true;
        dontUseTrials=[];
        return
    end
end

% find details for this unit
runit=regexp(unitname,'_');
unitOnlyName=unitname(1:runit(1)-1);
r=regexp(foldername,sep);
a=load([foldername(1:r(end)) 'unit_details' sep unitOnlyName '_unitdets.mat']);
dontUseTrials=a.unitdets.dontUseTrials;
undets=[a.unitdets.inStructure a.unitdets.isFS a.unitdets.isTAN a.unitdets.isSPN a.unitdets.isLowFRThin];

criteria=getCriteriaForUnitsToPlot();
criteria(criteria==-100)=undets(criteria==-100); % -100 is a wildcard
% [inStructure isFS isTAN isSPN isLowFRThin]
out=isequal(undets,criteria);

end