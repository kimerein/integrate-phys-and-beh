function [cuetuningindex,isCueTuned]=getCueTunedUnits(cueR,uncuedReachR)

ds=1;
smoo=80;
cueStartCutoff=-0.25;

% match which units and get responses unit by unit
r=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cueR,uncuedReachR,'meanAcrossUnits',ds,true);
cueresp=r.response1; uncuedreachresp=r.response2;
if ~isempty(smoo)
    for i=1:size(cueresp.unitbyunit,1)
        temp=cueresp.unitbyunit(i,:);
        cueresp.unitbyunit(i,:)=smoothdata(temp,'gaussian',smoo);
    end
    for i=1:size(uncuedreachresp.unitbyunit,1)
        temp=uncuedreachresp.unitbyunit(i,:);
        uncuedreachresp.unitbyunit(i,:)=smoothdata(temp,'gaussian',smoo);
    end
end

% get cue response vs response during uncued reach
t=cueresp.unittimes;
beforecue=max(cueresp.unitbyunit(:,t>-0.9 & t<cueStartCutoff),[],'omitnan'); 
duringcue=max(cueresp.unitbyunit(:,t>cueStartCutoff & t<0.8),[],'omitnan'); 
cuer=duringcue-beforecue;

t=cueresp.unittimes;
beforereach=max(uncuedreachresp.unitbyunit(:,t<-1),[],'omitnan'); 
duringreach=max(uncuedreachresp.unitbyunit(:,t>-1 & t<1),[],'omitnan'); 
reachr=duringreach-beforereach;

% define cue tuned units using cue tuning index
cuetuningindex=(cuer-reachr)/(cuer+reachr);
isCueTuned=cuetuningindex>0.5;

end