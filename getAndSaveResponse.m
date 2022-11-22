function [D1tagged_cueResponse,D1orD2taggingExpt,putAlignPeakAt]=getAndSaveResponse(dd_more,getThese,settings,putAlignPeakAtInput)

[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns,D1orD2taggingExpt,~,~,~,fromWhichSess,fromWhichUnit,fromWhichSess_forTrials,fromWhichTrial,isEventInThisTrial]=alignToCompanion(dd_more,false,[],[],[],getThese,settings,[]);
[~,ma]=nanmax(nanmean(aligncomp_y,1));
temp=nanmean(aligncomp_x,1);
putAlignPeakAt=temp(ma);
D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D1tagged_cueResponse.aligncomp_x=aligncomp_x;
D1tagged_cueResponse.aligncomp_y=aligncomp_y;
D1tagged_cueResponse.excluded=excluded;
D1tagged_cueResponse.ns=ns;
D1tagged_cueResponse.fromWhichSess=fromWhichSess;
D1tagged_cueResponse.fromWhichUnit=fromWhichUnit;
D1tagged_cueResponse.fromWhichTrial=fromWhichTrial;
D1tagged_cueResponse.isEventInThisTrial=isEventInThisTrial;
D1tagged_cueResponse.fromWhichSess_forTrials=fromWhichSess_forTrials;
if ~isempty(putAlignPeakAtInput)
    D1tagged_cueResponse=putAlignmentAtTime(putAlignPeakAtInput,D1tagged_cueResponse);
else
    D1tagged_cueResponse=putAlignmentAtTime(putAlignPeakAt,D1tagged_cueResponse);
end

end

function D1tagged_cueResponse=putAlignmentAtTime(putAlignPeakAt,D1tagged_cueResponse)

for i=1:size(D1tagged_cueResponse.aligncomp_y,1)
    [~,ma]=nanmax(D1tagged_cueResponse.aligncomp_y(i,:));
    temp=D1tagged_cueResponse.aligncomp_x(i,:);
    currAlignPeak=temp(ma);
    shiftByTime=putAlignPeakAt-currAlignPeak;
    D1tagged_cueResponse.unitbyunit_x(i,:)=D1tagged_cueResponse.unitbyunit_x(i,:)+shiftByTime;
    D1tagged_cueResponse.aligncomp_x(i,:)=D1tagged_cueResponse.aligncomp_x(i,:)+shiftByTime;
end

end