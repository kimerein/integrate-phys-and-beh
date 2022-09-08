function scriptToOrganizeD1vD2unitResponses(dd)

% load dd

% choose type of response and time window to analyze
responseType='uncued_failure';
timeWindow=[-0.22 0.5]; %[-0.25 0.3]; % relative to alignment companion onset, in seconds
responseBaseline=[]; %[-1.05 -0.25];
cueWindow=[0 0.5];
beforeCueBaseline=[-1.05 0.2];

% get cue response of each tagged type
dd_more=cell(1,length(dd));
for i=1:length(dd)
    dd_more{i}=[dd{i} '\cue'];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D1tagged_cueResponse.aligncomp_x=aligncomp_x;
D1tagged_cueResponse.aligncomp_y=aligncomp_y;
D1tagged_cueResponse.excluded=excluded;
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
D2tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D2tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D2tagged_cueResponse.aligncomp_x=aligncomp_x;
D2tagged_cueResponse.aligncomp_y=aligncomp_y;
D2tagged_cueResponse.excluded=excluded;

% get response of each tagged type
for i=1:length(dd)
    dd_more{i}=[dd{i} '\' responseType];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
activityD1tagged.unitbyunit_x=unitbyunit_x;
activityD1tagged.unitbyunit_y=unitbyunit_y;
activityD1tagged.aligncomp_x=aligncomp_x;
activityD1tagged.aligncomp_y=aligncomp_y;
activityD1tagged.excluded=excluded;
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
activityD2tagged.unitbyunit_x=unitbyunit_x;
activityD2tagged.unitbyunit_y=unitbyunit_y;
activityD2tagged.aligncomp_x=aligncomp_x;
activityD2tagged.aligncomp_y=aligncomp_y;
activityD2tagged.excluded=excluded;

[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);

scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindow, responseBaseline,cueWindow,beforeCueBaseline);

end

function [D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged)

excluded1=D1tagged_cueResponse.excluded; 
excluded2=activityD1tagged.excluded;
eitherExcluded=excluded1==1 | excluded2==1;
D1tagged_cueResponse.unitbyunit_x=D1tagged_cueResponse.unitbyunit_x(eitherExcluded==0,:);
D1tagged_cueResponse.unitbyunit_y=D1tagged_cueResponse.unitbyunit_y(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_x=D1tagged_cueResponse.aligncomp_x(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_y=D1tagged_cueResponse.aligncomp_y(eitherExcluded==0,:);
D1tagged_cueResponse.excluded=D1tagged_cueResponse.excluded(eitherExcluded==0);

activityD1tagged.unitbyunit_x=activityD1tagged.unitbyunit_x(eitherExcluded==0,:);
activityD1tagged.unitbyunit_y=activityD1tagged.unitbyunit_y(eitherExcluded==0,:);
activityD1tagged.aligncomp_x=activityD1tagged.aligncomp_x(eitherExcluded==0,:);
activityD1tagged.aligncomp_y=activityD1tagged.aligncomp_y(eitherExcluded==0,:);
activityD1tagged.excluded=activityD1tagged.excluded(eitherExcluded==0);

end
