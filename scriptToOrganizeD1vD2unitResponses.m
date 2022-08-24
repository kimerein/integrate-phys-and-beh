% load dd

% choose type of response and time window to analyze
responseType='cued_success';
timeWindow=[-0.3 1.5]; % relative to alignment companion onset, in seconds

% get cue response of each tagged type
dd_more=cell(1,length(dd));
for i=1:length(dd)
    dd_more{i}=[dd{i} '\cue'];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D1tagged_cueResponse.aligncomp_x=aligncomp_x;
D1tagged_cueResponse.aligncomp_y=aligncomp_y;
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
D2tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D2tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D2tagged_cueResponse.aligncomp_x=aligncomp_x;
D2tagged_cueResponse.aligncomp_y=aligncomp_y;

% get response of each tagged type
for i=1:length(dd)
    dd_more{i}=[dd{i} '\' responseType];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
activityD1tagged.unitbyunit_x=unitbyunit_x;
activityD1tagged.unitbyunit_y=unitbyunit_y;
activityD1tagged.aligncomp_x=aligncomp_x;
activityD1tagged.aligncomp_y=aligncomp_y;
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
activityD2tagged.unitbyunit_x=unitbyunit_x;
activityD2tagged.unitbyunit_y=unitbyunit_y;
activityD2tagged.aligncomp_x=aligncomp_x;
activityD2tagged.aligncomp_y=aligncomp_y;

scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindow, [0 1.5], [-1 0]);
