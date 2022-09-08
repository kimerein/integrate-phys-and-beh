function scriptToOrganizeD1vD2unitResponses(dd)

% load dd

pvalcutoff=1;

% choose type of response and time window to analyze
responseType='cued_success';
timeWindow=[-0.22 0.5]; %[-0.25 0.3]; % relative to alignment companion onset, in seconds
responseBaseline=[]; %[-1.05 -0.25];
cueWindow=[0 0.5];
beforeCueBaseline=[-1.05 0.2];

% get cue response of each tagged type
dd_more=cell(1,length(dd));
for i=1:length(dd)
    dd_more{i}=[dd{i} '\cue'];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
[~,ma]=nanmax(nanmean(aligncomp_y,1));
temp=nanmean(aligncomp_x,1);
putAlignPeakAt=temp(ma);
D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D1tagged_cueResponse.aligncomp_x=aligncomp_x;
D1tagged_cueResponse.aligncomp_y=aligncomp_y;
D1tagged_cueResponse.excluded=excluded;
D1tagged_cueResponse.ns=ns;
D1tagged_cueResponse=putAlignmentAtTime(putAlignPeakAt,D1tagged_cueResponse);
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
D2tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D2tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D2tagged_cueResponse.aligncomp_x=aligncomp_x;
D2tagged_cueResponse.aligncomp_y=aligncomp_y;
D2tagged_cueResponse.excluded=excluded;
D2tagged_cueResponse.ns=ns;
D2tagged_cueResponse=putAlignmentAtTime(putAlignPeakAt,D2tagged_cueResponse);

% get response of each tagged type
for i=1:length(dd)
    dd_more{i}=[dd{i} '\' responseType];
end
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns]=alignToCompanion(dd_more,true,[],[],[],'D1tagged');
activityD1tagged.unitbyunit_x=unitbyunit_x;
activityD1tagged.unitbyunit_y=unitbyunit_y;
activityD1tagged.aligncomp_x=aligncomp_x;
activityD1tagged.aligncomp_y=aligncomp_y;
activityD1tagged.excluded=excluded;
activityD1tagged.ns=ns;
activityD1tagged=putAlignmentAtTime(putAlignPeakAt,activityD1tagged);
[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns]=alignToCompanion(dd_more,true,[],[],[],'A2atagged');
activityD2tagged.unitbyunit_x=unitbyunit_x;
activityD2tagged.unitbyunit_y=unitbyunit_y;
activityD2tagged.aligncomp_x=aligncomp_x;
activityD2tagged.aligncomp_y=aligncomp_y;
activityD2tagged.excluded=excluded;
activityD2tagged.ns=ns;
activityD2tagged=putAlignmentAtTime(putAlignPeakAt,activityD2tagged);

[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);

figure(); 
histogram([activityD1tagged.ns; activityD2tagged.ns],50);
title('Histogram of trial counts across units');

% exclude units with too few trials
trial_n_cutoff=5; % at least this many trials, else exclude unit
activityD1tagged.excluded=activityD1tagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD1tagged.excluded==1)) ' D1 units because too few trials']);
activityD2tagged.excluded=activityD2tagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD2tagged.excluded==1)) ' D2 units because too few trials']);
[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);

scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindow, responseBaseline, cueWindow, beforeCueBaseline, pvalcutoff);

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

function [D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged)

excluded1=D1tagged_cueResponse.excluded; 
excluded2=activityD1tagged.excluded;
eitherExcluded=excluded1==1 | excluded2==1;
D1tagged_cueResponse.unitbyunit_x=D1tagged_cueResponse.unitbyunit_x(eitherExcluded==0,:);
D1tagged_cueResponse.unitbyunit_y=D1tagged_cueResponse.unitbyunit_y(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_x=D1tagged_cueResponse.aligncomp_x(eitherExcluded==0,:);
D1tagged_cueResponse.aligncomp_y=D1tagged_cueResponse.aligncomp_y(eitherExcluded==0,:);
D1tagged_cueResponse.excluded=D1tagged_cueResponse.excluded(eitherExcluded==0);
D1tagged_cueResponse.ns=D1tagged_cueResponse.ns(eitherExcluded==0);

activityD1tagged.unitbyunit_x=activityD1tagged.unitbyunit_x(eitherExcluded==0,:);
activityD1tagged.unitbyunit_y=activityD1tagged.unitbyunit_y(eitherExcluded==0,:);
activityD1tagged.aligncomp_x=activityD1tagged.aligncomp_x(eitherExcluded==0,:);
activityD1tagged.aligncomp_y=activityD1tagged.aligncomp_y(eitherExcluded==0,:);
activityD1tagged.excluded=activityD1tagged.excluded(eitherExcluded==0);
activityD1tagged.ns=activityD1tagged.ns(eitherExcluded==0);

end
