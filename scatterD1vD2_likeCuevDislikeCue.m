function scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindowToPlot, cueWindow, beforeCueBaseline)

% time windows are in seconds relative to onset of alignment companion

% find cue-responsive cells
[pvals_D1,cueResponseIncrease_D1]=findCellsWithResponse(D1tagged_cueResponse, cueWindow, beforeCueBaseline);
[pvals_D2,cueResponseIncrease_D2]=findCellsWithResponse(D2tagged_cueResponse, cueWindow, beforeCueBaseline);

% plot responses D1 v D2 and cued v uncued
responses_D1=getResponses(activityD1tagged, timeWindowToPlot);
responses_D2=getResponses(activityD2tagged, timeWindowToPlot);

figure();
temp=responses_D1(pvals_D1<0.05 & cueResponseIncrease_D1==true);
c=[0.8500 0.3250 0.0980];
s=scatter(zeros(size(temp)),temp,[],'filled','Color',c,'MarkerFaceAlpha',0.5);
hold on; 
temp=responses_D1(pvals_D1<0.05 & cueResponseIncrease_D1==false);
c=[0.4940 0.1840 0.5560];
s=scatter(temp,zeros(size(temp)),[],'filled','Color',c,'MarkerFaceAlpha',0.5);
% temp=responses_D1(pvals_D1>=0.1);
% c='k';
% s=scatter(zeros(size(temp)),temp,[],'filled','Color',c);
% s.AlphaData=0.2;

temp=responses_D2(pvals_D2<0.05 & cueResponseIncrease_D2==true);
c=[0.3010 0.7450 0.9330];
s=scatter(zeros(size(temp)),-temp,[],'filled','Color',c,'MarkerFaceAlpha',0.5);
hold on; 
temp=responses_D2(pvals_D2<0.05 & cueResponseIncrease_D2==false);
c=[0.4660 0.6740 0.1880];
s=scatter(-temp,zeros(size(temp)),[],'filled','Color',c,'MarkerFaceAlpha',0.5);
% temp=responses_D2(pvals_D2>=0.1);
% c='k';
% s=scatter(zeros(size(temp)),temp,[],'filled','Color',c);
% s.AlphaData=0.2;

end

function responses=getResponses(activityD1tagged, timeWindowToPlot)

% Get which cells respond to the cue or not
% Take baseline fluctuations and compare to fluctuations of activity level
% during cue window to get p-val for this cell
% Example format:
% activityD1tagged.unitbyunit_x=unitbyunit_x;
% activityD1tagged.unitbyunit_y=unitbyunit_y;
% activityD1tagged.aligncomp_x=aligncomp_x;
% activityD1tagged.aligncomp_y=aligncomp_y;
responses=nan(1,size(activityD1tagged.unitbyunit_y,1));
for i=1:size(activityD1tagged.unitbyunit_y,1)
    % Find alignment companion onset
    x=nanmean(activityD1tagged.aligncomp_x,1);
    y=nanmean(activityD1tagged.aligncomp_y,1);
    [~,ma]=nanmax(y);
    timeOfAlignCompOnset=x(ma);
    timeOfWindow=timeOfAlignCompOnset+timeWindowToPlot;
    % Get indices for baseline and response window
    x=activityD1tagged.unitbyunit_x(i,:);
    y=activityD1tagged.unitbyunit_y(i,:);
    [~,indsForWindow_start]=nanmin(abs(timeOfWindow(1)-x));
    [~,indsForWindow_end]=nanmin(abs(timeOfWindow(2)-x));
    responses(i)=nanmean(y(indsForWindow_start:indsForWindow_end));
end

end

function [pvals,cueResponseIncrease]=findCellsWithResponse(D1tagged_cueResponse, cueWindow, beforeCueBaseline)

% Get which cells respond to the cue or not
% Take baseline fluctuations and compare to fluctuations of activity level
% during cue window to get p-val for this cell
% Example format:
% D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
% D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
% D1tagged_cueResponse.aligncomp_x=aligncomp_x;
% D1tagged_cueResponse.aligncomp_y=aligncomp_y;
cueResponseIncrease=nan(1,size(D1tagged_cueResponse.unitbyunit_y,1));
pvals=nan(1,size(D1tagged_cueResponse.unitbyunit_y,1));
for i=1:size(D1tagged_cueResponse.unitbyunit_y,1)
    % Find alignment companion onset
    x=nanmean(D1tagged_cueResponse.aligncomp_x,1);
    y=nanmean(D1tagged_cueResponse.aligncomp_y,1);
    timeOfAlignCompOnset=x(find(y>0.05,1,'first'));
    timeOfBaseline=timeOfAlignCompOnset+beforeCueBaseline; % before cue baseline should be negative
    timeOfWindow=timeOfAlignCompOnset+cueWindow;
    % Get indices for baseline and response window
    x=D1tagged_cueResponse.unitbyunit_x(i,:);
    y=D1tagged_cueResponse.unitbyunit_y(i,:);
    [~,indsForBaseline_start]=nanmin(abs(timeOfBaseline(1)-x));
    [~,indsForBaseline_end]=nanmin(abs(timeOfBaseline(2)-x));
    [~,indsForWindow_start]=nanmin(abs(timeOfWindow(1)-x));
    [~,indsForWindow_end]=nanmin(abs(timeOfWindow(2)-x));
    % compare baseline to response during window
    pvals(i)=ranksum(y(indsForBaseline_start:indsForBaseline_end),y(indsForWindow_start:indsForWindow_end));
    cueResponseIncrease(i)=nanmean(y(indsForWindow_start:indsForWindow_end))>nanmean(y(indsForBaseline_start:indsForBaseline_end));
end

end