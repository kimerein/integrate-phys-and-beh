function [D1tagged_cueResponse,D2tagged_cueResponse,activityD1tagged,activityD2tagged]=scriptToOrganizeD1vD2unitResponses(dd,settings)

% load dd

pvalcutoff=settings.pvalcutoff;
responseType=settings.responseType;
timeWindow=settings.timeWindow;
responseBaseline=settings.responseBaseline;
cueWindow=settings.cueWindow;
beforeCueBaseline=settings.beforeCueBaseline;

% pvalcutoff=[-0.001 0.2]; %[0.5 1.001];

% choose type of response and time window to analyze
% responseType='uncued_failure';
% timeWindow=[-0.22 3.32]; %[-0.22 0.5]; %[-0.25 0.3]; % relative to alignment companion onset, in seconds
% responseBaseline=[]; %[-1.05 -0.25];
% cueWindow=[0 0.5]; %[5 16]; %[0 3]; %[0 0.5];
% beforeCueBaseline=[-1.05 -0.2];

if ~settings.skipCueAlignment
    % get cue response of each tagged type
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} '\cue'];
    end
    settings.studyOptoTag=true;
    [D1tagged_cueResponse,D1orD2taggingExpt,putAlignPeakAt,firstVal_D1tagged,optoFR_D1tagged,indsForAnalysisPerSess]=getAndSaveResponse(dd_more,'D1tagged',settings,[],[]);
    [D1untagged_cueResponse,~,~,firstVal_D1untagged,optoFR_D1untagged]=getAndSaveResponse(dd_more(D1orD2taggingExpt==1),'__',settings,putAlignPeakAt,indsForAnalysisPerSess);
    [D2tagged_cueResponse,~,~,firstVal_D2tagged,optoFR_D2tagged,indsForAnalysisPerSess]=getAndSaveResponse(dd_more,'A2atagged',settings,putAlignPeakAt,[]);
    [D2untagged_cueResponse,~,~,firstVal_D2untagged,optoFR_D2untagged]=getAndSaveResponse(dd_more(D1orD2taggingExpt==2),'__',settings,putAlignPeakAt,indsForAnalysisPerSess);
    figure();
    scatter(firstVal_D1tagged,optoFR_D1tagged,[],'k');
    hold on;
    scatter(firstVal_D1untagged,optoFR_D1untagged,[],'r');
    scatter(firstVal_D2tagged,optoFR_D2tagged,[],'c');
    scatter(firstVal_D2untagged,optoFR_D2untagged,[],'g');
    legend({'D1tagged','D1untagged','D2tagged','D2untagged'});
    D1tagged_cueResponse.excluded((firstVal_D1tagged<3 & optoFR_D1tagged<5) | firstVal_D1tagged<0)=1;
    D2tagged_cueResponse.excluded((firstVal_D2tagged<3 & optoFR_D2tagged<5) | firstVal_D2tagged<0)=1;
    figure();
    scatter(firstVal_D1tagged(D1tagged_cueResponse.excluded==0),optoFR_D1tagged(D1tagged_cueResponse.excluded==0),[],'k');
    hold on;
    scatter(firstVal_D1untagged,optoFR_D1untagged,[],'r');
    scatter(firstVal_D2tagged(D2tagged_cueResponse.excluded==0),optoFR_D2tagged(D2tagged_cueResponse.excluded==0),[],'c');
    scatter(firstVal_D2untagged,optoFR_D2untagged,[],'g');
    legend({'D1tagged','D1untagged','D2tagged','D2untagged'});
    title('After cut');
end

% get response of each tagged type
for i=1:length(dd)
    dd_more{i}=[dd{i} '\' responseType];
end
settings.studyOptoTag=false;
activityD1tagged=getAndSaveResponse(dd_more,'D1tagged',settings,putAlignPeakAt,[]);
activityD1untagged=getAndSaveResponse(dd_more(D1orD2taggingExpt==1),'__',settings,putAlignPeakAt,[]);
activityD2tagged=getAndSaveResponse(dd_more,'A2atagged',settings,putAlignPeakAt,[]);
activityD2untagged=getAndSaveResponse(dd_more(D1orD2taggingExpt==2),'__',settings,putAlignPeakAt,[]);

[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D1untagged_cueResponse,activityD1untagged]=cutExcluded(D1untagged_cueResponse,activityD1untagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);
[D2untagged_cueResponse,activityD2untagged]=cutExcluded(D2untagged_cueResponse,activityD2untagged);

figure(); 
histogram([activityD1tagged.ns; activityD2tagged.ns; activityD1untagged.ns; activityD2untagged.ns],50);
title('Histogram of trial counts across units');

% exclude units with too few trials
trial_n_cutoff=1; % at least this many trials, else exclude unit
activityD1tagged.excluded=activityD1tagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD1tagged.excluded==1)) ' D1 units because too few trials']);
activityD1untagged.excluded=activityD1untagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD1untagged.excluded==1)) ' untagged during D1 units because too few trials']);
activityD2tagged.excluded=activityD2tagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD2tagged.excluded==1)) ' D2 units because too few trials']);
activityD2untagged.excluded=activityD2untagged.ns'<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD2untagged.excluded==1)) ' untagged during D2 units because too few trials']);
[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D1untagged_cueResponse,activityD1untagged]=cutExcluded(D1untagged_cueResponse,activityD1untagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);
[D2untagged_cueResponse,activityD2untagged]=cutExcluded(D2untagged_cueResponse,activityD2untagged);

figure(); plot(nanmean(activityD1tagged.unitbyunit_x,1),activityD1tagged.unitbyunit_y'); hold on; plot(nanmean(activityD1tagged.aligncomp_x,1),nanmean(activityD1tagged.aligncomp_y,1),'Color','b'); title('D1');
figure(); plot(nanmean(activityD2tagged.unitbyunit_x,1),activityD2tagged.unitbyunit_y'); hold on; plot(nanmean(activityD2tagged.aligncomp_x,1),nanmean(activityD2tagged.aligncomp_y,1),'Color','b'); title('D2');

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
figure(); plot(timesD1,activityD1tagged.unitbyunit_y'); title('Check D1 times');
temp=nanmean(activityD1untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1untagged.aligncomp_y,1));
timesD1un=nanmean(activityD1untagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2tagged.aligncomp_y,1));
timesD2=nanmean(activityD2tagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2untagged.aligncomp_y,1));
timesD2un=nanmean(activityD2untagged.unitbyunit_x,1)-temp(f);
takewin1=[2.25 3.56]-3.2476; % relative to peak of alignment companion
takewin2=[3.56 5.56]-3.2476;
% takewin1=[2.25 3]-3.2476; % relative to peak of alignment companion
% takewin2=[3 3.88]-3.2476;
figure(); 
scatter([nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin1(2)),2)./nanmax(activityD2untagged.unitbyunit_y(:,1:200),[],2)],[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>takewin2(1) & timesD2un<=takewin2(2)),2)./nanmax(activityD2untagged.unitbyunit_y(:,1:200),[],2)],[],'k');
title('D1');
figure(); 
scatter([nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin1(2)),2)./nanmax(activityD1untagged.unitbyunit_y(:,1:200),[],2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)./nanmax(activityD2tagged.unitbyunit_y(:,1:200),[],2)],[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>takewin2(1) & timesD1un<=takewin2(2)),2)./nanmax(activityD1untagged.unitbyunit_y(:,1:200),[],2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>takewin2(1) & timesD2<=takewin2(2)),2)./nanmax(activityD2tagged.unitbyunit_y(:,1:200),[],2)],[],'r');
title('D2'); 
figure(); 
scatter(nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2),nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2),[],'k');
title('D1 definite'); 
figure(); 
scatter(nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)./nanmax(activityD2tagged.unitbyunit_y(:,1:200),[],2),nanmean(activityD2tagged.unitbyunit_y(:,timesD2>takewin2(1) & timesD2<=takewin2(2)),2)./nanmax(activityD2tagged.unitbyunit_y(:,1:200),[],2),[],'r');
title('D2 definite'); 

figure();
scatter([nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin1(2)),2)],[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>takewin2(1) & timesD2un<=takewin2(2)),2)],[],'k');
hold on;
scatter([nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin1(2)),2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)],[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>takewin2(1) & timesD1un<=takewin2(2)),2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>takewin2(1) & timesD2<=takewin2(2)),2)],[],'r');
title('D1 and D2 w inferred');
figure(); 
scatter(nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2),nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2),[],'k');
hold on;
scatter(nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2),nanmean(activityD2tagged.unitbyunit_y(:,timesD2>takewin2(1) & timesD2<=takewin2(2)),2),[],'r');
title('D1 and D2 definite'); 

D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)];
D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin2(1) & timesD1<=takewin2(2)),2)];
D2temp=[nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)];
D2temp_2ndwin=[nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin2(1) & timesD2<=takewin2(2)),2)];
% D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin1(2)),2)];
% D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2); nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>takewin2(1) & timesD2un<=takewin2(2)),2)];
% D2temp=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin1(2)),2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)];
% D2temp_2ndwin=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>takewin2(1) & timesD1un<=takewin2(2)),2); nanmean(activityD2tagged.unitbyunit_y(:,timesD2>takewin2(1) & timesD2<=takewin2(2)),2)];
D1temp(D1temp<0.0001)=0;
D1temp_2ndwin(D1temp_2ndwin<0.0001)=0;
D2temp(D2temp<0.0001)=0;
D2temp_2ndwin(D2temp_2ndwin<0.0001)=0;
[n,x]=hist(D1temp_2ndwin./D1temp,0:0.1:10);
[n_D1,x_D1]=cityscape_hist(n,x);
[n,x]=hist(D2temp_2ndwin./D2temp,0:0.1:10);
[n_D2,x_D2]=cityscape_hist(n,x);
figure(); plot(x_D1,n_D1./nansum(n_D1),'Color','k'); hold on; plot(x_D2,n_D2./nansum(n_D2),'Color','r'); legend({'D1','D2'});

% Zscore each time point across the population, align at time window 1 mean, then plot transition to time window 2
temp=activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin2(2));
wholePop=[activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin2(2)); activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin2(2)); activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin2(2)); activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin2(2))];
subD1Times=timesD1(timesD1>=takewin1(1) & timesD1<=takewin2(2));
subD2Times=timesD2(timesD2>=takewin1(1) & timesD2<=takewin2(2));
Zscored_D1tagged=(temp-repmat(mean(wholePop,1,'omitnan'),size(temp,1),1))./repmat(std(wholePop,0,1,'omitnan'),size(temp,1),1);
Zscored_D1tagged=Zscored_D1tagged-repmat(mean(Zscored_D1tagged(:,subD1Times>=takewin1(1) & subD1Times<=takewin1(2)),2,'omitnan'),1,size(temp,2));
temp=activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin2(2));
Zscored_D2tagged=(temp-repmat(mean(wholePop,1,'omitnan'),size(temp,1),1))./repmat(std(wholePop,0,1,'omitnan'),size(temp,1),1);
Zscored_D2tagged=Zscored_D2tagged-repmat(mean(Zscored_D2tagged(:,subD2Times>=takewin1(1) & subD2Times<=takewin1(2)),2,'omitnan'),1,size(temp,2));

useUnits=nanmean([D1temp D1temp_2ndwin],2)>0.2;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
figure(); 
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
    currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,subD1Times>=-1.05 & subD1Times<-0.9),2),1,size(temp,2)))';
%     currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
    if nanmean(temp(useUnits,subD1Times<0))>0.3
        c='k';
    else
        c='r';
    end
    plot(timesD1(timesD1>=takewin1(1) & timesD1<=takewin2(2)),currtoplot,'Color',c); hold on;
end
title('D1');
useUnits=nanmean([D2temp D2temp_2ndwin],2)>0.2;
temp=Zscored_D2tagged;
for i=1:size(temp,1)
temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
figure(); plot(timesD2(timesD2>=takewin1(1) & timesD2<=takewin2(2)),(temp(useUnits,:)-repmat(nanmean(temp(useUnits,subD2Times>=-1.05 & subD2Times<-0.9),2),1,size(temp,2)))'); title('D2');

scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindow, responseBaseline, cueWindow, beforeCueBaseline, pvalcutoff);

end


function [D1tagged_cueResponse,D1orD2taggingExpt,putAlignPeakAt,firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess]=getAndSaveResponse(dd_more,getThese,settings,putAlignPeakAtInput,indsForAnalysisPerSess)

[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns,D1orD2taggingExpt,firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess]=alignToCompanion(dd_more,false,[],[],[],getThese,settings,indsForAnalysisPerSess);
[~,ma]=nanmax(nanmean(aligncomp_y,1));
temp=nanmean(aligncomp_x,1);
putAlignPeakAt=temp(ma);
D1tagged_cueResponse.unitbyunit_x=unitbyunit_x;
D1tagged_cueResponse.unitbyunit_y=unitbyunit_y;
D1tagged_cueResponse.aligncomp_x=aligncomp_x;
D1tagged_cueResponse.aligncomp_y=aligncomp_y;
D1tagged_cueResponse.excluded=excluded;
D1tagged_cueResponse.ns=ns;
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
