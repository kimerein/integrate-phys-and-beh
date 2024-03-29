function [D1tagged_cueResponse,D2tagged_cueResponse,activityD1tagged,activityD2tagged]=scriptToOrganizeD1vD2unitResponses(dd,settings)

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
        if ismac()
            dd_more{i}=[dd{i} '/cue'];
        else
            dd_more{i}=[dd{i} '\cue'];
        end
    end
    
    [D1tagged_cueResponse,D1orD2taggingExpt,putAlignPeakAt]=getAndSaveResponse(dd_more,'D1tagged',settings,[],[]);
    [D1untagged_cueResponse]=getAndSaveResponse(dd_more,'__',settings,putAlignPeakAt,[]);
    D1untagged_cueResponse=takeOnlyUnitsFromSess(D1untagged_cueResponse,unique(D1tagged_cueResponse.fromWhichSess));
    [D2tagged_cueResponse]=getAndSaveResponse(dd_more,'A2atagged',settings,putAlignPeakAt,[]);
    [D2untagged_cueResponse]=getAndSaveResponse(dd_more,'__',settings,putAlignPeakAt,[]);
    D2untagged_cueResponse=takeOnlyUnitsFromSess(D2untagged_cueResponse,unique(D2tagged_cueResponse.fromWhichSess));

end

% get response of each tagged type
for i=1:length(dd)
    if ismac()
        dd_more{i}=[dd{i} '/' responseType];
    else
        dd_more{i}=[dd{i} '\' responseType];
    end
end
settings.studyOptoTag=false;
activityD1tagged=getAndSaveResponse(dd_more,'D1tagged',settings,putAlignPeakAt,[]);
activityD1untagged=getAndSaveResponse(dd_more,'__',settings,putAlignPeakAt,[]);
activityD1untagged=takeOnlyUnitsFromSess(activityD1untagged,unique(activityD1tagged.fromWhichSess));
activityD2tagged=getAndSaveResponse(dd_more,'A2atagged',settings,putAlignPeakAt,[]);
activityD2untagged=getAndSaveResponse(dd_more,'__',settings,putAlignPeakAt,[]);
activityD2untagged=takeOnlyUnitsFromSess(activityD2untagged,unique(activityD2tagged.fromWhichSess));

% [D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
% [D1untagged_cueResponse,activityD1untagged]=cutExcluded(D1untagged_cueResponse,activityD1untagged);
% [D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);
% [D2untagged_cueResponse,activityD2untagged]=cutExcluded(D2untagged_cueResponse,activityD2untagged);

figure(); 
histogram([activityD1tagged.ns; activityD2tagged.ns; activityD1untagged.ns; activityD2untagged.ns],50);
title('Histogram of trial counts across units');

% exclude units with too few trials
D1tagged_cueResponse.excluded=zeros(size(D1tagged_cueResponse.ns));
D1untagged_cueResponse.excluded=zeros(size(D1untagged_cueResponse.ns));
D2tagged_cueResponse.excluded=zeros(size(D2tagged_cueResponse.ns));
D2untagged_cueResponse.excluded=zeros(size(D2untagged_cueResponse.ns));
trial_n_cutoff=1; % at least this many trials, else exclude unit
activityD1tagged.excluded=activityD1tagged.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD1tagged.excluded==1)) ' D1 units because too few trials']);
activityD1untagged.excluded=activityD1untagged.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD1untagged.excluded==1)) ' untagged during D1 units because too few trials']);
activityD2tagged.excluded=activityD2tagged.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD2tagged.excluded==1)) ' D2 units because too few trials']);
activityD2untagged.excluded=activityD2untagged.ns<trial_n_cutoff;
disp(['excluding ' num2str(nansum(activityD2untagged.excluded==1)) ' untagged during D2 units because too few trials']);
[D1tagged_cueResponse,activityD1tagged]=cutExcluded(D1tagged_cueResponse,activityD1tagged);
[D1untagged_cueResponse,activityD1untagged]=cutExcluded(D1untagged_cueResponse,activityD1untagged);
[D2tagged_cueResponse,activityD2tagged]=cutExcluded(D2tagged_cueResponse,activityD2tagged);
[D2untagged_cueResponse,activityD2untagged]=cutExcluded(D2untagged_cueResponse,activityD2untagged);

makeSomePlots(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged);
pause; close all;
morePlots(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged);
pause; close all;
sessionBySessionPlot(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged);

scatterD1vD2_likeCuevDislikeCue(D1tagged_cueResponse, D2tagged_cueResponse, activityD1tagged, activityD2tagged, timeWindow, responseBaseline, cueWindow, beforeCueBaseline, pvalcutoff);

% Now UNTAGGED
pause;
scatterD1vD2_likeCuevDislikeCue(D2untagged_cueResponse, D1untagged_cueResponse, activityD2untagged, activityD1untagged, timeWindow, responseBaseline, cueWindow, beforeCueBaseline, pvalcutoff);

end

function D1tagged_cueResponse=takeOnlyUnitsFromSess(D1tagged_cueResponse,takeTheseSess)

D1tagged_cueResponse.unitbyunit_x=D1tagged_cueResponse.unitbyunit_x(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess),:);
D1tagged_cueResponse.unitbyunit_y=D1tagged_cueResponse.unitbyunit_y(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess),:);
D1tagged_cueResponse.aligncomp_x=D1tagged_cueResponse.aligncomp_x(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess),:);
D1tagged_cueResponse.aligncomp_y=D1tagged_cueResponse.aligncomp_y(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess),:);
D1tagged_cueResponse.excluded=D1tagged_cueResponse.excluded(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess));
D1tagged_cueResponse.ns=D1tagged_cueResponse.ns(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess));
D1tagged_cueResponse.fromWhichSess=D1tagged_cueResponse.fromWhichSess(ismember(D1tagged_cueResponse.fromWhichSess,takeTheseSess));

end

function makeSomePlots(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged)

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
% For success
takewin1=[2.25 3.56]-3.2476; % relative to peak of alignment companion
takewin2=[3.56 7.56]-3.2476;
% For failure
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
% temp=activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin2(2));
% wholePop=[activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin2(2)); activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin2(2)); activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin2(2)); activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin2(2))];
% subD1Times=timesD1(timesD1>=takewin1(1) & timesD1<=takewin2(2));
% subD2Times=timesD2(timesD2>=takewin1(1) & timesD2<=takewin2(2));
% Zscored_D1tagged=(temp-repmat(mean(wholePop,1,'omitnan'),size(temp,1),1))./repmat(std(wholePop,0,1,'omitnan'),size(temp,1),1);
% Zscored_D1tagged=Zscored_D1tagged-repmat(mean(Zscored_D1tagged(:,subD1Times>=takewin1(1) & subD1Times<=takewin1(2)),2,'omitnan'),1,size(temp,2));
% temp=activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin2(2));
% Zscored_D2tagged=(temp-repmat(mean(wholePop,1,'omitnan'),size(temp,1),1))./repmat(std(wholePop,0,1,'omitnan'),size(temp,1),1);
% Zscored_D2tagged=Zscored_D2tagged-repmat(mean(Zscored_D2tagged(:,subD2Times>=takewin1(1) & subD2Times<=takewin1(2)),2,'omitnan'),1,size(temp,2));

% Zscore within each unit, then plot
% transition to time window 2
ZeroAt=takewin1;
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D1tagged=temp-repmat(mean(temp(:,timesD1>=ZeroAt(1) & timesD1<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
% Zscored_D1tagged=temp;
temp=activityD2tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D2tagged=temp-repmat(mean(temp(:,timesD2>=ZeroAt(1) & timesD2<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
% Zscored_D2tagged=temp;

useUnits=nanmean([D1temp D1temp_2ndwin],2)>0.1;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD1<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD1>=0.2 & timesD1<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
    if nanmean(temp(useUnits,timesD1<-1))>popmedian
        c='k';
    else
        c='r';
    end
    plot(timesD1,currtoplot,'Color',c); hold on;
end
title('D1');

useUnits=nanmean([D2temp D2temp_2ndwin],2)>0.1;
temp=Zscored_D2tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD2<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD2>=0.2 & timesD2<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
    if nanmean(temp(useUnits,timesD2<-1))>popmedian
        c='k';
    else
        c='r';
    end
    plot(timesD2,currtoplot,'Color',c); hold on;
end
title('D2');

% Population vector after the outcome
% Which pre-outcome time points correlate best
timeBinsStep=0.25;
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D1tagged=temp;
temp=activityD2tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D2tagged=temp;

temp=activityD1untagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D1untagged=temp;
temp=activityD2untagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D2untagged=temp;

popD1=mean(Zscored_D1tagged(:,timesD1>0.3 & timesD1<4),2,'omitnan');
popD1(isnan(popD1))=0;
timeBins=timesD1(1):timeBinsStep:timesD1(end);
Rs=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD1=mean(Zscored_D1tagged(:,timesD1>timeBins(i) & timesD1<timeBins(i+1)),2,'omitnan');
    currpopD1(isnan(currpopD1))=0;
    R=corrcoef(popD1,currpopD1);
    Rs(i)=R(1,2);
end
figure(); 
plot(timeBins(1:end-1),Rs,'Color','k'); title('D1 v D2 corrcoef with post-outcome driven pop vec');

popD2=mean(Zscored_D2tagged(:,timesD2>0.3 & timesD2<4),2,'omitnan');
popD2(isnan(popD2))=0;
timeBins=timesD2(1):timeBinsStep:timesD2(end);
Rs=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD2=mean(Zscored_D2tagged(:,timesD2>timeBins(i) & timesD2<timeBins(i+1)),2,'omitnan');
    currpopD2(isnan(currpopD2))=0;
    R=corrcoef(popD2,currpopD2);
    Rs(i)=R(1,2);
end
hold on; 
plot(timeBins(1:end-1),Rs,'Color','r'); 

timeBins=timesD1(1):timeBinsStep:timesD1(end);
Rs=nan(length(timeBins)-1,length(timeBins)-1);
for i=1:length(timeBins)-1
    popD1=mean(Zscored_D1tagged(:,timesD1>timeBins(i) & timesD1<timeBins(i+1)),2,'omitnan');
    popD1(isnan(popD1))=0;
    for j=1:length(timeBins)-1
        currpopD1=mean(Zscored_D1tagged(:,timesD1>timeBins(j) & timesD1<timeBins(j+1)),2,'omitnan');
        currpopD1(isnan(currpopD1))=0;
        R=corrcoef(popD1,currpopD1);
        Rs(i,j)=R(1,2);
    end
end
Rs_D1=Rs;
figure(); 
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),Rs(timeBins<9.5,timeBins<9.5)); title('D1 corrcoef matrix');

timeBins=timesD2(1):timeBinsStep:timesD2(end);
RsD2=nan(length(timeBins)-1,length(timeBins)-1);
for i=1:length(timeBins)-1
    popD2=mean(Zscored_D2tagged(:,timesD2>timeBins(i) & timesD2<timeBins(i+1)),2,'omitnan');
    popD2(isnan(popD2))=0;
    for j=1:length(timeBins)-1
        currpopD2=mean(Zscored_D2tagged(:,timesD2>timeBins(j) & timesD2<timeBins(j+1)),2,'omitnan');
        currpopD2(isnan(currpopD2))=0;
        R=corrcoef(popD2,currpopD2);
        RsD2(i,j)=R(1,2);
    end
end
Rs_D2=RsD2;
figure(); 
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),RsD2(timeBins<9.5,timeBins<9.5)); title('D2 corrcoef matrix');

figure();
l=min([size(Rs,1) size(RsD2,1) find(timeBins<9.5,1,'last')]);
imagesc([Rs(1:l, 1:l) RsD2(1:l, 1:l)]);
title('D1 then D2 corrcoef matrix');

timeBins=timesD1un(1):timeBinsStep:timesD1un(end);
Rs=nan(length(timeBins)-1,length(timeBins)-1);
for i=1:length(timeBins)-1
    popD1=mean(Zscored_D1untagged(:,timesD1un>timeBins(i) & timesD1un<timeBins(i+1)),2,'omitnan');
    popD1(isnan(popD1))=0;
    for j=1:length(timeBins)-1
        currpopD1=mean(Zscored_D1untagged(:,timesD1un>timeBins(j) & timesD1un<timeBins(j+1)),2,'omitnan');
        currpopD1(isnan(currpopD1))=0;
        R=corrcoef(popD1,currpopD1);
        Rs(i,j)=R(1,2);
    end
end
RsD1un=Rs;
figure(); 
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),Rs(timeBins<9.5,timeBins<9.5)); title('D1 untagged corrcoef matrix');

timeBins=timesD2un(1):timeBinsStep:timesD2un(end);
RsD2=nan(length(timeBins)-1,length(timeBins)-1);
for i=1:length(timeBins)-1
    popD2=mean(Zscored_D2untagged(:,timesD2un>timeBins(i) & timesD2un<timeBins(i+1)),2,'omitnan');
    popD2(isnan(popD2))=0;
    for j=1:length(timeBins)-1
        currpopD2=mean(Zscored_D2untagged(:,timesD2un>timeBins(j) & timesD2un<timeBins(j+1)),2,'omitnan');
        currpopD2(isnan(currpopD2))=0;
        R=corrcoef(popD2,currpopD2);
        RsD2(i,j)=R(1,2);
    end
end
RsD2un=RsD2;
figure(); 
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),RsD2(timeBins<9.5,timeBins<9.5)); title('D2 untagged corrcoef matrix');

figure();
l=min([size(Rs_D1,1) size(Rs_D2,1) find(timeBins<9.5,1,'last')]);
imagesc([Rs_D1(1:l, 1:l) Rs_D2(1:l, 1:l) RsD1un(1:l, 1:l) RsD2un(1:l, 1:l)]);
title('D1 then D2 then D1un then D2un autocorr matrices');

end

function sessionBySessionPlot(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged)

projTimeWindow1=[-1 0]; %[-2 -0.4]; %[2.25 3.56]-3.2476;
projTimeWindow2=[1.5 5]; %[3.56 7.26]-3.2476;
projTimeWindow3=[-2.21859 -1.46859];
projTimeWindow4=[0.781415 1.78141];
winBeforeCue=[6 7]; %[-3 -2.5];

uSess=unique(activityD1tagged.fromWhichSess);
D1sess=true;
figure();
offset=0;
D1PreOutcome=nan(1,length(uSess));
D1PreOutcome_winBeforeCue=nan(1,length(uSess));
D1unPreOutcome_winBeforeCue=nan(1,length(uSess));
D1Random=nan(1,length(uSess));
D1unPreOutcome=nan(1,length(uSess));
D1unRandom=nan(1,length(uSess));
D1_tw1=[]; D1_tw2=[]; D1_tw3=[]; D1_tw4=[];
D1un_tw1=[]; D1un_tw2=[]; D1un_tw3=[]; D1un_tw4=[];
for i=1:length(uSess)
    [D1ontoCD, D1unontoCD]=sessionBySessionPlot_subFunc(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged,uSess(i),D1sess,offset);
    offset=offset+5;
    % for D1, use projection 2 onto projection 1
    % for D1 un, use projection 4 onto projection 3
    D1PreOutcome(i)=nanmean(D1ontoCD.projTimeWindow2.R(D1ontoCD.projTimeWindow2.t>projTimeWindow1(1) & D1ontoCD.projTimeWindow2.t<projTimeWindow1(2)))-nanmean(D1ontoCD.projTimeWindow2.R(D1ontoCD.projTimeWindow2.t>projTimeWindow2(1) & D1ontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D1PreOutcome_winBeforeCue(i)=nanmean(D1ontoCD.projTimeWindow2.R(D1ontoCD.projTimeWindow2.t>winBeforeCue(1) & D1ontoCD.projTimeWindow2.t<winBeforeCue(2)))-nanmean(D1ontoCD.projTimeWindow2.R(D1ontoCD.projTimeWindow2.t>projTimeWindow2(1) & D1ontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D1Random(i)=D1ontoCD.projTimeWindow2.mfr;
    D1_tw1=[D1_tw1; D1ontoCD.projTimeWindow1.R]; D1_tw2=[D1_tw2; D1ontoCD.projTimeWindow2.R]; D1_tw3=[D1_tw3; D1ontoCD.projTimeWindow3.R]; D1_tw4=[D1_tw4; D1ontoCD.projTimeWindow4.R];
    D1unPreOutcome(i)=nanmean(D1unontoCD.projTimeWindow2.R(D1unontoCD.projTimeWindow2.t>projTimeWindow1(1) & D1unontoCD.projTimeWindow2.t<projTimeWindow1(2)))-nanmean(D1unontoCD.projTimeWindow2.R(D1unontoCD.projTimeWindow2.t>projTimeWindow2(1) & D1unontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D1unPreOutcome_winBeforeCue(i)=nanmean(D1unontoCD.projTimeWindow2.R(D1unontoCD.projTimeWindow2.t>winBeforeCue(1) & D1unontoCD.projTimeWindow2.t<winBeforeCue(2)))-nanmean(D1unontoCD.projTimeWindow2.R(D1unontoCD.projTimeWindow2.t>projTimeWindow2(1) & D1unontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1; 
    D1unRandom(i)=D1unontoCD.projTimeWindow2.mfr;
    D1un_tw1=[D1un_tw1; D1unontoCD.projTimeWindow1.R]; D1un_tw2=[D1un_tw2; D1unontoCD.projTimeWindow2.R]; D1un_tw3=[D1un_tw3; D1unontoCD.projTimeWindow3.R]; D1un_tw4=[D1un_tw4; D1unontoCD.projTimeWindow4.R];
end
title('D1 experiments');
uSess=unique(activityD2tagged.fromWhichSess);
D1sess=false;
offset=0;
figure();
D2PreOutcome=nan(1,length(uSess));
D2PreOutcome_winBeforeCue=nan(1,length(uSess));
D2unPreOutcome_winBeforeCue=nan(1,length(uSess));
D2Random=nan(1,length(uSess));
D2unPreOutcome=nan(1,length(uSess));
D2unRandom=nan(1,length(uSess));
D2_tw1=[]; D2_tw2=[]; D2_tw3=[]; D2_tw4=[];
D2un_tw1=[]; D2un_tw2=[]; D2un_tw3=[]; D2un_tw4=[];
for i=1:length(uSess)
    [D2ontoCD, D2unontoCD]=sessionBySessionPlot_subFunc(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged,uSess(i),D1sess,offset);
    offset=offset+5;
    % for D2 un, use projection 2 onto projection 1
    % for D2, use projection 4 onto projection 3
    D2PreOutcome(i)=nanmean(D2ontoCD.projTimeWindow2.R(D2ontoCD.projTimeWindow2.t>projTimeWindow1(1) & D2ontoCD.projTimeWindow2.t<projTimeWindow1(2)))-nanmean(D2ontoCD.projTimeWindow2.R(D2ontoCD.projTimeWindow2.t>projTimeWindow2(1) & D2ontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D2PreOutcome_winBeforeCue(i)=nanmean(D2ontoCD.projTimeWindow2.R(D2ontoCD.projTimeWindow2.t>winBeforeCue(1) & D2ontoCD.projTimeWindow2.t<winBeforeCue(2)))-nanmean(D2ontoCD.projTimeWindow2.R(D2ontoCD.projTimeWindow2.t>projTimeWindow2(1) & D2ontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D2Random(i)=D2ontoCD.projTimeWindow2.mfr;
    D2_tw1=[D2_tw1; D2ontoCD.projTimeWindow1.R]; D2_tw2=[D2_tw2; D2ontoCD.projTimeWindow2.R]; D2_tw3=[D2_tw3; D2ontoCD.projTimeWindow3.R]; D2_tw4=[D2_tw4; D2ontoCD.projTimeWindow4.R];
    D2unPreOutcome(i)=nanmean(D2unontoCD.projTimeWindow2.R(D2unontoCD.projTimeWindow2.t>projTimeWindow1(1) & D2unontoCD.projTimeWindow2.t<projTimeWindow1(2)))-nanmean(D2unontoCD.projTimeWindow2.R(D2unontoCD.projTimeWindow2.t>projTimeWindow2(1) & D2unontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D2unPreOutcome_winBeforeCue(i)=nanmean(D2unontoCD.projTimeWindow2.R(D2unontoCD.projTimeWindow2.t>winBeforeCue(1) & D2unontoCD.projTimeWindow2.t<winBeforeCue(2)))-nanmean(D2unontoCD.projTimeWindow2.R(D2unontoCD.projTimeWindow2.t>projTimeWindow2(1) & D2unontoCD.projTimeWindow2.t<projTimeWindow2(2)))+1;
    D2unRandom(i)=D2unontoCD.projTimeWindow2.mfr;
    D2un_tw1=[D2un_tw1; D2unontoCD.projTimeWindow1.R]; D2un_tw2=[D2un_tw2; D2unontoCD.projTimeWindow2.R]; D2un_tw3=[D2un_tw3; D2unontoCD.projTimeWindow3.R]; D2un_tw4=[D2un_tw4; D2unontoCD.projTimeWindow4.R];
end
title('D2 experiments');

% figure(); 
% % scatter(-D1PreOutcome_winBeforeCue,-D1PreOutcome,[],'k');
% scatter(D1Random,-D1PreOutcome,[],'k');
% xlabel('Uncue'); ylabel('Cue');
% % hold on; 
% % scatter(D2unRandom,-D2unPreOutcome,'Color','g');
% % legend({'Tagged D1','Untagged during D2'});
% title('D1 space');
% figure(); 
% % scatter(-D2PreOutcome_winBeforeCue,-D2PreOutcome,[],'r');
% scatter(D2Random,-D2PreOutcome,[],'r');
% xlabel('Uncue'); ylabel('Cue');
% % hold on; 
% % scatter(D1unRandom,-D1unPreOutcome,'Color','g');
% % legend({'Tagged D2','Untagged during D1'});
% title('D2 space');

figure(); 
scatter(-D1PreOutcome_winBeforeCue,-D1PreOutcome,[],'k');
% scatter(D1Random,-D1PreOutcome,[],'k');
hold on;
% scatter(D1unRandom,-D1unPreOutcome,[],'g');
scatter(-D1unPreOutcome_winBeforeCue,-D1unPreOutcome,[],'g');
scatter(-D2PreOutcome_winBeforeCue,-D2PreOutcome,[],'r');
% scatter(D2Random,-D2PreOutcome,[],'r');
% scatter(D2unRandom,-D2unPreOutcome,[],'c');
scatter(-D2unPreOutcome_winBeforeCue,-D2unPreOutcome,[],'c');
xlabel('Uncue'); ylabel('Cue');
line([-1 1],[0 0]); line([0 0],[-1 1]); %daspect([1 1 1]);
% line([0 nanmean([-D1PreOutcome_winBeforeCue D2PreOutcome_winBeforeCue])],[0 nanmean([-D1PreOutcome D2PreOutcome])],'LineWidth',6,'Color','m');
title('D1 and D2 space');

figure();
% scatter(D1Random-D1unRandom,-D1PreOutcome-D1unPreOutcome,[],'black');
scatter((-D1PreOutcome_winBeforeCue)-(-D1unPreOutcome_winBeforeCue),-D1PreOutcome+D1unPreOutcome,[],'black');
hold on;
% scatter(D2unRandom-D2Random,-D2unPreOutcome-D2PreOutcome,[],'red');
scatter((-D2unPreOutcome_winBeforeCue)-(-D2PreOutcome_winBeforeCue),-D2unPreOutcome+D2PreOutcome,[],'red');
xlabel('Uncue'); ylabel('Cue');
line([-1 1],[0 0]); line([0 0],[-1 1]);
xmean=nanmean([(-D1PreOutcome_winBeforeCue)-(-D1unPreOutcome_winBeforeCue) (-D2unPreOutcome_winBeforeCue)-(-D2PreOutcome_winBeforeCue)]);
ymean=nanmean([-D1PreOutcome+D1unPreOutcome -D2unPreOutcome+D2PreOutcome]);
line([0 xmean],[0 ymean],'LineWidth',6,'Color','m');
title('D1 minus D2 space');

figure();
temp=[[-D1PreOutcome; -D1unPreOutcome]'; [-D2unPreOutcome; -D2PreOutcome]'];
tempx=zeros(size(temp));
tempx(:,2)=1;
line(tempx', temp');


figure(); plot(D1ontoCD.projTimeWindow1.t,nanmean(D1_tw1,1));
hold on; plot(D1ontoCD.projTimeWindow1.t,nanmean(D1_tw1,1)-nanstd(D1_tw1,[],1)./sqrt(size(D1_tw1,1)));
hold on; plot(D1ontoCD.projTimeWindow1.t,nanmean(D1_tw1,1)+nanstd(D1_tw1,[],1)./sqrt(size(D1_tw1,1)));
title('D1 after cue');
figure(); plot(D1ontoCD.projTimeWindow4.t,nanmean(D1_tw4,1));
hold on; plot(D1ontoCD.projTimeWindow4.t,nanmean(D1_tw4,1)-nanstd(D1_tw4,[],1)./sqrt(size(D1_tw4,1)));
hold on; plot(D1ontoCD.projTimeWindow4.t,nanmean(D1_tw4,1)+nanstd(D1_tw4,[],1)./sqrt(size(D1_tw4,1)));
title('D1 before cue');
% figure(); plot(D1unontoCD.projTimeWindow4.t,nanmean(D1un_tw4,1));
% hold on; plot(D1unontoCD.projTimeWindow4.t,nanmean(D1un_tw4,1)-nanstd(D1un_tw4,[],1)./sqrt(size(D1un_tw4,1)));
% hold on; plot(D1unontoCD.projTimeWindow4.t,nanmean(D1un_tw4,1)+nanstd(D1un_tw4,[],1)./sqrt(size(D1un_tw4,1)));
figure(); plot(D2ontoCD.projTimeWindow4.t,nanmean(D2_tw4,1));
hold on; plot(D2ontoCD.projTimeWindow4.t,nanmean(D2_tw4,1)-nanstd(D2_tw4,[],1)./sqrt(size(D2_tw4,1)));
hold on; plot(D2ontoCD.projTimeWindow4.t,nanmean(D2_tw4,1)+nanstd(D2_tw4,[],1)./sqrt(size(D2_tw4,1)));
title('D2 before cue');
figure(); plot(D2ontoCD.projTimeWindow2.t,nanmean(D2_tw2,1));
hold on; plot(D2ontoCD.projTimeWindow2.t,nanmean(D2_tw2,1)-nanstd(D2_tw2,[],1)./sqrt(size(D2_tw2,1)));
hold on; plot(D2ontoCD.projTimeWindow2.t,nanmean(D2_tw2,1)+nanstd(D2_tw2,[],1)./sqrt(size(D2_tw2,1)));
title('D2 after cue');
% figure(); plot(D2unontoCD.projTimeWindow2.t,nanmean(D2un_tw2,1));
% hold on; plot(D2unontoCD.projTimeWindow2.t,nanmean(D2un_tw2,1)-nanstd(D2un_tw2,[],1)./sqrt(size(D2un_tw2,1)));
% hold on; plot(D2unontoCD.projTimeWindow2.t,nanmean(D2un_tw2,1)+nanstd(D2un_tw2,[],1)./sqrt(size(D2un_tw2,1)));

end

function [D1ontoCD, D1unontoCD]=sessionBySessionPlot_subFunc(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged,whichSess,isD1sess,offset)

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD1untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1untagged.aligncomp_y,1));
timesD1un=nanmean(activityD1untagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2tagged.aligncomp_y,1));
timesD2=nanmean(activityD2tagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2untagged.aligncomp_y,1));
timesD2un=nanmean(activityD2untagged.unitbyunit_x,1)-temp(f);
% For success
takewin1=[2.25 3.56]-3.2476; % relative to peak of alignment companion
takewin2=[3.56 7.26]-3.2476;
% For failure, these windows interesting -- from autocorrelation matrix
% takewin1=[-2.21859 -1.46859]; % relative to peak of alignment companion
% takewin2=[0.781415 1.78141];

% Get units from this sess
if isD1sess==true
    useUnits_D1=ismember(activityD1tagged.fromWhichSess,whichSess);
    useUnits_D1un=ismember(activityD1untagged.fromWhichSess,whichSess);
else
    useUnits_D1=ismember(activityD2tagged.fromWhichSess,whichSess);
    useUnits_D1un=ismember(activityD2untagged.fromWhichSess,whichSess);
    activityD1tagged=activityD2tagged;
    activityD1untagged=activityD2untagged;
    timesD1=timesD2;
    timesD1un=timesD2un;
end

D1untemp=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin1(2)),2)];
D1untemp_2ndwin=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin2(1) & timesD1un<=takewin2(2)),2)];
D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)];
D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin2(1) & timesD1<=takewin2(2)),2)];
D1untemp(D1untemp<0.0001)=0;
D1untemp_2ndwin(D1untemp_2ndwin<0.0001)=0;
D1temp(D1temp<0.0001)=0;
D1temp_2ndwin(D1temp_2ndwin<0.0001)=0;

turnsOn_D1=(D1temp_2ndwin./D1temp)>0.5;
turnsOn_D1un=(D1untemp_2ndwin./D1untemp)>0.5;
ZeroAt=takewin1;
temp=activityD1tagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD1<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD1<9.5),0,2,'omitnan');
% Zscored_D1untagged=temp-repmat(mean(temp(:,timesD1un>=ZeroAt(1) & timesD1un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D1tagged=temp;
temp=activityD1untagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD1un<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD1un<9.5),0,2,'omitnan');
% Zscored_D1untagged=temp-repmat(mean(temp(:,timesD1un>=ZeroAt(1) & timesD1un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D1untagged=temp;

useUnits=useUnits_D1 & nanmean([D1temp D1temp_2ndwin],2)>0;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
% figure(); 
popmedian=median(mean(temp(useUnits,timesD1<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD1un>=0.2 & timesD1un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD1un<-1))>popmedian
    if turnsOn_D1(useUnits)
        c='k';
    else
        c='r';
    end
    plot(timesD1,currtoplot+offset,'Color',c); hold on;
end
% if isD1sess
%     addToTit='& this is D1 sess';
% else
%     addToTit='actually this is D2 sess';
% end
% title(['D1 ' addToTit]);

useUnits=useUnits_D1un & nanmean([D1untemp D1untemp_2ndwin],2)>0;
temp=Zscored_D1untagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
% figure(); 
popmedian=median(mean(temp(useUnits,timesD1un<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD2un>=0.2 & timesD2un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD2un<-1))>popmedian
    if turnsOn_D1un(useUnits)
        c='c';
    else
        c='y';
    end
    plot(timesD1un,currtoplot+offset,'Color',c); hold on;
end
% title(['D1 untagged ' addToTit]);

% Get projections

% For success
projTimeWindow1=[-1 0]; %[-2 -0.22]; %[2.25 3.56]-3.2476;
projTimeWindow2=[0 5]; %[3.56 7.26]-3.2476;
projTimeWindow3=[-2.21859 -1.46859];
projTimeWindow4=[0.781415 1.78141];
timecutoff=4;
temp=activityD1tagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD1<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD1<timecutoff),0,2,'omitnan');
Zscored_D1tagged=temp;
temp=activityD1untagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD1un<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD1un<timecutoff),0,2,'omitnan');
Zscored_D1untagged=temp;
timeBinsStep=0.1;
takeTimeRange=-4:0.25:9.5;
gausssmooth=1;
useUnits=useUnits_D1 & nanmean([D1temp D1temp_2ndwin],2)>0;
atLeastNUnits=1;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1tagged(useUnits,timesD1<9.5), timesD1(timesD1<9.5), projTimeWindow1, timeBinsStep, atLeastNUnits);
D1ontoCD.projTimeWindow1.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1ontoCD.projTimeWindow1.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1ontoCD.projTimeWindow1.vec=vec;
D1ontoCD.projTimeWindow1.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1tagged(useUnits,timesD1<9.5), timesD1(timesD1<9.5), projTimeWindow2, timeBinsStep, atLeastNUnits);
D1ontoCD.projTimeWindow2.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1ontoCD.projTimeWindow2.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1ontoCD.projTimeWindow2.vec=vec;
D1ontoCD.projTimeWindow2.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1tagged(useUnits,timesD1<9.5), timesD1(timesD1<9.5), projTimeWindow3, timeBinsStep, atLeastNUnits);
D1ontoCD.projTimeWindow3.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1ontoCD.projTimeWindow3.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1ontoCD.projTimeWindow3.vec=vec;
D1ontoCD.projTimeWindow3.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1tagged(useUnits,timesD1<9.5), timesD1(timesD1<9.5), projTimeWindow4, timeBinsStep, atLeastNUnits);
D1ontoCD.projTimeWindow4.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1ontoCD.projTimeWindow4.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1ontoCD.projTimeWindow4.vec=vec;
D1ontoCD.projTimeWindow4.mfr=mfr;

useUnits=useUnits_D1un & nanmean([D1untemp D1untemp_2ndwin],2)>0;
atLeastNUnits=1;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1untagged(useUnits,timesD1un<9.5), timesD1un(timesD1un<9.5), projTimeWindow1, timeBinsStep, atLeastNUnits);
D1unontoCD.projTimeWindow1.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1unontoCD.projTimeWindow1.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1unontoCD.projTimeWindow1.vec=vec;
D1unontoCD.projTimeWindow1.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1untagged(useUnits,timesD1un<9.5), timesD1un(timesD1un<9.5), projTimeWindow2, timeBinsStep, atLeastNUnits);
D1unontoCD.projTimeWindow2.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1unontoCD.projTimeWindow2.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1unontoCD.projTimeWindow2.vec=vec;
D1unontoCD.projTimeWindow2.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1untagged(useUnits,timesD1un<9.5), timesD1un(timesD1un<9.5), projTimeWindow3, timeBinsStep, atLeastNUnits);
D1unontoCD.projTimeWindow3.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1unontoCD.projTimeWindow3.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1unontoCD.projTimeWindow3.vec=vec;
D1unontoCD.projTimeWindow3.mfr=mfr;
[Rs,tb,vec,mfr]=getProjection(Zscored_D1untagged(useUnits,timesD1un<9.5), timesD1un(timesD1un<9.5), projTimeWindow4, timeBinsStep, atLeastNUnits);
D1unontoCD.projTimeWindow4.t=takeVecValsAtTimes(tb,tb,takeTimeRange);
D1unontoCD.projTimeWindow4.R=normVecMinus1to1(takeVecValsAtTimes(smoothdata(Rs,'gaussian',gausssmooth),tb,takeTimeRange),takeTimeRange<4);
D1unontoCD.projTimeWindow4.vec=vec;
D1unontoCD.projTimeWindow4.mfr=mfr;

end

function temp=normVecMinus1to1(temp,useindfornorm)

% temp=temp-nanmin(temp(useindfornorm));
% temp=temp./nanmax(temp(useindfornorm));
% temp=(temp.*2)-1;

end

function newVec=takeVecValsAtTimes(vec,vecTimes,newTimes)

newVec=nan(1,length(newTimes));
for i=1:length(newTimes)
    [~,mi]=nanmin(abs(newTimes(i)-vecTimes));
    if mi>length(vec)
        mi=length(vec);
    end
    newVec(i)=vec(mi);
end

end

function mfr=getAverageOfRandomProjections(vec,allunitsalltimes)

% just make this the mean firing rate across units
mfr=mean(vec-mean(allunitsalltimes,2,'omitnan'),1,'omitnan');

end

function [Rs,timeBins,popD1,mfr]=getProjection(ZscoredData, times, timeWindow, timeBinsStep, atLeastNUnits)

subtractOffShuffleCorr=false;

for i=1:size(ZscoredData,1)
    ZscoredData(i,:)=smoothdata(ZscoredData(i,:),'gaussian',20);
end
% ZscoredData=ZscoredData-repmat(nanmean(ZscoredData(:,times<8),2),1,size(ZscoredData,2));
% ZscoredData=ZscoredData-repmat(nanmin(ZscoredData(:,times<8),[],2),1,size(ZscoredData,2));
% ZscoredData=ZscoredData./repmat(nanmax(ZscoredData(:,times<8),[],2),1,size(ZscoredData,2));

popD1=mean(ZscoredData(:,times>timeWindow(1) & times<timeWindow(2)),2,'omitnan');
popD1(isnan(popD1))=0;
popD1=popD1-nanmin(popD1);
popD1=popD1./nanmax(popD1);
% average pop vector in this time window is popD1
mfr=getAverageOfRandomProjections(popD1,ZscoredData);
timeBins=times(1):timeBinsStep:times(end);
Rs=nan(1,length(timeBins)-1);
Ps=nan(1,length(timeBins)-1);
shuffle_Rs=nan(1,length(timeBins)-1);
if size(ZscoredData,1)<atLeastNUnits % must be at least units
    return
end
for i=1:length(timeBins)-1
    currpopD1=mean(ZscoredData(:,times>timeBins(i) & times<timeBins(i+1)),2,'omitnan');
    currpopD1(isnan(currpopD1))=0;
    currpopD1=currpopD1-nanmin(currpopD1);
    currpopD1=currpopD1./nanmax(currpopD1);
%     R=corrcoef(popD1-nanmean(popD1),currpopD1-nanmean(currpopD1));
    [R,P,RL,RU]=corrcoef(popD1,currpopD1,'Alpha',0.05);
    if length(popD1)==1
        Rs(i)=nan; % not enough units
        Ps(i)=nan;
    else
        Rs(i)=R(1,2);
        Ps(i)=P(1,2);
    end
%     Rs(i)=dot(popD1,currpopD1);
    if subtractOffShuffleCorr
        subShuffleRs=nan(1,50);
        for j=1:50
            R=corrcoef(popD1,currpopD1(randperm(length(currpopD1))));
            if length(popD1)==1
                subShuffleRs(j)=nan; % not enough units
            else
                subShuffleRs(j)=R(1,2);
            end
        end
        shuffle_Rs(i)=nanmean(subShuffleRs);
    end
end

if subtractOffShuffleCorr
    Rs=Rs-shuffle_Rs;
end

end

function morePlots(activityD1tagged,activityD1untagged,activityD2tagged,activityD2untagged)

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD1untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1untagged.aligncomp_y,1));
timesD1un=nanmean(activityD1untagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2tagged.aligncomp_y,1));
timesD2=nanmean(activityD2tagged.unitbyunit_x,1)-temp(f);
temp=nanmean(activityD2untagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD2untagged.aligncomp_y,1));
timesD2un=nanmean(activityD2untagged.unitbyunit_x,1)-temp(f);
% For success
takewin1=[2.25 3.56]-3.2476; % relative to peak of alignment companion
takewin2=[3.56 7.56]-3.2476;
% For failure, these windows interesting -- from autocorrelation matrix
% takewin1=[-2.21859 -1.46859]; % relative to peak of alignment companion
% takewin2=[0.781415 1.78141];

D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)];
D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin2(1) & timesD1<=takewin2(2)),2)];
D2temp=[nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin1(1) & timesD2<=takewin1(2)),2)];
D2temp_2ndwin=[nanmean(activityD2tagged.unitbyunit_y(:,timesD2>=takewin2(1) & timesD2<=takewin2(2)),2)];
D1temp(D1temp<0.0001)=0;
D1temp_2ndwin(D1temp_2ndwin<0.0001)=0;
D2temp(D2temp<0.0001)=0;
D2temp_2ndwin(D2temp_2ndwin<0.0001)=0;

D1untemp=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin1(1) & timesD1un<=takewin1(2)),2)];
D1untemp_2ndwin=[nanmean(activityD1untagged.unitbyunit_y(:,timesD1un>=takewin2(1) & timesD1un<=takewin2(2)),2)];
D2untemp=[nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>=takewin1(1) & timesD2un<=takewin1(2)),2)];
D2untemp_2ndwin=[nanmean(activityD2untagged.unitbyunit_y(:,timesD2un>=takewin2(1) & timesD2un<=takewin2(2)),2)];
D1untemp(D1untemp<0.0001)=0;
D1untemp_2ndwin(D1untemp_2ndwin<0.0001)=0;
D2untemp(D2untemp<0.0001)=0;
D2untemp_2ndwin(D2untemp_2ndwin<0.0001)=0;

[n,x]=hist(D1untemp_2ndwin./D1untemp,0:0.1:10);
[n_D1,x_D1]=cityscape_hist(n,x);
[n,x]=hist(D2untemp_2ndwin./D2untemp,0:0.1:10);
[n_D2,x_D2]=cityscape_hist(n,x);
figure(); plot(x_D1,n_D1./nansum(n_D1),'Color','k'); hold on; plot(x_D2,n_D2./nansum(n_D2),'Color','r'); legend({'D1','D2'});

turnsOn_D1=(D1temp_2ndwin./D1temp)>0.5;
turnsOn_D2=(D2temp_2ndwin./D2temp)>0.5;
ZeroAt=takewin1;
temp=activityD1tagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD1<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD1<9.5),0,2,'omitnan');
% Zscored_D1untagged=temp-repmat(mean(temp(:,timesD1un>=ZeroAt(1) & timesD1un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D1tagged=temp;
temp=activityD2tagged.unitbyunit_y;
% temp=temp-mean(temp(:,timesD2<9.5),2,'omitnan');
temp=temp./std(temp(:,timesD2<9.5),0,2,'omitnan');
% Zscored_D2untagged=temp-repmat(mean(temp(:,timesD2un>=ZeroAt(1) & timesD2un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D2tagged=temp;

useUnits=nanmean([D1temp D1temp_2ndwin],2)>0.25;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD1<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD1un>=0.2 & timesD1un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD1un<-1))>popmedian
    if turnsOn_D1(useUnits)
        c='k';
    else
        c='r';
    end
    plot(timesD1,currtoplot,'Color',c); hold on;
end
title('D1');

useUnits=nanmean([D2temp D2temp_2ndwin],2)>0.25;
temp=Zscored_D2tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD2<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD2un>=0.2 & timesD2un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD2un<-1))>popmedian
    if turnsOn_D2(useUnits)
        c='k';
    else
        c='r';
    end
    plot(timesD2,currtoplot,'Color',c); hold on;
end
title('D2 tagged');

% unit color coding scheme
turnsOn_D1un=(D1untemp_2ndwin./D1untemp)>0.5;
turnsOn_D2un=(D2untemp_2ndwin./D2untemp)>0.5;

% Zscore within each unit, then plot
% transition to time window 2
ZeroAt=takewin1;
temp=activityD1untagged.unitbyunit_y;
% temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
% Zscored_D1untagged=temp-repmat(mean(temp(:,timesD1un>=ZeroAt(1) & timesD1un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D1untagged=temp;
temp=activityD2untagged.unitbyunit_y;
% temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
% Zscored_D2untagged=temp-repmat(mean(temp(:,timesD2un>=ZeroAt(1) & timesD2un<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D2untagged=temp;

useUnits=nanmean([D1untemp D1untemp_2ndwin],2)>0.25;
temp=Zscored_D1untagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD1un<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD1un>=0.2 & timesD1un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD1un<-1))>popmedian
    if turnsOn_D1un(useUnits)
        c='k';
    else
        c='r';
    end
    plot(timesD1un,currtoplot,'Color',c); hold on;
end
title('D1 untagged');

useUnits=nanmean([D2untemp D2untemp_2ndwin],2)>0.25;
temp=Zscored_D2untagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',20);
end
useUnitss=find(useUnits);
figure(); 
popmedian=median(mean(temp(useUnits,timesD2un<-1),2,'omitnan'),1,'omitnan');
for i=1:length(useUnitss)
    useUnits=useUnitss(i);
%     currtoplot=(temp(useUnits,:)-repmat(nanmean(temp(useUnits,timesD2un>=0.2 & timesD2un<0.3),2),1,size(temp,2)))';
    currtoplot=temp(useUnits,:)';
%     currtoplot=currtoplot./nanmax(abs(currtoplot));
%     if nanmean(temp(useUnits,timesD2un<-1))>popmedian
    if turnsOn_D2un(useUnits)
        c='k';
    else
        c='r';
    end
    plot(timesD2un,currtoplot,'Color',c); hold on;
end
title('D2 untagged');

% Population vector after the outcome
% Which pre-outcome time points correlate best
timeBinsStep=0.25;
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D1tagged=temp;
temp=activityD2tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D2tagged=temp;

temp=activityD1untagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D1untagged=temp;
temp=activityD2untagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
Zscored_D2untagged=temp;

popD1=mean(Zscored_D1untagged(:,timesD1un>0.3 & timesD1un<4),2,'omitnan');
popD1(isnan(popD1))=0;
timeBins=timesD1un(1):timeBinsStep:timesD1un(end);
Rs=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD1=mean(Zscored_D1untagged(:,timesD1un>timeBins(i) & timesD1un<timeBins(i+1)),2,'omitnan');
    currpopD1(isnan(currpopD1))=0;
    R=corrcoef(popD1,currpopD1);
    Rs(i)=R(1,2);
end
figure(); 
plot(timeBins(1:end-1),Rs,'Color','k'); title('D1 un v D2 untagged corrcoef with post-outcome driven pop vec');

popD2=mean(Zscored_D2untagged(:,timesD2un>0.3 & timesD2un<4),2,'omitnan');
popD2(isnan(popD2))=0;
timeBins=timesD2un(1):timeBinsStep:timesD2un(end);
Rs=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD2=mean(Zscored_D2untagged(:,timesD2un>timeBins(i) & timesD2un<timeBins(i+1)),2,'omitnan');
    currpopD2(isnan(currpopD2))=0;
    R=corrcoef(popD2,currpopD2);
    Rs(i)=R(1,2);
end
hold on; 
plot(timeBins(1:end-1),Rs,'Color','r'); 

end


function [D1tagged_cueResponse,D1orD2taggingExpt,putAlignPeakAt,firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess]=getAndSaveResponse(dd_more,getThese,settings,putAlignPeakAtInput,indsForAnalysisPerSess)

[~,~,~,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns,D1orD2taggingExpt,firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess,fromWhichSess]=alignToCompanion(dd_more,false,[],[],[],getThese,settings,indsForAnalysisPerSess);
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
D1tagged_cueResponse.fromWhichSess=D1tagged_cueResponse.fromWhichSess(eitherExcluded==0);

activityD1tagged.unitbyunit_x=activityD1tagged.unitbyunit_x(eitherExcluded==0,:);
activityD1tagged.unitbyunit_y=activityD1tagged.unitbyunit_y(eitherExcluded==0,:);
activityD1tagged.aligncomp_x=activityD1tagged.aligncomp_x(eitherExcluded==0,:);
activityD1tagged.aligncomp_y=activityD1tagged.aligncomp_y(eitherExcluded==0,:);
activityD1tagged.excluded=activityD1tagged.excluded(eitherExcluded==0);
activityD1tagged.ns=activityD1tagged.ns(eitherExcluded==0);
activityD1tagged.fromWhichSess=activityD1tagged.fromWhichSess(eitherExcluded==0);

end
