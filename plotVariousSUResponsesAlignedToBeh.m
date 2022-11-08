function out=plotVariousSUResponsesAlignedToBeh(varargin)

if length(varargin)<2
    error('plotVariousSUResponsesAlignedToBeh takes at least 2 input arguments');
end

Response=varargin{2};

out=[];

switch varargin{1}
    case 'meanAcrossUnits'
        downSampFac=varargin{3};
        out=meanAcrossUnits(Response,downSampFac);
    case 'scatterInTimeWindows'
        timeWindow1=varargin{3};
        timeWindow2=varargin{4};
        [out.timeWindow1_FR,out.timeWindow2_FR]=scatterInTimeWindows(Response,timeWindow1,timeWindow2);
    case 'modulationIndex'
        timeWindow1=varargin{3};
        timeWindow2=varargin{4};
        out.modIndex=plotModulationIndex(Response,timeWindow1,timeWindow2);
    case 'ZscoredAndSmoothed'
        timeWindow1=varargin{3};
        timeWindow2=varargin{4};
        ZscoredAndSmoothed(Response,timeWindow1,timeWindow2);
    case 'populationVectorCorrelation'
        timeBinsStep=varargin{3};
        sliceAtTimeWindow=varargin{4};
        populationVectorCorrelation(Response,timeBinsStep,sliceAtTimeWindow);
    case 'trialVectorCorrelation'
        timeBinsStep=varargin{3};
        sliceAtTimeWindow=varargin{4};
        trialVectorCorrelation_wrapper(Response,timeBinsStep,sliceAtTimeWindow);
        % for a single unit
        % trialVectorCorrelation(filterResponseToOneSU(Response,1),timeBinsStep,sliceAtTimeWindow);
    case 'scatterResponseVsResponse'
        Response2=varargin{3};
        typeOfPlot=varargin{4};
        [Response,Response2]=matchUnitsAcrossResponses(Response,Response2);
        outResponse3=[];
        if any(strcmp([varargin],'ColorLabel'))
            Response3=varargin{find(strcmp([varargin],'ColorLabel'))+1};
            Response3=matchUnitsAcrossResponses(Response3,Response);
            outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
            outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
            if length(varargin)==find(strcmp([varargin],'ColorLabel'))+1
                outResponse3=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response3,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
            else
                outResponse3=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response3,varargin{find(strcmp([varargin],'ColorLabel'))+2:end});
            end
        else
            outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{5:end});
            outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{5:end});
        end
        switch typeOfPlot
            case 'meanAcrossUnits'
                overlayPlots(outResponse1,outResponse2);
            case 'scatterInTimeWindows'
                if ~isempty(outResponse3)
                    col=outResponse3.timeWindow2_FR./outResponse3.timeWindow1_FR;
                else
                    col=[];
                end
                scatterTwoResponses(outResponse1.timeWindow1_FR,outResponse2.timeWindow1_FR,col); title('Time window 1 firing rates');
                scatterTwoResponses(outResponse1.timeWindow2_FR,outResponse2.timeWindow2_FR,col); title('Time window 2 firing rates');
                scatterTwoResponses(outResponse1.timeWindow2_FR./outResponse1.timeWindow1_FR,outResponse2.timeWindow2_FR./outResponse2.timeWindow1_FR,col); title('Time window 2 to 1 ratio');
            case 'modulationIndex'
                if ~isempty(outResponse3)
                    col=outResponse3.modIndex;
                else
                    col=[];
                end
                scatterTwoResponses(outResponse1.modIndex,outResponse2.modIndex,col); title('Modulation index');
            case 'populationVectorCorrelation'
            case 'trialVectorCorrelation'
            otherwise
                disp(['do not recognize scatterResponseVsResponse for plot type ' typeOfPlot]);
        end
    otherwise
        error('Do not recognize first argument value in plotVariousSUResponsesAlignedToBeh.m');
end

end

function overlayPlots(r1,r2)

figure();
plot(r1.t,r1.me,'Color','k'); hold on;
plot(r1.t,r1.plusSe,'Color','k');
plot(r1.t,r1.minusSe,'Color','k');

plot(r2.t,r2.me,'Color','r'); hold on;
plot(r2.t,r2.plusSe,'Color','r');
plot(r2.t,r2.minusSe,'Color','r');

end

function scatterTwoResponses(r1,r2,c)

figure();
if isempty(c)
    scatter(r1,r2);
else
    scatter(r1,r2,[],c);
end
colormap('jet');
xlabel('Response 1');
ylabel('Response 2');

end

function [Response,Response2]=matchUnitsAcrossResponses(Response,Response2)

Response1_indsIntoExcluded=find(Response.excluded==0);
Response2_indsIntoExcluded=find(Response2.excluded==0);
useTheseUnitsResponse1=find(ismember(Response1_indsIntoExcluded,Response2_indsIntoExcluded));
useTheseUnitsResponse2=find(ismember(Response2_indsIntoExcluded,Response1_indsIntoExcluded));
Response=filterResponseToOneSU(Response,useTheseUnitsResponse1);
Response2=filterResponseToOneSU(Response2,useTheseUnitsResponse2);

end

function out=filterResponseToOneSU(Response,whichUnitToUse)

f=fieldnames(Response);
fieldLikeResponseSize=size(Response.unitbyunit_y,1);
if length(Response.fromWhichUnit)==fieldLikeResponseSize
    % took all trials
    whichUnit=Response.fromWhichUnit;
elseif length(Response.fromWhichSess)==fieldLikeResponseSize
    % took unit by unit
    whichUnit=1:fieldLikeResponseSize;
else
    error('do not recognize structure of Response in plotVariousSUResponsesAlignedToBeh.m');
end

for i=1:length(f)
    temp=Response.(f{i});
    if size(temp,1)==fieldLikeResponseSize
        % filter this field
        if length(size(temp))>1
            % 2D
            out.(f{i})=temp(ismember(whichUnit,whichUnitToUse),:);
        else
            % 1D
            out.(f{i})=temp(ismember(whichUnit,whichUnitToUse));
        end
    else
        % don't filter this field
        out.(f{i})=temp;
    end
end

end

function out=meanAcrossUnits(activityD1tagged,downSampFac)

if size(activityD1tagged.unitbyunit_y,1)>100
    % skip unit by unit plot because too crowded
else
    figure();
    plot(nanmean(activityD1tagged.unitbyunit_x,1),activityD1tagged.unitbyunit_y');
    hold on; plot(nanmean(activityD1tagged.aligncomp_x,1),nanmean(activityD1tagged.aligncomp_y,1),'Color','b');
    xlabel('Time (sec)');
    ylabel('Firing rate');
    title('With raw times');
end

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end
if size(activityD1tagged.unitbyunit_y,1)>100
    % skip unit by unit plot because too crowded
else
    figure(); plot(timesD1,activityD1tagged.unitbyunit_y');
    title('Align companion at time=0');
    xlabel('Time (sec)');
    ylabel('Firing rate');
end

[out.me,out.plusSe,out.minusSe,out.t]=plotMeanAndSE(downSampMatrix(activityD1tagged.unitbyunit_y,downSampFac),'k',downSampAv(timesD1,downSampFac));
title('Mean and se across units');
xlabel('Time (sec)');
ylabel('Firing rate');

end

function [me,plusSe,minusSe,t]=plotMeanAndSE(data1,color1,timeBins)

plotAsCityscape=false;

data1_mean=mean(data1,1,'omitnan');
data1_se=std(data1,0,1,'omitnan')./sqrt(size(data1,1));

figure();
if plotAsCityscape==true
    [n,x]=cityscape_hist(data1_mean,timeBins);
    plot(x,n,'Color',color1); hold on;
    t=x;
    me=n;
    [n,x]=cityscape_hist(data1_mean+data1_se,timeBins);
    plot(x,n,'Color',color1);
    plusSe=n;
    [n,x]=cityscape_hist(data1_mean-data1_se,timeBins);
    plot(x,n,'Color',color1);
    minusSe=n;
else
    plot(timeBins,data1_mean,'Color',color1); hold on;
    plot(timeBins,data1_mean+data1_se,'Color',color1);
    plot(timeBins,data1_mean-data1_se,'Color',color1);
    t=timeBins;
    me=data1_mean;
    plusSe=data1_mean+data1_se;
    minusSe=data1_mean-data1_se;
end

end

function [D1temp,D1temp_2ndwin]=scatterInTimeWindows(activityD1tagged,takewin1,takewin2)

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

figure(); 
scatter(nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2),nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2)./nanmax(activityD1tagged.unitbyunit_y(:,1:200),[],2),[],'k');
xlabel('Peak-normed FR in time window 1');
ylabel('Peak-normed FR in time window 2');

figure();
scatter(nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2),nanmean(activityD1tagged.unitbyunit_y(:,timesD1>takewin2(1) & timesD1<=takewin2(2)),2),[],'k');
xlabel('FR in time window 1');
ylabel('FR in time window 2');


D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)];
D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin2(1) & timesD1<=takewin2(2)),2)];
D1temp(D1temp<0.0001)=0;
D1temp_2ndwin(D1temp_2ndwin<0.0001)=0;
[n,x]=hist(D1temp_2ndwin./D1temp,0:0.1:10);
[n_D1,x_D1]=cityscape_hist(n,x);
figure(); 
plot(x_D1,n_D1,'Color','k'); 
xlabel('Ratio of FR in win 2 to FR in win 1');
ylabel('Count');
figure(); 
plot(x_D1,n_D1./nansum(n_D1),'Color','k'); 
xlabel('Ratio of FR in win 2 to FR in win 1');
ylabel('Count integral-normed');

end

function ZscoredAndSmoothed(activityD1tagged,takewin1,takewin2)

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

D1temp=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin1(1) & timesD1<=takewin1(2)),2)];
D1temp_2ndwin=[nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=takewin2(1) & timesD1<=takewin2(2)),2)];
D1temp(D1temp<0.0001)=0;
D1temp_2ndwin(D1temp_2ndwin<0.0001)=0;

% Zscore within each unit, then plot transition to time window 2
ZeroAt=takewin1;
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp,0,2,'omitnan');
% Zscored_D1tagged=temp-repmat(mean(temp(:,timesD1>=ZeroAt(1) & timesD1<=ZeroAt(2)),2,'omitnan'),1,size(temp,2));
Zscored_D1tagged=temp;

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
xlabel('Time (sec)');
ylabel('Zscored FR');
title('Zscored and Gaussian-smoothed');

end

function populationVectorCorrelation(activityD1tagged,timeBinsStep,sliceAtTimeWindow)

% Zscore within each unit
% Then show cov(A,B)/(std(A)*std(B))
% where A and B are vectors of the population activity across all units at
% time points a and b

% timesD1 is times relative to alignment companion peak
temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

% Population vector after the outcome (or whatever is in alignment
% companion)
% Which pre-outcome (or pre-alignment companion) time points correlate best
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan'); % subtract off each unit's own mean
temp=temp./std(temp,0,2,'omitnan'); % Z-score each unit's own activity
Zscored_D1tagged=temp;

% Matrix
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
figure(); 
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),Rs(timeBins<9.5,timeBins<9.5)); 
title('corrcoef matrix');

% Plot a slice through this matrix
popD1=mean(Zscored_D1tagged(:,timesD1>sliceAtTimeWindow(1) & timesD1<sliceAtTimeWindow(2)),2,'omitnan');
popD1(isnan(popD1))=0;
%timeBinsStep=0.25;
timeBins=timesD1(1):timeBinsStep:timesD1(end);
Rs=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD1=mean(Zscored_D1tagged(:,timesD1>timeBins(i) & timesD1<timeBins(i+1)),2,'omitnan');
    currpopD1(isnan(currpopD1))=0;
    R=corrcoef(popD1,currpopD1);
    Rs(i)=R(1,2);
end
figure(); 
plot(timeBins(1:end-1),Rs,'Color','k'); 
title(['Slice at time window ' num2str(sliceAtTimeWindow(1)) ' to ' num2str(sliceAtTimeWindow(2)) ' seconds']);

end

function C=addMatricesIgnoreNans(A,B)

tmp = cat(3,A,B); C = nansum(tmp,3);

end

function trialVectorCorrelation_wrapper(Response,timeBinsStep,sliceAtTimeWindow)

suppressPlots=true;

u=unique(Response.fromWhichUnit);
allRs=[];
allRs_slice=[];
% timesD1 is times relative to alignment companion peak
temp=nanmean(Response.aligncomp_x,1);
[~,f]=nanmax(nanmean(Response.aligncomp_y,1));
timesD1=nanmean(Response.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end
disp(['averaging this many units: ' num2str(length(u))]);
for i=1:length(u)
    if mod(i,50)==0
        disp(['on unit ' num2str(i)]);
    end
    [Rs,Rs_slice]=trialVectorCorrelation(filterResponseToOneSU(Response,u(i)),timeBinsStep,sliceAtTimeWindow,suppressPlots,timesD1);
    if all(isnan(Rs)) | isempty(Rs)
        continue
    end
    if  isempty(allRs)
        allRs=Rs;
        allRs_slice=Rs_slice;
    else
        allRs=addMatricesIgnoreNans(allRs,Rs);
        allRs_slice=addMatricesIgnoreNans(allRs_slice,Rs_slice);
    end
end
allRs=allRs./length(u);
allRs_slice=allRs_slice./length(u);

timeBins=timesD1(1):timeBinsStep:timesD1(end);
figure();
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),allRs(timeBins<9.5,timeBins<9.5));
title('corrcoef matrix average across units');
fillVal=mean(allRs,'all','omitnan');
for i=1:size(allRs,1)
    for j=i-2:i+2
        if j>=1 && j<=size(allRs,1)
            allRs(i,j)=fillVal;
        end
    end
end
figure();
imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),allRs(timeBins<9.5,timeBins<9.5));
title('corrcoef matrix average across units -- filling in near diagnonal');
figure();
plot(timeBins(1:end-1),nanmean(allRs_slice,1),'Color','k');
title(['Slice at time window ' num2str(sliceAtTimeWindow(1)) ' to ' num2str(sliceAtTimeWindow(2)) ' seconds - av across units']);

end

function [Rs,Rs_slice]=trialVectorCorrelation(activityD1tagged,timeBinsStep,sliceAtTimeWindow,suppressPlots,timesD1)

% Zscore using mean within each trial, std across all trials
% Then show cov(A,B)/(std(A)*std(B))
% where A and B are vectors of the single unit's activity across trials at
% time points a and b

Rs=[];
Rs_slice=[];

s=settingsForStriatumUnitPlots();
if s.keepAllSingleTrials~=true
    disp('Trial vector correlation requires single trial data ... cannot continue');
    return
end

% timesD1 is times relative to alignment companion peak
if isempty(timesD1)
    temp=nanmean(activityD1tagged.aligncomp_x,1);
    [~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
    timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
    if any(isnan(timesD1))
        timesD1=fillmissing(timesD1,'linear');
    end
end

% Population vector after the outcome (or whatever is in alignment
% companion)
% Which pre-outcome (or pre-alignment companion) time points correlate best
temp=activityD1tagged.unitbyunit_y;
temp=temp-mean(temp,'all','omitnan'); % subtract off mean across all trials, all timepoints
temp=temp./std(temp,0,'all','omitnan'); % Z-score using std across all trials, all timepoints
Zscored_D1tagged=temp;

% Matrix
timeBins=timesD1(1):timeBinsStep:timesD1(end);
Rs=nan(length(timeBins)-1,length(timeBins)-1);
for i=1:length(timeBins)-1
    popD1=mean(Zscored_D1tagged(:,timesD1>timeBins(i) & timesD1<timeBins(i+1)),2,'omitnan');
    popD1(isnan(popD1))=0;
    for j=1:length(timeBins)-1
        currpopD1=mean(Zscored_D1tagged(:,timesD1>timeBins(j) & timesD1<timeBins(j+1)),2,'omitnan');
        currpopD1(isnan(currpopD1))=0;
        R=corrcoef(popD1,currpopD1);
        if size(R,1)<2
            continue
        else
            Rs(i,j)=R(1,2);
        end
    end
end
for i=1:length(timeBins)-1
    Rs(i,i)=mean(Rs,'all','omitnan');
end
if suppressPlots==false
    figure();
    imagesc(timeBins(timeBins<9.5),timeBins(timeBins<9.5),Rs(timeBins<9.5,timeBins<9.5));
    title('corrcoef matrix');
end

% Plot a slice through this matrix
popD1=mean(Zscored_D1tagged(:,timesD1>sliceAtTimeWindow(1) & timesD1<sliceAtTimeWindow(2)),2,'omitnan');
popD1(isnan(popD1))=0;
%timeBinsStep=0.25;
timeBins=timesD1(1):timeBinsStep:timesD1(end);
Rs_slice=nan(1,length(timeBins)-1);
for i=1:length(timeBins)-1
    currpopD1=mean(Zscored_D1tagged(:,timesD1>timeBins(i) & timesD1<timeBins(i+1)),2,'omitnan');
    currpopD1(isnan(currpopD1))=0;
    R=corrcoef(popD1,currpopD1);
    if size(R,1)<2
        continue
    else
        Rs_slice(i)=R(1,2);
    end
end
if suppressPlots==false
    figure();
    plot(timeBins(1:end-1),Rs_slice,'Color','k');
    title(['Slice at time window ' num2str(sliceAtTimeWindow(1)) ' to ' num2str(sliceAtTimeWindow(2)) ' seconds']);
end

end

function modIndex=plotModulationIndex(activityD1tagged,inWindow1,inWindow2)

% inWindow in seconds relative to alignmentCompanion peak

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

inWindow1=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=inWindow1(1) & timesD1<=inWindow1(2)),2);
inWindow2=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>inWindow2(1) & timesD1<=inWindow2(2)),2);

modIndex=(inWindow2-inWindow1)./(inWindow2+inWindow1);

modBins=[-1-0.1:0.1:1]+0.05;
[N,edges]=histcounts(modIndex,modBins);
binCenters=edges(1:end-1)+0.05;
figure();
bar(binCenters,N,'LineWidth',0.1);
xlabel('Modulation index');
ylabel('Count');

end