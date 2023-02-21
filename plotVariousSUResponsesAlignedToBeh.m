function out=plotVariousSUResponsesAlignedToBeh(varargin)

if length(varargin)<2
    error('plotVariousSUResponsesAlignedToBeh takes at least 2 input arguments');
end

Response=varargin{2};

out=[];

switch varargin{1}
    case 'modRatioPSTH'
        col=[];
        Response2=varargin{3};
        Response3=varargin{4};
        [Response,Response2,Response3]=matchUnitsAcrossResponses(Response,Response2,Response3);
        outResponse1=plotVariousSUResponsesAlignedToBeh('modulationIndex',Response,varargin{5:end});
        outResponse2=plotVariousSUResponsesAlignedToBeh('modulationIndex',Response2,varargin{5:end});
        outResponse3=plotVariousSUResponsesAlignedToBeh('modulationIndex',Response3,[-1 -0.5],[-0.4 1]);
        % col=outResponse2.modIndex./outResponse1.modIndex;
        % figure(); scatter(outResponse1.modIndex,outResponse2.modIndex,[],'k');
        div=false;
        if div==true
            disp(['modRatio is y./x']);
            col=outResponse2.modIndex./outResponse1.modIndex;
            colRelativeTo=1;
        else
            disp(['modRatio is y-x']);
            col=outResponse2.modIndex-outResponse1.modIndex;
            colRelativeTo=0;
        end
        modRatioPSTH(Response,Response2,Response3,col,colRelativeTo,varargin{5:end},outResponse1.modIndex,outResponse2.modIndex,outResponse3.modIndex);
    case 'coloredPSTH'
        % do not need to process anything
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
    case 'getResponsive'
        out.isResponsive=isResponsive(Response,varargin{3},varargin{4});
    case 'matchUnitsAcrossResponses'
        Response2=varargin{3};
        Response3=varargin{4};
        Response4=varargin{5};
        Response5=varargin{6};
        [out.Response1,out.Response2,out.Response3,out.Response4,out.Response5]=matchUnitsAcrossResponses(Response,Response2,Response3,Response4,Response5);
    case 'scatterResponseVsResponse'
        grayOutNonResponse=true; %true;
        Response2=varargin{3};
        typeOfPlot=varargin{4};
        [Response,Response2]=matchUnitsAcrossResponses(Response,Response2);
        outResponse3=[];
        isResp=[];
        if any(strcmp([varargin],'ColorLabel'))
            Response3=varargin{find(strcmp([varargin],'ColorLabel'))+1};
            [Response,Response2,Response3]=matchUnitsAcrossResponses(Response,Response2,Response3);
            outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
            outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
            if length(Response.fromWhichUnit)==size(Response.unitbyunit_y,1)
                % took trial by trial
                ubyu_Response3=Response3;
            else
                ubyu_Response3=collapseToUnitByUnit(Response3);
            end
            alwaysUseModForColor=true;
            if alwaysUseModForColor==true
                tOp='modulationIndex';
            else
                tOp=typeOfPlot;
            end
            ZscoreR3=false;
            if ZscoreR3==true
                timesD1=getTimesD1(ubyu_Response3);
                meanForZscore=mean([Response.unitbyunit_y(:,1:200) Response2.unitbyunit_y(:,1:200) ubyu_Response3.unitbyunit_y(:,ismember(1:200,find(timesD1>-1)))],2,'omitnan');
                sdForZscore=std([Response.unitbyunit_y(:,1:200) Response2.unitbyunit_y(:,1:200) ubyu_Response3.unitbyunit_y(:,ismember(1:200,find(timesD1>-1)))],0,2,'omitnan');
                ubyu_Response3=ZscoreDontSmooth(ubyu_Response3,meanForZscore,sdForZscore);
            end
            if length(varargin)==find(strcmp([varargin],'ColorLabel'))+1
                outResponse3=plotVariousSUResponsesAlignedToBeh(tOp,ubyu_Response3,varargin{5:find(strcmp([varargin],'ColorLabel'))-1});
                isResp=isResponsive(Response3,varargin{5:find(strcmp([varargin],'ColorLabel'))-1}); % find responsive for color plot
            else
                outResponse3=plotVariousSUResponsesAlignedToBeh(tOp,ubyu_Response3,varargin{find(strcmp([varargin],'ColorLabel'))+2:end});
                isResp=isResponsive(Response3,varargin{find(strcmp([varargin],'ColorLabel'))+2:end}); % find responsive for color plot
            end
        else
            if any(strcmp(varargin,'timeWindowsForResponse1'))
                % syntax
                % plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',cue_Response,cuedSuccess_Response,'scatterInTimeWindows','timeWindowsForResponse1',[-1 -0.5],[0 1],'timeWindowsForResponse2',[-3 0.25],[0.5 4]);
                outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{6:7});
                outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{9:10});
            else
                outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{5:end});
                outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{5:end});
            end
        end
        switch typeOfPlot
            case 'coloredPSTH'
                col=[];
                if ~isempty(outResponse3) 
                    if alwaysUseModForColor==true
                        temp=outResponse3.modIndex;
                        if ~isempty(isResp) && grayOutNonResponse==true
                            temp(isResp.isSig==0)=nan;
                        end
                        col=temp;
                    end
                end
                coloredPSTH(Response,Response2,col);
            case 'meanAcrossUnits'
                overlayPlots(outResponse1,outResponse2);
            case 'scatterInTimeWindows'
                if ~isempty(outResponse3)
                    if alwaysUseModForColor==true
                        temp=outResponse3.modIndex;
                    else
%                         temp=outResponse3.timeWindow2_FR./outResponse3.timeWindow1_FR;
%                         temp=outResponse3.timeWindow2_FR-outResponse3.timeWindow1_FR;
                        temp=log(outResponse3.timeWindow2_FR);
%                         temp=outResponse3.timeWindow2_FR;
%                         temp(temp>1)=1;

%                         temp=outResponse1.timeWindow1_FR-outResponse2.timeWindow1_FR;
%                         temp=temp-mean(temp,'all','omitnan');
%                         temp=temp./std(temp,1,'all','omitnan');
%                         temp(temp>1)=1; temp(temp<-1)=-1;
                    end
                    if ~isempty(isResp) && grayOutNonResponse==true
                        temp(isResp.isSig==false)=nan;
                    end
                    col=temp;
                else
                    col=[];
                end
                scatterTwoResponses(outResponse1.timeWindow1_FR,outResponse2.timeWindow1_FR,col); title('Time window 1 firing rates');
                scatterTwoResponses(outResponse1.timeWindow2_FR,outResponse2.timeWindow2_FR,col); title('Time window 2 firing rates');
                scatterTwoResponses(outResponse1.timeWindow2_FR./outResponse1.timeWindow1_FR,outResponse2.timeWindow2_FR./outResponse2.timeWindow1_FR,col); title('Time window 2 to 1 ratio');
            case 'modulationIndex'
                if ~isempty(outResponse3)
                    temp=outResponse3.modIndex;
                    if ~isempty(isResp) && grayOutNonResponse==true
                        temp(isResp.isSig==false)=nan;
                    end
                    col=temp;
                else
                    col=[];
                end
                scatterTwoResponses(outResponse1.modIndex,outResponse2.modIndex,col); title('Modulation index');
                out.response2minus1=outResponse2.modIndex-outResponse1.modIndex;
            case 'populationVectorCorrelation'
            case 'trialVectorCorrelation'
            case 'scatter3D'
            otherwise
                disp(['do not recognize scatterResponseVsResponse for plot type ' typeOfPlot]);
        end
    case 'scatter3D'
        Response2=varargin{3};
        Response3=varargin{4};
        [Response,Response2,Response3]=matchUnitsAcrossResponses(Response,Response2,Response3);
        typeOfPlot='modulationIndex';
        outResponse1=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response,varargin{5:end});
        outResponse2=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response2,varargin{5:end});
        outResponse3=plotVariousSUResponsesAlignedToBeh(typeOfPlot,Response3,varargin{5:end});
        figure();
        scatter3(outResponse1.modIndex,outResponse2.modIndex,outResponse3.modIndex);
        xlabel('Response 1');
        ylabel('Response 2');
        zlabel('Response 3');
        out.modIndex1=outResponse1.modIndex;
        out.modIndex2=outResponse2.modIndex;
        out.modIndex3=outResponse3.modIndex;
        out.meanFR1=nanmean(Response.unitbyunit_y,2);
        out.meanFR2=nanmean(Response2.unitbyunit_y,2);
        out.meanFR3=nanmean(Response3.unitbyunit_y,2);
    otherwise
        error('Do not recognize first argument value in plotVariousSUResponsesAlignedToBeh.m');
end

end

function out=collapseToUnitByUnit(Response)

fieldLikeResponseSize=size(Response.unitbyunit_y,1);
out=Response;
if length(Response.fromWhichUnit)==fieldLikeResponseSize
    % took all trials
    whichUnit=Response.fromWhichUnit;
    u=unique(whichUnit);
    out.unitbyunit_x=nan(length(u),size(Response.unitbyunit_x,2));
    out.unitbyunit_y=nan(length(u),size(Response.unitbyunit_y,2));
    out.aligncomp_x=nan(length(u),size(Response.aligncomp_x,2));
    out.aligncomp_y=nan(length(u),size(Response.aligncomp_y,2));
    for i=1:length(u)
        currU=u(i);
        out.unitbyunit_x(i,:)=mean(Response.unitbyunit_x(whichUnit==currU,:),1,'omitnan');
        out.unitbyunit_y(i,:)=mean(Response.unitbyunit_y(whichUnit==currU,:),1,'omitnan');
        out.aligncomp_x(i,:)=mean(Response.aligncomp_x(whichUnit==currU,:),1,'omitnan');
        out.aligncomp_y(i,:)=mean(Response.aligncomp_y(whichUnit==currU,:),1,'omitnan');
    end
elseif length(Response.fromWhichSess)==fieldLikeResponseSize
    % already took unit by unit
end

end

function timesD1=getTimesD1(Response1)

temp=nanmean(Response1.aligncomp_x,1);
[~,f]=nanmax(nanmean(Response1.aligncomp_y,1));
timesD1=nanmean(Response1.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

end

function Response1=ZscoreAndSmooth(varargin)

if length(varargin)==1
    Response1=varargin{1};
    meforZ=[];
    sdforZ=[];
elseif length(varargin)==3
    Response1=varargin{1};
    meforZ=varargin{2};
    sdforZ=varargin{3};
end

% Zscore within each unit, then plot transition to time window 2
temp=Response1.unitbyunit_y;
if isempty(meforZ)
    temp=temp-repmat(mean(temp,2,'omitnan'),1,size(temp,2));
    temp=temp./repmat(std(temp(:,1:200),0,2,'omitnan'),1,size(temp,2));
else
    temp=temp-repmat(meforZ,1,size(temp,2));
    temp=temp./repmat(sdforZ,1,size(temp,2));
end
Zscored_D1tagged=temp;

temp=Zscored_D1tagged;
for i=1:size(temp,1)
%     temp(i,:)=smoothdata(temp(i,:),'gaussian',2);
%     temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
    temp(i,:)=smoothdata(temp(i,:),'gaussian',7);
end
Response1.unitbyunit_y=temp;

end

function [Response,newCol,sinew,cmap,colBins,m1,m2]=binByCol(Response,col,binStep,mod1,mod2)

% colBins=min(col,[],'all','omitnan'):binStep:max(col,[],'all','omitnan')+0.01;
colBins=prctile(col,[0:binStep:100]);
temp=Response.unitbyunit_y;
out=nan(length(colBins)-1,size(temp,2));
newCol=nan(length(colBins)-1,1);
m1=nan(length(colBins)-1,1);
m2=nan(length(colBins)-1,1);
for i=1:length(colBins)-1
    out(i,:)=mean(temp(col>=colBins(i) & col<colBins(i+1),:),1,'omitnan');
    newCol(i)=mean(col(col>=colBins(i) & col<colBins(i+1)),'all','omitnan');
    m1(i)=mean(mod1(col>=colBins(i) & col<colBins(i+1)),'all','omitnan');
    m2(i)=mean(mod2(col>=colBins(i) & col<colBins(i+1)),'all','omitnan');
end
out=out(~isnan(newCol),:);
newCol=newCol(~isnan(newCol));
Response.unitbyunit_y=out;
[~,sinewp]=sort(newCol(~isnan(newCol)));
sinew=nan(size(sinewp));
for i=1:length(sinewp)
    sinew(i)=find(sinewp==i);
end
cmap=jet(length(newCol(~isnan(newCol))));

end

function plotModScatter(varargin)

if length(varargin)==5
    mod1=varargin{1};
    mod2=varargin{2};
    cmap=varargin{3};
    si=varargin{4};
    col=varargin{5};
elseif length(varargin)==6
    mod1=varargin{1};
    mod2=varargin{2};
    cmap=varargin{3};
    si=varargin{4};
    col=varargin{5};
    takeHalf=varargin{6};
    if strcmp(takeHalf,'takeTopHalf')
        cmap=jet(2*length(col(~isnan(col))));
        cmap=cmap(length(col(~isnan(col)))+1:end,:);
    elseif strcmp(takeHalf,'takeBottomHalf')
        cmap=jet(2*length(col(~isnan(col))));
        cmap=cmap(1:length(col(~isnan(col))),:);
    end
end

figure();
scatter(mod1(~isnan(col)),mod2(~isnan(col)),[],cmap(si,:));

end

function plotResponseLineByLine(varargin)

if length(varargin)==5
    timesD1=varargin{1};
    Response1=varargin{2};
    cmap=varargin{3};
    si=varargin{4};
    col=varargin{5};
    takeHalf=[];
elseif length(varargin)==6
    timesD1=varargin{1};
    Response1=varargin{2};
    cmap=varargin{3};
    si=varargin{4};
    col=varargin{5};
    takeHalf=varargin{6};
    if strcmp(takeHalf,'takeTopHalf')
        cmap=jet(2*length(col(~isnan(col))));
        cmap=cmap(length(col(~isnan(col)))+1:end,:);
    elseif strcmp(takeHalf,'takeBottomHalf')
        cmap=jet(2*length(col(~isnan(col))));
        cmap=cmap(1:length(col(~isnan(col))),:);
    end
else
    error('incorrect number of args to plotResponseLineByLine in plotVariousSUResponsesAlignedToBeh.m');
end

skipGray=true;
LineW=1;
figure();
incColorInd=1;
if ~isempty(takeHalf)
    if strcmp(takeHalf,'takeTopHalf')
        takeNotNan=~isnan(col);
        col=col(takeNotNan);
        Response1.unitbyunit_y(takeNotNan,:);
        [~,sortedRind]=sort(col,'descend');
        Response1.unitbyunit_y=Response1.unitbyunit_y(sortedRind,:);
        si=si(sortedRind);
    elseif strcmp(takeHalf,'takeBottomHalf')
        takeNotNan=~isnan(col);
        col=col(takeNotNan);
        Response1.unitbyunit_y(takeNotNan,:);
        [~,sortedRind]=sort(col,'ascend');
        Response1.unitbyunit_y=Response1.unitbyunit_y(sortedRind,:);
        si=si(sortedRind);
    end
end
offsetLines=true;
plotY=0;
if offsetLines==true
    takeNotNan=~isnan(col);
    col=col(takeNotNan);
    Response1.unitbyunit_y(takeNotNan,:);
    [~,sortedRind]=sort(col,'ascend');
    Response1.unitbyunit_y=Response1.unitbyunit_y(sortedRind,:);
    si=si(sortedRind);
end
for i=1:size(Response1.unitbyunit_y,1)
    if isnan(col(i))
        if skipGray==false
            lh=plot(plotY+timesD1,Response1.unitbyunit_y(i,:),'Color',[0.9 0.9 0.9]); hold on;
            lh.Color=[0.9 0.9 0.9 0.2];
            if offsetLines==true
                plotY=plotY+max(Response1.unitbyunit_y(i,1:200),[],2,'omitnan')+1;
            end
        end
    else
        alph=1;
        if strcmp(takeHalf,'takeTopHalf') 
            if si(incColorInd)<size(cmap,1)/2
                alph=1;
                LineW=2;
            else
                alph=0.5;
            end
        end
        if strcmp(takeHalf,'takeBottomHalf')
            if si(incColorInd)>size(cmap,1)/2
                alph=0.9;
                LineW=2;
            else
                alph=0.5;
            end
        end
        lh=plot(timesD1,plotY+Response1.unitbyunit_y(i,:),'Color',cmap(si(incColorInd),:),'LineWidth',LineW); hold on;
        lh.Color=[cmap(si(incColorInd),:),alph];
        incColorInd=incColorInd+1;
        if offsetLines==true
            plotY=plotY+max(Response1.unitbyunit_y(i,1:200),[],2,'omitnan')+1;
        end
    end
end

end

function modIndex=getModIndex(activityD1tagged,timesD1,inWindow1,inWindow2)

inWindow1=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=inWindow1(1) & timesD1<=inWindow1(2)),2);
inWindow2=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>inWindow2(1) & timesD1<=inWindow2(2)),2);
inWindow1(inWindow1<0.0001)=0;
inWindow2(inWindow2<0.0001)=0;

modIndex=(inWindow2-inWindow1)./(inWindow2+inWindow1);

end

function plotRainbowMaps(Response1,cmap,si,col,colRelativeTo,responseN,binStep,modIndex1,modIndex2)

timesD1=getTimesD1(Response1);
plotResponseLineByLine(timesD1,Response1,cmap,si,col); title(['Response ' num2str(responseN) ' -- all lines']);
plotModScatter(modIndex1,modIndex2,cmap,si,col);

% [R,newcol,sinew,newcmap,~,m1,m2]=binByCol(Response1,col,binStep,modIndex1,modIndex2); 
% plotResponseLineByLine(timesD1,R,newcmap,sinew,newcol);  title(['All colormap ' num2str(responseN) ' -- only mod ratio BINNED']);
% plotModScatter(m1,m2,newcmap,sinew,newcol);

tempcol=col; tempcol(col<colRelativeTo)=nan; [sinew,newcmap]=getNewSi(tempcol);
plotResponseLineByLine(timesD1,Response1,newcmap,sinew,tempcol,'takeTopHalf'); title(['Response ' num2str(responseN) ' -- only mod ratio >']);
plotModScatter(modIndex1(~isnan(tempcol)),modIndex2(~isnan(tempcol)),newcmap,sinew,tempcol(~isnan(tempcol)),'takeTopHalf');

[R,newcol,sinew,newcmap,~,m1,m2]=binByCol(Response1,tempcol,binStep,modIndex1,modIndex2); 
plotResponseLineByLine(timesD1,R,newcmap,sinew,newcol,'takeTopHalf');  title(['Response ' num2str(responseN) ' -- only mod ratio > BINNED']);
plotModScatter(m1,m2,newcmap,sinew,newcol,'takeTopHalf');

tempcol=col; tempcol(col>colRelativeTo)=nan; [sinew,newcmap]=getNewSi(tempcol);
plotResponseLineByLine(timesD1,Response1,newcmap,sinew,tempcol,'takeBottomHalf'); title(['Response ' num2str(responseN) ' -- only mod ratio <']);
plotModScatter(modIndex1(~isnan(tempcol)),modIndex2(~isnan(tempcol)),newcmap,sinew,tempcol(~isnan(tempcol)),'takeBottomHalf');

[R,newcol,sinew,newcmap,~,m1,m2]=binByCol(Response1,tempcol,binStep,modIndex1,modIndex2); 
plotResponseLineByLine(timesD1,R,newcmap,sinew,newcol,'takeBottomHalf');  title(['Response ' num2str(responseN) ' -- only mod ratio < BINNED']);
plotModScatter(m1,m2,newcmap,sinew,newcol,'takeBottomHalf');

end

function [si,cmap]=getNewSi(col)

maxOutCmap=false;
if maxOutCmap==true
    fullrange=[-2 2]+2;
    currcolrange=[min(col,[],'all','omitnan') max(col,[],'all','omitnan')]+2;
    fracThroughCmapRange=[currcolrange(1)./fullrange(2) currcolrange(2)./fullrange(2)];
    padl=[fracThroughCmapRange(1)*length(col(~isnan(col)))/range(currcolrange) fracThroughCmapRange(2)*length(col(~isnan(col)))/range(currcolrange)];
    cmap=jet(ceil(padl(1)+length(col(~isnan(col)))+padl(2)));
    cmap=cmap(ceil(fracThroughCmapRange(1)*length(col(~isnan(col)))/range(currcolrange))+1:ceil(padl(1)+length(col(~isnan(col)))),:);
else
    cmap=jet(length(col(~isnan(col))));
end

[~,sip]=sort(col(~isnan(col)));
si=nan(size(sip));
for i=1:length(sip)
    si(i)=find(sip==i);
end

end

function modRatioPSTH(Response1,Response2,Response3,col,colRelativeTo,win1,win2,modIndex1,modIndex2,modIndex3)

modIndexAfterSmooth=false;
% binStep=15;
% binStep=7;
binStep=30;

[si,cmap]=getNewSi(col);

timesD1=getTimesD1(Response3);
meanForZscore=mean([Response1.unitbyunit_y(:,1:200) Response2.unitbyunit_y(:,1:200) Response3.unitbyunit_y(:,ismember(1:200,find(timesD1>-1)))],2,'omitnan');
sdForZscore=std([Response1.unitbyunit_y(:,1:200) Response2.unitbyunit_y(:,1:200) Response3.unitbyunit_y(:,ismember(1:200,find(timesD1>-1)))],0,2,'omitnan');

timesD1=getTimesD1(Response1);
% Response1.unitbyunit_y(any(Response1.unitbyunit_y>15,2),:)=nan;
Response1=ZscoreAndSmooth(Response1,meanForZscore,sdForZscore);
timesD1=getTimesD1(Response2);
% Response2.unitbyunit_y(any(Response2.unitbyunit_y>15,2),:)=nan;
Response2=ZscoreAndSmooth(Response2,meanForZscore,sdForZscore);
timesD1=getTimesD1(Response3);
% for Response 3, zero out optotag window
% and realign after cue
Response3=ZscoreAndSmooth(Response3,meanForZscore,sdForZscore);
Response3.unitbyunit_y(:,timesD1<-1)=nan;
% Response3.unitbyunit_y(any(Response3.unitbyunit_y>15,2),:)=nan;
% Response3.unitbyunit_y=Response3.unitbyunit_y-repmat(mean(Response3.unitbyunit_y(:,timesD1>-1 & timesD1<-0.5),2,'omitnan'),1,size(Response3.unitbyunit_y,2));

if modIndexAfterSmooth==true
    timesD1=getTimesD1(Response1);
    mod1=getModIndex(Response1,timesD1,win1,win2);
    timesD1=getTimesD1(Response2);
    mod2=getModIndex(Response2,timesD1,win1,win2);
    timesD1=getTimesD1(Response3);
    mod3=getModIndex(Response3,timesD1,win1,win2);
    if colRelativeTo==1
        % divide
         modRat=mod2./mod1;
    elseif colRelativeTo==0
        % subtract
        modRat=mod2-mod1;
    end
    col=modRat;
end

% Plot scatters vs Response3
plotJustRainbowScatter(modIndex1,modIndex3,cmap,si,col); title('Mod index 3 on y axis, Mod index 1 on x axis');
plotJustRainbowScatter(modIndex2,modIndex3,cmap,si,col); title('Mod index 3 on y axis, Mod index 2 on x axis');

% Response 1
plotRainbowMaps(Response1,cmap,si,col,colRelativeTo,1,binStep,modIndex1,modIndex2);
% Response 2
plotRainbowMaps(Response2,cmap,si,col,colRelativeTo,2,binStep,modIndex1,modIndex2);
% Response 3
plotRainbowMaps(Response3,cmap,si,col,colRelativeTo,3,binStep,modIndex1,modIndex2);

end

function Response1=ZscoreDontSmooth(varargin)

if length(varargin)==1
    Response1=varargin{1};
    meforZ=[];
    sdforZ=[];
elseif length(varargin)==3
    Response1=varargin{1};
    meforZ=varargin{2};
    sdforZ=varargin{3};
end

% Zscore within each unit, then plot transition to time window 2
temp=Response1.unitbyunit_y;
if isempty(meforZ)
    temp=temp-repmat(mean(temp,2,'omitnan'),1,size(temp,2));
    temp=temp./repmat(std(temp(:,1:200),0,2,'omitnan'),1,size(temp,2));
else
    temp=temp-repmat(meforZ,1,size(temp,2));
    temp=temp./repmat(sdforZ,1,size(temp,2));
end
Response1.unitbyunit_y=temp;

end

function inWindow1=getRateInWindow(activityD1tagged,timesD1,inWin1)

inWindow1=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=inWin1(1) & timesD1<=inWin1(2)),2);

end

function plotJustRainbowScatter(modIndex1,modIndex2,cmap,si,col)

plotModScatter(modIndex1,modIndex2,cmap,si,col);

end

function coloredPSTH(Response1,Response2,col)

temp=nanmean(Response1.aligncomp_x,1);
[~,f]=nanmax(nanmean(Response1.aligncomp_y,1));
timesD1=nanmean(Response1.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end
% Zscore within each unit, then plot transition to time window 2
temp=Response1.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp(:,1:200),0,2,'omitnan');
Zscored_D1tagged=temp;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
Response1.unitbyunit_y=temp;
cmap=jet(length(col(~isnan(col))));
[~,si]=sort(col(~isnan(col)));
figure();
incColorInd=1;
offsetLines=true;
lineOff=0;
for i=1:size(Response1.unitbyunit_y,1)
    if isnan(col(i))
        lh=plot(timesD1,Response1.unitbyunit_y(i,:)+lineOff,'Color',[0.9 0.9 0.9]); hold on;
        lh.Color=[0.9 0.9 0.9 0.2];
    else
        lh=plot(timesD1,Response1.unitbyunit_y(i,:)+lineOff,'Color',cmap(si(incColorInd),:)); hold on;
        lh.Color=[cmap(si(incColorInd),:),0.5];
        incColorInd=incColorInd+1;
    end
    if offsetLines==true
        if all(isnan(Response1.unitbyunit_y(i,1:200)))
        else
            lineOff=lineOff+nanmax(Response1.unitbyunit_y(i,1:200));
        end
    end
end
title('Response 1');

temp=nanmean(Response2.aligncomp_x,1);
[~,f]=nanmax(nanmean(Response2.aligncomp_y,1));
timesD1=nanmean(Response2.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end
% Zscore within each unit, then plot transition to time window 2
temp=Response2.unitbyunit_y;
temp=temp-mean(temp,2,'omitnan');
temp=temp./std(temp(:,1:200),0,2,'omitnan');
Zscored_D1tagged=temp;
temp=Zscored_D1tagged;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
end
Response2.unitbyunit_y=temp;
cmap=jet(length(col(~isnan(col))));
[~,si]=sort(col(~isnan(col)));
figure();
incColorInd=1;
lineOff=0;
for i=1:size(Response2.unitbyunit_y,1)
    if isnan(col(i))
        lh=plot(timesD1,Response2.unitbyunit_y(i,:)+lineOff,'Color',[0.9 0.9 0.9]); hold on;
        lh.Color=[0.9 0.9 0.9 0.2];
    else
        lh=plot(timesD1,Response2.unitbyunit_y(i,:)+lineOff,'Color',cmap(si(incColorInd),:)); hold on;
        lh.Color=[cmap(si(incColorInd),:),0.5];
        incColorInd=incColorInd+1;
    end
    if offsetLines==true
        if all(isnan(Response2.unitbyunit_y(i,1:200)))
        else
            lineOff=lineOff+nanmax(Response2.unitbyunit_y(i,1:200));
        end
    end
end
title('Response 2');

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

cmap='jet'; %'hsv';

if length(c)~=length(r1)
    disp('lengths of Response and color vec do not match in scatterTwoResponses in plotVariousSUResponsesAlignedToBeh.m');
    disp('ignoring the color vec');
    c=[];
end

% if nans in c, gray out these
grayout=isnan(c);
if any(grayout)
    figure();
    scatter(r1(grayout==true),r2(grayout==true),[],[0.9 0.9 0.9],'LineWidth',1); hold on;
    scatter(r1(grayout==false),r2(grayout==false),[],c(grayout==false),'LineWidth',1);
else
    figure();
    if isempty(c)
        scatter(r1,r2);
    else
        scatter(r1,r2,[],c,'LineWidth',1);
    end
end
colormap(cmap);
xlabel('Response 1');
ylabel('Response 2');

end

function [Response,Response2,Response3,Response4,Response5]=matchUnitsAcrossResponses(varargin)

Response=[]; Response2=[]; Response3=[]; Response4=[]; Response5=[];
if length(varargin)==2
    Response=varargin{1};
    Response2=varargin{2};
    Response3=[];
    Response1_indsIntoExcluded=find(Response.excluded==0);
    Response2_indsIntoExcluded=find(Response2.excluded==0);
    useTheseUnitsResponse1=find(ismember(Response1_indsIntoExcluded,Response2_indsIntoExcluded));
    Response.excluded(Response1_indsIntoExcluded(~ismember(Response1_indsIntoExcluded,Response2_indsIntoExcluded)))=1;
    useTheseUnitsResponse2=find(ismember(Response2_indsIntoExcluded,Response1_indsIntoExcluded));
    Response2.excluded(Response2_indsIntoExcluded(~ismember(Response2_indsIntoExcluded,Response1_indsIntoExcluded)))=1;
    Response=filterResponseToOneSU(Response,useTheseUnitsResponse1);
    Response2=filterResponseToOneSU(Response2,useTheseUnitsResponse2);
elseif length(varargin)==3
    Response=varargin{1};
    Response2=varargin{2};
    Response3=varargin{3};
    Response1_indsIntoExcluded=find(Response.excluded==0);
    Response2_indsIntoExcluded=find(Response2.excluded==0);
    Response3_indsIntoExcluded=find(Response3.excluded==0);

    u=unique([Response1_indsIntoExcluded; Response2_indsIntoExcluded; Response3_indsIntoExcluded]);
    isInAll=zeros(size(u));
    for i=1:length(u)
        if ismember(u(i),Response1_indsIntoExcluded) && ismember(u(i),Response2_indsIntoExcluded) && ismember(u(i),Response3_indsIntoExcluded)
            isInAll(i)=1;
        else
            isInAll(i)=0;
        end
    end
    inAll=u(isInAll==1);

    temp=ismember(Response1_indsIntoExcluded,inAll);
    useTheseUnitsResponse1=find(temp);
    Response.excluded(Response1_indsIntoExcluded(~temp))=1;

    temp=ismember(Response2_indsIntoExcluded,inAll);
    useTheseUnitsResponse2=find(temp);
    Response2.excluded(Response2_indsIntoExcluded(~temp))=1;

    temp=ismember(Response3_indsIntoExcluded,inAll);
    useTheseUnitsResponse3=find(temp);
    Response3.excluded(Response3_indsIntoExcluded(~temp))=1;

    Response=filterResponseToOneSU(Response,useTheseUnitsResponse1);
    Response2=filterResponseToOneSU(Response2,useTheseUnitsResponse2);
    Response3=filterResponseToOneSU(Response3,useTheseUnitsResponse3);
elseif length(varargin)==5
    Response=varargin{1};
    Response2=varargin{2};
    Response3=varargin{3};
    Response4=varargin{4};
    Response5=varargin{5};
    Response1_indsIntoExcluded=find(Response.excluded==0);
    Response2_indsIntoExcluded=find(Response2.excluded==0);
    if ~isempty(Response3)
        Response3_indsIntoExcluded=find(Response3.excluded==0);
    else
        Response3_indsIntoExcluded=[];
    end
    if ~isempty(Response4)
        Response4_indsIntoExcluded=find(Response4.excluded==0);
    else
        Response4_indsIntoExcluded=[];
    end
    if ~isempty(Response5)
        Response5_indsIntoExcluded=find(Response5.excluded==0);
    else
        Response5_indsIntoExcluded=[];
    end

    u=unique([Response1_indsIntoExcluded; Response2_indsIntoExcluded; Response3_indsIntoExcluded; Response4_indsIntoExcluded; Response5_indsIntoExcluded]);
    isInAll=zeros(size(u));
    for i=1:length(u)
        if ismember(u(i),Response1_indsIntoExcluded) && ismember(u(i),Response2_indsIntoExcluded) && ismember(u(i),Response3_indsIntoExcluded) && ismember(u(i),Response4_indsIntoExcluded) && ismember(u(i),Response5_indsIntoExcluded)
            isInAll(i)=1;
        else
            isInAll(i)=0;
        end
    end
    inAll=u(isInAll==1);

    temp=ismember(Response1_indsIntoExcluded,inAll);
    useTheseUnitsResponse1=find(temp);
    Response.excluded(Response1_indsIntoExcluded(~temp))=1;

    temp=ismember(Response2_indsIntoExcluded,inAll);
    useTheseUnitsResponse2=find(temp);
    Response2.excluded(Response2_indsIntoExcluded(~temp))=1;

    if ~isempty(Response3)
        temp=ismember(Response3_indsIntoExcluded,inAll);
        useTheseUnitsResponse3=find(temp);
        Response3.excluded(Response3_indsIntoExcluded(~temp))=1;
    else
        useTheseUnitsResponse3=[];
    end

    if ~isempty(Response4)
        temp=ismember(Response4_indsIntoExcluded,inAll);
        useTheseUnitsResponse4=find(temp);
        Response4.excluded(Response4_indsIntoExcluded(~temp))=1;
    else
        useTheseUnitsResponse4=[];
    end

    if ~isempty(Response5)
        temp=ismember(Response5_indsIntoExcluded,inAll);
        useTheseUnitsResponse5=find(temp);
        Response5.excluded(Response5_indsIntoExcluded(~temp))=1;
    else
        useTheseUnitsResponse5=[];
    end

    Response=filterResponseToOneSU(Response,useTheseUnitsResponse1);
    Response2=filterResponseToOneSU(Response2,useTheseUnitsResponse2);
    Response3=filterResponseToOneSU(Response3,useTheseUnitsResponse3);
    Response4=filterResponseToOneSU(Response4,useTheseUnitsResponse4);
    Response5=filterResponseToOneSU(Response5,useTheseUnitsResponse5);

    fieldLikeResponseSize=size(Response.unitbyunit_y,1);
    if length(Response.fromWhichUnit)==fieldLikeResponseSize
        % took all trials
        Response=checkThatAllUnitsRepresented(Response,useTheseUnitsResponse1);
        Response2=checkThatAllUnitsRepresented(Response2,useTheseUnitsResponse2);
        Response3=checkThatAllUnitsRepresented(Response3,useTheseUnitsResponse3);
        Response4=checkThatAllUnitsRepresented(Response4,useTheseUnitsResponse4);
        Response5=checkThatAllUnitsRepresented(Response5,useTheseUnitsResponse5);
    end
end

end

function out=checkThatAllUnitsRepresented(Response,whichUnitToUse)

if isempty(Response)
    out=[];
    return
end

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

out=Response;
if any(~ismember(whichUnitToUse,whichUnit))
    % Response is missing this unit
    % add a row for this unit(s)
    addBack=whichUnitToUse(~ismember(whichUnitToUse,whichUnit));
    for i=1:length(f)
        temp=Response.(f{i});
        if size(temp,1)==fieldLikeResponseSize
            % add back to this field
            if all(size(temp)>1)
                % 2D
                out.(f{i})=[temp; nan(length(addBack),size(temp,2))];
            else
                % 1D
                if strcmp(f{i},'fromWhichUnit')
                    out.(f{i})=[temp; addBack];
                else
                    out.(f{i})=[temp; nan(length(addBack),1)];
                end
            end
        end
    end

    % sort according to unit order
    [~,si]=sort(out.fromWhichUnit);
    for i=1:length(f)
%         temp=Response.(f{i});
        temp=out.(f{i});
        if size(temp,1)==fieldLikeResponseSize
            if all(size(temp)>1)
                % 2D
                out.(f{i})=temp(si,:);
            else
                % 1D
                out.(f{i})=temp(si);
            end
        end
    end
end

end

function out=filterResponseToOneSU(Response,whichUnitToUse)

if isempty(Response)
    out=[];
    return
end

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
        if all(size(temp)>1)
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

rmOutliers=false;
if rmOutliers==true
    [activityD1tagged.unitbyunit_y,Tfrm]=rmoutliers(activityD1tagged.unitbyunit_y,"mean","ThresholdFactor",15);
    disp(['removed ' num2str(sum(Tfrm)) ' outliers']);
end
temp=activityD1tagged.unitbyunit_y;
for i=1:size(temp,1)
%     temp(i,:)=smoothdata(temp(i,:),'gaussian',2);
%     temp(i,:)=smoothdata(temp(i,:),'gaussian',10);
    temp(i,:)=smoothdata(temp(i,:),'gaussian',7);
end
activityD1tagged.unitbyunit_y=temp;

if size(activityD1tagged.unitbyunit_y,1)>1000
    % skip unit by unit plot because too crowded
else
    figure();
    plot(nanmean(activityD1tagged.unitbyunit_x,1),activityD1tagged.unitbyunit_y');
    hold on; plot(nanmean(activityD1tagged.aligncomp_x,1),nanmean(activityD1tagged.aligncomp_y,1).*max(activityD1tagged.unitbyunit_y,[],'all','omitnan'),'Color','b');
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

figure();
unitoffset=0;
for i=1:size(activityD1tagged.unitbyunit_y,1)
    plot(timesD1,activityD1tagged.unitbyunit_y(i,:)+unitoffset,'Color','k'); hold on;
    uptooo=200; 
    if uptooo>size(activityD1tagged.unitbyunit_y,2)
        uptooo=size(activityD1tagged.unitbyunit_y,2);
    end
    if isnan(max(activityD1tagged.unitbyunit_y(i,1:uptooo),[],'all','omitnan'))
        continue
    end
    unitoffset=unitoffset+max(activityD1tagged.unitbyunit_y(i,1:uptooo),[],'all','omitnan');
end
title('Align companion at time=0');
xlabel('Time (sec)');
ylabel('Firing rate');

[out.me,out.plusSe,out.minusSe,out.t]=plotMeanAndSE(downSampMatrix(activityD1tagged.unitbyunit_y,downSampFac),'k',downSampAv(timesD1,downSampFac));
title('Mean and se across units');
xlabel('Time (sec)');
ylabel('Firing rate');

end

function [me,plusSe,minusSe,t]=plotMeanAndSE(data1,color1,timeBins)

plotAsCityscape=false;

data1=data1(any(~isnan(data1),2),:);

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
temp=temp./std(temp(:,1:200),0,2,'omitnan');
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
    if max(currtoplot(1:200),[],'all','omitnan')>4
        pause;
    end
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
temp=temp./std(temp(:,1:200),0,2,'omitnan'); % Z-score each unit's own activity
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
fillVal=mean(Rs,'all','omitnan');
for i=1:size(Rs,1)
    for j=i-0:i+0
        if j>=1 && j<=size(Rs,1)
            Rs(i,j)=fillVal;
        end
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
%         Rs(isnan(Rs))=0;
%         allRs=allRs+Rs;
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
figure();
imagesc(timeBins(timeBins<4),timeBins(timeBins<4),allRs(timeBins<4,timeBins<4));
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

% get rid of trials where event did not occur (all nan)
activityD1tagged.unitbyunit_y=activityD1tagged.unitbyunit_y(any(~isnan(activityD1tagged.unitbyunit_y),2),:);

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
temp=temp./std(temp(:,1:200),0,'all','omitnan'); % Z-score using std across all trials, all timepoints
Zscored_D1tagged=temp;

Zscored_D1tagged=Zscored_D1tagged(~all(isnan(Zscored_D1tagged),2),:);

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

indicateWhichAreSig=true;

% inWindow in seconds relative to alignmentCompanion peak

temp=nanmean(activityD1tagged.aligncomp_x,1);
[~,f]=nanmax(nanmean(activityD1tagged.aligncomp_y,1));
timesD1=nanmean(activityD1tagged.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

inWindow1=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>=inWindow1(1) & timesD1<=inWindow1(2)),2);
inWindow2=nanmean(activityD1tagged.unitbyunit_y(:,timesD1>inWindow2(1) & timesD1<=inWindow2(2)),2);
inWindow1(inWindow1<0.0001)=0;
inWindow2(inWindow2<0.0001)=0;

modIndex=(inWindow2-inWindow1)./(inWindow2+inWindow1);
if any(modIndex>1,'all') || any(modIndex<-1,'all')
    pause;
end

modBins=[-1-0.1:0.1:1]+0.05;
[N,edges]=histcounts(modIndex,modBins);
% binCenters=edges(1:end-1)+0.05;
[n,x]=cityscape_hist(N,edges);
figure();
plot(x,n,'Color',[0.8 0.8 0.8]); hold on;
% bar(binCenters,N,'LineWidth',0.1);
xlabel('Modulation index');
ylabel('Count');

if indicateWhichAreSig==true
    if isfield(activityD1tagged,'isSig')
        modBins=[-1-0.1:0.1:1]+0.05;
        [N,edges]=histcounts(modIndex(activityD1tagged.isSig==1),modBins);
        [n,x]=cityscape_hist(N,edges);
        plot(x,n,'Color','k');
    end
end

end

function out=isResponsive(Response,takewin1,takewin2)

% binomial or ranksum stats
nTimesStdevAtBaseline=1.2; % for sustained test
ds_for_sustained=10; % downsamp bins
% sigThreshold=0.01; 
sigThreshold=0.2;

temp=nanmean(Response.aligncomp_x,1);
[~,f]=nanmax(nanmean(Response.aligncomp_y,1));
timesD1=nanmean(Response.unitbyunit_x,1)-temp(f);
if any(isnan(timesD1))
    timesD1=fillmissing(timesD1,'linear');
end

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

u=unique(whichUnit);

indsForWin1=find(timesD1>=takewin1(1) & timesD1<=takewin1(2));
indsForWin2=find(timesD1>=takewin2(1) & timesD1<=takewin2(2));

binsize=250; % in ms
binsize=binsize/1000; % in sec

nindsperbin=ceil(binsize/(timesD1(2)-timesD1(1)));
binsForWin2=indsForWin2(1):nindsperbin:indsForWin2(end);
isSig=nan(size(u));
pvals=nan(size(u));
positiveMod=nan(size(u));
sustained=nan(size(u));
inWindow1Mean=nan(size(u));
inWindow2Mean=nan(size(u));
allWindowMean=nan(size(u));
allWindowVar=nan(size(u));
for i=1:length(u)
    currU=u(i);
    temp=Response.unitbyunit_y;
    unittemp=temp(whichUnit==currU,:);
    if length(Response.fromWhichUnit)==fieldLikeResponseSize
        % took all trials
        % look for stat across trials
        unittemp=unittemp(~all(isnan(unittemp),2),:);
        inWindow1=nanmean(unittemp(:,indsForWin1),2);
        inWindow2=nanmean(unittemp(:,indsForWin2),2);
        if all(isnan(inWindow1)) || all(isnan(inWindow2))
            isSig(i)=false;
            continue
        end
        rank_p=ranksum(inWindow1,inWindow2);
        bino_p=compareCounts(inWindow1>0,inWindow2>0);
        p=min([rank_p bino_p]);
        pvals(i)=p;
        isSig(i)=p<sigThreshold;
        positiveMod(i)=mean(inWindow2,'all','omitnan')>mean(inWindow1,'all','omitnan');
        inWindow1Mean(i)=mean(inWindow1,'all','omitnan');
        inWindow2Mean(i)=mean(inWindow2,'all','omitnan');
        allWindowMean(i)=mean(unittemp,'all','omitnan');
        allWindowVar(i)=var(unittemp,0,'all','omitnan');
        % also check whether one timepoint in win 2 is modulated 
%         pval_timepoints_win2=nan(1,length(binsForWin2)-1);
%         countsDuringWin1=countEventsInBin(unittemp,indsForWin1);
%         for j=1:length(binsForWin2)-1
%             currbin=binsForWin2(j):binsForWin2(j+1);
%             countsDuringWin2=countEventsInBin(unittemp,currbin);
%             pval_timepoints_win2(j)=compareCounts(countsDuringWin1,countsDuringWin2);
%         end
%         timepoints_pval_thresh=0;
%         if any(pval_timepoints_win2<timepoints_pval_thresh)
%             pvals(i)=min(min(pval_timepoints_win2,[],'all','omitnan'),pvals(i));
%             isSig(i)=true;
%         end
    elseif length(Response.fromWhichSess)==fieldLikeResponseSize
        % took unit by unit
        % compare stat in time window 1 to 2
        if size(unittemp,1)>1
            error('unittemp size unexpected in plotVariousSUResponsesAlignedToBeh.m');
        end
        inWindow1=unittemp(indsForWin1);
        inWindow2=unittemp(indsForWin2);
        if all(isnan(inWindow1)) || all(isnan(inWindow2))
            isSig(i)=false;
            continue
        end
        rank_p=ranksum(inWindow1,inWindow2);
        bino_p=compareCounts(inWindow1,inWindow2);
        p=min([rank_p bino_p]);
        pvals(i)=p;
        isSig(i)=p<sigThreshold;
        positiveMod(i)=mean(inWindow2,'all','omitnan')>mean(inWindow1,'all','omitnan');
        inWindow1Mean(i)=mean(inWindow1,'all','omitnan');
        inWindow2Mean(i)=mean(inWindow2,'all','omitnan');
        allWindowMean(i)=mean(unittemp,'all','omitnan');
        allWindowVar(i)=var(unittemp,0,'all','omitnan');
    end
    te=-2*ones(1,size(unittemp,2));
    te(indsForWin1)=0; % specific win1 inds
    te(indsForWin2)=1; % specific win2 inds
    [increase,decrease]=isSustained(nanmean(unittemp,1),te,nTimesStdevAtBaseline,ds_for_sustained);
    if increase || decrease
        sustained(i)=true;
    else
        sustained(i)=false;
    end
end

out.isSig=isSig;
out.pvals=pvals;
out.positiveMod=positiveMod;
out.sustained=sustained;
out.inWindow1Mean=inWindow1Mean;
out.inWindow2Mean=inWindow2Mean;
out.allWindowMean=allWindowMean;
out.allWindowVar=allWindowVar;

% out.isSig=(isSig + (positiveMod==0))>1.5;

end

function pval=compareCounts(base,test)

% binomial stats
% what is the probability of observing the number of spikes during opto 
% if the real probability of a spike were unchanged from baseline
% one-sided test that opto spiking is greater than baseline

p=sum(base>0)/length(base);
k=sum(test>0);
n=length(test);
pval=0;
for i=k:n
    pval=pval+binopdf(i,n,p);
end

end

function out=countEventsInBin(dataMatrix,binInInds)

dataMatrix(dataMatrix>0)=1;
out=sum(dataMatrix(:,binInInds),2,'omitnan');

end

function [increase,decrease]=isSustained(avAlignedToOpto,optoOnInUnitTimes,nTimesStdevAtBaseline,ds)

avAlignedToOpto=downSampAv(avAlignedToOpto,ds);
optoOnInUnitTimes=downSampAv(optoOnInUnitTimes,ds);

meanAtBaseline=mean(avAlignedToOpto(optoOnInUnitTimes==0),'all','omitnan');
stdAtBaseline=std(avAlignedToOpto(optoOnInUnitTimes==0),0,'all','omitnan');
meanDuringOpto=mean(avAlignedToOpto(optoOnInUnitTimes==1),'all','omitnan');
increase=meanDuringOpto>meanAtBaseline+(stdAtBaseline*nTimesStdevAtBaseline);
decrease=meanDuringOpto<meanAtBaseline-(stdAtBaseline*nTimesStdevAtBaseline);

end