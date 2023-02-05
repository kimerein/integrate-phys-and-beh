function [cueCD,uncueCD,lateUncueCD,unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,excluded,ns,D1orD2taggingExpt,firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess,fromWhichSess,fromWhichUnit,fromWhichSess_forTrials,fromWhichTrial,isEventInThisTrial,D1taggedCells,A2ataggedCells]=alignToCompanion(datadir,getCDs,cueCD,uncueCD,lateUncueCD,onlyTakeTheseUnits,settings,indsForAnalysisPerSess)

firstValNtimesBaseVar=[]; optoFRoverBaseline=[];

% excludeHigherFR=settings.excludeHigherFR;
cutAtTime=settings.cutAtTime;
ds=settings.ds;
normalizePSTHs=settings.normalizePSTHs;
suppressPlots=settings.suppressPlots;
% excludeAboveFR=settings.excludeAboveFR; % in spikes / sec, exclude units with average firing rate above this
padsize=settings.padsize;
testForAlignment=settings.testForAlignment;
unitbaseline=settings.unitbaseline;
% beforeOptoBaseline=settings.beforeOptoBaseline;
keepAllSingleTrials=settings.keepAllSingleTrials;

% excludeHigherFR=false;
% cutAtTime=3; % stop plotting this many seconds after max of alignment companion
% ds=6; % downsample bin size
% onlyTakeTheseUnits='D1tagged'; % if is not empty, will only take units with this string in filename, e.g., 'D1tagged'
% or make cutAtTime empty to plot all time points
% getCDs=true;
% normalizePSTHs=false;

maxUnitsPerSess=settings.maxUnitsPerSess;
if keepAllSingleTrials==true
    maxUnitsPerSess=maxUnitsPerSess*30;
end
if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
unitbyunit_x=[];
unitbyunit_y=[];
ns=nan(maxUnitsPerSess*length(dd),1);
aligncomp_x=[];
aligncomp_y=[];
excluded=zeros(maxUnitsPerSess*length(dd),1);
D1orD2taggingExpt=nan(length(dd),1); % will be 1 for D1, 2 for A2a tagging session
fromWhichSess=nan(maxUnitsPerSess*length(dd),1);
fromWhichSess_forTrials=nan(maxUnitsPerSess*length(dd),1);
fromWhichUnit=nan(maxUnitsPerSess*length(dd),1);
fromWhichTrial=nan(maxUnitsPerSess*length(dd),1);
isEventInThisTrial=nan(maxUnitsPerSess*length(dd),1);
D1taggedCells=zeros(maxUnitsPerSess*length(dd),1);
A2ataggedCells=zeros(maxUnitsPerSess*length(dd),1);
excluded_count=1;
units_count=1;
trials_count=1;
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    disp(['reading in from ' datadir]);
    ls=dir(datadir);
    for i=3:length(ls)
        a=[];
        if ismac()
            a=load([ls(i).folder '/' ls(i).name]);
        else
            a=load([ls(i).folder '\' ls(i).name]);
        end
        
        if ~isempty(regexp(ls(i).name, 'D1tagged'))
            D1orD2taggingExpt(j)=1;
            D1taggedCells(excluded_count)=1;
        elseif ~isempty(regexp(ls(i).name, 'A2atagged'))
            D1orD2taggingExpt(j)=2;
            A2ataggedCells(excluded_count)=1;
        end

        if ~isempty(onlyTakeTheseUnits)
            if isempty(regexp(ls(i).name, onlyTakeTheseUnits))
                % doesn't contain string, continue and skip this unit
                excluded(excluded_count)=1;
                excluded_count=excluded_count+1;
                continue
            end
        end

        [unitTest,dontUseTrials]=doUnitTest(ls(i).folder, ls(i).name);
        if unitTest
        else
            excluded(excluded_count)=1;
            excluded_count=excluded_count+1;
            continue % failed unit test
        end

        if isempty(a.dataout)
            disp([ls(i).folder '\' ls(i).name ' is empty ... skipping']);
            excluded(excluded_count)=1;
            excluded_count=excluded_count+1;
            continue
        end

        if length(dontUseTrials)~=size(a.dataout.y,1)
            disp('dontUseTrials length does not match size of a.dataout.y in line 91 of alignToCompanion.m, so skipping this unit');
            excluded(excluded_count)=1;
            excluded_count=excluded_count+1;
            continue
        end
        excluded_count=excluded_count+1;

        % throw out trials where unit dead or out of range
        % keep track of which trials included
        if any(dontUseTrials==1)
            a.dataout.y=a.dataout.y(dontUseTrials==0,:);
            a.alignComp.y=a.alignComp.y(dontUseTrials==0,:);
            usingTrialsInds=find(dontUseTrials==0);
        else
            usingTrialsInds=1:size(a.dataout.y,1);
        end

        % was there an event on each trial
        eventHappens=any(~isnan(a.alignComp.y),2);
        
%         if excludeHigherFR
%             if nanmean(nanmean(a.dataout.y(:,1:beforeOptoBaseline),1),2)>excludeAboveFR
%                 disp(['excluding unit ' ls(i).name ' based on high FR']);
%                 excluded=[excluded 1];
%                 continue
%             else
%                 excluded=[excluded 0];
%             end
%         else
%             if nanmean(nanmean(a.dataout.y(:,1:beforeOptoBaseline),1),2)>excludeAboveFR
%                 excluded=[excluded 1];
%             else
%                 excluded=[excluded 0];
%             end
%         end
        
        % realign first!
        timestep_for_aligncomp=mode(diff(a.alignComp.x));
        timestep_for_unit=mode(diff(a.dataout.x));
        [~,ma]=max(a.alignComp.y,[],2,'omitnan');
        ma=mode(ma(any(~isnan(a.alignComp.y),2)));
        if isnan(ma)
            ma=1; % no behavior event to align to
            if any(a.dataout.y(~isnan(a.dataout.y(1:end)))~=0)
                error('in alignToCompanion.m, max of alignment companion is nan but there are behavior alignments');
            end
        end
        [~,mi]=nanmin(abs(a.dataout.x-a.alignComp.x(ma)));
        if testForAlignment==true
            a.dataout.y(:,mi)=a.dataout.y(:,mi)+100;
        end
        % ensure that each unit data has the same length baseline before mi
        if mi<unitbaseline
            a.dataout.x=[nan(1,unitbaseline-mi) a.dataout.x];
            a.dataout.y=[nan(size(a.dataout.y,1),unitbaseline-mi) a.dataout.y];
        elseif unitbaseline<mi
            a.dataout.x=a.dataout.x(mi-unitbaseline:end);
            a.dataout.y=a.dataout.y(:,mi-unitbaseline:end);
        end
        curr_n=nansum(any(a.dataout.y>0,2));
        if ~isempty(unitbyunit_x)
            upTo2=size(aligncomp_x,2);
        else
            upTo2=length(a.alignComp.x)+(padsize-ma);
        end
        if upTo2-(padsize-ma)>length(a.alignComp.x)
            a.alignComp.x=[a.alignComp.x nan(1,upTo2-(padsize-ma)-length(a.alignComp.x))];
            a.alignComp.y=[a.alignComp.y nan(size(a.alignComp.y,1),upTo2-(padsize-ma)-size(a.alignComp.y,2))];
        end
        if isempty(aligncomp_x)
            % initialize
            aligncomp_x=nan(maxUnitsPerSess*length(dd),length([nan(1,padsize-ma) a.alignComp.x(1:upTo2-(padsize-ma))]));
            aligncomp_y=nan(maxUnitsPerSess*length(dd),length([nan(1,padsize-ma) nanmean(a.alignComp.y(:,1:upTo2-(padsize-ma)),1)]));
        end
        if keepAllSingleTrials==false
            aligncomp_x(units_count,:)=[nan(1,padsize-ma) a.alignComp.x(1:upTo2-(padsize-ma))];
            aligncomp_y(units_count,:)=[nan(1,padsize-ma) nanmean(a.alignComp.y(:,1:upTo2-(padsize-ma)),1)];
        else
            aligncomp_x(trials_count:trials_count+size(a.dataout.y,1)-1,:)=repmat([nan(1,padsize-ma) a.alignComp.x(1:upTo2-(padsize-ma))],size(a.dataout.y,1),1);
            aligncomp_y(trials_count:trials_count+size(a.dataout.y,1)-1,:)=repmat([nan(1,padsize-ma) nanmean(a.alignComp.y(:,1:upTo2-(padsize-ma)),1)],size(a.dataout.y,1),1);
        end
        if ~isempty(unitbyunit_x)
            upTo=size(unitbyunit_x,2);
        else
            upTo=length(a.dataout.x);
        end
        if upTo>size(a.dataout.x,2)
            % pad dataout
            a.dataout.x=[a.dataout.x nan(size(a.dataout.x,1),upTo-size(a.dataout.x,2))];
            a.dataout.y=[a.dataout.y nan(size(a.dataout.y,1),upTo-size(a.dataout.y,2))];
            % truncate
            % unitbyunit_x=unitbyunit_x(:,1:size(a.dataout.x,2));
            % unitbyunit_y=unitbyunit_y(:,1:size(a.dataout.x,2));
            % upTo=size(unitbyunit_x,2);
        end
        if isempty(unitbyunit_x)
            % initialize
            unitbyunit_x=nan(maxUnitsPerSess*length(dd),length([a.dataout.x(1:upTo)]));
            unitbyunit_y=nan(maxUnitsPerSess*length(dd),length([nanmean(a.dataout.y(:,1:upTo),1)]));
        end
        if isempty(a.dataout.y)
            if isempty(unitbyunit_x)
                disp('first unit must have trials for matrix memory allocation');
                error(disp([ls(i).name ' has no trials and cannot initialize']));
            end
            disp([ls(i).name ' has no trials']);
            unitbyunit_x(units_count,:)=[nan(size(unitbyunit_x(1,:)))];
            unitbyunit_y(units_count,:)=[nan(size(unitbyunit_x(1,:)))];
            fromWhichUnit(trials_count)=units_count;
            fromWhichTrial(trials_count)=nan;
            isEventInThisTrial(trials_count)=nan;
            fromWhichSess_forTrials(trials_count)=j*ones(size(curr_n));
            trials_count=trials_count+1;
        else
            if keepAllSingleTrials==false
                unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
            else
                unitbyunit_x(trials_count:trials_count+size(a.dataout.y,1)-1,:)=repmat(a.dataout.x(1:upTo),size(a.dataout.y,1),1);
                unitbyunit_y(trials_count:trials_count+size(a.dataout.y,1)-1,:)=[a.dataout.y(:,1:upTo)];
            end
            fromWhichUnit(trials_count:trials_count+size(a.dataout.y,1)-1)=units_count;
            fromWhichTrial(trials_count:trials_count+size(a.dataout.y,1)-1)=usingTrialsInds;
            isEventInThisTrial(trials_count:trials_count+size(a.dataout.y,1)-1)=eventHappens;
            fromWhichSess_forTrials(trials_count:trials_count+size(a.dataout.y,1)-1)=j*ones(size(curr_n));
            trials_count=trials_count+size(a.dataout.y,1);
        end
        ns(units_count)=curr_n;
        fromWhichSess(units_count)=j*ones(size(curr_n));

        disp(['Added ' ls(i).name]);
        units_count=units_count+1;
    end
end
excluded_count=excluded_count-1;
units_count=units_count-1;
if isempty(a.dataout)
    lastsetoftrialsn=1;
else
    lastsetoftrialsn=size(a.dataout.y,1);
end
trials_count=trials_count-lastsetoftrialsn;
ns=ns(1:units_count);
fromWhichSess=fromWhichSess(1:units_count);
excluded=excluded(1:excluded_count);
D1taggedCells=D1taggedCells(1:excluded_count);
A2ataggedCells=A2ataggedCells(1:excluded_count);
fromWhichUnit=fromWhichUnit(1:trials_count);
fromWhichTrial=fromWhichTrial(1:trials_count);
isEventInThisTrial=isEventInThisTrial(1:trials_count);
fromWhichSess_forTrials=fromWhichSess_forTrials(1:trials_count);
if keepAllSingleTrials==false
    unitbyunit_x=unitbyunit_x(1:units_count,:);
    unitbyunit_y=unitbyunit_y(1:units_count,:);
    aligncomp_x=aligncomp_x(1:units_count,:);
    aligncomp_y=aligncomp_y(1:units_count,:);
else
    unitbyunit_x=unitbyunit_x(1:trials_count,:);
    unitbyunit_y=unitbyunit_y(1:trials_count,:);
    aligncomp_x=aligncomp_x(1:trials_count,:);
    aligncomp_y=aligncomp_y(1:trials_count,:);

end


% if settings.studyOptoTag==true
%     [firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess]=studyOptoTagging(unitbyunit_y, unitbyunit_x, settings, fromWhichSess, indsForAnalysisPerSess);
% else
%     firstValNtimesBaseVar=[];
%     optoFRoverBaseline=[];
%     indsForAnalysisPerSess=[];
% end

if ds~=1
    unitbyunit_x=downSampMatrix(unitbyunit_x,ds);
    unitbyunit_y=downSampMatrix(unitbyunit_y,ds);
    aligncomp_x=downSampMatrix(aligncomp_x,ds);
    aligncomp_y=downSampMatrix(aligncomp_y,ds);
end

if normalizePSTHs==true
    fi=find(nanmean(unitbyunit_x,1)>min(nanmean(unitbyunit_x,1))+cutAtTime*2,1,'first');
    mi=min(unitbyunit_y(:,1:fi),[],2,'omitnan');
    unitbyunit_y=unitbyunit_y-repmat(mi,1,size(unitbyunit_y,2));
    ma=max(unitbyunit_y(:,1:fi),[],2,'omitnan');
    unitbyunit_y=unitbyunit_y./repmat(ma,1,size(unitbyunit_y,2));
end

if ~suppressPlots
    figure();
    plot(nanmean(unitbyunit_x,1),unitbyunit_y');
    hold on;
    plot(nanmean(aligncomp_x,1),aligncomp_y','Color','b');

    figure();
    plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1),'Color','k');
    hold on;
    plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)-nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');
    plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1)+nanstd(unitbyunit_y,[],1)./sqrt(size(unitbyunit_y,1)),'Color','k');

    plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');
    hold on;
    plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)-nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');
    plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1)+nanstd(aligncomp_y,[],1)./sqrt(size(aligncomp_y,1)),'Color','b');

    if getCDs==true
        %     % for cue
        cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 1.25]); % cue duration is 250 ms
        %     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 0.25]); % cue duration is 250 ms
        uncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-1.75 -0.25]);
        lateUncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[4 9.5]);
        % for reach
        %     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-1 0]); % around reach
        % %     cueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[-0.25 0.25]); % cue duration is 250 ms
        %     uncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[3 9]); % after reach
        %     lateUncueCD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,[4 9.5]);
    else
        % use CDs passed in
    end

    % If want to cut off the plotting at a particular time after the
    % aligncomp_y max, use cutAtTime
    % else
    % cutAtTime=[];
    plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD,cutAtTime);
end

end

% function [firstValNtimesBaseVar,optoFRoverBaseline,indsForAnalysisPerSess]=studyOptoTagging(unitbyunit_y, unitbyunit_x, settings, fromWhichSess, indsForAnalysisPerSess)
% 
% % find opto onset for these experiments
% uniqueSess=unique(fromWhichSess);
% firstValNtimesBaseVar=[];
% optoFRoverBaseline=[];
% if isempty(indsForAnalysisPerSess)
%     indsForAnalysisPerSess=nan(length(uniqueSess),2);
%     usePassedInInds=false;
% else
%     usePassedInInds=true;
% end
% for i=1:length(uniqueSess)
%     takeTheseCells=fromWhichSess==uniqueSess(i);
%     if nansum(takeTheseCells)==0
%         continue
%     end
%     if isnan(indsForAnalysisPerSess(i,1))
%         usePassedInIndsLocal=false;
%     else
%         if usePassedInInds==true
%             usePassedInIndsLocal=true;
%         else
%             usePassedInIndsLocal=false;
%         end
%     end
%     if takeTheseCells(end)>size(unitbyunit_y,1)
%         takeTheseCells=i:size(unitbyunit_y,1);
%     end 
%     temp=nanmean(unitbyunit_y(takeTheseCells,1:settings.unitbaseline),1);
%     tempvar=var(temp(1:settings.beforeOptoBaseline));
%     tempmean=mean(temp(1:settings.beforeOptoBaseline));
%     if usePassedInIndsLocal
%         baselineVariance=var(unitbyunit_y(takeTheseCells,1:indsForAnalysisPerSess(i,1)-1),0,2);
%         baselineMean=mean(unitbyunit_y(takeTheseCells,1:indsForAnalysisPerSess(i,1)-1),2,'omitnan');
%         firstOptoVal=unitbyunit_y(takeTheseCells,indsForAnalysisPerSess(i,1));
%         timestep=nanmean(diff(nanmean(unitbyunit_x,1)));
%         baselineVariance(baselineVariance<0.001)=0;
%         baselineMean(baselineMean<0.001)=0;
%         firstOptoVal(firstOptoVal<0.001)=0;
%         firstValRelativeToBaselineInVar=(firstOptoVal-baselineMean)./baselineVariance;
%         avDuringOpto=mean(unitbyunit_y(takeTheseCells,indsForAnalysisPerSess(i,1):indsForAnalysisPerSess(i,2)),2,'omitnan');
%     else
%         ma=nanmax((temp-tempmean)./tempvar);
%         f1=find(((temp-tempmean)./tempvar)>0.3*ma,1,'first');
%         if isempty(f1)
%             [~,f1]=nanmax((temp-tempmean)./tempvar);
%         end
%         baselineVariance=var(unitbyunit_y(takeTheseCells,1:f1-1),0,2);
%         baselineMean=mean(unitbyunit_y(takeTheseCells,1:f1-1),2,'omitnan');
%         firstOptoVal=unitbyunit_y(takeTheseCells,f1);
%         timestep=nanmean(diff(nanmean(unitbyunit_x,1)));
%         baselineVariance(baselineVariance<0.001)=0;
%         baselineMean(baselineMean<0.001)=0;
%         firstOptoVal(firstOptoVal<0.001)=0;
%         firstValRelativeToBaselineInVar=(firstOptoVal-baselineMean)./baselineVariance;
%         avDuringOpto=mean(unitbyunit_y(takeTheseCells,f1:f1+floor(settings.optoTagDuration/timestep)),2,'omitnan');
%         indsForAnalysisPerSess(i,1)=f1; indsForAnalysisPerSess(i,2)=f1+floor(settings.optoTagDuration/timestep);
%     end
%     firstValNtimesBaseVar=[firstValNtimesBaseVar; firstValRelativeToBaselineInVar];
%     optoFRoverBaseline=[optoFRoverBaseline; avDuringOpto./baselineMean];
% end
% 
% end

function CD=findCD(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,windowWrtAlignCompMax)

% Find cue coding dimension
% Find max of alignment companion
[~,ma]=nanmax(nanmean(aligncomp_y,1));
aligncomp_times=nanmean(aligncomp_x,1);
subtractOff=aligncomp_times(ma);
aligncomp_times=aligncomp_times-subtractOff;
unit_times=nanmean(unitbyunit_x,1);
unit_times=unit_times-subtractOff;
[~,beginIndsForCD]=nanmin(abs(unit_times-windowWrtAlignCompMax(1)));
[~,endIndsForCD]=nanmin(abs(unit_times-windowWrtAlignCompMax(2)));
indsForCD=beginIndsForCD:endIndsForCD;
disp(['Calculating CD for time window ' num2str(unit_times(indsForCD(1))) ' to ' num2str(unit_times(indsForCD(end)))]);

% Find population vector during window
CD=nanmean(unitbyunit_y(:,indsForCD),2);

end

function plotEvolutionOfPopulationVector(unitbyunit_x,unitbyunit_y,aligncomp_x,aligncomp_y,cueCD,uncueCD,cutAtTime)

[~,ma]=nanmax(nanmean(aligncomp_y,1));
aligncomp_times=nanmean(aligncomp_x,1);
subtractOff=aligncomp_times(ma);
aligncomp_times=aligncomp_times-subtractOff;
unit_times=nanmean(unitbyunit_x,1);
unit_times=unit_times-subtractOff;

if ~isempty(cutAtTime)
    % find ind into aligncomp for cutoff
    [~,cutoffind]=nanmin(abs(aligncomp_times-cutAtTime));
    % find ind into unitbyunit for cutoff
    [~,cutoffind2]=nanmin(abs(unit_times-cutAtTime));
    aligncomp_x=aligncomp_x(:,1:cutoffind);
    aligncomp_times=aligncomp_times(1:cutoffind);
    aligncomp_y=aligncomp_y(:,1:cutoffind);
    unitbyunit_x=unitbyunit_x(:,1:cutoffind2);
    unit_times=unit_times(1:cutoffind2);
    unitbyunit_y=unitbyunit_y(:,1:cutoffind2);    
end

evolveProjectOntoCue=nan(1,size(unitbyunit_y,2));
evolveProjectOntoUncue=nan(1,size(unitbyunit_y,2));
for i=1:size(unitbyunit_y,2)
    evolveProjectOntoCue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),cueCD);
    evolveProjectOntoUncue(i)=projectTimePointOntoCD(unitbyunit_y(:,i),uncueCD);
end

cmap=colormap('cool');
stepintocmap=size(cmap,1)/nansum(~isnan(evolveProjectOntoCue));
indsintocmap=0:stepintocmap:(nansum(~isnan(evolveProjectOntoCue))-1)*stepintocmap;
indsintocmap(1)=0.0001;

figure();
k=1;
for i=1:length(evolveProjectOntoCue)
    if isnan(evolveProjectOntoCue(i)) | isnan(evolveProjectOntoUncue(i))
        continue
    end
    scatter(evolveProjectOntoUncue(i),evolveProjectOntoCue(i),[],cmap(ceil(indsintocmap(k)),:));
    k=k+1;
    hold on;
end
daspect([1 1 1]);
xlabel('Projection onto uncued CD');
ylabel('Projection onto cued CD');
linemi=min([evolveProjectOntoUncue evolveProjectOntoCue],[],2,'omitnan');
linema=max([evolveProjectOntoUncue evolveProjectOntoCue],[],2,'omitnan');
line([linemi linema],[linemi linema],'Color',[0.8 0.8 0.8]);

figure(); plot(nanmean(aligncomp_x,1),nanmean(aligncomp_y,1),'Color','b');
hold on; 
k=1;
for i=1:length(evolveProjectOntoCue)
    if isnan(evolveProjectOntoCue(i)) | isnan(evolveProjectOntoUncue(i))
        continue
    end
    scatter(nanmean(unitbyunit_x(:,i),1),nanmean(unitbyunit_y(:,i),1),[],cmap(ceil(indsintocmap(k)),:));
    k=k+1;
    hold on;
end
plot(nanmean(unitbyunit_x,1),nanmean(unitbyunit_y,1),'Color','k');

end

function proj=projectTimePointOntoCD(activity,CD)

% Project activity onto CD
% Find unit vector of CD
unit_CD=CD./norm(CD);
proj=dot(activity,unit_CD);

end