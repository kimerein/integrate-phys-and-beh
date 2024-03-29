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
excluded=nan(maxUnitsPerSess*length(dd),1);
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
    currTrainingSet=[]; % different training set from each session
    if isempty(ls)
        continue
    end
    % COMMONtrainingSet.mat is always in cue folder
    rcue=regexp(ls(1).folder,'\');
    if exist([ls(1).folder(1:rcue(end)) 'cue' sep 'COMMONtrainingSet.mat'],'file')
        a=load([ls(1).folder(1:rcue(end)) 'cue' sep 'COMMONtrainingSet.mat']);
        currTrainingSet=a.currTrainingSet;
    end
    for i=3:length(ls)
        a=[];
        if contains(ls(i).name,'trainingSet') || contains(ls(i).name,'testSet') % ignore these files
            continue
        end

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

        if ~isfield(settings,'isPhotometry')
            settings.isPhotometry=false;
        end

        if settings.isPhotometry==false
            [unitTest,dontUseTrials]=doUnitTest(ls(i).folder, ls(i).name);
            if unitTest
            else
                if settings.isPhotometry==true
                else
                    excluded(excluded_count)=1;
                    excluded_count=excluded_count+1;
                    continue % failed unit test
                end
            end
        end

        if isempty(a.dataout)
            disp([ls(i).folder '\' ls(i).name ' is empty ... skipping']);
            excluded(excluded_count)=1;
            excluded_count=excluded_count+1;
            continue
        end

        if settings.isPhotometry==true
            dontUseTrials=zeros(size(a.dataout.y,1),1);
        end

        if length(dontUseTrials)~=size(a.dataout.y,1)
            % as a default, take all trials
            % but I might want to modify this behavior later
            % sometimes dontUseTrials is empty if figure opening failed,
            % or other
            disp('dontUseTrials length does not match size of a.dataout.y in line 91 of alignToCompanion.m, but taking all trials');
            dontUseTrials=zeros(size(a.dataout.y,1),1);
        end
        excluded(excluded_count)=0; % if got here, use this unit
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
            % only trials where opto on
            if settings.onlyTrialsWhereOptoDuringCue==true
                rlastslash=regexp(ls(i).folder,sep);
                if exist([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'opto_was_on.txt'],'file')
                    % then continue
                    tempOpto=load([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'physiology_tbt.mat']);
                    minTrialInds=floor(settings.minTrialLength./mode(diff(nanmean(tempOpto.physiology_tbt.cuetimes_wrt_trial_start,1))));
                    isOverlapping=any(tempOpto.physiology_tbt.opto(:,1:minTrialInds)==1 & tempOpto.physiology_tbt.cue(:,1:minTrialInds)==1,2);
                    a.dataout.y(~isOverlapping(usingTrialsInds)==1,:)=nan;
                    a.alignComp.y(~isOverlapping(usingTrialsInds)==1,:)=nan;
                else
                    % can't use any trials
                    a.dataout.y(1:size(a.dataout.y,1),:)=nan;
                    a.alignComp.y(1:size(a.dataout.y,1),:)=nan;
                end
            end

            % discard opto trials if opto was on during cue
            if settings.discardTrialsWhereOptoDuringCue==true
                rlastslash=regexp(ls(i).folder,sep);
                if exist([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'opto_was_on.txt'],'file')
                    % then continue
                    tempOpto=load([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'physiology_tbt.mat']);
                    minTrialInds=floor(settings.minTrialLength./mode(diff(nanmean(tempOpto.physiology_tbt.cuetimes_wrt_trial_start,1))));
                    isOverlapping=any(tempOpto.physiology_tbt.opto(:,1:minTrialInds)==1 & tempOpto.physiology_tbt.cue(:,1:minTrialInds)==1,2);
                    a.dataout.y(isOverlapping(usingTrialsInds)==1,:)=nan;
                    a.alignComp.y(isOverlapping(usingTrialsInds)==1,:)=nan;
                end
            end

            % discard trials if any opto
            if settings.discardTrialsIfAnyOpto==true
                rlastslash=regexp(ls(i).folder,sep);
                if exist([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'opto_was_on.txt'],'file')
                    % then continue
                    tempOpto=load([ls(i).folder(1:rlastslash(end-1)) 'tbt' sep 'physiology_tbt.mat']);
                    minTrialInds=floor(settings.minTrialLength./mode(diff(nanmean(tempOpto.physiology_tbt.cuetimes_wrt_trial_start,1))));
                    isOverlapping=any(tempOpto.physiology_tbt.opto(:,1:minTrialInds)==1,2);
                    a.dataout.y(isOverlapping(usingTrialsInds)==1,:)=nan;
                    a.alignComp.y(isOverlapping(usingTrialsInds)==1,:)=nan;
                end
            end
            
            % discard drops
            if settings.discardDrops==true
                rlastslash=regexp(ls(i).folder,sep);
                lastunderscore=regexp(ls(i).name,'_');
                tempDrop=load([ls(i).folder(1:rlastslash(end)) settings.dropFolderName sep ls(i).name(1:lastunderscore(end)) settings.dropFileName '.mat']);
                if ~isempty(tempDrop.dataout)
                    tempDrop.dataout.y=tempDrop.dataout.y(usingTrialsInds,:);
                    tempDrop.alignComp.y=tempDrop.alignComp.y(usingTrialsInds,:);
                    dropHappens=any(~isnan(tempDrop.alignComp.y),2);
                    % nan out the drops
                    a.dataout.y(dropHappens==1,:)=nan;
                    a.alignComp.y(dropHappens==1,:)=nan;
                end
            end

            if keepAllSingleTrials==false
                if settings.useSameTrainingSetForAllNeurons==true && settings.useTestSet==true
                    a.dataout.y(ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                    unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                    unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
                elseif settings.useSameTrainingSetForAllNeurons==true && settings.useTrainingSet==true
                        a.dataout.y(~ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                        unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                        unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
                elseif settings.makeTrainingSet==true
                    f=find(eventHappens==1);
                    temp=randperm(nansum(eventHappens==1));
                    tempindsend=ceil(nansum(eventHappens==1)*settings.fracForTrainingSet);
                    if tempindsend>length(temp)
                        tempindsend=length(temp);
                    end
                    trainingSet=temp(1:tempindsend);
                    trainingSetTrials=f(trainingSet);
                    tempindsstart=tempindsend+1;
                    if tempindsstart>length(temp)
                        % no test set
                        testSet=[];
                        testSetTrials=[];
                    else
                        testSet=temp(tempindsstart:end);
                        testSetTrials=f(testSet);
                    end
                    save([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_testSet.mat'],'testSet','testSetTrials');
                    save([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_trainingSet.mat'],'trainingSet','trainingSetTrials');
                    unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                    unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(trainingSetTrials,1:upTo),1)];
                elseif settings.useTestSet==true && settings.useSameTrainingSetForAllNeurons==false
                    if ~isempty(settings.useTheseTestSets)
                        b.testSet=[]; b.testSetTrials=[];
                        for itthroughtest=1:length(settings.useTheseTestSets)
                            if contains(settings.useTheseTestSets{itthroughtest},'Intersect')
                                % only take trials from this test set that
                                % are also in current eventHappens
                                % i.e., intersection of test set and
                                % current set
                                doingIntersect=true;
                                tsstr=settings.useTheseTestSets{itthroughtest};
                                settings.useTheseTestSets{itthroughtest}=tsstr(1:regexp(settings.useTheseTestSets{itthroughtest},'Intersect')-1);
                            else
                                doingIntersect=false;
                            end 
                            rlastslash=regexp(ls(i).folder,sep);
                            lastunderscore=regexp(ls(i).name,'_');
                            btemp=load([ls(i).folder(1:rlastslash(end)) settings.useTheseTestSets{itthroughtest} sep ls(i).name(1:lastunderscore(end)) settings.useTheseTestFilenames{itthroughtest} '_testSet.mat']);
                            if ~isempty(btemp.testSet)
                                if doingIntersect==true
                                    % find intersection of test set and
                                    % current eventHappens trials
                                    f=find(eventHappens==1);
                                    whichInBoth=find(ismember(f,btemp.testSetTrials));
                                    whichInBoth2=find(ismember(btemp.testSetTrials,f));
                                    b.testSet=[b.testSet btemp.testSet(whichInBoth2)]; b.testSetTrials=[b.testSetTrials; f(whichInBoth)];
                                else
                                    b.testSet=[b.testSet btemp.testSet]; b.testSetTrials=[b.testSetTrials; btemp.testSetTrials];
                                end
                            end
                        end
                        unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                        if isempty(b.testSet) % only one trial
                            unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                            unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
                        else
                            unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(b.testSetTrials,1:upTo),1)];
                        end
                    else
                        b=load([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_testSet.mat']);
                        unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                        if isempty(b.testSet) % only one trial
                            unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                            unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
                        else
                            unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(b.testSetTrials,1:upTo),1)];
                        end
                    end
                elseif settings.useTrainingSet==true && settings.useSameTrainingSetForAllNeurons==false
                    b=load([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_trainingSet.mat']);
                    unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                    unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(b.trainingSetTrials,1:upTo),1)];
                else
                    unitbyunit_x(units_count,:)=[a.dataout.x(1:upTo)];
                    unitbyunit_y(units_count,:)=[nanmean(a.dataout.y(:,1:upTo),1)];
                end
            else
                if settings.useTestSet==true && settings.useSameTrainingSetForAllNeurons==true
                    a.dataout.y(ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                    a.alignComp.y(ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                elseif settings.useSameTrainingSetForAllNeurons==true && settings.useTrainingSet==true
                    if isempty(currTrainingSet)
                        b=load([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_trainingSet.mat']);
                        a.dataout.y(~ismember(1:size(a.dataout.y,1),b.trainingSetTrials),:)=nan;
                        a.alignComp.y(~ismember(1:size(a.dataout.y,1),b.trainingSetTrials),:)=nan;
                        currTrainingSet=b.trainingSetTrials;
                        save([ls(i).folder sep 'COMMONtrainingSet.mat'],'currTrainingSet');
                    else
                        a.dataout.y(~ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                        a.alignComp.y(~ismember(1:size(a.dataout.y,1),currTrainingSet),:)=nan;
                    end
                elseif settings.useTrainingSet==true && settings.useSameTrainingSetForAllNeurons==false
                    b=load([ls(i).folder sep ls(i).name(1:regexp(ls(i).name,'.mat','once')-1) '_trainingSet.mat']);
                    a.dataout.y(~ismember(1:size(a.dataout.y,1),b.trainingSetTrials),:)=nan;
                    a.alignComp.y(~ismember(1:size(a.dataout.y,1),b.trainingSetTrials),:)=nan;
                elseif settings.useTestSet==true && settings.useSameTrainingSetForAllNeurons==false
                    error('Have not yet implemented settings.useTestSet==true && settings.useSameTrainingSetForAllNeurons==false && keepAllSingleTrials==true');
                end
                unitbyunit_x(trials_count:trials_count+size(a.dataout.y,1)-1,:)=repmat(a.dataout.x(1:upTo),size(a.dataout.y,1),1);
                unitbyunit_y(trials_count:trials_count+size(a.dataout.y,1)-1,:)=[a.dataout.y(:,1:upTo)];
            end
            eventHappens=any(~isnan(a.alignComp.y),2);
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
trials_count=trials_count-1;
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