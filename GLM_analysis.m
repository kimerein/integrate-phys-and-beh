function GLM_analysis(whichSess,dd)
% wrapper for first pass at GLM, calls B's poissModel

if length(whichSess)>1
    dd=dd(whichSess);
end

whichUnitsToGrab='_'; plotUnitCriteria=[-100 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
setForUn=settingsForStriatumUnitPlots;
if setForUn.keepAllSingleTrials~=true
    error('need trial by trial data for GLM analysis');
end
response_to_plot='cue';
if length(whichSess)>1
    dd_more=cell(1,length(dd));
    dd_m=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
        dd_m{i}=[dd{i}];
    end
    ResponseCued=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
    [~,ma]=max(mean(ResponseCued.aligncomp_y,1,'omitnan'),[],2,'omitnan');
    temp=mean(ResponseCued.aligncomp_x,1,'omitnan');
    grabOtherBehaviorEvents(dd_m,mean(ResponseCued.unitbyunit_x,1,'omitnan')-temp(ma));
else
    ResponseCued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
    [~,ma]=max(mean(ResponseCued.aligncomp_y,1,'omitnan'),[],2,'omitnan');
    temp=mean(ResponseCued.aligncomp_x,1,'omitnan');
    grabOtherBehaviorEvents(dd{whichSess},mean(ResponseCued.unitbyunit_x,1,'omitnan')-temp(ma));
end
ResponseCued.unitbyunit_x=downSampMatrix(ResponseCued.unitbyunit_x,downSampBy);
ResponseCued.unitbyunit_y=downSampMatrix(ResponseCued.unitbyunit_y,downSampBy);
ResponseCued=makeUnitsUnique(ResponseCued);




% temp=Response.unitbyunit_y; temp(isnan(temp))=0;
% behEvents=zeros(size(temp));
% behEvents(:,50)=1;
% poissModel(temp,0:0.04:(size(temp,2)-1)*0.04,behEvents);

end

function grabOtherBehaviorEvents(datadir,unitTimes)

getEventsFromPhysTbt={'opto','distractor'};
getEventsFromBehTbt={'all_reachBatch','isFidgeting','success_fromPerchOrWheel',...
    'drop_fromPerchOrWheel','misses_and_pelletMissing','misses_and_pelletMissing_and_drop',...
    'success batch when pellet dislodged','drop batch when pellet dislodged','isChewing'};

if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    f=regexp(datadir,sep);
    tbtin=[datadir(1:f(end)) 'tbt'];
    disp(['reading in beh events from ' tbtin]);
    ls=dir(tbtin);
    for i=3:length(ls)
        a=[];
        if strcmp(ls(i).name,'beh2_tbt.mat')
            a=load([ls(i).folder sep ls(i).name]);
            beh2_tbt=a.beh2_tbt;
        elseif strcmp(ls(i).name,'physiology_tbt.mat')
            a=load([ls(i).folder sep ls(i).name]);
            phys_tbt=a.physiology_tbt;
        else
            continue
        end
    end
    beh2_tbt=addFidgeting(datadir(1:f(end)-1), beh2_tbt);
    if ~all(size(beh2_tbt.isFidgeting)==size(beh2_tbt.cue))
        % problem with isFidgeting size
        disp(['problem with fidgets size for ' datadir(1:f(end))]);
        beh2_tbt.isFidgeting=nan(size(beh2_tbt.cue));
    end
    % get behavior events
    fromPhystbtTimes=phys_tbt.('cuetimes_wrt_trial_start');
    cuePhystbt=phys_tbt.('cue');
    behtimes=getEventsOfType('times_wrt_trial_start',beh2_tbt);
    % how to map to unit times
    [indsIntoBeh_step1,indsIntoBeh_step2]=mapToPhysTimes(behtimes,fromPhystbtTimes,cuePhystbt,unitTimes);
    for i=1:length(getEventsFromBehTbt)
        temp=getEventsOfType(getEventsFromBehTbt{i},beh2_tbt);
        % map to unit times
        tempinunittimes=mapToUnitTimes(temp,true,indsIntoBeh_step1,indsIntoBeh_step2,fromPhystbtTimes,unitTimes);
    end
    for i=1:length(getEventsFromPhysTbt)
        temp=phys_tbt.(getEventsFromPhysTbt{i});
        % map to unit times
        tempinunittimes=mapToUnitTimes(temp,false,[],indsIntoBeh_step2,fromPhystbtTimes,unitTimes);
    end
end
end

function step2=mapToUnitTimes(behVals,isBeh,indsStep1,indsStep2,fromPhystbtTimes,unitTimes)

if isBeh==true
    % behavior tbt
    step1=zeros(size(fromPhystbtTimes));
    step2=zeros(size(fromPhystbtTimes,1),size(unitTimes,2));
    for i=1:size(step1,1)
        % gives inds into fromPhystbtTimes but of length of behtimes
        for j=1:length(behVals(i,:))
            step1(i,indsStep1(i,j))=step1(i,indsStep1(i,j))+behVals(i,j);
        end
        % gives inds into unitTimes but of length of fromPhystbtTimes
        for j=1:length(step1(i,:))
            step2(i,indsStep2(i,j))=step2(i,indsStep2(i,j))+step1(i,j);
        end
    end
    step2=step2>0.5;
else
    % physiology tbt
    % skip step 1
    step2=zeros(size(fromPhystbtTimes,1),size(unitTimes,2));
    for i=1:size(step2,1)
        % gives inds into unitTimes but of length of fromPhystbtTimes
        for j=1:length(behVals(i,:))
            step2(i,indsStep2(i,j))=step2(i,indsStep2(i,j))+behVals(i,j);
        end
    end
    step2=step2>0.5;
end

end

function [indsIntoBeh_step1,indsIntoBeh_step2]=mapToPhysTimes(behtimes,fromPhystbtTimes,cuePhystbt,unitTimes)

if size(behtimes,1)~=size(fromPhystbtTimes,1)
    error('behtimes and fromPhystbtTimes must have same number rows in GLM_analysis.m');
end

% behtimes and fromPhystbtTimes should be the same times
% just different sized arrays
% unitTimes is the matrix from Response.unitbyunit_x, i.e., aligned to cue
indsIntoBeh_step1=nan(size(behtimes));
indsIntoBeh_step2=nan(size(fromPhystbtTimes));
% every behavior timepoint goes to a phys timepoint
for i=1:size(fromPhystbtTimes,1)
    % do it this way because cue detection in physiology is more precise
    % find closest to behtimes in fromPhystbtTimes  
    % gives inds into fromPhystbtTimes but of length of behtimes
    closest=findClosest(behtimes(i,:),fromPhystbtTimes(i,:));
    % also realigned at cue
    fcue=find(cuePhystbt(i,:)>0.5,1,'first');
    cuetime=fromPhystbtTimes(i,fcue);
    % find closest in unitTimes
    % gives inds into unitTimes but of length of fromPhystbtTimes
    nextstep=findClosest(fromPhystbtTimes(i,:)-cuetime,unitTimes);
    indsIntoBeh_step1(i,:)=closest;
    indsIntoBeh_step2(i,:)=nextstep;
end

end

function closest=findClosest(a,in_b)

closest=nan(size(a));
for i=1:length(a)
    [~,closest(i)]=min(abs(a(i)-in_b));
end

end

function beh2_tbt=addFidgeting(direc, beh2_tbt)
% processed_data directory must contain
mustcontain={'humanchecked.txt','fixed_tbt_success_v_drop.txt','fixed_miss_v_grab.txt'};

% go into O2 output directory
% load autoReachSettings.mat and fidget.mat for each video

% first find processed_data folders and associated names
direc=[direc sep 'O2 output'];
ls=dir(direc);
vidnames={};
foldnames={};
k=1;
for i=3:length(ls)
    r=regexp(ls(i).name,'processed_data');
    if ~isempty(r)
        % check whether this data was humanchecked
        dontuse=false;
        for check=1:length(mustcontain)
            if dircontains([direc sep ls(i).name],mustcontain{check})==false
                dontuse=true;
                break
            end
        end
        if dontuse==true
            continue
        end
        foldnames{k}=ls(i).name;
        f=regexp(ls(i).name,'_processed_data');
        vidnames{k}=ls(i).name(1:f-1);
        k=k+1;
    else
        continue
    end
end
    
% then order names
[vidnames,si]=sort(vidnames);
foldnames=foldnames(si);

% for each, load autoReachSettings and fidget
for i=1:length(vidnames)
    a=load([direc sep vidnames{i} '_autoReachSettings.mat']);
    settings=a.settings;
    a=load([direc sep vidnames{i} '_fidget.mat']);
    fidget=a.fidget;
    beh2_tbt=addBackFidgets(beh2_tbt,settings,fidget,i);
end

end

function out=dircontains(d,fn)

ls=dir(d);
out=false;
for i=3:length(ls)
    if ~isempty(regexp(ls(i).name,fn))
        out=true;
        break
    end
end

end

function out=getEventsOfType(evT,behavior_tbt)

useCombo=false;
switch evT
    case 'cue_followedby_success'
        useReach='combo';
        [~,f]=nanmax(nanmean(behavior_tbt.cueZone_onVoff,1));
        temp=behavior_tbt.cueZone_onVoff;
        temp(~any(behavior_tbt.reachBatch_success_reachStarts(:,f:end)>0.5,2),:)=0;
        useCombo=temp;
    case 'cue_followedby_late_success'
        useReach='combo';
        [~,f]=nanmax(nanmean(behavior_tbt.cueZone_onVoff,1));
        temp=behavior_tbt.cueZone_onVoff;
        temp(~any(behavior_tbt.reachBatch_success_reachStarts(:,f+75:end)>0.5,2),:)=0;
        useCombo=temp;
    case 'cue_noReach'
        useReach='combo';
        [~,f]=nanmax(nanmean(behavior_tbt.cueZone_onVoff,1));
        temp=behavior_tbt.cueZone_onVoff;
        temp(any(behavior_tbt.reachStarts(:,1:f+115)>0.5,2),:)=0;
        useCombo=temp;
    case 'success_fromPerchOrWheel'
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_success_reachStarts+behavior_tbt.reachBatch_success_reachStarts_pawOnWheel;
    case 'drop_fromPerchOrWheel'
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
    case 'miss_fromPerchOrWheel'
        useReach='combo';
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel;
    case 'misses_and_pelletMissing'
        useReach='combo';
        excludeAfterSuccess=true;
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel+behavior_tbt.pelletmissingreach_reachStarts;
        if excludeAfterSuccess==true
            for i=1:size(useCombo,1)
                f=find(useCombo(i,:)>0.5,1,'first');
                if isempty(f)
                    continue
                end
                for j=1:length(f)
                    currf=f(j);
                    starter=currf-excludeWithinInds;
                    if starter<1
                        starter=1;
                    end
                    if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                        % success before this, exclude
                        useCombo(i,f)=0;
                    end
                end
            end
        end
    case 'failure_noSuccessBeforeAndNoReachingAfter'
        useReach='combo';
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        excludeAfterTimeWindow=[0.5 4]; % seconds after
        excludeAfterInds=[floor(excludeAfterTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))))];
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel+behavior_tbt.pelletmissingreach_reachStarts+behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
        for i=1:size(useCombo,1)
            f=find(useCombo(i,:)>0.05,1,'first');
            if isempty(f)
                continue
            end
            % exclude if this reach preceded by a success within 6 sec
            for j=1:length(f)
                currf=f(j);
                starter=currf-excludeWithinInds;
                if starter<1
                    starter=1;
                end
                if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                    % success before this, exclude
                    upto=currf+10;
                    if upto>length(useCombo(i,:))
                        upto=length(useCombo(i,:));
                    end
                    useCombo(i,starter:upto)=0;
                end
            end
            % exclude if this reach followed by any additional reach more than
            % 0.5 seconds later
            for j=1:length(f)
                currf=f(j);
                starter=currf+excludeAfterInds(1);
                if starter>size(behavior_tbt.success_reachStarts,2)
                    starter=size(behavior_tbt.success_reachStarts,2);
                end
                ender=currf+excludeAfterInds(2);
                if ender>size(behavior_tbt.success_reachStarts,2)
                    ender=size(behavior_tbt.success_reachStarts,2);
                end
                if any(behavior_tbt.all_reachBatch(i,starter:ender)>0.5)
                    % reach after this, exclude
                    upto=currf+10;
                    if upto>length(useCombo(i,:))
                        upto=length(useCombo(i,:));
                    end
                    fromto=currf-10;
                    if fromto<1
                        fromto=1;
                    end
                    useCombo(i,fromto:upto)=0;
                end
            end
        end
    case 'misses_and_pelletMissing_and_drop'
        useReach='combo';
        excludeAfterSuccess=true;
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.reachBatch_miss_reachStarts_pawOnWheel+behavior_tbt.pelletmissingreach_reachStarts+behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
        if excludeAfterSuccess==true
            for i=1:size(useCombo,1)
                f=find(useCombo(i,:)>0.05,1,'first');
                if isempty(f)
                    continue
                end
                for j=1:length(f)
                    currf=f(j);
                    starter=currf-excludeWithinInds;
                    if starter<1
                        starter=1;
                    end
                    if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                        % success before this, exclude
                        upto=currf+10; 
                        if upto>length(useCombo(i,:))
                            upto=length(useCombo(i,:));
                        end
                        useCombo(i,starter:upto)=0;
                    end
                end
            end
        end
    case 'misses_and_pelletMissing_and_drop_noPawOnWheel'
        useReach='combo';
        excludeAfterSuccess=true;
        excludeWithinTimeWindow=6; % in sec
        excludeWithinInds=floor(excludeWithinTimeWindow/mode(diff(nanmean(behavior_tbt.times,1))));
        useCombo=behavior_tbt.reachBatch_miss_reachStarts+behavior_tbt.pelletmissingreach_reachStarts+behavior_tbt.reachBatch_drop_reachStarts;
        if excludeAfterSuccess==true
            for i=1:size(useCombo,1)
                f=find(useCombo(i,:)>0.05,1,'first');
                if isempty(f)
                    continue
                end
                for j=1:length(f)
                    currf=f(j);
                    starter=currf-excludeWithinInds;
                    if starter<1
                        starter=1;
                    end
                    if any(behavior_tbt.success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts(i,starter:currf)>0.5 | behavior_tbt.reachBatch_success_reachStarts_pawOnWheel(i,starter:currf)>0.5)
                        % success before this, exclude
                        upto=currf+10; 
                        if upto>length(useCombo(i,:))
                            upto=length(useCombo(i,:));
                        end
                        useCombo(i,starter:upto)=0;
                    end
                end
            end
        end
    case 'reachBatch_success_reachStarts'
        useReach=alignTo;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
    case 'reachBatch_drop_reachStarts'
        useReach=alignTo;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
    case 'success batch when pellet dislodged'
        useReach='combo';
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'success batch av reach and pellet dislodged'
        useReach='combo';
        behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'drop batch when pellet dislodged'
        useReach='combo';
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'drop batch av reach and pellet dislodged'
        useReach='combo';
        behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,'reachBatch_drop_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    otherwise 
        useReach=evT;
end
if strcmp(useReach,'combo')
    out=useCombo;
else
    out=behavior_tbt.(useReach);
end

end