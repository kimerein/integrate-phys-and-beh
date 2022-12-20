function GLM_analysis(whichSess,dd,downSampBy)
% wrapper for first pass at GLM, calls B's poissModel

% if want to put in more than one session, will need to debug this code

cutAllTrialsAtThisTime=10.5; % trial length in seconds

if length(whichSess)>1
    dd=dd(whichSess);
end

whichUnitsToGrab='_'; 
% will include unit if unitdets match the following
% [inStructure isFS isTAN isSPN isLowFRThin]
plotUnitCriteria=[-100 -100 -100 -100 -100]; 
% plotUnitCriteria=[-100 0 0 1 0]; % -100 is a wildcard, else 0 (false) and 1 (true)
getCriteriaForUnitsToPlot(plotUnitCriteria);
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
    [phystbtout,behtbtout,fromwhichday]=grabOtherBehaviorEvents(dd_m,mean(ResponseCued.unitbyunit_x,1,'omitnan')-temp(ma));
else
    ResponseCued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
    [~,ma]=max(mean(ResponseCued.aligncomp_y,1,'omitnan'),[],2,'omitnan');
    temp=mean(ResponseCued.aligncomp_x,1,'omitnan');
    [phystbtout,behtbtout,fromwhichday]=grabOtherBehaviorEvents(dd{whichSess},mean(ResponseCued.unitbyunit_x,1,'omitnan')-temp(ma));
end
% phystbtout.cue=shiftCueBack(phystbtout.cue);
ResponseCued=makeUnitsUnique(ResponseCued);

% cut all trials at this time
[phystbtout,behtbtout,ResponseCued]=cutAllTrialsToLength(phystbtout,behtbtout,ResponseCued,cutAllTrialsAtThisTime);

% figure(); plot(nanmean(ResponseCued.unitbyunit_y(ResponseCued.fromWhichUnit==7,:),1));
% hold on; plot(nanmean(phystbtout.cue,1),'Color','b');
% zero out artifact in behtbtout at the beginning of trial
% fie=fieldnames(behtbtout);
% for i=1:length(fie)
%     temp=behtbtout.(fie{i});
%     temp(:,7)=0;
%     behtbtout.(fie{i})=temp;
% end

u=unique(ResponseCued.fromWhichUnit);
neuron_data_matrix=[];
neuron_disappears=zeros(length(u),1);
gotBehEvents=false;
for i=1:length(u)
    % get one neuron's activity pattern across all trials
    whichNeuron=u(i);
    whichTrialsForThisCell=ResponseCued.fromWhichTrial(ResponseCued.fromWhichUnit==whichNeuron);
    dataMat=ResponseCued.unitbyunit_y(ResponseCued.fromWhichUnit==whichNeuron,:);
    whichSessForThisCell=ResponseCued.fromWhichSess_forTrials(ResponseCued.fromWhichUnit==whichNeuron);
    if nansum(ResponseCued.fromWhichUnit==i)<0.8*size(phystbtout.cue,1) 
        % throw out any neuron that is present for less than 80% of the
        % full session
        neuron_disappears(i)=1;
        continue
    elseif nansum(ResponseCued.fromWhichUnit==i)~=size(phystbtout.cue,1) 
        % fill in missing trials with zero spiking
        % find missing trials
        alltri=1:size(phystbtout.cue,1); 
        misstri=alltri(~ismember(alltri,whichTrialsForThisCell));
        newDataMat=nan(length(alltri),size(dataMat,2));
        newDataMat(whichTrialsForThisCell,:)=dataMat;
        newDataMat(misstri,:)=zeros(size(newDataMat(misstri,:)));
        dataMat=newDataMat;
    end
    if gotBehEvents==false
%         behEvents=assembleBehEvents(phystbtout,behtbtout,fromwhichday,whichTrialsForThisCell,whichSessForThisCell);
        behEvents=assembleBehEvents(phystbtout,behtbtout,fromwhichday,1:size(phystbtout.cue,1),mode(whichSessForThisCell)*ones(size(1:size(phystbtout.cue,1))));
        gotBehEvents=true;
    end
    dataMat=dataMat';
    dataMat=dataMat(1:end);
    if i==1
        neuron_data_matrix=nan(length(u),length(dataMat));
    end
    neuron_data_matrix(i,:)=dataMat;
end
neuron_data_matrix=neuron_data_matrix(neuron_disappears==0,:);

timestep=mode(abs(diff(mean(ResponseCued.unitbyunit_x,1,'omitnan'))));
neuron_data_matrix=downSampMatrix(neuron_data_matrix,downSampBy);
timepoints=downSampMatrix(0:timestep:(size(behEvents,2)-1)*timestep,downSampBy);
behEvents=downSampMatrix(behEvents,downSampBy);
behEvents(behEvents~=0)=1;
disp('behEvents rows are');
disp(phystbtout);
disp(behtbtout);
do_glm(neuron_data_matrix,timepoints,behEvents);

end

function [phystbtout,behtbtout,ResponseCued]=cutAllTrialsToLength(phystbtout,behtbtout,ResponseCued,cutAllTrialsAtThisTime)

timepoints=mean(ResponseCued.unitbyunit_x,1,'omitnan');
timepoints=timepoints-min(timepoints,[],'all','omitnan');
figure(); plot(timepoints); title('Trial timepoints'); ylabel('sec');
f=find(timepoints>cutAllTrialsAtThisTime,1,'first');
if isempty(f)
    return
end
ResponseCued.unitbyunit_x=ResponseCued.unitbyunit_x(:,1:f);
ResponseCued.unitbyunit_y=ResponseCued.unitbyunit_y(:,1:f);
fnames=fieldnames(phystbtout);
for i=1:length(fnames)
    temp=phystbtout.(fnames{i});
    phystbtout.(fnames{i})=temp(:,1:f);
end
fnames=fieldnames(behtbtout);
for i=1:length(fnames)
    temp=behtbtout.(fnames{i});
    behtbtout.(fnames{i})=temp(:,1:f);
end

end

function cue=shiftCueBack(cue)

% have to do this because small capacitance to IR LED only on WHISPER rig
cue=[cue(:,2:end) zeros(size(cue,1),1)];

end

function behEvents=assembleBehEvents(phystbtout,behtbtout,fromwhichday,whichTrialsForThisCell,whichSessForThisCell)

behEvents=[];
f=fieldnames(phystbtout);
countEvType=1;
for i=1:length(f)
    temp=phystbtout.(f{i}); % structure is row by row trials, col by col timepoints
    % want to concatenate trials where unit present
    % first take day with this unit
    % take trials for which this unit present
    temp=temp(fromwhichday==mode(whichSessForThisCell),:);
    temp=temp(ismember(1:size(temp,1),whichTrialsForThisCell),:);
    temp=temp';
    temp=temp(1:end);
    if i==1
        % initialize
        % concatenate all trials
        behEvents=nan(length(fieldnames(phystbtout))+length(fieldnames(behtbtout)),length(temp));
    end 
    behEvents(countEvType,:)=temp;
    countEvType=countEvType+1;
end
f=fieldnames(behtbtout);
for i=1:length(f)
    temp=behtbtout.(f{i}); 
    temp=temp(fromwhichday==mode(whichSessForThisCell),:);
    temp=temp(ismember(1:size(temp,1),whichTrialsForThisCell),:);
    temp=temp';
    temp=temp(1:end);
    behEvents(countEvType,:)=temp;
    countEvType=countEvType+1;
end

end

function do_glm(dataMatrix,tRange,events)

% each row of events is a different type of behavior event
neuronFiring=dataMatrix; % rows are neurons, cols are timepoints
neuronFiring(neuronFiring<0)=0;
nNeurons=size(dataMatrix,1);
nEventTypes=size(events,1);
bins=length(tRange);

% prepare time shift matrix of events for the glm
% maxShifts=100; % how many +/- bins to consider for glm
maxShifts=15;
shifts=-maxShifts:maxShifts;
nShifts=length(shifts);
allEvents=zeros(nEventTypes*length(shifts), bins);
for counter=1:nShifts
    e2=circshift(events, shifts(counter), 2);
    allEvents(nShifts*((1:nEventTypes)-1)+counter, :)=e2;
end

% setup for doing glms
testNeuron=1:nNeurons;
% a small gaussian filter used below to test how good the coefficients are

% do glms 
% we fit the spiking rate but look at the reconstruction of the 
% underlying time varying rates.  I guess we really model the lambdas, but
% same thing here

for neuronIndex=testNeuron
    figure('NumberTitle','off', 'Name', ['neuron ' num2str(neuronIndex)])
    set(gcf, 'Position', [  99         549        1386         317])

    % run model with time shifted events. This is what we would do
    mdl=fitglm(allEvents', neuronFiring(neuronIndex,:)', 'linear', 'Link', 'identity'); 
    % You can add 'Distribution', 'poisson' but I find that I get identical
    % results and without it, it runs much faster and converges better
    % Maybe it gets the right answer with a linear link  
    % because of the smearing and summing across multiple poisson
    % processes?

    subplot(1, 4, 1)
    pva=mdl.Coefficients.pValue(2:end);
    coef=mdl.Coefficients.Estimate(2:end);
    whichcoef=1:length(coef);
    scatter(whichcoef(pva<0.2),coef(pva<0.2),[],'k');
    hold on; 
    scatter(whichcoef(pva>=0.2),coef(pva>=0.2),[],[0.5 0.5 0.5]);
    cm=colormap('jet');
    colstep=floor(size(cm,1)/size(events,1));
    cols=1:colstep:size(cm,1);
    for i=1:size(events,1)
        line([(i-1)*nShifts+1 (i-1)*nShifts+nShifts],[0 0],'Color',cm(cols(i),:),'LineWidth',5);
        hold on;
    end
    title('model coeffs');

    % reconstruct the underlying rates, not the spiking!
    yy=predict(mdl, allEvents');

    subplot(1, 4, 2)
    plot(yy, 'DisplayName','recon rate')
    hold on
    plot(neuronFiring(neuronIndex, :), 'DisplayName','real neuron')
    legend
    title('reconstruction')

    subplot(1, 4, 3)
    % very smoothed
    plot(smoothdata(yy, 'gaussian', 30), 'DisplayName','recon rate')
    hold on
    plot(smoothdata(neuronFiring(neuronIndex, :), 'gaussian', 30), 'DisplayName','real neuron')
    legend
    title('very smoothed')

    subplot(1, 4, 4)
    pvalbins=0:0.01:1;
    for i=1:length(pvalbins)
        fracBelow(i)=nansum(mdl.Coefficients.pValue<pvalbins(i))/length(mdl.Coefficients.pValue);
    end
    scatter(pvalbins,fracBelow);
    xlabel('pval bins'); ylabel('fraction below this pval');
    line([0 1],[0 1]);
end

end

function Response=makeUnitsUnique(Response)

uSess=unique(Response.fromWhichSess_forTrials);
sessUOffset=0;
triOffset=0;
backupSess=Response.fromWhichSess_forTrials;
backupUs=Response.fromWhichUnit;
backupTrials=Response.fromWhichTrial;
for i=1:length(uSess)
    currSess=uSess(i);
    whichUs=unique(Response.fromWhichUnit(backupSess==currSess));
    for j=1:length(whichUs)
        currU=whichUs(j);
        Response.fromWhichUnit(backupSess==currSess & backupUs==currU)=currU+sessUOffset;
    end
    currTris=unique(backupTrials(backupSess==currSess));
    for j=1:length(currTris)
        currT=currTris(j);
        Response.fromWhichTrial(backupSess==currSess & backupTrials==currT)=currT+triOffset;
    end
    sessUOffset=sessUOffset+max(whichUs,[],'all','omitnan')+1;
    triOffset=triOffset+max(currTris,[],'all','omitnan')+1;
end

end

function [phystbtout,behtbtout,fromwhichday]=grabOtherBehaviorEvents(datadir,unitTimes)

getEventsFromPhysTbt={'cue','opto','distractor'};
getEventsFromBehTbt={'all_reachBatch','isFidgeting','success_fromPerchOrWheel',...
    'drop_fromPerchOrWheel','misses_and_pelletMissing','misses_and_pelletMissing_and_drop','isChewing'};

if iscell(datadir)
    dd=datadir;
else
    dd=1;
end
phystbtout.(getEventsFromPhysTbt{1})=[];
behtbtout.(getEventsFromBehTbt{1})=[];
fromwhichday=[];
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
        if ~isfield(behtbtout,getEventsFromBehTbt{i})
            behtbtout.(getEventsFromBehTbt{i})=tempinunittimes;
        elseif isempty(behtbtout.(getEventsFromBehTbt{i}))
            behtbtout.(getEventsFromBehTbt{i})=tempinunittimes;
        else
            behtbtout.(getEventsFromBehTbt{i})=[behtbtout.(getEventsFromBehTbt{i}); tempinunittimes];
        end
    end
    for i=1:length(getEventsFromPhysTbt)
        temp=phys_tbt.(getEventsFromPhysTbt{i});
        % map to unit times
        tempinunittimes=mapToUnitTimes(temp,false,[],indsIntoBeh_step2,fromPhystbtTimes,unitTimes);
        if ~isfield(phystbtout,getEventsFromPhysTbt{i})
            phystbtout.(getEventsFromPhysTbt{i})=tempinunittimes;
        elseif isempty(phystbtout.(getEventsFromPhysTbt{i}))
            phystbtout.(getEventsFromPhysTbt{i})=tempinunittimes;
        else
            phystbtout.(getEventsFromPhysTbt{i})=[phystbtout.(getEventsFromPhysTbt{i}); tempinunittimes];
        end
    end
    fromwhichday=[fromwhichday; j*ones(size(tempinunittimes,1),1)];
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

beh2_tbt=fixMissingVideos(beh2_tbt);
% for each, load autoReachSettings and fidget
for i=1:length(vidnames)
    a=load([direc sep vidnames{i} '_autoReachSettings.mat']);
    settings=a.settings;
    a=load([direc sep vidnames{i} '_fidget.mat']);
    fidget=a.fidget;
    beh2_tbt=addBackFidgets(beh2_tbt,settings,fidget,i);
end

end

function beh2_tbt=fixMissingVideos(beh2_tbt)

% there was a bug in physiology_triggered_on_reach.m, it seems,
% where second video not being populated
% if first and third videos, intervening trials must be second video
if isfield(beh2_tbt,'from_first_video') && isfield(beh2_tbt,'from_third_video')
    fromfirst=mean(beh2_tbt.from_first_video,2,'omitnan');
    fromthird=mean(beh2_tbt.from_third_video,2,'omitnan');
    if all(beh2_tbt.from_second_video==0,'all')
        beh2_tbt.from_second_video(~fromfirst & ~fromthird,:)=1;
        figure();
        imagesc([beh2_tbt.from_first_video beh2_tbt.from_second_video beh2_tbt.from_third_video]);
        title('Checking fix for missing values in field from_second_video in GLM_analysis');
    end
end

end

function tbt=fixDropOrSuccess(tbt)

% only one drop or success per trial -- take first
% just sometimes gets stuck at 1, interp failure
for i=1:size(tbt)
    temp=tbt(i,:);
    f=find(temp>0.5,1,'first');
    tempie=zeros(size(temp));
    tempie(f)=1;
    tbt(i,:)=tempie;
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

% behavior_tbt.reachBatch_drop_reachStarts=fixDropOrSuccess(behavior_tbt.reachBatch_drop_reachStarts);
% behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel=fixDropOrSuccess(behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel);
% behavior_tbt.reachBatch_success_reachStarts=fixDropOrSuccess(behavior_tbt.reachBatch_success_reachStarts);
% behavior_tbt.reachBatch_success_reachStarts_pawOnWheel=fixDropOrSuccess(behavior_tbt.reachBatch_success_reachStarts_pawOnWheel);

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
    case 'success_dislodged'
        useReach='combo';
        behavior_tbt.temp=behavior_tbt.reachBatch_success_reachStarts+behavior_tbt.reachBatch_success_reachStarts_pawOnWheel;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'temp','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
        behavior_tbt=rmfield(behavior_tbt,'temp');
    case 'success batch av reach and pellet dislodged'
        useReach='combo';
        behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,'reachBatch_success_reachStarts','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
    case 'drop_dislodged'
        useReach='combo';
        behavior_tbt.temp=behavior_tbt.reachBatch_drop_reachStarts+behavior_tbt.reachBatch_drop_reachStarts_pawOnWheel;
        behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,'temp','pelletPresent');
        useCombo=behavior_tbt.pelletDislodgedAfterSuccess;
        behavior_tbt=rmfield(behavior_tbt,'temp');
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

function behavior_tbt=findPelletDislodgedAvWithSuccess(behavior_tbt,successField,pelletPresent)

successes=behavior_tbt.(successField);
pellet=behavior_tbt.(pelletPresent);
behavior_tbt.pelletDislodgedAfterSuccess=nan(size(successes));
for i=1:size(successes,1)
    f=find(successes(i,:)>0.5);
    pelletGoneInds=nan(1,length(f));
    currpellet=pellet(i,:);
    for j=1:length(f)
        % find first index when pellet dislodged after success
        temp=f(j)-1+find(currpellet(f(j):end)<0.5,1,'first');
        if isempty(temp)
            continue
        end
        pelletGoneInds(j)=round(nanmean([temp f(j)],2));
    end
    behavior_tbt.pelletDislodgedAfterSuccess(i,:)=zeros(size(behavior_tbt.pelletDislodgedAfterSuccess(i,:)));
    behavior_tbt.pelletDislodgedAfterSuccess(i,pelletGoneInds)=1;
end

end

function behavior_tbt=findPelletDislodgedAfterSuccess(behavior_tbt,successField,pelletPresent)

successes=behavior_tbt.(successField);
pellet=behavior_tbt.(pelletPresent);
behavior_tbt.pelletDislodgedAfterSuccess=nan(size(successes));
for i=1:size(successes,1)
    f=find(successes(i,:)>0.5);
    pelletGoneInds=nan(1,length(f));
    currpellet=pellet(i,:);
    for j=1:length(f)
        % find first index when pellet dislodged after success
        temp=f(j)-1+find(currpellet(f(j):end)<0.5,1,'first');
        if isempty(temp)
            continue
        end
        pelletGoneInds(j)=temp;
    end
    behavior_tbt.pelletDislodgedAfterSuccess(i,:)=zeros(size(behavior_tbt.pelletDislodgedAfterSuccess(i,:)));
    if isnan(pelletGoneInds)
        behavior_tbt.pelletDislodgedAfterSuccess(i,end)=1;
    else
        behavior_tbt.pelletDislodgedAfterSuccess(i,pelletGoneInds)=1;
    end
end

end