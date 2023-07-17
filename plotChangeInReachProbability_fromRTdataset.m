function [out,sessidsperrow]=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,cueName,shuffleTrialOrder,settings)

epsilon_cue=settings.epsilon_cue; % in seconds
epsilon_uncue=settings.epsilon_uncue; % in seconds
epsilon_beforecue=settings.epsilon_beforecue; % in seconds
percentOfReachesFromSess_forInitCond=settings.percentOfReachesFromSess_forInitCond; % use this fraction of reaches from beginning of session to get initial conditions
percentOfReachesFromSess_forInitRate=settings.percentOfReachesFromSess_forInitRate; % ... for initial rates
maxTrialLength=settings.maxTrialLength; % in sec, wrt cue
minTrialLength=settings.minTrialLength; % wrt cue, in sec
acrossSess_window1=settings.acrossSess_window1;
acrossSess_window2=settings.acrossSess_window2;
acrossSess_window3=settings.acrossSess_window3;
initWindows=settings.initWindows;
addSatietyLines=settings.addSatietyLines;
useRateMethod=settings.useRateMethod; % Approach 1, 2 or 3 to plot and return
useWindowsForUncued=settings.useWindowsForUncued; % which windows to use to calculate uncued reaching, can be 2 or 3 or both
scatterPointSize=settings.scatterPointSize;

% find cue ind
if size(alltbt.(cueName),2)~=size(dataset.realDistributions.rawReaching_event_trial1InSeq{1},2)
    error('dataset size does not fit alltbt size');
end
[~,cueInd]=nanmax(nanmean(alltbt.(cueName),1));
% get timestep
timeStep=mode(diff(nanmean(alltbt.times)));

% three windows
% if t_n is time of reach on trial n
% t_cue is time of cue
% epsilon is size of time window including cue
% then these windows apply to reach time on trial n+1
% if t_n > t_cue -- if using reaction times for t_n, this must be true
window1 = @(t_cue,t_n,epsilon) [t_cue nanmin([t_n+epsilon maxTrialLength])];                             % cued window 
window2 = @(t_cue,t_n,epsilon) [nanmin([t_n+epsilon maxTrialLength]) maxTrialLength];                    % after cued window
window3 = @(t_cue,t_n,epsilon) [minTrialLength t_cue-epsilon];                  % before cue window     
% then need to divide number of reaches in window by duration of window to
% get rate

% if RT is nan, might have reached before cue, and Approaches 1 and 2 will
% give nan, but Approach 3 will give 0 as the cued reach rate
% count up number of reaches in each window
rts_trial_n=dataset.realDistributions.event_RT_trial1InSeq{1}; % this is the reaction time of trial n
rts_trial_nplus1=dataset.realDistributions.event_RT_trialiInSeq{1}; % this is the reaction time trial n+1 (or n+i, where i was specified when built dataset)
allreaches_trial_n=dataset.realDistributions.rawReaching_event_trial1InSeq{1}; % this is the "PSTH" of reaching for first trial in sequence
allreaches_trial_nplus1=dataset.realDistributions.rawReaching_event_trialiInSeq{1}; % this is the "PSTH" of reaching for last trial in sequence (i.e., trial n+i)

% get fixed windows from average RT at beginning of each session
matchesEventCond_trial_n=dataset.realDistributions.event_isSeq{1}==1;
matchesEventCond_trial_all=dataset.realDistributions.allTrialsSequence_isSeq{1}==1;
matchesEventCond_trial_nplus1=dataset.realDistributions.templateSequence2_end;
% check whether sessids are unique, else make them unique
u=unique(metadata.mouseid);
alreadyUsedSessIDs=[];
sessIDsAreUnique=true;
for i=1:length(u)
    currmouse=u(i);
    if any(ismember(alreadyUsedSessIDs,unique(metadata.sessid(metadata.mouseid==currmouse))))
        sessIDsAreUnique=false;
        break
    else
        alreadyUsedSessIDs=[alreadyUsedSessIDs; unique(metadata.sessid(metadata.mouseid==currmouse))];
    end
end
if sessIDsAreUnique==true
    % use sessids
else
    % make unique sessids
    u=unique(metadata.mouseid);
    j=0;
    for i=1:length(u)
        metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
        j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
    end
end
nth_sessions=metadata.sessid(matchesEventCond_trial_n);
nth_sessions_alltrials=metadata.sessid(matchesEventCond_trial_all);
mouseid=metadata.mouseid(matchesEventCond_trial_n);
% [metadata,fractionThroughSess]=howFarThroughSession(metadata,false,[]);
fractionThroughSess=alltbt.fractionThroughSess_adjusted;
fractionThroughSess=fractionThroughSess(matchesEventCond_trial_n); % how far through session was each trial in sequence

u=unique(nth_sessions_alltrials);
sessidsperrow=u;
n_for_init_cond=nan(length(u),1);
n_for_init_rate=nan(length(u),1);
averageRT_trial_n=nan(length(u),1);
if isempty(initWindows)
    init_fixed_window1=nan(length(u),2);
    init_fixed_window2=nan(length(u),2);
    init_fixed_window3=nan(length(u),2);
else
    init_fixed_window1=initWindows{1};
    init_fixed_window2=initWindows{2};
    init_fixed_window3=initWindows{3};
end
totalTrialsPerSess=nan(length(u),1);
for i=1:length(u)
    totalTrialsPerSess(i)=nansum(nth_sessions==u(i));
end
ratein_fixed_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 1
ratein_fixed_window2=nan(length(u),nanmax(totalTrialsPerSess));
ratein_fixed_window3=nan(length(u),nanmax(totalTrialsPerSess));
reachprobin_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 2
reachprobin_window2=nan(length(u),nanmax(totalTrialsPerSess));
ratein_window3=nan(length(u),nanmax(totalTrialsPerSess));
ratein_acrossSess_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 3
ratein_acrossSess_window2=nan(length(u),nanmax(totalTrialsPerSess));
ratein_acrossSess_window3=nan(length(u),nanmax(totalTrialsPerSess));
trial1_window1=nan(length(u),nanmax(totalTrialsPerSess));
trial1_window2=nan(length(u),nanmax(totalTrialsPerSess));
trial1_window3=nan(length(u),nanmax(totalTrialsPerSess));
fracsThroughSess=nan(length(u),nanmax(totalTrialsPerSess));
for i=1:length(u) % for each session
    n_for_init_cond(i)=ceil((percentOfReachesFromSess_forInitCond/100)*nansum(nth_sessions==u(i))); % for calculating where animal's behavior begins this session
    n_for_init_rate(i)=ceil((percentOfReachesFromSess_forInitRate/100)*nansum(nth_sessions==u(i))); % for calculating where animal's behavior begins this session
    temp=rts_trial_n(nth_sessions==u(i)); % reaction time first trial in sequence
    temp_rawReach=allreaches_trial_nplus1(nth_sessions==u(i),:); % PSTH of reaching last trial in sequence
    temp_rawReach_trial1=allreaches_trial_n(nth_sessions==u(i),:); % PSTH of reaching last trial in sequence
    subfracsForUsedTrials=fractionThroughSess(nth_sessions==u(i)); % how far through session were trials in sequence
    if shuffleTrialOrder==true
        p=randperm(length(temp));
        temp=temp(p);
        temp_rawReach=temp_rawReach(p,:);
    end
    averageRT_trial_n(i)=nanmean(temp(1:n_for_init_cond(i)));

    if isempty(initWindows)
        init_fixed_window1(i,:)=window1(0,averageRT_trial_n(i),epsilon_cue); % this is in terms of reaction time; thus, cue time is 0
        init_fixed_window2(i,:)=window2(0,averageRT_trial_n(i),epsilon_uncue); % this is in terms of reaction time; thus, cue time is 0
        init_fixed_window3(i,:)=window3(0,averageRT_trial_n(i),epsilon_beforecue); % this is in terms of reaction time; thus, cue time is 0
    end
        
    % 3 approaches
    % Approach 1: get reach rate in each window over course of session, where windows are specified
    % wrt animal's behavior at the beginning of the session (fixed windows)
    % Approach 2: get reach rate in each window over course of session, where
    % windows change with changing behavior over the course of session
    % (non-fixed windows)
    % Approach 3: get reach rate in each window over course of session,
    % where windows are fixed across all sessions
    for j=1:size(temp_rawReach,1)
        % APPROACH 1    
        % rows are different sessions, columns are different trials in each
        % session
        ratein_fixed_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach(j,:),cueInd,timeStep); % getReachRate returns number of reaches in time window, per unit of time
        ratein_fixed_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach(j,:),cueInd,timeStep);
        ratein_fixed_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach(j,:),cueInd,timeStep);
        % APPROACH 2
        if isnan(temp(j))
            % didn't reach at all after the cue on trial n, so any reach on
            % trial n+1 is a change
            reachprobin_window1(i,j)=getReachProbability(window1(0,maxTrialLength,epsilon_cue),temp_rawReach(j,:),cueInd,timeStep); % getReachProbability returns 1 or 0 if any reach in time window
            reachprobin_window2(i,j)=getReachProbability([maxTrialLength maxTrialLength],temp_rawReach(j,:),cueInd,timeStep);
            ratein_window3(i,j)=getReachRate(window3(0,[],epsilon_beforecue),temp_rawReach(j,:),cueInd,timeStep);
        else
            reachprobin_window1(i,j)=getReachProbability(window1(0,temp(j),epsilon_cue),temp_rawReach(j,:),cueInd,timeStep);
            reachprobin_window2(i,j)=getReachProbability(window2(0,temp(j),epsilon_uncue),temp_rawReach(j,:),cueInd,timeStep);
            ratein_window3(i,j)=getReachRate(window3(0,temp(j),epsilon_beforecue),temp_rawReach(j,:),cueInd,timeStep);
        end
        % APPROACH 3
        ratein_acrossSess_window1(i,j)=getReachRate(acrossSess_window1,temp_rawReach(j,:),cueInd,timeStep);
        ratein_acrossSess_window2(i,j)=getReachRate(acrossSess_window2,temp_rawReach(j,:),cueInd,timeStep);
        ratein_acrossSess_window3(i,j)=getReachRate(acrossSess_window3,temp_rawReach(j,:),cueInd,timeStep);
        
        
        % Also get reach rate of trial 1
        if useRateMethod==1
            trial1_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
            trial1_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
            trial1_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
        elseif useRateMethod==2
            % not a meaningful thing to do
        elseif useRateMethod==3
            trial1_window1(i,j)=getReachRate(acrossSess_window1,temp_rawReach_trial1(j,:),cueInd,timeStep);
            trial1_window2(i,j)=getReachRate(acrossSess_window2,temp_rawReach_trial1(j,:),cueInd,timeStep);
            trial1_window3(i,j)=getReachRate(acrossSess_window3,temp_rawReach_trial1(j,:),cueInd,timeStep);
        end
        
        
        % save when trial occured
        fracsThroughSess(i,j)=subfracsForUsedTrials(j);  
    end
end
out.init_fixed_window1=init_fixed_window1;
out.init_fixed_window2=init_fixed_window2;
out.init_fixed_window3=init_fixed_window3;
out.fracsThroughSess=fracsThroughSess;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot output
if isempty(ratein_fixed_window1)
    out=[];
    return
end
[moreout,meansForProportionality_x,meansForProportionality_y]=plotAverageAcrossSessions_trialByTrialWithinSession(useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,n_for_init_rate,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3,addSatietyLines);
Y_distance_from_proportionality(useRateMethod,settings,meansForProportionality_x,meansForProportionality_y,moreout,scatterPointSize,n_for_init_rate);
plotTrialByTrialWithinSession_noAveraging(settings.binThisManyTrials,useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3);
plotTrialByTrial_compareFirstToLastTrialInSequence(settings.nBinsForZones,useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3);
f=fieldnames(moreout);
for i=1:length(f)
    out.(f{i})=moreout.(f{i});
end

end

function out=plotTrialByTrialWithinSession_noAveraging(binThisManyTrials,useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3)

trialBins=[1:binThisManyTrials:size(ratein_fixed_window1,2) size(ratein_fixed_window1,2)];
if settings.suppressPlots==false
    figure(); % Approach 1
end
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/(length(trialBins)-1));
approach_alltrials_cued=nan(size(ratein_fixed_window1,1),length(trialBins)-1);
approach_alltrials_uncued=nan(size(ratein_fixed_window1,1),length(trialBins)-1);
approach_trial1_cued=nan(size(ratein_fixed_window1,1),length(trialBins)-1);
approach_trial1_uncued=nan(size(ratein_fixed_window1,1),length(trialBins)-1);
for i=1:length(trialBins)-1 % across trials binned
    takeTheseTrials=trialBins(i):trialBins(i+1);
    if useRateMethod==1
        rateinwindow1=ratein_fixed_window1;
        rateinwindow2=ratein_fixed_window2;
        rateinwindow3=ratein_fixed_window3;
    elseif useRateMethod==2
        rateinwindow1=reachprobin_window1;
        rateinwindow2=reachprobin_window2;
        rateinwindow3=ratein_window3;
    elseif useRateMethod==3
        rateinwindow1=ratein_acrossSess_window1;
        rateinwindow2=ratein_acrossSess_window2;
        rateinwindow3=ratein_acrossSess_window3;
    end
    if all(ismember([2 3],useWindowsForUncued))
        % average uncued reaching in these two windows
        temp_uncued=nanmean([rateinwindow2(:,takeTheseTrials) rateinwindow3(:,takeTheseTrials)],2);
        temp_uncued_trial1=nanmean([trial1_window2(:,takeTheseTrials) trial1_window3(:,takeTheseTrials)],2);
    elseif ismember(2,useWindowsForUncued)
        temp_uncued=nanmean([rateinwindow2(:,takeTheseTrials)],2);
        temp_uncued_trial1=nanmean([trial1_window2(:,takeTheseTrials)],2);
    elseif ismember(3,useWindowsForUncued)
        temp_uncued=nanmean([rateinwindow3(:,takeTheseTrials)],2);
        temp_uncued_trial1=nanmean([trial1_window3(:,takeTheseTrials)],2);
    end
    % rows are different sessions, columns are different trials in each
    % session
    temp_cued=nanmean(rateinwindow1(:,takeTheseTrials),2);
    temp_cued_trial1=nanmean(trial1_window1(:,takeTheseTrials),2);
    a=-0.01;
    b=0.01;
    r=(b-a).*rand(size(temp_cued))+a;
    s=[];
    if i==1 && settings.suppressPlots==false
        trial1start=scatter(temp_uncued+r,temp_cued+r,[],'k'); hold on; % first trial in session, last trial in sequence
    elseif settings.suppressPlots==false
        s=scatter(temp_uncued+r,temp_cued+r,[],cmap(k,:)); hold on;
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
    approach_alltrials_cued(:,i)=temp_cued;
    approach_alltrials_uncued(:,i)=temp_uncued;
    approach_trial1_cued(:,i)=temp_cued_trial1;
    approach_trial1_uncued(:,i)=temp_uncued_trial1;
end
if settings.suppressPlots==false
    xlabel('Uncued reach rate (1/sec)');
    ylabel('Cued reach rate (1/sec)');
    if useRateMethod==1
        title('Approach 1, window 2/3 vs window 1 NO AVERAGING ACROSS SESSIONS');
    elseif useRateMethod==2
        title('Approach 2, window 2 vs window 1 NO AVERAGING ACROSS SESSIONS');
        xlabel('Uncued reach rate (1/sec)');
        ylabel('Probability that reached faster after cue on second trial');
    elseif useRateMethod==3
        title('Approach 3, window 2/3 vs window 1 NO AVERAGING ACROSS SESSIONS');
    end
end
l=[];
if settings.addSessionLines==true
    for i=1:size(approach_alltrials_cued,1)
        l=line(approach_alltrials_uncued(i,:),approach_alltrials_cued(i,:),'Color',[0.5 0.5 0.5],'LineWidth',0.1);
    end
end
out.alltrials_cued=approach_alltrials_cued;
out.alltrials_uncued=approach_alltrials_uncued;
out.trial1_alltrials_uncued=approach_trial1_uncued;
out.trial1_alltrials_cued=approach_trial1_cued;
if settings.suppressPlots==false
    legend([trial1start s l],{'Each trial (end of sequence) in each session','Other trials in session, last trial in sequence','Lines connect trials within same session'});
end

end

function out=plotTrialByTrial_compareFirstToLastTrialInSequence(nBinsForZones,useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3)

if settings.suppressPlots==false
    figure(); % Approach 1
end
if useRateMethod==1
    rateinwindow1=ratein_fixed_window1;
    rateinwindow2=ratein_fixed_window2;
    rateinwindow3=ratein_fixed_window3;
elseif useRateMethod==2
    rateinwindow1=reachprobin_window1;
    rateinwindow2=reachprobin_window2;
    rateinwindow3=ratein_window3;
elseif useRateMethod==3
    rateinwindow1=ratein_acrossSess_window1;
    rateinwindow2=ratein_acrossSess_window2;
    rateinwindow3=ratein_acrossSess_window3;
end
% throw out sessions with not enough trials in session
atleastntrials=20;
rateinwindow1(nansum(~isnan(rateinwindow1),2)<atleastntrials)=nan;
rateinwindow2(nansum(~isnan(rateinwindow2),2)<atleastntrials)=nan;
rateinwindow3(nansum(~isnan(rateinwindow3),2)<atleastntrials)=nan;
if all(ismember([2 3],useWindowsForUncued))
    % average uncued reaching in these two windows
    temp_uncued=nanmean([rateinwindow2 rateinwindow3],2);
    temp_uncued_trial1=nanmean([trial1_window2 trial1_window3],2);
elseif ismember(2,useWindowsForUncued)
    temp_uncued=nanmean([rateinwindow2],2);
    temp_uncued_trial1=nanmean([trial1_window2],2);
elseif ismember(3,useWindowsForUncued)
    temp_uncued=nanmean([rateinwindow3],2);
    temp_uncued_trial1=nanmean([trial1_window3],2);
end
% rows are different sessions, columns are different trials in each
% session
temp_cued=nanmean(rateinwindow1,2);
temp_cued_trial1=nanmean(trial1_window1,2);
bigq=[];
if settings.suppressPlots==false
    bigq=quiver(nanmean(temp_uncued_trial1),nanmean(temp_cued_trial1),nanmean(temp_uncued-temp_uncued_trial1),nanmean(temp_cued-temp_cued_trial1),'Color','k','Linewidth',4); hold on;
end
% remove outliers using vector length
[~,isout]=rmoutliers(sqrt((temp_cued-temp_cued_trial1).^2+(temp_uncued-temp_uncued_trial1).^2),'percentiles',[10 90]);
% [~,isout]=rmoutliers(temp_cued,'percentiles',[5 95]);
% [~,isout2]=rmoutliers(temp_uncued,'percentiles',[5 95]);
% [~,isout3]=rmoutliers(temp_cued_trial1,'percentiles',[5 95]);
% [~,isout4]=rmoutliers(temp_uncued_trial1,'percentiles',[5 95]);
% isout=isout | isout2 | isout3 | isout4;
% isout=zeros(size(temp_cued));
temp_cued=temp_cued(~isout);
temp_uncued=temp_uncued(~isout);
temp_cued_trial1=temp_cued_trial1(~isout);
temp_uncued_trial1=temp_uncued_trial1(~isout);

approach_alltrials_cued=temp_cued;
approach_alltrials_uncued=temp_uncued;
approach_trial1_cued=temp_cued_trial1;
approach_trial1_uncued=temp_uncued_trial1;

rangex=[nanmin([temp_uncued; temp_uncued_trial1]) nanmax([temp_uncued; temp_uncued_trial1])];
binsForZones_x=rangex(1):range(rangex)/(nBinsForZones-1):rangex(2);
rangey=[nanmin([temp_cued; temp_cued_trial1]) nanmax([temp_cued; temp_cued_trial1])];
binsForZones_y=rangey(1):range(rangey)/(nBinsForZones-1):rangey(2);
if rangey(1)-rangey(2)<1e5
    return
end
% find which bins actually occupied
occupiedBins=zeros(length(binsForZones_x)-1,length(binsForZones_y)-1);
testToUse='nanmean([approach_trial1_cued approach_alltrials_cued],2)>=currbin_y(1) & nanmean([approach_trial1_cued approach_alltrials_cued],2)<=currbin_y(2) & nanmean([approach_trial1_uncued approach_alltrials_uncued],2)>=currbin_x(1) & nanmean([approach_trial1_uncued approach_alltrials_uncued],2)<=currbin_x(2)';
for i=1:length(binsForZones_y)-1
    for j=1:length(binsForZones_x)-1
        currbin_x=[binsForZones_x(j) binsForZones_x(j+1)];
        currbin_y=[binsForZones_y(i) binsForZones_y(i+1)];
        occupiedBins(j,i)=any(eval(testToUse));
    end
end
totalBins=nansum(occupiedBins(1:end));
cmap=colormap('cool');
kstep=ceil(size(cmap,1)/totalBins);
k=1;
js_for_colormap=nan(size(occupiedBins));
for i=1:length(binsForZones_y)-1
    for j=1:length(binsForZones_x)-1
        if occupiedBins(j,i)==1
            js_for_colormap(j,i)=k;
            k=k+kstep;
            if k>size(cmap,1)
                k=size(cmap,1);
            end
        end
    end
end
l=[];
la=[];
for i=1:length(binsForZones_y)-1
    for j=1:length(binsForZones_x)-1
        currbin_x=[binsForZones_x(j) binsForZones_x(j+1)];
        currbin_y=[binsForZones_y(i) binsForZones_y(i+1)];
        k=js_for_colormap(j,i);
        whichToUse=eval(testToUse);
        if nansum(whichToUse)>0 && isnan(k)
            error('Problem with calculating occupied bins in plotChangeInReachProbability_fromRTdataset.m');
        end
        if nansum(whichToUse)==0
            continue
        end
        temp_cued=approach_alltrials_cued(whichToUse);
        temp_uncued=approach_alltrials_uncued(whichToUse);
        temp_cued_trial1=approach_trial1_cued(whichToUse);
        temp_uncued_trial1=approach_trial1_uncued(whichToUse);
        if settings.suppressPlots==false
            lastTrialInSeq=scatter(temp_uncued,temp_cued,[],cmap(k,:)); hold on; % last trial in sequence across sessions
            l=line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
                [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:));
            line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:));
            
            firstTrialInSeq=scatter(temp_uncued_trial1,temp_cued_trial1,[],'k'); hold on; % first trial in sequence across sessions
            scatter(nanmean(temp_uncued_trial1),nanmean(temp_cued_trial1),[],cmap(k,:),'filled');
            la=line([nanmean(temp_uncued_trial1)-nanstd(temp_uncued_trial1)./sqrt(nansum(~isnan(temp_uncued_trial1))) nanmean(temp_uncued_trial1)+nanstd(temp_uncued_trial1)./sqrt(nansum(~isnan(temp_uncued_trial1)))],...
                [nanmean(temp_cued_trial1) nanmean(temp_cued_trial1)],'Color','k');
            line([nanmean(temp_uncued_trial1) nanmean(temp_uncued_trial1)],...
                [nanmean(temp_cued_trial1)-nanstd(temp_cued_trial1)./sqrt(nansum(~isnan(temp_cued_trial1))) nanmean(temp_cued_trial1)+nanstd(temp_cued_trial1)./sqrt(nansum(~isnan(temp_cued_trial1)))],'Color','k');
            if nansum(whichToUse)<3 % not enough sessions to feel confident
                continue
            end
            quiver(nanmean(temp_uncued_trial1),nanmean(temp_cued_trial1),nanmean(temp_uncued-temp_uncued_trial1),nanmean(temp_cued-temp_cued_trial1),'Color',cmap(k,:),'Linewidth',2);
        end
    end
end
if settings.suppressPlots==false
    xlabel('Uncued reach rate (1/sec)');
    ylabel('Cued reach rate (1/sec)');
    if useRateMethod==1
        title('Approach 1, window 2/3 vs window 1 trial 1 in sequence vs last trial in sequence');
    elseif useRateMethod==2
        title('Approach 2, window 2 vs window 1 trial 1 in sequence vs last trial in sequence');
        xlabel('Uncued reach rate (1/sec)');
        ylabel('Probability that reached faster after cue on second trial');
    elseif useRateMethod==3
        title('Approach 3, window 2/3 vs window 1 trial 1 in sequence vs last trial in sequence');
    end
end
out.alltrials_cued=approach_alltrials_cued;
out.alltrials_uncued=approach_alltrials_uncued;
out.trial1_alltrials_uncued=approach_trial1_uncued;
out.trial1_alltrials_cued=approach_trial1_cued;
if settings.suppressPlots==false
    line([rangex(1) rangex(2)],[rangex(1) rangex(2)],'Color',[0.5 0.5 0.5]);
    legend([bigq lastTrialInSeq l firstTrialInSeq la],{'Av over all sess','Last trial in sequence across session','Mean and se last trial in seq','First trial in sequence across session','Mean and se first trial in seq'});
end

end

function [out,meansForProportionality_x,meansForProportionality_y]=plotAverageAcrossSessions_trialByTrialWithinSession(useRateMethod,settings,ratein_fixed_window1,ratein_fixed_window2,ratein_fixed_window3,ratein_acrossSess_window1,ratein_acrossSess_window2,ratein_acrossSess_window3,n_for_init_rate,useWindowsForUncued,trial1_window1,trial1_window2,trial1_window3,reachprobin_window1,reachprobin_window2,ratein_window3,addSatietyLines)

meansForProportionality_x=[]; % return uncued reach rates for first few trials in each session
meansForProportionality_y=[]; % return cued reach rates for first few trials in each session
suppressBootstrap=true;
if useRateMethod==1 || useRateMethod==3
    if settings.suppressPlots==false
        figure(); % Approach 1
    end
    cmap=colormap('cool');
    k=1;
    if ~isempty(settings.stopPlottingTrialsAfterN)
        kstep=ceil(size(cmap,1)/settings.stopPlottingTrialsAfterN);
    else
        kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(ratein_fixed_window1,1))));
    end
    approach_cued=nan(1,size(ratein_fixed_window1,2));
    approach_uncued=nan(1,size(ratein_fixed_window1,2));
    approach_alltrials_cued=nan(size(ratein_fixed_window1));
    approach_alltrials_uncued=nan(size(ratein_fixed_window1));
    approach_trial1_cued=nan(size(ratein_fixed_window1));
    approach_trial1_uncued=nan(size(ratein_fixed_window1));
    approach_m=nan(1,size(ratein_fixed_window1,2));
    nIndsForfirstRates=nanmean(n_for_init_rate);
    meansForProportionality_x=nan(size(ratein_fixed_window1,1),floor(nIndsForfirstRates));
    meansForProportionality_y=nan(size(ratein_fixed_window1,1),floor(nIndsForfirstRates));
    currbincued=[];
    currbinuncued=[];
    currbincounter=0;
    for i=1:size(ratein_fixed_window1,2) % across trials
        if useRateMethod==1
            rateinwindow1=ratein_fixed_window1;
            rateinwindow2=ratein_fixed_window2;
            rateinwindow3=ratein_fixed_window3;
        elseif useRateMethod==3
            rateinwindow1=ratein_acrossSess_window1;
            rateinwindow2=ratein_acrossSess_window2;
            rateinwindow3=ratein_acrossSess_window3;
        end
        if all(ismember([2 3],useWindowsForUncued))
            % average uncued reaching in these two windows
            temp_uncued=nanmean([rateinwindow2(:,i) rateinwindow3(:,i)],2);
            temp_uncued_trial1=nanmean([trial1_window2(:,i) trial1_window3(:,i)],2);
        elseif ismember(2,useWindowsForUncued)
            temp_uncued=nanmean([rateinwindow2(:,i)],2);
            temp_uncued_trial1=nanmean([trial1_window2(:,i)],2);
        elseif ismember(3,useWindowsForUncued)
            temp_uncued=nanmean([rateinwindow3(:,i)],2);
            temp_uncued_trial1=nanmean([trial1_window3(:,i)],2);
        end
        % rows are different sessions, columns are different trials in each
        % session
        % so ACROSS ALL SESSIONS, take each trial in session
        temp_cued=rateinwindow1(:,i); % reach rate in trial n+i (last trial) of sequence
        temp_cued_trial1=trial1_window1(:,i); % reach rate in trial n (trial 1) of sequence
        if i==1 && settings.suppressPlots==false
            trial1start=scatter(nanmean(temp_uncued),nanmean(temp_cued),[],'k'); % first trial in SESSION, last trial in sequence
        end
        targets=[];
        if ~isempty(settings.stopPlottingTrialsAfterN)
            if i<=settings.stopPlottingTrialsAfterN
                goAhead=true;
            else
                goAhead=false;
            end
        else
            goAhead=true;
        end
        if settings.suppressPlots==false && goAhead
            if settings.binTrialsForAvAcrossSess==true
                currbincued=[currbincued temp_cued];
                currbinuncued=[currbinuncued temp_uncued];
                currbincounter=currbincounter+1;
                if currbincounter==settings.binThisManyTrials
                    % plot and reset
                    currbincued=nanmean(currbincued,2);
                    currbinuncued=nanmean(currbinuncued,2);
                    currbincounter=0;
                    backup_temp_cued=temp_cued;
                    backup_temp_uncued=temp_uncued;
                    temp_cued=currbincued;
                    temp_uncued=currbinuncued;
                    currbincued=[];
                    currbinuncued=[];
                    line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
                        [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:),'LineWidth',1); % later trials in SESSION, last trial in sequence
                    hold on;
                    if i==1
                        targets=line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                            [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                    else
                        line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                            [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                    end
                    if size(ratein_fixed_window1,1)==1
                        scatter(nanmean(temp_uncued),nanmean(temp_cued),[],cmap(k,:),'filled');
                    end
                    temp_cued=backup_temp_cued;
                    temp_uncued=backup_temp_uncued;
                end
            else
                line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
                    [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:),'LineWidth',1); % later trials in SESSION, last trial in sequence
                hold on;
                if i==1
                    targets=line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                        [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                else
                    line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                        [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                end
            end
        end
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
        approach_cued(i)=nanmean(temp_cued);
        approach_uncued(i)=nanmean(temp_uncued);
        approach_alltrials_cued(:,i)=temp_cued;
        approach_alltrials_uncued(:,i)=temp_uncued;
        approach_trial1_cued(:,i)=temp_cued_trial1;
        approach_trial1_uncued(:,i)=temp_uncued_trial1;
        if i<=nIndsForfirstRates
            meansForProportionality_x(:,i)=temp_uncued;
            meansForProportionality_y(:,i)=temp_cued;
        end
        delta_uncued=(nanmean(temp_uncued))/(nanmean(temp_uncued)+nanmean(temp_cued));
        delta_cued=(nanmean(temp_cued))/(nanmean(temp_cued)+nanmean(temp_uncued));
        m=delta_cued/delta_uncued; % this is just y over x ratio
        approach_m(i)=m; % for each trial in session, approach_m saves y/x (i.e., cued over uncued)
        if addSatietyLines==true && settings.suppressPlots==false
            length_line_segment=nanmean([nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))])/m;
            y_x_proportionality_lines=line([nanmean(temp_uncued)-length_line_segment/2 nanmean(temp_uncued)+length_line_segment/2],[nanmean(temp_cued)-(length_line_segment*m)/2 nanmean(temp_cued)+(length_line_segment*m)/2],'Color','k','LineWidth',0.25);
            hold on;
        else
            y_x_proportionality_lines=[];
        end
    end
    if settings.suppressPlots==false && settings.showFitLine==true
        lastToUse=settings.stopPlottingTrialsAfterN;
        if lastToUse>length(approach_uncued)
            lastToUse=length(approach_uncued);
        end
        X=[ones(length(approach_uncued(1:lastToUse)),1) approach_uncued(1:lastToUse)'];
        b=X\approach_cued(1:lastToUse)';
        plot(approach_uncued(1:lastToUse)',X*b,'Color','k');
%         plot(downSampAv(approach_uncued(1:lastToUse),10),downSampAv(approach_cued(1:lastToUse),10),'Color','k');
    end
    % tmp=cat(3,ratein_fixed_window2,ratein_fixed_window3);
    % C=nansum(tmp,3);
    % tempuncued=nanmean(nanmean(C/2,1),2);
    % tempcued=nanmean(nanmean(ratein_fixed_window1,1),2);
    tempuncued=nanmean(nanmean(meansForProportionality_x,1),2);
    tempcued=nanmean(nanmean(meansForProportionality_y,1),2);
    tempuncued(isnan(tempuncued) | isinf(tempuncued))=0;
    tempcued(isnan(tempcued) | isinf(tempcued))=0;
    % quiver(0,0,tempuncued,tempcued,'Color','k');
    if settings.suppressPlots==false
        trial1line=line([0 approach_uncued(1)],[0 (tempcued/tempuncued)*approach_uncued(1)],'Color','k','LineWidth',3);
        xlabel('Uncued reach rate (1/sec)');
        ylabel('Cued reach rate (1/sec)');
        if useRateMethod==1
            title('Approach 1, window 2/3 vs window 1');
        elseif useRateMethod==3
            title('Approach 3, window 2/3 vs window 1');
        end
    end
    out.cued=approach_cued;
    out.uncued=approach_uncued;
    out.alltrials_cued=approach_alltrials_cued;
    out.alltrials_uncued=approach_alltrials_uncued;
    out.m=approach_m;
    out.trial1_alltrials_uncued=approach_trial1_uncued;
    out.trial1_alltrials_cued=approach_trial1_cued;
    
    nIndsForfirstRates=floor(nanmean(n_for_init_rate));
    if ~isempty(settings.stopPlottingTrialsAfterN)
        upTo=settings.stopPlottingTrialsAfterN;
        if upTo>length(approach_cued)
            upTo=length(approach_cued);
        end
    else
        upTo=length(approach_cued);
    end
    if settings.suppressPlots==false
        quiver(nanmean(approach_uncued(1:nIndsForfirstRates)),nanmean(approach_cued(1:nIndsForfirstRates)),...
            nanmean(approach_uncued(upTo-nIndsForfirstRates:upTo))-nanmean(approach_uncued(1:nIndsForfirstRates)),...
            nanmean(approach_cued(upTo-nIndsForfirstRates:upTo))-nanmean(approach_cued(1:nIndsForfirstRates)),'Color','k');
    end
    quiver_y=nanmean(approach_cued(upTo-nIndsForfirstRates:upTo))-nanmean(approach_cued(1:nIndsForfirstRates));
    quiver_x=nanmean(approach_uncued(upTo-nIndsForfirstRates:upTo))-nanmean(approach_uncued(1:nIndsForfirstRates));
    quiver_m=quiver_y/quiver_x;
    ang=atand(quiver_m);
    if quiver_y<0 && quiver_x<0
        ang=180+abs(ang);
    end
    disp(['vector angle ' num2str(ang)]);
    proportional_m=nanmean(approach_cued(1:nIndsForfirstRates))/nanmean(approach_uncued(1:nIndsForfirstRates));
    ang2=atand(proportional_m);
    ang2=180+abs(ang2);
    disp(['proportional angle ' num2str(ang2)]);
    disp(['diff between angles ' num2str(ang2-ang)]);
    
    % Bootstrap means include all trials in session
    altogether_prob_cued=approach_alltrials_cued(1:end);
    altogether_prob_uncued=approach_alltrials_uncued(1:end);
    takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
    altogether_prob_cued=altogether_prob_cued(takeTrials==1);
    altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
    % Show bootstrapped 95% CI
    takeFracForBootstrap=0.66;
    takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
    nRuns=100;
    bootMeans=nan(2,nRuns);
    for i=1:nRuns
        takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
        sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
        sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
        bootMeans(1,i)=nanmean(sub_prob_uncued);
        bootMeans(2,i)=nanmean(sub_prob_cued);
    end
    s_alltrials=[];
    if settings.suppressPlots==false && suppressBootstrap~=true
        s_alltrials=scatter(bootMeans(1,:),bootMeans(2,:),20,'k','filled');
        s_alltrials.AlphaData = 0.5*ones(1,size(bootMeans,2));
        s_alltrials.MarkerFaceAlpha = 'flat';
        scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'k','filled');
    end
    out.boot_trialn_realmean_x=nanmean(altogether_prob_uncued);
    out.boot_trialn_realmean_y=nanmean(altogether_prob_cued);
    out.boot_trialn=bootMeans;
    % Bootstrap trial 1
    altogether_prob_cued=approach_trial1_cued(1:end);
    altogether_prob_uncued=approach_trial1_uncued(1:end);
    takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
    altogether_prob_cued=altogether_prob_cued(takeTrials==1);
    altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
    % Show bootstrapped 95% CI
    takeFracForBootstrap=0.66;
    takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
    nRuns=100;
    bootMeans=nan(2,nRuns);
    for i=1:nRuns
        takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
        sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
        sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
        bootMeans(1,i)=nanmean(sub_prob_uncued);
        bootMeans(2,i)=nanmean(sub_prob_cued);
    end
    s=[];
    if settings.suppressPlots==false && suppressBootstrap~=true
        s=scatter(bootMeans(1,:),bootMeans(2,:),20,'c','filled');
        s.AlphaData = 0.5*ones(1,size(bootMeans,2));
        s.MarkerFaceAlpha = 'flat';
        scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'c','filled');
    end
    out.boot_trial1_realmean_x=nanmean(altogether_prob_uncued);
    out.boot_trial1_realmean_y=nanmean(altogether_prob_cued);
    out.boot_trial1=bootMeans;
    if settings.suppressPlots==false
        legend([trial1start targets y_x_proportionality_lines trial1line s_alltrials s],{'First trial in session','Each trial in session','X-Y proportionality lines','Y over X first trial in session','Bootstrap all trials','Bootstrap trial 1'});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useRateMethod==2
    if settings.suppressPlots==false
        figure(); % Approach 2, proportionality lines won't make sense
    end
    cmap=colormap('cool');
    k=1;
    kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(reachprobin_window1,1))));
    approach2_cued=nan(1,size(ratein_fixed_window1,2));
    approach2_uncued=nan(1,size(ratein_fixed_window1,2));
    approach2_alltrials_cued=nan(size(ratein_fixed_window1));
    approach2_alltrials_uncued=nan(size(ratein_fixed_window1));
    for i=1:size(reachprobin_window1,2) % across trials
        if all(ismember([2 3],useWindowsForUncued))
            % average uncued reaching in these two windows
            temp_uncued=nanmean([reachprobin_window2(:,i) ratein_window3(:,i)],2);
        elseif ismember(2,useWindowsForUncued)
            temp_uncued=nanmean([reachprobin_window2(:,i)],2);
        elseif ismember(3,useWindowsForUncued)
            temp_uncued=nanmean([ratein_window3(:,i)],2);
        end
        temp_cued=reachprobin_window1(:,i);
        if i==1 && settings.suppressPlots==false
            trial1start=scatter(nanmean(temp_uncued),nanmean(temp_cued),[],'k');
        end
        if settings.suppressPlots==false
            line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
                [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:),'LineWidth',1);
            hold on;
            if i==1
                targets=line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                    [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
            else
                line([nanmean(temp_uncued) nanmean(temp_uncued)],...
                    [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
            end
        end
        approach2_cued(i)=nanmean(temp_cued);
        approach2_uncued(i)=nanmean(temp_uncued);
        approach2_alltrials_cued(:,i)=temp_cued;
        approach2_alltrials_uncued(:,i)=temp_uncued;
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    altogether_prob_cued=approach2_alltrials_cued(1:end);
    altogether_prob_uncued=approach2_alltrials_uncued(1:end);
    takeTrials=~isnan(altogether_prob_cued) & ~isnan(altogether_prob_uncued);
    altogether_prob_cued=altogether_prob_cued(takeTrials==1);
    altogether_prob_uncued=altogether_prob_uncued(takeTrials==1);
    % Show bootstrapped 95% CI
    takeFracForBootstrap=0.66;
    takeIndsForBootstrap=ceil(takeFracForBootstrap*length(altogether_prob_cued));
    nRuns=100;
    bootMeans=nan(2,nRuns);
    for i=1:nRuns
        takeTheseForBoot=randi(length(altogether_prob_cued),1,takeIndsForBootstrap); % with replacement
        sub_prob_cued=altogether_prob_cued(takeTheseForBoot);
        sub_prob_uncued=altogether_prob_uncued(takeTheseForBoot);
        bootMeans(1,i)=nanmean(sub_prob_uncued);
        bootMeans(2,i)=nanmean(sub_prob_cued);
    end
    if settings.suppressPlots==false
        s=scatter(bootMeans(1,:),bootMeans(2,:),20,'k','filled');
        s.AlphaData = 0.5*ones(1,size(bootMeans,2));
        s.MarkerFaceAlpha = 'flat';
        scatter(nanmean(altogether_prob_uncued),nanmean(altogether_prob_cued),50,'k','filled');
        %     tmp=cat(3,reachprobin_window2,ratein_window3);
        %     C=nansum(tmp,3);
        %     tempuncued=nanmean(nanmean(C/2,1),2);
        %     tempuncued(isnan(tempuncued) | isinf(tempuncued))=0;
        %     tempcued=nanmean(nanmean(reachprobin_window1,1),2);
        %     tempcued(isnan(tempcued) | isinf(tempcued))=0;
        %     quiver(0,0,tempuncued,tempcued,'Color','k');
        xlabel('Uncued reach rate (1/sec)');
        %     xlabel('Probability that he reached in uncued window (-2 to -1 sec before cue)');
        % ylabel('Cued reach rate (1/sec)');
        ylabel('Probability that reached faster after cue on second trial');
        title('Approach 2, window 2 vs window 1');
    end
    out.cued=approach2_cued;
    out.uncued=approach2_uncued;
    out.alltrials_cued=approach2_alltrials_cued;
    out.alltrials_uncued=approach2_alltrials_uncued;
    out.m=nan; % doesn't really make sense to return slopes for this approach   
    
    % Plot average and s.e. cued vs. uncued for all trials, all sessions
    temp_cued=approach2_cued(1:end);
    temp_uncued=approach2_uncued(1:end);
    [~,isout]=rmoutliers(temp_cued,'percentiles',[5 95]);
    [~,isout2]=rmoutliers(temp_uncued,'percentiles',[5 95]);
    isout=isout | isout2;
%     isout=zeros(size(temp_cued));
    temp_cued=temp_cued(~isout);
    temp_uncued=temp_uncued(~isout);
    if settings.suppressPlots==false
        l=line([mean(temp_uncued,'omitnan') mean(temp_uncued,'omitnan')],[mean(temp_cued,'omitnan')-std(temp_cued,[],2,'omitnan')./sqrt(sum(~isnan(temp_cued))) ...
                                                                        mean(temp_cued,'omitnan')+std(temp_cued,[],2,'omitnan')./sqrt(sum(~isnan(temp_cued)))],'Color','k','LineWidth',2);
        line([mean(temp_uncued,'omitnan')-std(temp_uncued,[],2,'omitnan')./sqrt(sum(~isnan(temp_uncued))) ...
          mean(temp_uncued,'omitnan')+std(temp_uncued,[],2,'omitnan')./sqrt(sum(~isnan(temp_uncued)))],[mean(temp_cued,'omitnan') mean(temp_cued,'omitnan')],'Color','k','LineWidth',2);                                                           
        legend([trial1start targets s l],{'First trial in session','Each trial in session','Bootstrap of probability faster cued & uncued reach rate','Mean and std err'});
    end
end

end

function Y_distance_from_proportionality(useRateMethod,settings,meansForProportionality_x,meansForProportionality_y,out,scatterPointSize,n_for_init_rate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Y distance from proportionality line -- only reasonable for Approach
% 1 or 3
nIndsForfirstRates=nanmean(n_for_init_rate);
if (useRateMethod==1 || useRateMethod==3) && settings.suppressPlots==false
    cmap=colormap('cool');
    tempuncued=nanmean(nanmean(meansForProportionality_x,2),1);
    tempcued=nanmean(nanmean(meansForProportionality_y,2),1);
    m=tempcued/tempuncued;
    expected_y=m*out.uncued;
    actual_y=out.cued; % average across all sessions, cued reaching rate for each trial in session, i.e., trial 1 to last trial in session
    actual_x=out.uncued; % average across all sessions, uncued reaching rate for each trial in session, i.e., trial 1 to last trial in session
    figure();
    k=1;
    if ~isempty(settings.stopPlottingTrialsAfterN)
        kstep=ceil(size(cmap,1)/length(settings.stopPlottingTrialsAfterN));
    else
        kstep=ceil(size(cmap,1)/length(actual_y));
    end
    suppressUnfilled=true;
    for i=1:length(actual_y)
        if suppressUnfilled==false
            if ~isempty(settings.stopPlottingTrialsAfterN)
                if i>settings.stopPlottingTrialsAfterN
                else
                    s=scatter(actual_x(i),actual_y(i)-expected_y(i),scatterPointSize,cmap(k,:),'LineWidth',0.8);
                end
            else
                s=scatter(actual_x(i),actual_y(i)-expected_y(i),scatterPointSize,cmap(k,:),'LineWidth',0.8);
            end            
        else
            s=[];
        end
        hold on;
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    line([0 nanmax(out.uncued)],[0 -m*nanmax(out.uncued)],'Color',[0.5 0.5 0.5]);
    downSampForFilled=3;
    binned_x=downSampAv(actual_x,downSampForFilled);
    ydiff=actual_y-expected_y;
    binned_ydiff=downSampAv(ydiff,downSampForFilled);
    if ~isempty(settings.stopPlottingTrialsAfterN)
        binnedStopPlotting=floor(settings.stopPlottingTrialsAfterN/downSampForFilled);
        if binnedStopPlotting<1
            binnedStopPlotting=1;
        end
    end
    k=1;
    if ~isempty(settings.stopPlottingTrialsAfterN)
        kstep=ceil(size(cmap,1)/binnedStopPlotting);
    else
        kstep=ceil(size(cmap,1)/length(binned_ydiff));
    end
    for i=1:length(binned_ydiff)
        if ~isempty(settings.stopPlottingTrialsAfterN)
            if i>binnedStopPlotting
            else
                sbin=scatter(binned_x(i),binned_ydiff(i),scatterPointSize,cmap(k,:),'filled','LineWidth',1);
            end
        else
            sbin=scatter(binned_x(i),binned_ydiff(i),scatterPointSize,cmap(k,:),'filled','LineWidth',1);
        end
        hold on;
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    % rmoutliers only for plotting
    [~,isout]=rmoutliers(actual_x);
    [~,isout2]=rmoutliers(ydiff);
    isout=isout | isout2;
    lims=[min([actual_x(~isout) ydiff(~isout)],[],'omitnan') max([actual_x(~isout) ydiff(~isout)],[],'omitnan')];
    if lims(2)==0
        isout=zeros(size(actual_x));
        lims=[min([actual_x(~isout) ydiff(~isout)],[],'omitnan') max([actual_x(~isout) ydiff(~isout)],[],'omitnan')];
    end
    if ~isempty(lims)
        if lims(2)==0
            lims(2)=0.01;
        end
        xlim([0 lims(2)]);
        ylim([lims(1) lims(2)]);
    end
    tempx=actual_x(~isout);
    tempy=ydiff(~isout);
    firstPart=[mean(tempx(1:floor(nIndsForfirstRates)),'omitnan') mean(tempy(1:floor(nIndsForfirstRates)),'omitnan')];
    secondPart=[mean(tempx(end-floor(nIndsForfirstRates):end),'omitnan') mean(tempy(end-floor(nIndsForfirstRates):end),'omitnan')];
    if all(~isnan([firstPart secondPart])) && all(~isinf([firstPart secondPart]))
        q=quiver(firstPart(1),firstPart(2),secondPart(1)-firstPart(1),secondPart(2)-firstPart(2),'Color','k','LineWidth',2);
    else
        q=[];
    end
    daspect([1 1 1]);
    ylabel('Cued reach rate distance from proportionality line');
    xlabel('Real uncued reach rate');
    disp('Slope of quiver');
    disp((secondPart(2)-firstPart(2))/(secondPart(1)-firstPart(1)));
    legend([s sbin q],{'Actual uncued reach rate versus difference between expected and actual cued rate','Down-sampled','First part of session to last part of session'});
end

end

function rr=getReachRate(window,reachData,cueInd,timeStep)

thresh=0.05; % for reaching, should be 1 if reaching, else 0

% window is wrt cue, given that cue time is 0 sec
% cueInd is cue timing wrt reachData vector
% timeStep is for reachData vector

window_inds(1)=cueInd+round(window(1)/timeStep); % convert from real time to indices
window_inds(2)=cueInd+round(window(2)/timeStep);

rr=nansum(reachData(window_inds(1):window_inds(2))>thresh)/(window(2)-window(1)); % divide by total duration of window to get rate
% prevent infinity
if (window(2)-window(1))==0
    rr=nan;
end

end
    
function rr=getReachProbability(window,reachData,cueInd,timeStep)

thresh=0.05; % for reaching, should be 1 if reaching, else 0

% window is wrt cue, given that cue time is 0 sec
% cueInd is cue timing wrt reachData vector
% timeStep is for reachData vector

window_inds(1)=cueInd+round(window(1)/timeStep); % convert from real time to indices
window_inds(2)=cueInd+round(window(2)/timeStep);

rr=any(reachData(window_inds(1):window_inds(2))>thresh);

end    
    
    
    
    
