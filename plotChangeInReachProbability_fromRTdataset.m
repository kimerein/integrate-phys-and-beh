function plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,cueName)

epsilon=0.25; % in seconds
percentOfReachesFromSess_forInitCond=20; % use this fraction of reaches from beginning of session to get initial conditions
maxTrialLength=9; % in sec, wrt cue
minTrialLength=-2; % wrt cue, in sec
acrossSess_window1=[0 1.5];
acrossSess_window2=[1.5 maxTrialLength];
acrossSess_window3=[minTrialLength -0.5];

% find cue ind
if size(alltbt.(cueName),2)~=size(dataset.realDistributions.rawReaching_event_trial1InSeq{1},2)
    error('dataset size does not fit alltbt size');
end
[~,cueInd]=nanmax(nanmean(alltbt.(cueName),1));
% get timestep
timeStep=mode(diff(nanmean(alltbt.times_wrt_trial_start)));

% three windows
% if t_n is time of reach on trial n
% t_cue is time of cue
% epsilon is size of time window including cue
% then these windows apply to reach time on trial n+1
% if t_n > t_cue -- if using reaction times for t_n, this must be true
% window1 = @(t_cue,t_n,epsilon) [t_cue t_n+epsilon];                                 % cued window 
% window2 = @(t_cue,t_n,epsilon) [t_n+epsilon ((t_n+epsilon)-t_cue)+t_n+epsilon];     % after cued window
% window3 = @(t_cue,t_n,epsilon) [t_cue-((t_n+epsilon)-t_cue) t_cue];                 % before cue window     
window1 = @(t_cue,t_n,epsilon) [t_cue nanmin([t_n+epsilon maxTrialLength])];                             % cued window 
window2 = @(t_cue,t_n,epsilon) [nanmin([t_n+epsilon maxTrialLength]) maxTrialLength];                    % after cued window
window3 = @(t_cue,t_n,epsilon) [minTrialLength t_cue-epsilon];                  % before cue window     
% then need to divide number of reaches in window by duration of window to
% get rate

% if RT is nan, might have reached before cue -- will need to check this
% count up number of reaches in each window
rts_trial_n=dataset.realDistributions.event_RT_trial1InSeq{1};
rts_trial_nplus1=dataset.realDistributions.event_RT_trialiInSeq{1};
allreaches_trial_n=dataset.realDistributions.rawReaching_event_trial1InSeq{1};
allreaches_trial_nplus1=dataset.realDistributions.rawReaching_event_trialiInSeq{1};

% get fixed windows from average RT at beginning of each session
matchesEventCond_trial_n=dataset.realDistributions.templateSequence2_cond;
matchesEventCond_trial_nplus1=dataset.realDistributions.templateSequence2_end;
nth_sessions=metadata.nth_session(matchesEventCond_trial_n);
u=unique(nth_sessions);
n_for_init_cond=nan(length(u),1);
averageRT_trial_n=nan(length(u),1);
init_fixed_window1=nan(length(u),2);
init_fixed_window2=nan(length(u),2);
init_fixed_window3=nan(length(u),2);
totalTrialsPerSess=nan(length(u),1);
for i=1:length(u)
    totalTrialsPerSess(i)=nansum(nth_sessions==u(i));
end
changein_fixed_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 1
changein_fixed_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_fixed_window3=nan(length(u),nanmax(totalTrialsPerSess));
changein_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 2
changein_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_window3=nan(length(u),nanmax(totalTrialsPerSess));
changein_acrossSess_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 3
changein_acrossSess_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_acrossSess_window3=nan(length(u),nanmax(totalTrialsPerSess));
for i=1:length(u)
    n_for_init_cond(i)=ceil((percentOfReachesFromSess_forInitCond/100)*nansum(nth_sessions==u(i)));
    temp=rts_trial_n(nth_sessions==u(i));
    temp_rawReach=allreaches_trial_nplus1(nth_sessions==u(i),:);
    averageRT_trial_n(i)=nanmean(temp(1:n_for_init_cond(i)));
    if isnan(averageRT_trial_n(i)) % there were no reaches after the cue ... no opportunity to learn
        continue
    end
    init_fixed_window1(i,:)=window1(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
    init_fixed_window2(i,:)=window2(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
    init_fixed_window3(i,:)=window3(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
    
    % two approaches
    % Approach 1: get reach rate in each window over course of session, where windows are specified
    % wrt animal's behavior at the beginning of the session (fixed windows)
    % Approach 2: get reach rate in each window over course of session, where
    % windows change with changing behavior over the course of session
    % (non-fixed windows)
    % Approach 3: get reach rate in each window over course of session,
    % where windows are fixed across all sessions
    for j=1:size(temp_rawReach,1)
        % APPROACH 1    
        changein_fixed_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach(j,:),cueInd,timeStep);
        changein_fixed_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach(j,:),cueInd,timeStep);
        changein_fixed_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach(j,:),cueInd,timeStep);
        % APPROACH 2
        changein_window1(i,j)=getReachRate(window1(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
        changein_window2(i,j)=getReachRate(window2(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
        changein_window3(i,j)=getReachRate(window3(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
        % APPROACH 3
        changein_acrossSess_window1(i,j)=getReachRate(acrossSess_window1,temp_rawReach(j,:),cueInd,timeStep);
        changein_acrossSess_window2(i,j)=getReachRate(acrossSess_window2,temp_rawReach(j,:),cueInd,timeStep);
        changein_acrossSess_window3(i,j)=getReachRate(acrossSess_window3,temp_rawReach(j,:),cueInd,timeStep);
    end
end

% Plot output
figure(); % Approach 1
k=1;
currColors={'k','r','c','b','m','y','g'};
for i=1:length(u)
    currColor=currColors{k};
    k=k+1;
    if k>length(currColors)
        k=1;
    end
    for j=2:size(changein_fixed_window1,2)
        quiver(changein_fixed_window2(i,j-1),changein_fixed_window1(i,j-1),changein_fixed_window2(i,j),changein_fixed_window1(i,j),currColor);
        hold on;
    end
end
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');
title('Approach 1, window 2 vs window 1');

end

function rr=getReachRate(window,reachData,cueInd,timeStep)

thresh=0.05; % for reaching, should be 1 if reaching, else 0

% window is wrt cue, given that cue time is 0 sec
% cueInd is cue timing wrt reachData vector
% timeStep is for reachData vector

window_inds(1)=cueInd+round(window(1)/timeStep); % convert from real time to indices
window_inds(2)=cueInd+round(window(2)/timeStep);

rr=nansum(reachData(window_inds(1):window_inds(2))>thresh)/(window(2)-window(1)); % divide by total duration of window to get rate

end
    
    
    
    
    
    
