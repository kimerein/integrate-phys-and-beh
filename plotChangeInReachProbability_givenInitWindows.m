function out=plotChangeInReachProbability_givenInitWindows(dataset,metadata,alltbt,cueName,shuffleTrialOrder,initWindows,sessidPerRow)

init_fixed_window1=initWindows.init_fixed_window1;
init_fixed_window2=initWindows.init_fixed_window2;
init_fixed_window3=initWindows.init_fixed_window3;

epsilon=0.25; % in seconds
percentOfReachesFromSess_forInitCond=20; % use this fraction of reaches from beginning of session to get initial conditions
maxTrialLength=9; % in sec, wrt cue
minTrialLength=-2; % wrt cue, in sec
% acrossSess_window1=[0 1.5];
% acrossSess_window2=[1.5 maxTrialLength];
% acrossSess_window3=[minTrialLength -0.5];
acrossSess_window1=[0 1.5];
acrossSess_window2=[minTrialLength -0.5];
acrossSess_window3=[minTrialLength -0.5];
addSatietyLines=false;

[metadata,fractionThroughSess]=howFarThroughSession(metadata);

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
matchesEventCond_trial_n=dataset.realDistributions.event_isSeq{1}==1;
matchesEventCond_trial_nplus1=dataset.realDistributions.templateSequence2_end;

% metadata needs to have unique sessids
nth_sessions=metadata.sessid(matchesEventCond_trial_n);
mouseid=metadata.mouseid(matchesEventCond_trial_n);
fractionThroughSess=fractionThroughSess(matchesEventCond_trial_n);

u=unique(sessidPerRow);
n_for_init_cond=nan(length(u),1);
averageRT_trial_n=nan(length(u),1);
totalTrialsPerSess=nan(length(u),1);
for i=1:length(u)
    totalTrialsPerSess(i)=nansum(nth_sessions==u(i));
end
changein_fixed_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 1
changein_fixed_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_fixed_window3=nan(length(u),nanmax(totalTrialsPerSess));
trial1_fixed_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 1
trial1_fixed_window2=nan(length(u),nanmax(totalTrialsPerSess));
trial1_fixed_window3=nan(length(u),nanmax(totalTrialsPerSess));
changein_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 2
changein_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_window3=nan(length(u),nanmax(totalTrialsPerSess));
changein_acrossSess_window1=nan(length(u),nanmax(totalTrialsPerSess)); % Approach 3
changein_acrossSess_window2=nan(length(u),nanmax(totalTrialsPerSess));
changein_acrossSess_window3=nan(length(u),nanmax(totalTrialsPerSess));
if size(init_fixed_window1,1)~=length(u)
    error('Problem: init windows passed in do not match number of sessions');
end
fracsThroughSess=nan(length(u),nanmax(totalTrialsPerSess));
for i=1:length(u)
    n_for_init_cond(i)=ceil((percentOfReachesFromSess_forInitCond/100)*nansum(nth_sessions==u(i)));
    temp=rts_trial_n(nth_sessions==u(i));
    temp_rawReach=allreaches_trial_nplus1(nth_sessions==u(i),:);
    temp_rawReach_trial1=allreaches_trial_n(nth_sessions==u(i),:);
    subfracsForUsedTrials=fractionThroughSess(nth_sessions==u(i));
    if shuffleTrialOrder==true
        p=randperm(length(temp));
        temp=temp(p);
        temp_rawReach=temp_rawReach(p,:);
    end
    averageRT_trial_n(i)=nanmean(temp(1:n_for_init_cond(i)));
    
%     if isnan(averageRT_trial_n(i)) % there were no reaches after the cue ... no opportunity to learn
%         continue
%     end
    
%     init_fixed_window1(i,:)=window1(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
%     init_fixed_window2(i,:)=window2(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
%     init_fixed_window3(i,:)=window3(0,averageRT_trial_n(i),epsilon); % this is in terms of reaction time; thus, cue time is 0
    
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
%         changein_fixed_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach(j,:),cueInd,timeStep);
%         changein_fixed_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach(j,:),cueInd,timeStep);
%         changein_fixed_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach(j,:),cueInd,timeStep);
        
        % Also get reach rate of trial 1
        trial1_fixed_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
        trial1_fixed_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
        trial1_fixed_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach_trial1(j,:),cueInd,timeStep);
%         trial1_fixed_window1(i,j)=getReachRate(acrossSess_window1,temp_rawReach_trial1(j,:),cueInd,timeStep);
%         trial1_fixed_window2(i,j)=getReachRate(acrossSess_window2,temp_rawReach_trial1(j,:),cueInd,timeStep);
%         trial1_fixed_window3(i,j)=getReachRate(acrossSess_window3,temp_rawReach_trial1(j,:),cueInd,timeStep);
        
        % save when trial occured
        fracsThroughSess(i,j)=subfracsForUsedTrials(j);  
        
        % APPROACH 2
%         changein_window1(i,j)=getReachRate(window1(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
%         changein_window2(i,j)=getReachRate(window2(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
%         changein_window3(i,j)=getReachRate(window3(0,temp(j),epsilon),temp_rawReach(j,:),cueInd,timeStep);
%         % APPROACH 3
        changein_acrossSess_window1(i,j)=getReachRate(init_fixed_window1(i,:),temp_rawReach(j,:),cueInd,timeStep);
        changein_acrossSess_window2(i,j)=getReachRate(init_fixed_window2(i,:),temp_rawReach(j,:),cueInd,timeStep);
        changein_acrossSess_window3(i,j)=getReachRate(init_fixed_window3(i,:),temp_rawReach(j,:),cueInd,timeStep);
    end
end
out.init_fixed_window1=init_fixed_window1;
out.init_fixed_window2=init_fixed_window2;
out.init_fixed_window3=init_fixed_window3;
out.fracsThroughSess=fracsThroughSess;
 
% Plot output
figure(); % Approach 1
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(changein_fixed_window1,1))));
out.cued=nan(1,size(changein_fixed_window1,2));
out.uncued=nan(1,size(changein_fixed_window1,2));
out.alltrials_cued=nan(size(changein_fixed_window1));
out.alltrials_uncued=nan(size(changein_fixed_window1));
out.trial1_cued=nan(1,size(trial1_fixed_window1,2));
out.trial1_uncued=nan(1,size(trial1_fixed_window1,2));
out.trial1_alltrials_cued=nan(size(trial1_fixed_window1));
out.trial1_alltrials_uncued=nan(size(trial1_fixed_window1));
if addSatietyLines==true
    out.m=nan(1,size(changein_fixed_window1,2));
end
for i=1:size(changein_fixed_window1,2) % across trials
    if nanmean(changein_fixed_window1(:,i))==0 || nanmean([changein_fixed_window2(:,i); changein_fixed_window3(:,i)])==0
        continue
    end
    temp_uncued=nanmean([changein_acrossSess_window2(:,i) changein_acrossSess_window3(:,i)],2);
    temp_cued=changein_acrossSess_window1(:,i);
    temp_trial1_uncued=nanmean([trial1_fixed_window2(:,i) trial1_fixed_window3(:,i)],2);
    temp_trial1_cued=trial1_fixed_window1(:,i);
    if i==1
        scatter(nanmean(temp_uncued),nanmean(temp_cued),[],'k');
    end
    line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
         [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:),'LineWidth',1);
    hold on;
    line([nanmean(temp_uncued) nanmean(temp_uncued)],...
         [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
    scatter(nanmean(temp_trial1_uncued),nanmean(temp_trial1_cued),[],'k','filled');
    line([nanmean(temp_trial1_uncued) nanmean(temp_uncued)],[nanmean(temp_trial1_cued) nanmean(temp_cued)],'Color','k','LineWidth',0.3);
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
    out.cued(i)=nanmean(temp_cued);
    out.uncued(i)=nanmean(temp_uncued);
    out.alltrials_cued(:,i)=temp_cued;
    out.alltrials_uncued(:,i)=temp_uncued;
    out.trial1_cued(i)=nanmean(temp_trial1_cued);
    out.trial1_uncued(i)=nanmean(temp_trial1_uncued);
    out.trial1_alltrials_cued(:,i)=temp_trial1_cued;
    out.trial1_alltrials_uncued(:,i)=temp_trial1_uncued;
    
    if addSatietyLines==true
        delta_uncued=(nanmean(temp_uncued)^2)/(nanmean(temp_uncued)+nanmean(temp_cued));
        delta_cued=(nanmean(temp_cued)^2)/(nanmean(temp_cued)+nanmean(temp_uncued));
        m=delta_cued/delta_uncued;
        out.m(i)=m;
        length_line_segment=nanmean([nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))])/m;
        line([nanmean(temp_uncued)-length_line_segment/2 nanmean(temp_uncued)+length_line_segment/2],[nanmean(temp_cued)-(length_line_segment*m)/2 nanmean(temp_cued)+(length_line_segment*m)/2],'Color','k','LineWidth',0.25);
    end
end
tmp=cat(3,changein_fixed_window2,changein_fixed_window3); 
C=nansum(tmp,3);
tempuncued=nanmean(nanmean(C/2,1),2);
tempuncued(isnan(tempuncued) | isinf(tempuncued))=0;
tempcued=nanmean(nanmean(changein_fixed_window1,1),2);
tempcued(isnan(tempcued) | isinf(tempcued))=0;
quiver(0,0,tempuncued,tempcued,'Color','k');
xlabel('Uncued reach rate (1/sec)');
ylabel('Cued reach rate (1/sec)');
title('Approach 1, window 2 vs window 1');

temp=out.alltrials_cued;
f=find(isnan(temp));
f2=find(isnan(fracsThroughSess));
if any(~ismember(f2,f))
    error('Fracs through session not properly assigned');
end

% figure(); % Approach 2
% cmap=colormap('cool');
% k=1;
% kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(changein_window1,1))));
% for i=1:size(changein_window1,2) % across trials
%     if nanmean(changein_window1(:,i))==0 || nanmean([changein_window2(:,i); changein_window3(:,i)])==0
%         continue
%     end
%     temp_uncued=nanmean([changein_window2(:,i) changein_window3(:,i)],2);
%     temp_cued=changein_window1(:,i);
%     if i==1
%         scatter(nanmean(temp_uncued),nanmean(temp_cued),[],'k');
%     end
%     line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued))) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)))],...
%          [nanmean(temp_cued) nanmean(temp_cued)],'Color',cmap(k,:),'LineWidth',1);
%     hold on;
%     line([nanmean(temp_uncued) nanmean(temp_uncued)],...
%          [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
%     k=k+kstep;
%     if k>size(cmap,1)
%         k=size(cmap,1);
%     end
%     out.approach2_alltrials_cued(:,i)=temp_cued;
%     out.approach2_alltrials_uncued(:,i)=temp_uncued;
% end
% tmp=cat(3,changein_window2,changein_window3); 
% C=nansum(tmp,3);
% tempuncued=nanmean(nanmean(C/2,1),2);
% tempuncued(isnan(tempuncued) | isinf(tempuncued))=0;
% tempcued=nanmean(nanmean(changein_window1,1),2);
% tempcued(isnan(tempcued) | isinf(tempcued))=0;
% quiver(0,0,tempuncued,tempcued,'Color','k');
% xlabel('Uncued reach rate (1/sec)');
% ylabel('Cued reach rate (1/sec)');
% title('Approach 2, window 2 vs window 1');

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