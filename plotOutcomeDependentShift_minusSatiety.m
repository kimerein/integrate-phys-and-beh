function plotOutcomeDependentShift_minusSatiety(alltbt,metadata,trialTypes)

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.isLongITI_1back=[0; trialTypes.isLongITI(1:end-1)];
temp=trialTypes.optoGroup~=1;
trialTypes.isNotGroup1_1forward=[temp(2:end); 0];
trialTypes.isNotGroup1_1back=[0; temp(1:end-1)];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.optoGroup_1back=[0; trialTypes.optoGroup(1:end-1)];

% Note the opto group trials appear to come in pairs
% opto group 1 then 2 then 3

% trialTypes=reachingPeaksBeforeLateOpto(alltbt,metadata,trialTypes);

trial1Type='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==0 & trialTypes.led==0';
% trial1Type='(trialTypes.after_cue_drop_1forward==1 | trialTypes.after_cue_success_1forward==1) & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0) & trialTypes.led==0';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==0) & trialTypes.isLongITI_1back==1';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==0)';
% nbins=50;
% trial1Type='(trialTypes.after_cue_drop_1forward==1 | trialTypes.after_cue_success_1forward==1) & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==1 & trialTypes.isNotGroup1_1forward==1) & trialTypes.reachingPeaksBeforeLateOpto==1';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==1 & trialTypes.isNotGroup1_1back==1) & trialTypes.isLongITI_1back==1 & trialTypes.reachingPeaksBeforeLateOpto==1';
% trial1Type='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.optoGroup_1forward==2 & trialTypes.reachingPeaksBeforeLateOpto==0 & trialTypes.led==0';
% trial1Type='(trialTypes.after_cue_drop_1forward==1 | trialTypes.after_cue_success_1forward==1) & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==1 & (trialTypes.optoGroup_1forward==3 | trialTypes.optoGroup_1forward==2)) & trialTypes.isLongITI_1forward==1'; % & trialTypes.reachingPeaksBeforeLateOpto==1';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==1)'; %& trialTypes.reachingPeaksBeforeLateOpto==1';
nbins=50;

% trial1Type='(trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1)';
trial2Type='(trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1) & trialTypes.led==0';
% nbins=50;

% trial1Type='(trialTypes.after_cue_no_pellet_1forward==1 | trialTypes.after_cue_miss_1forward==1) & (trialTypes.led_1forward==0)';
% trial2Type='(trialTypes.after_cue_no_pellet_1back==1 | trialTypes.after_cue_miss_1back==1) & (trialTypes.led_1back==0)';
% nbins=50;
% trial1Type='(trialTypes.after_cue_no_pellet_1forward==1 | trialTypes.after_cue_miss_1forward==1) & (trialTypes.led_1forward==1 & trialTypes.isNotGroup1_1forward==1)';
% trial2Type='(trialTypes.after_cue_no_pellet_1back==1 | trialTypes.after_cue_miss_1back==1) & (trialTypes.led_1back==1 & trialTypes.isNotGroup1_1back==1)';
% nbins=50;

% trial1Type='(trialTypes.after_cue_drop_1forward==1 | trialTypes.after_cue_success_1forward==1) & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0)';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==0) & trialTypes.isLongITI_1back==1';
% nbins=50;
% trial1Type='(trialTypes.after_cue_drop_1forward==1 | trialTypes.after_cue_success_1forward==1) & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==1 & trialTypes.isNotGroup1_1forward==1) ';
% trial2Type='(trialTypes.after_cue_drop_1back==1 | trialTypes.after_cue_success_1back==1) & trialTypes.touched_pellet_1back==1 & (trialTypes.led_1back==1 & trialTypes.isNotGroup1_1back==1) & trialTypes.isLongITI_1back==1';
% nbins=50;

saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\RT pairs data sets\random\']; 

% make unique sessids for each mouse
u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

% Use all trials to find initial windows for calculating cued reaching and 
% slopes capturing satiety change
tbt_filter.name='dunno';
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=false; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; % in indices, must pass something in, but is not used
skipCorrected=true;
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % this function builds the RT pairs dataset
[rrs_alltoall,sessidsperrow]=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',false);
% Get m's for removing satiety component
[m_binned,refsToMbinned,cued_binned,uncued_binned]=plotDeviationFromSatiety(rrs_alltoall,nbins);

% Get trial to trial change using previously defined init windows
test.nInSequence=[3]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1=trial1Type;
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
trial2=trial2Type;
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
fakeCueInd=50; % in indices, must pass something in, but is not used
skipCorrected=true;
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % this function builds the RT pairs dataset
rrs_trialtotrial=plotChangeInReachProbability_givenInitWindows(dataset,metadata,alltbt,'cueZone_onVoff',false,rrs_alltoall,sessidsperrow);

% Plot trial to trial change minus satiety component
plotTrialToTrialChangeMinusSatiety(rrs_trialtotrial,m_binned,cued_binned,uncued_binned);