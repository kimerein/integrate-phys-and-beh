% runReactionTimeAnalysis.m

% script for running a frequently used subset of analyses

%% load in data

exptDataDir='Z:\Kim\for_orchestra\combineReachData\O2 output\alltbt02Mar2022145751'; % directory containing experimental data

% only load these fields of alltbt
disp('loading alltbt');
whichFieldsToLoad={'cue','all_reachBatch','cueZone_onVoff','dprimes','isChewing','maybeDrop_reachStarts','maybeDrop_reachStarts_pawOnWheel','isHold','optoOn','optoZone','pawOnWheel','pelletPresent','reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts','pelletmissingreach_reachStarts','reachStarts_pelletPresent','times','times_wrt_trial_start','timesFromSessionStart'};
alltbt=loadStructFieldByField([exptDataDir '\alltbt'],whichFieldsToLoad); % load alltbt
disp('loading out');
trialTypes=loadStructFieldByField([exptDataDir '\out']); % load out
disp('loading metadata');
metadata=loadStructFieldByField([exptDataDir '\metadata']); % load metadata
a=load([exptDataDir '\reachExptAnalysis_settings.mat']); % load reach expt analysis settings 
reachExptSettings=a.settings;

% Optional
% Back-up full, unfiltered alltbt in workspace
backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

%% choose additional settings for reaction time analysis

settings=RTanalysis_settings('display settings','clear');

% find trials with long ITIs
trialTypes=getLongITIs(alltbt,trialTypes,settings);

trialTypes=getTimingOfOpto(alltbt,'optoOn',trialTypes,settings.multipleOptoTimes);

trialTypes.mouseid=metadata.mouseid;

% Optional: Fix sessids to match nth_sessions
u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

% Optional: dprimes for each mouse, each session
[alltbt,trialTypes,metadata]=get_dprime_per_mouse(alltbt,trialTypes,metadata);

% Optional
% Back-up full, unfiltered alltbt in workspace
backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

%% perform any filtering on alltbt
% for example, filter by d-prime

temp=datestr(datetime('now'));
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
saveDir=['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\RT pairs data sets\' temp]; % where to save details of alltbt filtering and RT pairs data set

% filter settings
tbt_filter.sortField='dprimes';
tbt_filter.range_values=[-100 100];
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;

% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

%% build relevant data sets

% settings for paired RT data set
% test.nInSequence=[2 3 4 7 11]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% test.nInSequence=[2 3 4 5 6 7 8 9 10 11]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% requirement for first trial in pair
% trial1='any(alltbt.all_reachBatch>0.5,2) & trialTypes.touched_pellet==1 & trialTypes.led==0 & trialTypes.paw_during_wheel==0'; % e.g., mouse reached, touched pellet, no LED, no paw_during_wheel
% trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==0 & trialTypes.consumed_pellet==0 & (trialTypes.led==0)';
% trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1 & (trialTypes.led==1) & trialTypes.isLongITI==1'; % e.g., mouse reached, touched pellet, no LED
% trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1 & (trialTypes.led==0) & trialTypes.isLongITI==1'; % e.g., mouse reached, touched pellet, no LED
% trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1 & (trialTypes.led==1)'; % e.g., mouse reached, touched pellet
% trial1='trialTypes.isLongITI==1';
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
% trial1='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup~=1';
% trial1='trialTypes.touch_in_cued_window_1forward==1';
% trial1='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & (trialTypes.led==0)'; 
% trial1='trialTypes.led==0';
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
% generally, second trial in pair should take all trial types
% trial2='(trialTypes.optoGroup~=1) & trialTypes.led==1';
% trial2='trialTypes.led==1';
% trial2='trialTypes.optoGroup~=1';
trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
% trial2='any(alltbt.all_reachBatch(:,94:end)>0.5,2)';
% trial2='trialTypes.touch_in_cued_window==1';
% trial2='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1 & (trialTypes.led==0)';
% trial2='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==0 & trialTypes.consumed_pellet==0 & (trialTypes.led==0)';
% trial2='any(alltbt.all_reachBatch(:,94:end)>0.5,2) & trialTypes.touched_pellet==1';
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];

saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);

% save settings for paired RT data set
save([saveDir2 '\test_settings.mat'],'test');

% build paired RT data set
fakeCueInd=50; % in indices
skipCorrected=true;
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); % actually, this function just builds the RT pairs dataset
% function also outputs PCA-based model, but I'm going to deprecate that
% analysis approach because may be confusing to some people
%
% RT pair data sets saved to saveDir in buildReachingRTModel.m

%% fit model terms

% Model Assumptions:
% 
% 1. Model instantaneous reaction time distribution as a gamma distribution
% with shape parameter a_instant and scale parameter b_instant
% 
% 2. Markov assumption, i.e., pdf(reaction times at time t+1 | RT at time t, RT at time t-1, etc.)=
% pdf(reaction times at time t+1 | RT at time t)
%
% 3. Then model the cued component of reaction time 
% distribution across time (e.g., sum of a_instant over the 
% course of a session) as itself a gamma-distributed random variable 
% with the same scale parameter b_instant and a shape parameter a_sum 
% that is the sum of the shape parameters for the constituent summands
%
% 4. Model the non-cued component of reaction time distribution as a
% uniform Poisson-distributed random variable with rate lambda
%
% Effects of behavior event on rate of non-cued reaching does not vary
% with time / over the course of learning ... ?
%
% 5. Interaction of non-cued reaching rate and RPE effects:
%
% What about trials with 2 wheel turns in between?

% load reference data set

% putTogetherTermsForRTModel.m

settings.trialDuration=9.5;
settings.n_sems=1;
settings.baselineWindow=[1 2];
settings.cueName='cueZone_onVoff';
settings.useAllDatasetTrialsAsRef=false;
settings.n_trials_away=1;
settings.histo_nbins=[-9:0.1:9];
settings.bins=0:0.1:9;
settings.nTrialsAheadForRT1=-1;
settings.trialStartCountAt1=2;
settings.useReachType='all_reachBatch';


%% session-by-session change in RTs

plotWithinSessionRTshift(dataset.realDistributions,alltbt,trialTypes,metadata,settings,[0 3]);

%% plot data sets

returnThis=plotBehaviorEventFx(dataset,alltbt,ref)

% Dim 2 conversion is 0.7885 sec in dim 2 is 1 sec real
% Dim 1 conversion is 0.611 sec in dim 1 is 1 sec real

plot_rawReaching=false;
plot_rawReaching_cdf=false;
plot_rt_pdf=false;
plot_rt_cdf=false;
plot_delta_rt_pdf=false;
plot_delta_rt_pdf_2D=false;
plot_delta_rt_cdf_2D=false;
plot_delta_rt_cdf=true;
plot_delta_rt_asFunc_rt=false;
plot_dim1_delta_asFunc_rt=false;
plot_dim2_delta_asFunc_rt=false;
plot_3D_dim1_dim2_asFunc_rt=false;
plot_delta_rt_asFunc_rt_removeMeanRegression=false;
plot_earth_mover_wander=false;

% histo_nbins=200; % number of bins for reaction time histogram
histo_nbins=[-4*12.4245:0.2510:4*12.4245];
% histo_nbins=[-4*12.4245-2*0.2510:4*0.2510:4*12.4245];
% histo_nbins=[-4*12.4245-0.75*0.2510:1.5*0.2510:4*12.4245];
backup_histo_nbins=histo_nbins;
scatterJitter=0.5;
% alpha=0.01;
alpha=0.1;
%spotSize=1;
spotSize=10;
nBinsFor2Dhist=100;






%% plot

% plotBehaviorEventFx.m

% fitRTMtermToOtherReaches.m

% fitLocalCorrelationStructureOfReaches.m

% fitReachingRTModel.m

% pelletTouchAlignedReaching.m

