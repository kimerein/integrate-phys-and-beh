% script_for_reaching_rate_analysis.m
% function script_for_reaching_rate_analysis()

% script for running a frequently used subset of analyses

% DON'T FORGET: FEB_3 (mouse_id 4), FEB_4 (mouse_id 5) AND MITCH_NONE (mouse_id 19) WERE CONTROLS

%% load in data

exptDataDir='Z:\MICROSCOPE\Kim\for_orchestra\combineReachData\O2 output\alltbt01Aug2023121014\'; % directory containing experimental data
behaviorLogDir='C:\Users\sabatini\Downloads\Combo Behavior Log - Slimmed down w old mice added.csv'; % directory containing behavior log, download from Google spreadsheet as .tsv, change extension to .csv
mouseDBdir='Z:\MICROSCOPE\Kim\for_orchestra\combineReachData\O2 output\alltbt01Aug2023121014\mouse_database.mat'; % directory containing mouse database, constructed during prepToCombineReachData_short.m

if ismac==true
    sprtr='/';
else
    sprtr='\';
end

% only load these fields of alltbt
disp('loading alltbt');
whichFieldsToLoad={'cue','all_reachBatch','cueZone_onVoff','dprimes','isChewing','isHold','optoOn','optoZone','pawOnWheel','pelletPresent','reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts','pelletmissingreach_reachStarts','reachStarts_pelletPresent','times','timesFromSessionStart','movie_distractor'};
alltbt=loadStructFieldByField([exptDataDir sprtr 'alltbt'],whichFieldsToLoad); % load alltbt
disp('loading out');
trialTypes=loadStructFieldByField([exptDataDir sprtr 'out']); % load out
disp('loading metadata');
metadata=loadStructFieldByField([exptDataDir sprtr 'metadata']); % load metadata
a=load([exptDataDir sprtr 'reachExptAnalysis_settings.mat']); % load reach expt analysis settings 
reachExptSettings=a.settings;

% Set relevant fields to single
alltbt=setToSingle(alltbt,{'all_reachBatch','cueZone_onVoff','isChewing','isHold','optoOn','pawOnWheel','pelletPresent','pelletmissingreach_reachStarts','reachBatch_all_pawOnWheel','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','reachBatch_success_reachStarts','reachStarts','reachStarts_pelletPresent'},0.1); % last argument is low threshold for conversion, below this will be 0, above this will be 1

% Use behavior log table to fix nth_session, where possible
a=load(mouseDBdir); mouse_database=a.mouse_database;
metadata=getNthSession(behaviorLogDir,mouse_database,metadata,true,true); % last two args are alsoFixOptoOnHere, then excludeTrainingRig

% Optional: get day 1 for learning curves
trialTypes.mouseid=metadata.mouseid;
[~,~,~,isreachout_permouse,permouse_mouseid]=get_dprime_per_mouse(alltbt,trialTypes,metadata,false,settingsForDprimes(alltbt,'cueZone_onVoff',false)); % last arg is filler, dprimes will be recalculated later
[day1,metadata]=defineDay1(alltbt,trialTypes,metadata,isreachout_permouse,permouse_mouseid);
alltbt.sess_wrt_day1=metadata.sess_wrt_day1; trialTypes.sess_wrt_day1=metadata.sess_wrt_day1;

% Optional
% Back-up full, unfiltered alltbt in workspace
backup.alltbt=alltbt;
backup.trialTypes=trialTypes;
backup.metadata=metadata;

% Optional: correct any LED trials for blinded control mice
% [alltbt,metadata,trialTypes]=turnOffLED(alltbt,metadata,trialTypes,[4 5 19]);

% Optional: discard preemptive
[alltbt,trialTypes,metadata]=discardPreemptive(alltbt,trialTypes,metadata);

% fix weird bug where reach batch sometimes get stuck at 1 (in less than 0.1% of trials), possibly an
% interp problem somewhere?? not sure
alltbt=fixReachesStuckAtOne(alltbt);

% fix weird opto bug, where opto output from Arduino was flickering
% alltbt.optoOn(nansum(alltbt.optoOn,2)>100,:)=0;

%% choose additional settings for reaction time analysis

settings=RTanalysis_settings('display settings','clear');


% find trials with long ITIs
trialTypes=getLongITIs(alltbt,trialTypes,settings);


trialTypes=getTimingOfOpto(alltbt,'optoOn',trialTypes,settings.multipleOptoTimes);
if ~isfield(trialTypes,'optoGroup')
    trialTypes.optoGroup=zeros(size(trialTypes.led));
end

% Optional: Fix sessids to match nth_sessions
u=unique(metadata.mouseid);
j=0;
for i=1:length(u)
    metadata.sessid(metadata.mouseid==u(i))=metadata.nth_session(metadata.mouseid==u(i))+j;
    j=j+nanmax(metadata.sessid(metadata.mouseid==u(i)));
end

% Optional: discard no pellet reaches within certain time after success
% bcz may represent chewing arm movements
% this does not really affect results, only needed w small subset of mice
% alltbt=ignoreReachesAfterSuccess(alltbt,metadata,9); % last arg is seconds after success, i.e., how long to chew pellet

% get miss or no pellet reach
alltbt.missOrNoPellet=alltbt.pelletmissingreach_reachStarts+alltbt.reachBatch_miss_reachStarts;
alltbt.missOrNoPellet=single(alltbt.missOrNoPellet>0.5);
% get any fail
alltbt.anyFail=alltbt.pelletmissingreach_reachStarts+alltbt.reachBatch_miss_reachStarts+alltbt.reachBatch_drop_reachStarts;
alltbt.anyFail=single(alltbt.anyFail>0.5);

% Optional: Exclude trials with paw on wheel
% [alltbt,metadata,trialTypes]=excludePawOnWheel(alltbt,metadata,trialTypes,'cueZone_onVoff');

% Optional: Exclude sessions where mouse was cheating
[alltbt,metadata,trialTypes]=excludePreemptiveSess(alltbt,metadata,trialTypes,[3.255-2 3.255-1],[3.255-1 3.255-0],2);

% Optional: dprimes for each mouse, each session
settingsDp=settingsForDprimes(alltbt,'cueZone_onVoff',true); % Check settings in settingsForDprimes
[alltbt,trialTypes,metadata]=get_dprime_per_mouse(alltbt,trialTypes,metadata,false,settingsDp); % last arg is whether to get rates instead
alltbt.dprimes(isinf(alltbt.dprimes))=3; 
% Get cued vs uncued reach rates
settingsRR=settingsForReachRates(alltbt,'cueZone_onVoff',false);
[~,~,metadata]=get_dprime_per_mouse(alltbt,trialTypes,metadata,true,settingsRR); % last arg is whether to get rates instead

% Get dprimes after distractor
% [~,~,met]=get_DistractorDprime_per_mouse(alltbt,trialTypes,metadata); alltbt.distract_dprimes(isinf(alltbt.distract_dprimes))=3;
% alltbt.distract_dprimes=met.distract_dprimes; metadata.distract_dprimes=met.distract_dprimes; trialTypes.distract_dprimes=met.distract_dprimes;

% Optional: how far through session is each trial
excludeNonReachingBeginAndEnd=true;
[metadata,fractionThroughSess]=howFarThroughSession(metadata,excludeNonReachingBeginAndEnd,trialTypes);
if excludeNonReachingBeginAndEnd
    alltbt.fractionThroughSess_adjusted=metadata.fractionThroughSess_adjusted;
    trialTypes.fractionThroughSess_adjusted=metadata.fractionThroughSess_adjusted;
else
    alltbt.fractionThroughSess=metadata.fractionThroughSess;
    trialTypes.fractionThroughSess=metadata.fractionThroughSess;
end

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
saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set

% filter settings
alltbt.mouseid=metadata.mouseid;
alltbt.sessid=metadata.sessid;
trialTypes.sessid=metadata.sessid;
% tbt_filter.sortField='sessid';
% tbt_filter.sortField='dprimes';
% tbt_filter.sortField='day1formice';
% tbt_filter.sortField='distractor_immediate_after_cue';
% tbt_filter.sortField='higherThanBefore';
% tbt_filter.sortField='fractionThroughSess_adjusted';
% tbt_filter.sortField='opto_enhanced_reach';
% tbt_filter.sortField='mouseLearned';
% tbt_filter.sortField='initiallyLowLEDsess';
tbt_filter.sortField='takemice';
% tbt_filter.range_values=[1 6 7 8 10 14 18];
% tbt_filter.range_values=[1 2 6 9 10 11 12 18];
% tbt_filter.range_values=[-100 0.9]; %[0.75 100];
% tbt_filter.range_values=[2 3 6 7 8 9];
% tbt_filter.range_values=[-100 0.75] ;%0.471];
% tbt_filter.range_values=[0 100]; % beginner: d<0.25, intermediate: 0.25<=d<0.75, expert: d>=0.75
% tbt_filter.range_values=[0.75 100]; % maybe 2,6,7,12
tbt_filter.range_values=[0.5 1.5]; % maybe 2,6,7,12
% tbt_filter.range_values=[2 3 4 5 6 7 8 9 10 11 12 14 15 17 18 19]; % which mice start at non-learning 
% tbt_filter.range_values=[1 2 4 5 6 7 8 9 10 11 12 17 18 19];
% tbt_filter.range_values=[1     2     3     6     7     8     9    10    11    12    14    15    17    18];
tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
temp=tbt_filter.name;
temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))=''; 
temp=temp(~isspace(temp));
tbt_filter.name=temp;
tbt_filter.clock_progress=true;

% filter alltbt
[alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);

%% check for opto-enhanced reaching
alltbt.sessid=metadata.sessid;
alltbt=checkForOptoEnhancedReach(alltbt,metadata,trialTypes,'all_reachBatch','trialTypes.optoGroup==2','cueZone_onVoff',[-0.25 0.5],20);
% alltbt=checkForOptoEnhancedReach(alltbt,metadata,trialTypes,'all_reachBatch','trialTypes.led==1','cueZone_onVoff',[-0.25 0.5],20);
trialTypes.opto_enhanced_reach=alltbt.opto_enhanced_reach;

%% find sessions where mouse learned
part1_fracThroughSess=[0 0.25];
part2_fracThroughSess=[0.25 1];
learningThresh=0.1;
alltbt=findSessWhereMouseLearned(alltbt,metadata,trialTypes,part1_fracThroughSess,part2_fracThroughSess,learningThresh);
learned1=alltbt.mouseLearned;
part1_fracThroughSess=[0 0.5];
part2_fracThroughSess=[0.5 1];
learningThresh=0.1;
alltbt=findSessWhereMouseLearned(alltbt,metadata,trialTypes,part1_fracThroughSess,part2_fracThroughSess,learningThresh);
learned2=alltbt.mouseLearned;
part1_fracThroughSess=[0 0.75];
part2_fracThroughSess=[0.75 1];
learningThresh=0.1; % used 0.1 for most of Fig 3
alltbt=findSessWhereMouseLearned(alltbt,metadata,trialTypes,part1_fracThroughSess,part2_fracThroughSess,learningThresh);
learned3=alltbt.mouseLearned;
alltbt.mouseLearned=learned1 | learned2 | learned3;
trialTypes.mouseLearned=alltbt.mouseLearned;

%% learning curves
% DURING INITIAL LEARNING CONTROL OR SILENCING
[learningC,days,reachrate_cued,reachrate_uncued,dayNdprime,day1dprime,quiverTips,biases]=learningCurves(alltbt,trialTypes,metadata,'sess_wrt_day1',[1],[15:20],false);
% for i=1:size(reachrate_cued,1)
%     reachrate_cued(i,days==1)=nanmean(reachrate_cued(i,days<=1),2);
%     reachrate_uncued(i,days==1)=nanmean(reachrate_uncued(i,days<=1),2);
% end
% plotCuedAndUncuedReachingOverDays(reachrate_cued,reachrate_uncued,days,1:20);

% ALIGN RECOVERY TO FIRST SESSION
% metadata.sess_wrt_day1=metadata.nth_session; 
% ui=unique(metadata.mouseid);
% for i=1:length(ui)
%     metadata.sess_wrt_day1(metadata.mouseid==ui(i))=metadata.sess_wrt_day1(metadata.mouseid==ui(i))-min(metadata.nth_session(metadata.mouseid==ui(i)),[],1,'omitnan')+1;
% end
% [learningC,days,reachrate_cued,reachrate_uncued,dayNdprime,day1dprime,quiverTips]=learningCurves(alltbt,trialTypes,metadata,'sess_wrt_day1','first','last',false);

%% If want to remove trials where distractor turns on immediately after cue
% Optional: discard trials where distractor turns on immediately after cue
alltbt.distractor_immediate_after_cue=any(alltbt.movie_distractor(:,93:123)>0.5,2); % over 1 sec after cue
trialTypes.distractor_immediate_after_cue=alltbt.distractor_immediate_after_cue; metadata.distractor_immediate_after_cue=alltbt.distractor_immediate_after_cue;
trialTypes.sess_wrt_day1=metadata.sess_wrt_day1; alltbt.sess_wrt_day1=metadata.sess_wrt_day1;
% use filtTbt

%% build relevant data sets

cuedreachtimewindow=0.4; % in seconds

% cuedreachinds=94:94+floor(cuedreachtimewindow./0.035);
% trialTypes.did_cued_reach=(any(alltbt.reachBatch_miss_reachStarts(:,cuedreachinds)>0.5,2) | any(alltbt.pelletmissingreach_reachStarts(:,cuedreachinds)>0.5,2)) & trialTypes.touched_pellet==0; % cued failure and no uncued success
% % trialTypes.did_cued_reach=any(alltbt.reachBatch_miss_reachStarts(:,94+128:94+214)>0.5,2) | any(alltbt.pelletmissingreach_reachStarts(:,94+128:94+214)>0.5,2); % delayed failure, 4.5 sec to 7.5 sec
% trialTypes.did_cued_reach_1forward=[trialTypes.did_cued_reach(2:end); 0];

% settings for paired RT data set
test.nInSequence=[3]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
% requirement for first trial in pair
% trial1='trialTypes.led==0';
trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.isLongITI_1back=[1; trialTypes.isLongITI(1:end-1)];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.led_11back=[ones(11,1); trialTypes.led(1:end-11)];
trialTypes.led_7back=[ones(7,1); trialTypes.led(1:end-7)];
trialTypes.led_6back=[ones(6,1); trialTypes.led(1:end-6)];

% % trial1='trialTypes.led==0';
% % trial1='trialTypes.led==0'; % | trialTypes.led_1back==0 | trialTypes.led_2back==0';
% % memory
% %this %trial1='trialTypes.led~=1'; 
% % trial1='trialTypes.isLongITI==1';
% trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
% trial1='trialTypes.optoGroup~=1 & trialTypes.did_cued_reach_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.isLongITI_1forward==1';
% % trial1='trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1'; % & trialTypes.isLongITI_1forward==1'];
% % trial1='trialTypes.touch_in_cued_window_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1';
% % trial1='trialTypes.cued_reach_1forward==1 & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0) & trialTypes.optoGroup~=1  & trialTypes.optoGroup_1forward~=1';
% % trial1='trialTypes.cued_reach_1forward==0  & trialTypes.touched_pellet_1forward==1 & (trialTypes.led_1forward==0) & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
% % trial1='trialTypes.cued_reach_1forward==1 & trialTypes.consumed_pellet_1forward==0 & trialTypes.led_1forward==0 & trialTypes.optoGroup_1forward~=1 & trialTypes.optoGroup~=1 & trialTypes.isLongITI_1forward==1';
% % trial1='trialTypes.optoGroup~=1 & trialTypes.consumed_pellet_1back==1 & trialTypes.after_cue_success_1forward==1 & trialTypes.consumed_pellet_1forward==1 & trialTypes.led_1forward==1 & trialTypes.optoGroup_1forward~=1';
% trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
% % memory
% % this %trial2='trialTypes.led==1 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3 & trialTypes.led_1forward==0';
% % trial2='trialTypes.optoGroup~=1 & trialTypes.led==0 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
% trial2='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
% trial2='trialTypes.optoGroup~=1 & trialTypes.led==1';
% % trial2='trialTypes.led==1 & trialTypes.led_1forward==0 & trialTypes.led_6back==0'; 
% % trial2='trialTypes.led==0 & trialTypes.led_1forward==1 & trialTypes.led_7back==1'; % & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';

% No opto inh
[trialTypes,trial1,trial2,~,~]=whichTrialTypesToUse(alltbt,trialTypes,metadata,'cued failure',[0 1],'missOrNoPellet'); % no LED
% Opto inh
% [trialTypes,~,~,trial1,trial2]=whichTrialTypesToUse(alltbt,trialTypes,metadata,'cued failure',[0 1],'missOrNoPellet'); % LED

% trial1=['(trialTypes.optoGroup~=1 & trialTypes.led==0 & trialTypes.consumed_pellet_1back==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0) | ' ...
%         '(trialTypes.optoGroup~=1 & trialTypes.led==0 & trialTypes.isLongITI_1forward==1 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
% trial2='trialTypes.optoGroup~=1 & trialTypes.led==0 & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
% 
% trial1='trialTypes.optoGroup==3';
% trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';

test.templateSequence2_cond=eval(trial1);
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
% test.onlyConsiderReachType='reachBatch_success_reachStarts';
test.onlyConsiderReachType=[];
if ~isempty(test.onlyConsiderReachType)
    alltbt.backup_all_reachBatch=alltbt.all_reachBatch;
    alltbt.all_reachBatch=alltbt.(test.onlyConsiderReachType);
end

saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);

% save settings for paired RT data set
save([saveDir2 '\test_settings.mat'],'test');

% build paired RT data set
fakeCueInd=50; % in indices, this is not relevant if not using PCA-based RT model
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,correctedDistributions]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 
if ~isempty(test.onlyConsiderReachType)
    alltbt.all_reachBatch=alltbt.backup_all_reachBatch;
end

%% plot events in session

plotBehaviorSessEvents_wrapper(alltbt,metadata,trialTypes,[37.5 38.5]); % last arg selects sessid, e.g., 1501; can enter "optoOn" at thresh prompt

%% plot features in data set

% last argument chooses type of plot
% see function plotBehaviorEventFx.m for options
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
% returnThisCDF=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching_cdf'); 

%% plot trial to trial change in reach CDF

plotChangeInReachCDF(dataset.realDistributions,alltbt,6); %9.5);

%% measure how reach rate changes over course of session

shuffleTrialOrder=false; % if want to randomly permute trial order to test for ordering effects

reachratesettings.epsilon_cue=0; % in seconds
reachratesettings.epsilon_uncue=2; % in seconds
reachratesettings.epsilon_beforecue=1; % in seconds
reachratesettings.percentOfReachesFromSess_forInitCond=20; % use this fraction of trials from beginning of session to get initial windows
reachratesettings.percentOfReachesFromSess_forInitRate=20; % use this fraction of trials from beginning of session to get initial reach rates
reachratesettings.maxTrialLength=9.5; % in sec, wrt cue
reachratesettings.minTrialLength=-2; % wrt cue, in sec
reachratesettings.suppressPlots=false;
 % sec wrt cue onset
reachratesettings.acrossSess_window1=[0.05 1]; % cued window [0.05 1]
% reachratesettings.acrossSess_window1=[4 7];
% note that after mouse gets a pellet, reaching is suppressed
reachratesettings.acrossSess_window2=[7 reachratesettings.maxTrialLength]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[reachratesettings.minTrialLength -1];
reachratesettings.scatterPointSize=50; % size for points in scatter plot
reachratesettings.addSatietyLines=false; % whether to add proportionality lines to figure
% reachratesettings.stopPlottingTrialsAfterN=500; % will stop plotting after this nth trial in session, also only use this many trials for regression fit -- see next line, also controls colormap
% reachratesettings.stopPlottingTrialsAfterN=175; % will stop plotting
reachratesettings.stopPlottingTrialsAfterN=200; % will stop plotting
% after this nth trial in session, also only use this many trials for
% regression fit -- see next line, also controls colormap
reachratesettings.showFitLine=true; % whether to show linear fit to change across trials
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
reachratesettings.initWindows=[]; % empty if want to calculate from dataset
reachratesettings.addSessionLines=false; % for no averaging across sessions plot, whether to connect trial bins within same session with lines
reachratesettings.binTrialsForAvAcrossSess=true; % whether to bin multiple trials for first figure, will bin into binThisManyTrials
reachratesettings.binThisManyTrials=30; %6; % how many trials to bin within each session
reachratesettings.nBinsForZones=40; % will be nBinsForZones squared total bins, this is # bins for each x and y axis
reachratesettings.useRateMethod=3; % 1, 2 or 3 (see explanation below)
% There are 3 approaches available for determing reach rates
% Code will calculate all three but only return useRateMethod
%
% Approach 1: get reach rate in each window over course of session, where windows are specified
% wrt animal's behavior at the beginning of the session (fixed windows)
% i.e., take the average reaction time over the first X% of trials, where X
% is percentOfReachesFromSess_forInitCond, then windows are wrt this av
% reaction time, such that the following eqs hold given that t_n is this
% average reaction time
% window1 = @(t_cue,t_n,epsilon_cue) [t_cue nanmin([t_n+epsilon_cue maxTrialLength])];                % cued window 
% window2 = @(t_cue,t_n,epsilon_uncue) [nanmin([t_n+epsilon_uncue maxTrialLength]) maxTrialLength];   % after cued window
% window3 = @(t_cue,t_n,epsilon_beforecue) [minTrialLength t_cue-epsilon_beforecue];                  % before cue window     
%
% Approach 2: get reach rate in each window over course of session, where
% windows change with changing behavior over the course of session
% (non-fixed windows)
% i.e., using eqs for window1, window2, window3 above, calculate window
% using the reaction time on the nth trial, then apply this window for
% calculating the reach rate on the n+1th trial
%
% Approach 3: get reach rate in each window over course of session,
% where windows are fixed across all sessions
% i.e., use windows specified above as acrossSess_window1, acrossSess_window2, acrossSess_window3 

reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 

%% Plot dprimes
plotVersusFrac=false;
[dprimes,fracs_over_sess]=plotDprimesFromReachRates(reachrates,false,plotVersusFrac);

%% Plot hallmarks of satiety

plotHallmarksOfSatiety(reachrates,dataset,alltbt,metadata,trialTypes);

%% memory effect
% get bestLearningSoFar.m to find days when new higher dprime achieved,
% using whole session to calc dprime so as not to select for any particular 
% shape of within-session change
% [alltbt,metadata,trialTypes,higherThanBefore_sessIDs]=bestLearningSoFar(alltbt,metadata,trialTypes,true);
% also filter for mouse learned
% nInSeq=5; % WHAT I USED FOR OPTO GRC TALK
% useFractionThroughSession=[0.4 0.6];
% plotCDFUpTo=3;
% memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession,[],plotCDFUpTo);
nInSeq=2; 
useFractionThroughSession=[0.5 1];
withinCueTimeWindow=0.4; % in sec
beginOfSess=0.25; % up to this fraction of session
plotCDFUpTo=9.5;
memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession,[],plotCDFUpTo,withinCueTimeWindow,beginOfSess);
% ANOTHER WAY
trialTypes.lowLEDsequence=findSeqsWithN_of_condition(trialTypes, 'led', 55, 150, false);
[~,ui]=unique(trialTypes.sessid);
trialTypes.initiallyLowLED=zeros(size(trialTypes.lowLEDsequence));
trialTypes.initiallyLowLED(ui)=trialTypes.lowLEDsequence(ui);
disp(nansum(trialTypes.initiallyLowLED))
theseAreLowLEDsess=unique(trialTypes.sessid(trialTypes.initiallyLowLED==1));
trialTypes.initiallyLowLEDsess=ismember(trialTypes.sessid,theseAreLowLEDsess); alltbt.initiallyLowLEDsess=trialTypes.initiallyLowLEDsess;
% filter tbt according to initiallyNoLED, then
nInSeq=2; useFractionThroughSession=[0.2 0.4]; plotCDFUpTo=9.5;
memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession,[],plotCDFUpTo);

%% led ongoing reach motor effects
nInSeq=2; 
useFractionThroughSession=[0 1];
plotCDFUpTo=3;
memoryEffect(alltbt,metadata,trialTypes,nInSeq,useFractionThroughSession,'reachBatch_miss_reachStarts',plotCDFUpTo);

%% p(reach preceded by cue) and p(reach followed by cue) for a sequence
% first run "build relevant data sets" to find sequences
usePositions=-10:20;
returnout=plotByPosition(dataset, alltbt, trialTypes, metadata, usePositions, 'all_reachBatch', 0.5);

%% shift in reach rate between trial pair

backup_test=test;

% get initial reach rates within each session using all trials
test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
test.trial1=trial1;
test.templateSequence2_cond=eval(trial1);
test.trial2=trial2;
test.templateSequence2_end=eval(trial2);
test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');
skipCorrected=true;
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,50,saveDir,test,skipCorrected); 
reachratesettings.initWindows=[]; 
reachratesettings.suppressPlots=true;
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
initWindows{1}=reachrates.init_fixed_window1;
initWindows{2}=reachrates.init_fixed_window2;
initWindows{3}=reachrates.init_fixed_window3;
reachratesettings.initWindows=initWindows; 
% then run trial pairs
test=backup_test;
dataset=buildReachingRTModel(alltbt,trialTypes,metadata,50,saveDir,test,skipCorrected); 
reachrates=plotChangeInReachProbability_fromRTdataset(dataset,metadata,alltbt,'cueZone_onVoff',shuffleTrialOrder,reachratesettings); 
plotPairedChangeMinusSatiety(reachrates);

%% plot outcome-dependent shifts
% scriptToMakeOutcomeFigure4(alltbt,trialTypes,metadata,cuedreachtimewindow,optoDuration,sad,savePlotsDir,tbt_filter);
reachratesettings.acrossSess_window1=[0 cuedreachtimewindow]; % time window wrt cue onset to classify reach as cued
reachratesettings.acrossSess_window2=[7 9.5]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=[-2 0]; % time window wrt cue onset to classify reach as uncued
reachratesettings.useWindowsForUncued=[3]; % to use window2 or window3 or both for the uncued reach rate
optoDuration=1; % in sec, for these data
timeWindowToCountAsEventReach=[3.5 8]; % set as nan if want to use default in outcomeDependentShift_wrapper.m
% test.trial1=nan; % set as nan if want to use default in outcomeDependentShift_wrapper.m, else will inherit from this function
[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,metadata,'uncued failure accumulate',timeWindowToCountAsEventReach,'all_reachBatch');
test.trial1=trial1; test.trial2=trial2; test.trial1_LED=trial1_LED; test.trial2_LED=trial2_LED;
% for accumulate
test.nInSequence=4; test.fillInBetweenWithAnything=false;
whichToPlot='false alarm'; % whichToPlot can be 'success','delayed success','drop','cued touch','cued touch and switch color','failed cued reach','false alarm','no reach','basic','wildcard','backward success'
[~,~,returnout]=outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,saveDir,[],[],reachratesettings,timeWindowToCountAsEventReach,test,whichToPlot); 

%% plot outcome-dependent shifts AND separate by dprime

alltbt=backup.alltbt; trialTypes=backup.trialTypes; metadata=backup.metadata;
outcomeDependentShift_acrossDprimes(alltbt,trialTypes,metadata);

%% plot change in behavior within session and across days for a single mouse (or average mouse)

% if want to do average mouse, use nth_session instead of sessid
doAverageMouse=true;
if doAverageMouse==true
    alltbt=mouse_alltbt; metadata=mouse_metadata; trialTypes=mouse_trialTypes;
%     alltbt=backup.alltbt; trialTypes=backup.trialTypes; metadata=backup.metadata;
    alltbt.sessid=metadata.sessid;
    backup_sessid=alltbt.sessid;
%     
    % for mice where didn't analyze each day, fix nth_session
    alltbt.sessid(backup_sessid==1)=1;
    alltbt.sessid(backup_sessid==2)=3;
    alltbt.sessid(backup_sessid==3)=7;
    alltbt.sessid(backup_sessid==4)=9;
    
    alltbt.sessid(backup_sessid==5)=2-1;
    alltbt.sessid(backup_sessid==6)=4-1;
    alltbt.sessid(backup_sessid==7)=8-1;
    alltbt.sessid(backup_sessid==8)=10-1;
    alltbt.sessid(backup_sessid==9)=12-1;
    alltbt.sessid(backup_sessid==10)=14-1; 
    alltbt.sessid(backup_sessid==11)=19; % 16
    
    alltbt.sessid(backup_sessid==16)=2-1;
    alltbt.sessid(backup_sessid==17)=4-1;
    alltbt.sessid(backup_sessid==18)=7; %6-1;
    alltbt.sessid(backup_sessid==19)=20-1;
    alltbt.sessid(backup_sessid==20)=21;
    alltbt.sessid(backup_sessid==21)=22; % 22 
    
    alltbt.sessid(backup_sessid==37)=3;
    alltbt.sessid(backup_sessid==38)=4+1;
    
    alltbt.sessid(backup_sessid==75)=1;
    alltbt.sessid(backup_sessid==76)=21;
    alltbt.sessid(backup_sessid==77)=22; 
    alltbt.sessid(backup_sessid==78)=23;
    
    alltbt.sessid(backup_sessid==153)=1;
    alltbt.sessid(backup_sessid==154)=5;
    alltbt.sessid(backup_sessid==155)=7;
    alltbt.sessid(backup_sessid==156)=9;
    alltbt.sessid(backup_sessid==157)=13;
    alltbt.sessid(backup_sessid==158)=16-1; % 15 watched video there is some cued reaching but also to distractor and box was open, not sure
    alltbt.sessid(backup_sessid==159)=17;
    alltbt.sessid(backup_sessid==160)=19; % 19 so much preemptive reaching
    alltbt.sessid(backup_sessid==161)=19; %21; 
    alltbt.sessid(backup_sessid==162)=23;
    
    alltbt.sessid(backup_sessid==315)=6-1;
    alltbt.sessid(backup_sessid==316)=10-1;
    alltbt.sessid(backup_sessid==317)=12-1;
    alltbt.sessid(backup_sessid==318)=14-1;
    alltbt.sessid(backup_sessid==319)=16-1;
    alltbt.sessid(backup_sessid==320)=20-1;
    alltbt.sessid(backup_sessid==321)=19; %22-1;
    alltbt.sessid(backup_sessid==322)=24-1; %
    alltbt.sessid(backup_sessid==323)=29; % 27 but it's the only one, average
    alltbt.sessid(backup_sessid==324)=29;
    
    alltbt.sessid(backup_sessid==638)=1;
    alltbt.sessid(backup_sessid==639)=3;
    alltbt.sessid(backup_sessid==640)=7;
    alltbt.sessid(backup_sessid==641)=9;
    alltbt.sessid(backup_sessid==642)=11;
    alltbt.sessid(backup_sessid==643)=13;
    alltbt.sessid(backup_sessid==644)=19; % 21; 
    alltbt.sessid(backup_sessid==645)=23;
    alltbt.sessid(backup_sessid==646)=25; % 25
    alltbt.sessid(backup_sessid==647)=26+1;
    alltbt.sessid(backup_sessid==648)=29;
    alltbt.sessid(backup_sessid==649)=29; % 30 
    
    alltbt.sessid(backup_sessid==1287)=2-1;
    alltbt.sessid(backup_sessid==1288)=4-1;
    alltbt.sessid(backup_sessid==1289)=6-1;
    alltbt.sessid(backup_sessid==1290)=9-2;
    alltbt.sessid(backup_sessid==1291)=10-1;
    alltbt.sessid(backup_sessid==1292)=12-1;
    alltbt.sessid(backup_sessid==1293)=14-1; 
    alltbt.sessid(backup_sessid==1294)=16-1;
    alltbt.sessid(backup_sessid==1295)=18-1;
    alltbt.sessid(backup_sessid==1296)=20-1;
    alltbt.sessid(backup_sessid==1297)=19; %22-1;
    alltbt.sessid(backup_sessid==1298)=29; % 31
    metadata.sessid=alltbt.sessid;
    trialTypes.sessid=alltbt.sessid;
    
    
%     alltbt.sessid=metadata.nth_session;
%     metadata.sessid=metadata.nth_session;
%     trialTypes.sessid=metadata.nth_session;
% %     start_sessid=nanmin(metadata.nth_session);
% %     end_sessid=nanmax(metadata.nth_session); % but may be quite noisy
    start_sessid=1;
%     end_sessid=26;
%     end_sessid=12;
    end_sessid=29;
%     end_sessid=18;
else
    start_sessid=33;
    end_sessid=71;
end
% outlierTest='tempcued(j)>1 && i<10'; % make empty if don't want to remove any outlier points
outlierTest=[]; % make empty if don't want to remove any outlier points
% learningFigSaveDir='Z:\Kim\example mouse learning\';
learningFigSaveDir=[];
plotDprimes=true;
plotWithinSession_and_dayByDay(alltbt,metadata,trialTypes,start_sessid,end_sessid,learningFigSaveDir,'.png',outlierTest,plotDprimes);
if doAverageMouse==true
    alltbt.sessid=backup_sessid;
    metadata.sessid=backup_sessid;
    trialTypes.sessid=backup_sessid;
end

%% Mouse by mouse change in dprime

[dfirst,dsecond]=getMouseByMouseChangeInDprime(alltbt,metadata,trialTypes,1,19,learningFigSaveDir,'.png',outlierTest,plotDprimes);

%% Plot successes drops etc and dprime
[~,ui]=unique(metadata.nth_session);
dprimes=metadata.dprimes(ui);
nthsess=unique(metadata.nth_session);
success=nan(1,length(nthsess));
drop=nan(1,length(nthsess));
miss=nan(1,length(nthsess));
for i=1:length(nthsess)
    currsess=nthsess(i);
    success(i)=nansum(trialTypes.success_in_cued_window(metadata.nth_session==currsess)==1);
    drop(i)=nansum(trialTypes.touch_in_cued_window(metadata.nth_session==currsess)==1 & trialTypes.consumed_pellet(metadata.nth_session==currsess)==0);
    miss(i)=nansum(trialTypes.cued_reach(metadata.nth_session==currsess)==1 & trialTypes.touched_pellet(metadata.nth_session==currsess)==0);
end
figure(); 
plot(nthsess,dprimes,'Color','k'); title('dprimes');
figure(); 
plot(nthsess,success,'Color','g'); hold on; title('outcomes');
plot(nthsess,drop,'Color','r');
plot(nthsess,miss,'Color','c');

%%
figure();
whichTrials=6:19;
for j=1:length(whichTrials)-1
i=whichTrials(j);
quiver(fraction_trials_reach_in_uncued(i)/1.5,fraction_trials_reach_in_cued_window(i)/1.5,(fraction_trials_reach_in_uncued(i+1)/1.5)-(fraction_trials_reach_in_uncued(i)/1.5),(fraction_trials_reach_in_cued_window(i+1)/1.5)-(fraction_trials_reach_in_cued_window(i)/1.5),'Color',cmap(mapintocmap(j),:));
hold on;
end