function variedTimingOfOptoWrtReach(alltbt,trialTypes,metadata,cuedreachtimewindow,acrossSess_window3,useWindowsForUncued,timeWindowToCountAsEventReach,whichEventName,sad,savePlotsDir,addOn,tbt_filter)
% sad is current saveDir

% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[0 optoDuration],'backwards cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,['win3_minus3to025_optoDurforcueBACKWARD' addtoall],tbt_filter);

whichReachName='reachBatch_success_reachStarts';
isAccum=false;
nInSequence=3;
whichToPlotName='success';

if isempty(acrossSess_window3)
    acrossSess_window3=[-3 -0.25];
end
if isempty(cuedreachtimewindow)
    cuedreachtimewindow=[-0.1 0.25];
end
if isempty(useWindowsForUncued)
    useWindowsForUncued=3;
end

cueonsetadjust=0; %-0.25; % cue duration
cuedreachtimewindow=cuedreachtimewindow+cueonsetadjust;
acrossSess_window3=acrossSess_window3+cueonsetadjust;
timeWindowToCountAsEventReach=timeWindowToCountAsEventReach+cueonsetadjust;

reachratesettings.acrossSess_window1=cuedreachtimewindow;
reachratesettings.acrossSess_window2=[7+cueonsetadjust 9.5+cueonsetadjust]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=acrossSess_window3; %[-2 0]; % time window wrt cue onset to classify reach as uncued
reachratesettings.useWindowsForUncued=useWindowsForUncued; %[3]; % to use window2 or window3 or both for the uncued reach rate

[trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,metadata,whichEventName,timeWindowToCountAsEventReach,whichReachName);
test.trial1=trial1; test.trial2=trial2; test.trial1_LED=trial1_LED; test.trial2_LED=trial2_LED;
if isAccum==true
    % for accumulate
    test.nInSequence=nInSequence; test.fillInBetweenWithAnything=false;
else
    test.nInSequence=3; test.fillInBetweenWithAnything=true;
end
whichToPlot=whichToPlotName; % whichToPlot can be 'success','delayed success','drop','cued touch','cued touch and switch color','failed cued reach','false alarm','no reach','basic','wildcard','backward success'
close all;
% reach pdf and cdf -- CONTROL
plotReachPDFandCDF(alltbt,trialTypes,metadata,test,tbt_filter,sad,fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn]),'control');
% scatter plots for outcome change
close all;
[~,~,returnout]=outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,sad,[],[],reachratesettings,timeWindowToCountAsEventReach,test,whichToPlot); 
% save plots and returnout
if ~exist(fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn]),'dir')
    mkdir(fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn]));
end
save(fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn],'returnout.mat'),'returnout');
saveAllOpenFigures(fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn]));
close all;
% reach pdf and cdf -- LED
test.trial1=test.trial1_LED;
test.trial2=test.trial2_LED;
plotReachPDFandCDF(alltbt,trialTypes,metadata,test,tbt_filter,sad,fullfile(savePlotsDir,[whichEventName '_seqLength' num2str(nInSequence) '_' addOn]),'pDMStinh');
close all;

end

function saveAllOpenFigures(FolderName)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [FigName 'fig' num2str(iFig) '.fig']));
end

end


function plotReachPDFandCDF(alltbt,trialTypes,metadata,test,tbt_filter,saveDir,savePlotsDir,addOn)

if ~isempty(regexp(test.trial1,'SPLIT'))
    test.templateSequence2_cond=test.trial1;
else
    test.templateSequence2_cond=eval(test.trial1);
end
test.templateSequence2_end=eval(test.trial2);
test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
saveDir2=[saveDir '\' test.event_name];
mkdir(saveDir2);
save([saveDir2 '\test_settings.mat'],'test');

% build paired RT data set
fakeCueInd=50; % in indices, this is not relevant if not using PCA-based RT model
skipCorrected=true;
% this function builds the dataset using the trial type sequences specified above
[dataset,~]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir,test,skipCorrected); 

% last argument chooses type of plot
% see function plotBehaviorEventFx.m for options
[returnThis,returnThisRef]=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching');
returnThisCDF=plotBehaviorEventFx(dataset.realDistributions,alltbt,[],'plot_rawReaching_cdf'); 

% save figures and results
if ~exist(fullfile(savePlotsDir,['pdfcdf' addOn]),'dir')
    mkdir(fullfile(savePlotsDir,['pdfcdf' addOn]));
end
save(fullfile(savePlotsDir,['pdfcdf' addOn],'returnThis.mat'),'returnThis');
save(fullfile(savePlotsDir,['pdfcdf' addOn],'returnThisRef.mat'),'returnThisRef');
save(fullfile(savePlotsDir,['pdfcdf' addOn],'returnThisCDF.mat'),'returnThisCDF');
saveAllOpenFigures(fullfile(savePlotsDir,['pdfcdf' addOn]));
close all;

end

function [trialTypes,trial1,trial2,trial1_LED,trial2_LED]=whichTrialTypesToUse(alltbt,trialTypes,metadata,whichEventType,timeWindow,whichReachInTimeWindow)

[~,cueindma]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
timestep=mode(diff(nanmean(alltbt.times,1)));
if isempty(timeWindow)
    trialTypes.reachedInTimeWindow=ones(size(trialTypes.led));
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 1];
    trialTypes.reachedInTimeWindow_1back=[1; trialTypes.reachedInTimeWindow(1:end-1)];
else
    timeWindowInds(1)=floor(timeWindow(1)/timestep);
    timeWindowInds(2)=floor(timeWindow(2)/timestep);
    temp=alltbt.(whichReachInTimeWindow);
    trialTypes.reachedInTimeWindow=any(temp(:,cueindma+timeWindowInds(1):cueindma+timeWindowInds(2))>0.05,2);
    trialTypes.reachedInTimeWindow_1forward=[trialTypes.reachedInTimeWindow(2:end); 0];
    trialTypes.reachedInTimeWindow_1back=[0; trialTypes.reachedInTimeWindow(1:end-1)];
end

trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
trialTypes.isLongITI_2forward=[trialTypes.led(1+2:end); zeros(2,1)];
trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
trialTypes.optoGroup_1back=[0; trialTypes.optoGroup(1:end-1)];
trialTypes.isLongITI_1back=[0; trialTypes.isLongITI(1:end-1)];
trialTypes.optoGroup_2back=[0; 0; trialTypes.optoGroup(1:end-2)];
trialTypes.noReach=~any(alltbt.all_reachBatch>0.05,2);
trialTypes.noReach_1forward=[trialTypes.noReach(2:end); 0];
trialTypes.noReach_1back=[0; trialTypes.noReach(1:end-1)];
trialTypes.reachedBeforeCue=any(alltbt.all_reachBatch(:,1:cueindma-1)>0.05,2);
trialTypes.reachedAfterCue=any(alltbt.all_reachBatch(:,cueindma:end)>0.05,2);
trialTypes.reachToPelletBeforeCue=any(alltbt.reachStarts_pelletPresent(:,1:cueindma-1)>0.05,2);
trialTypes.reachedBeforeCue_1forward=[trialTypes.reachedBeforeCue(2:end); 0];
trialTypes.reachedBeforeCue_1back=[0; trialTypes.reachedBeforeCue(1:end-1)];
trialTypes.reachToPelletBeforeCue_1forward=[trialTypes.reachToPelletBeforeCue(2:end); 0];
trialTypes.reachToPelletBeforeCue_1back=[0; trialTypes.reachToPelletBeforeCue(1:end-1)];
trialTypes.reachedAfterCue_1forward=[trialTypes.reachedAfterCue(2:end); 0];
trialTypes.led_5forward=[trialTypes.led(1+5:end); zeros(5,1)];
trialTypes.led_6forward=[trialTypes.led(1+6:end); zeros(6,1)];
trialTypes.led_7forward=[trialTypes.led(1+7:end); zeros(7,1)];
trialTypes.led_8forward=[trialTypes.led(1+8:end); zeros(8,1)];
trialTypes.led_9forward=[trialTypes.led(1+9:end); zeros(9,1)];
trialTypes.led_10forward=[trialTypes.led(1+10:end); zeros(10,1)];
trialTypes.led_11forward=[trialTypes.led(1+11:end); zeros(11,1)];
trialTypes.led_12forward=[trialTypes.led(1+12:end); zeros(12,1)];
trialTypes.reachToPelletBeforeCue=any(alltbt.reachStarts_pelletPresent(:,1:cueindma-1)>0.05,2);
trialTypes.reachToPelletBeforeCue_1forward=[trialTypes.reachToPelletBeforeCue(2:end); 0];

% FOR SESSIONS WITH INTERLEAVED OPTO
% linkerForNoLED_accumulate=[' & ((trialTypes.led_1forward==1 & trialTypes.led_2forward==1) | (trialTypes.led_2forward==1 & trialTypes.led_3forward==1) | (trialTypes.led_3forward==1 & trialTypes.led_4forward==1) | (trialTypes.led_4forward==1 & trialTypes.led_5forward==1) | (trialTypes.led_5forward==1 & trialTypes.led_6forward==1)' ...
%                          ' | (trialTypes.led_6forward==1 & trialTypes.led_7forward==1) | (trialTypes.led_7forward==1 & trialTypes.led_8forward==1) | (trialTypes.led_8forward==1 & trialTypes.led_9forward==1) | (trialTypes.led_9forward==1 & trialTypes.led_10forward==1) | (trialTypes.led_10forward==1 & trialTypes.led_11forward==1) | (trialTypes.led_11forward==1 & trialTypes.led_12forward==1))'];
% linkerForLED_accumulate=[' & ((trialTypes.led_1forward==0 & trialTypes.led_2forward==0) | (trialTypes.led_2forward==0 & trialTypes.led_3forward==0) | (trialTypes.led_3forward==0 & trialTypes.led_4forward==0) | (trialTypes.led_4forward==0 & trialTypes.led_5forward==0) | (trialTypes.led_5forward==0 & trialTypes.led_6forward==0)' ...
%                        ' | (trialTypes.led_6forward==0 & trialTypes.led_7forward==0) | (trialTypes.led_7forward==0 & trialTypes.led_8forward==0) | (trialTypes.led_8forward==0 & trialTypes.led_9forward==0) | (trialTypes.led_9forward==0 & trialTypes.led_10forward==0) | (trialTypes.led_10forward==0 & trialTypes.led_11forward==0) | (trialTypes.led_11forward==0 & trialTypes.led_12forward==0))'];
% linkerForNoLED=' & (trialTypes.led_1forward==1 | trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_1back==1)';
% linkerForNoLEDBACKWARDS=' & (trialTypes.led_2forward==1 | trialTypes.led_3forward==1 | trialTypes.led_4forward==1 | trialTypes.led_5forward==1)';
% FOR NO OPTO SESSIONS
linkerForNoLED_accumulate='';
linkerForLED_accumulate='';
linkerForNoLED='';
linkerForNoLEDBACKWARDS='';

% FOR VARIED TIMING
linkerForVariedTimingForward=' & trialTypes.optoGroup_1forward==2';
linkerForVariedTimingSame=' & trialTypes.optoGroup==2';
linkerForVariedTimingUncued='';
% ELSE
% linkerForVariedTimingForward='';
% linkerForVariedTimingSame='';
% linkerForVariedTimingUncued='';

switch whichEventType
    case 'cued success'
%         trial1=['(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)']; 
%         trial2=['trialTypes.optoGroup~=1' linkerForNoLED];
%         trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ')']; 
%         trial2_LED='trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';
        % FOR VARIED TIMING EXPERIMENT
        trial1=['(trialTypes.optoGroup~=1 & trialTypes.optoGroup~=2 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.led_1forward==0)']; 
        trial2=['trialTypes.optoGroup~=1 & trialTypes.optoGroup~=2' linkerForNoLED];
        trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.optoGroup~=2 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.led_1forward==1' linkerForVariedTimingForward ')']; 
        trial2_LED='trialTypes.optoGroup~=1 & trialTypes.optoGroup~=2 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)';

    case 'delayed success'
%         trial1=['(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.cued_reach_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)' linkerForVariedTimingUncued]; 
%         trial2=['trialTypes.optoGroup~=1' linkerForNoLED linkerForVariedTimingUncued];
%         trial1_LED=['(trialTypes.optoGroup~=1 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.cued_reach_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)' linkerForVariedTimingForward linkerForVariedTimingUncued]; 
%         trial2_LED=['trialTypes.optoGroup~=1 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)' linkerForVariedTimingUncued];
        % FOR VARIED TIMING EXPERIMENT
        trial1=['(trialTypes.led==0 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.cued_reach_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==0)' linkerForVariedTimingUncued]; 
        trial2=['trialTypes.led==0' linkerForNoLED linkerForVariedTimingUncued];
        trial1_LED=['(trialTypes.led==0 & trialTypes.reachedBeforeCue_1forward==0 & trialTypes.cued_reach_1forward==0 & trialTypes.reachedInTimeWindow_1forward==1 & trialTypes.optoGroup_1forward~=1 & trialTypes.led_1forward==1)' linkerForVariedTimingForward linkerForVariedTimingUncued]; 
        trial2_LED=['trialTypes.led==0 & (trialTypes.led_1forward==0 | trialTypes.led_2forward==0 | trialTypes.led_3forward==0 | trialTypes.led_4forward==0 | trialTypes.led_1back==0)' linkerForVariedTimingUncued];

end

end