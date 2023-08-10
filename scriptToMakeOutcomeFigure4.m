% scriptToMakeOutcomeFigure4.m
function scriptToMakeOutcomeFigure4(alltbt,trialTypes,metadata,cuedreachtimewindow,optoDuration,sad,savePlotsDir,tbt_filter,addtoall)

% CONTROLS
% Current trial CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[],'no no opto','all_reachBatch',false,3,'success',sad,savePlotsDir,['win3_minus3to025' addtoall],tbt_filter);


% Cued success CONTROL
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,['win2and3_minus2to0BACKWARD' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[0 optoDuration],'backwards cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,['win3_minus3to025_optoDurforcueBACKWARD' addtoall],tbt_filter);


% Cued failure CONTROL
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards cued failure','missOrNoPellet',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0BACKWARD' addtoall],tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'backwards cued failure','missOrNoPellet',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcueBACKWARD' addtoall],tbt_filter);


% All cued failures including drops CONTROL
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards all cued failures','anyFail',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0BACKWARD' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[0 optoDuration],'backwards all cued failures','anyFail',false,3,'failed cued reach',sad,savePlotsDir,['win3_minus3to025_optoDurforcueBACKWARD' addtoall],tbt_filter);


% Cued drops CONTROL
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards cued drop','reachBatch_drop_reachStarts',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0BACKWARD' addtoall],tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'backwards cued drop','reachBatch_drop_reachStarts',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcueBACKWARD' addtoall],tbt_filter);


% Uncued drops CONTROL
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'backwards uncued drop','reachBatch_drop_reachStarts',false,3,'false alarm',sad,savePlotsDir,['win2and3_minus2to0BACKWARD' addtoall],tbt_filter);


% Uncued success CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[3.5 7],'backwards delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,['win3_minus3to025BACKWARD' addtoall],tbt_filter);


% Uncued failure CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[3.5 7],'backwards uncued failure','anyFail',false,3,'false alarm',sad,savePlotsDir,['win3_minus3to025BACKWARD' addtoall],tbt_filter);


% EXPERIMENTS
% Cued success
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[0 optoDuration],'cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,['win3_minus3to025_optoDurforcue' addtoall],tbt_filter);


% Cued failure
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued failure','missOrNoPellet',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure','missOrNoPellet',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);


% All cued failures
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'all cued failures','anyFail',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[0 optoDuration],'all cued failures','anyFail',false,3,'failed cued reach',sad,savePlotsDir,['win3_minus3to025_optoDurforcue' addtoall],tbt_filter);


% Cued drop
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued drop','reachBatch_drop_reachStarts',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued drop','reachBatch_drop_reachStarts',false,3,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);


% Uncued drop
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued drop','reachBatch_drop_reachStarts',false,3,'false alarm',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);


% Uncued success
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[3.5 7],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,['win3_minus3to025' addtoall],tbt_filter);


% Uncued failure
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-3 -0.25],[3],[3.5 7],'uncued failure','anyFail',false,3,'false alarm',sad,savePlotsDir,['win3_minus3to025' addtoall],tbt_filter);
return

% Cued success accumulate
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued success accumulate','reachBatch_success_reachStarts',true,4,'success',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued success accumulate','reachBatch_success_reachStarts',true,4,'success',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);
X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued success accumulate','reachBatch_success_reachStarts',true,5,'success',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);


% Cued failure accumulate 
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued failure accumulate','missOrNoPellet',true,4,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure accumulate','missOrNoPellet',true,4,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);
X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure accumulate','missOrNoPellet',true,5,'failed cued reach',sad,savePlotsDir,['win2and3_minus2to0_optoDurforcue' addtoall],tbt_filter);


% Uncued success accumulate 
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,4,'delayed success',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,5,'delayed success',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);


% Uncued failure accumulate
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','missOrNoPellet',true,4,'false alarm',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','missOrNoPellet',true,5,'false alarm',sad,savePlotsDir,['win2and3_minus2to0' addtoall],tbt_filter);


% Uncued success TWO HALF SECS BEFORE OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[-1 -0.5],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,['win2and3_minus2to01SECBEFOREOPTO' addtoall],tbt_filter);
% Uncued success RIGHT AFTER OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[optoDuration optoDuration+0.5],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,['win2and3_minus2to0HALFSECAFTEROPTO' addtoall],tbt_filter);
% Uncued success RIGHT BEFORE OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[-0.5 0],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,['win2and3_minus2to0HALFSECBEFOREOPTO' addtoall],tbt_filter);


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

function plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,acrossSess_window3,useWindowsForUncued,timeWindowToCountAsEventReach,whichEventName,whichReachName,isAccum,nInSequence,whichToPlotName,sad,savePlotsDir,addOn,tbt_filter)

cueonsetadjust=0; %-0.25; % cue duration
cuedreachtimewindow=cuedreachtimewindow+cueonsetadjust;
acrossSess_window3=acrossSess_window3+cueonsetadjust;
timeWindowToCountAsEventReach=timeWindowToCountAsEventReach+cueonsetadjust;

reachratesettings.acrossSess_window1=[0+cueonsetadjust-0.1 cuedreachtimewindow]; % time window wrt cue onset to classify reach as cued
reachratesettings.acrossSess_window2=[7+cueonsetadjust 9.5+cueonsetadjust]; % beware reach suppression after a success
reachratesettings.acrossSess_window3=acrossSess_window3; %[-2 0]; % time window wrt cue onset to classify reach as uncued
reachratesettings.useWindowsForUncued=useWindowsForUncued; %[3]; % to use window2 or window3 or both for the uncued reach rate
% optoDuration=1; % in sec, for these data
% timeWindowToCountAsEventReach=[3.5 8]; % set as nan if want to use default in outcomeDependentShift_wrapper.m
% test.trial1=nan; % set as nan if want to use default in outcomeDependentShift_wrapper.m, else will inherit from this function
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