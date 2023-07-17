% scriptToMakeOutcomeFigure4.m
function scriptToMakeOutcomeFigure4(alltbt,trialTypes,metadata,cuedreachtimewindow,optoDuration,sad,savePlotsDir,tbt_filter)

% CONTROLS
% Cued success CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,'win2and3_minus2to0BACKWARD',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'backwards cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcueBACKWARD',tbt_filter);


% Cued failure CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'backwards cued failure','all_reachBatch',false,3,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0BACKWARD',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'backwards cued failure','all_reachBatch',false,3,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcueBACKWARD',tbt_filter);


% Uncued failure CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'backwards uncued failure','all_reachBatch',false,3,'false alarm',sad,savePlotsDir,'win2and3_minus2to0BACKWARD',tbt_filter);


% Uncued failure CONTROL
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'backwards uncued failure','all_reachBatch',false,3,'false alarm',sad,savePlotsDir,'win2and3_minus2to0BACKWARD',tbt_filter);

return

% EXPERIMENTS
% Cued success
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued success','reachBatch_success_reachStarts',false,3,'success',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);


% Cued failure
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued failure','all_reachBatch',false,3,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure','all_reachBatch',false,3,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);


% Uncued success
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);


% Uncued failure
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure','all_reachBatch',false,3,'false alarm',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);

% return

% Cued success accumulate
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued success accumulate','reachBatch_success_reachStarts',true,4,'success',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued success accumulate','reachBatch_success_reachStarts',true,4,'success',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);
X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued success accumulate','reachBatch_success_reachStarts',true,5,'success',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);


% Cued failure accumulate 
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[0 cuedreachtimewindow],'cued failure accumulate','all_reachBatch',true,4,'failed cued reach',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 cuedreachtimewindow],'cued failure accumulate','all_reachBatch',true,4,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[0 cuedreachtimewindow],'cued failure accumulate','all_reachBatch',true,4,'failed cued reach',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[0 cuedreachtimewindow],'cued failure accumulate','all_reachBatch',true,4,'failed cued reach',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);

plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure accumulate','all_reachBatch',true,4,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);

X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[0 optoDuration],'cued failure accumulate','all_reachBatch',true,5,'failed cued reach',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[0 optoDuration],'cued failure accumulate','all_reachBatch',true,5,'failed cued reach',sad,savePlotsDir,'win2and3_minus2to0_optoDurforcue',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[0 optoDuration],'cued failure accumulate','all_reachBatch',true,5,'failed cued reach',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[0 optoDuration],'cued failure accumulate','all_reachBatch',true,5,'failed cued reach',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);


% Uncued success accumulate 
X=2;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 4 sequence, "2forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,4,'delayed success',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,4,'delayed success',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,4,'delayed success',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,4,'delayed success',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);

X=3;
trialTypes.isLongITI_Xforward=[trialTypes.led(1+X:end); zeros(X,1)]; % for 5 sequence, "3forward"
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,5,'delayed success',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,5,'delayed success',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,5,'delayed success',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'delayed success accumulate','reachBatch_success_reachStarts',true,5,'delayed success',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);


% Uncued failure accumulate
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);

% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window3_minus2to0',tbt_filter);
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'win2and3_minus2to0',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window3_minus2tominus1',tbt_filter);
% plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window2and3_minus2tominus1',tbt_filter);


% Uncued success TWO HALF SECS BEFORE OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[-1 -0.5],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,'win2and3_minus2to01SECBEFOREOPTO',tbt_filter);
% Uncued success RIGHT AFTER OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[optoDuration optoDuration+0.5],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,'win2and3_minus2to0HALFSECAFTEROPTO',tbt_filter);
% Uncued success RIGHT BEFORE OPTO
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[-0.5 0],'delayed success','reachBatch_success_reachStarts',false,3,'delayed success',sad,savePlotsDir,'win2and3_minus2to0HALFSECBEFOREOPTO',tbt_filter);


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

reachratesettings.acrossSess_window1=[0 cuedreachtimewindow]; % time window wrt cue onset to classify reach as cued
reachratesettings.acrossSess_window2=[7 9.5]; % beware reach suppression after a success
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
for iFig = 1:length(FigList)-1
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [FigName 'fig' num2str(iFig) '.fig']));
end

end