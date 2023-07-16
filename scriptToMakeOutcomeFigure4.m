% scriptToMakeOutcomeFigure4.m
function scriptToMakeOutcomeFigure4(alltbt,trialTypes,metadata,cuedreachtimewindow,optoDuration,sad,savePlotsDir)

% Cued success

% Cued failure

% Uncued success

% Uncued failure

% Cued success accumulate

% Cued failure accumulate

% Uncued success accumulate

% Uncued failure accumulate
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window3_minus2to0');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window2and3_minus2to0');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window3_minus2tominus1');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,4,'false alarm',sad,savePlotsDir,'window2and3_minus2tominus1');

plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window3_minus2to0');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 0],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window2and3_minus2to0');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window3_minus2tominus1');
plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,[-2 -1],[2 3],[3.5 8],'uncued failure accumulate','all_reachBatch',true,5,'false alarm',sad,savePlotsDir,'window2and3_minus2tominus1');

end

function plotSaveOutShift(alltbt,trialTypes,metadata,cuedreachtimewindow,acrossSess_window3,useWindowsForUncued,timeWindowToCountAsEventReach,whichEventName,whichReachName,isAccum,nInSequence,whichToPlotName,sad,savePlotsDir,addOn)

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
[~,~,returnout]=outcomeDependentShift_wrapper(alltbt,trialTypes,metadata,sad,[],[],reachratesettings,timeWindowToCountAsEventReach,test,whichToPlot); 
% save plots and returnout
mkdir(fullfile(savePlotsDir,[whichReachName '_seqLength' num2str(nInSequence) '_' addOn]));
save(fullfile(savePlotsDir,[whichReachName '_seqLength' num2str(nInSequence) '_' addOn],'returnout.mat'),'returnout');
saveAllOpenFigures(fullfile(savePlotsDir,[whichReachName '_seqLength' num2str(nInSequence) '_' addOn]));
close all;

end

function saveAllOpenFigures(FolderName)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)-1
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [FigName num2str(iFig)], '.fig'));
end

end