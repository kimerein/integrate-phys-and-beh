function [tensor, allLabels, timepoints_for_tensor]=getTensorsForOneSession(whichSess, downSampBy, takeNPointsAfterEvent, takeNPointsBeforeEvent, dd)

% whichSess=81; % 81 has a fair number of cells
% downSampBy=4; % downsamp 60 ms bins by this much
% takeNPointsAfterEvent=15;
% takeNPointsBeforeEvent=0;
% response_to_plot1='cue'; response_to_plot2='uncued_success';
response_to_plot1='cued_success'; response_to_plot2='uncued_success';
putTogetherTensor(whichSess,downSampBy,takeNPointsAfterEvent,takeNPointsBeforeEvent,response_to_plot1,response_to_plot2,dd,'C:\Users\sabatini\Documents\MATLAB\cued_success_then_uncued_success');
response_to_plot1='cued_failure'; response_to_plot2='uncued_failure';
putTogetherTensor(whichSess,downSampBy,takeNPointsAfterEvent,takeNPointsBeforeEvent,response_to_plot1,response_to_plot2,dd,'C:\Users\sabatini\Documents\MATLAB\cued_failure_then_uncued_failure');
load('cued_success_then_uncued_success_tensor.mat')
load('cued_success_then_uncued_success_labels.mat')
tensor=tens; allLabels=labels;
load('cued_failure_then_uncued_failure_tensor.mat')
load('cued_failure_then_uncued_failure_labels.mat')
tensor=cat(3, tensor, tens); allLabels=[allLabels; labels+2];
timepoints=0:(downSampBy*0.06):(size(tensor,2)-1)*(downSampBy*0.06);
timepoints_for_tensor=timepoints;
mkdir(['C:\Users\sabatini\Documents\MATLAB\session' num2str(whichSess)]);
save(['C:\Users\sabatini\Documents\MATLAB\session' num2str(whichSess) '\allLabels.mat'],'allLabels');
save(['C:\Users\sabatini\Documents\MATLAB\session' num2str(whichSess) '\tensor.mat'],'tensor');
save(['C:\Users\sabatini\Documents\MATLAB\session' num2str(whichSess) '\timepoints_for_tensor.mat'],'timepoints_for_tensor');

end