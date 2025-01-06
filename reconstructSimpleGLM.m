function reconstructSimpleGLM(allcoef,ts,behEvents,neuron_data_matrix,timepoints,ds)

% last paragraph of this function assumes first row behEvents is cue

% allcoef, ts and so forth from plotGLMcoef
before0=nansum(ts<0);
after0=nansum(ts>0);
rowByrowSum=zeros(size(allcoef,1),size(behEvents,2));
for i=1:size(allcoef,1)
    currsum=zeros(1,size(behEvents,2));
    currevs=behEvents(i,:);
    currcoef=allcoef(i,:);
    f=find(currevs==1);
    for j=1:length(f)
        inds=f(j)-before0:f(j)+after0;
        currsum(inds(inds>=1 & inds<=length(currsum)))=currsum(inds(inds>=1 & inds<=length(currsum)))+currcoef(inds>=1 & inds<=length(currsum));
    end
    rowByrowSum(i,:)=currsum;
end

figure();
pre=sum(rowByrowSum,1,'omitnan');
plot(timepoints,pre,'Color','r'); hold on;
plot(downSampAv(timepoints,ds),downSampAv(neuron_data_matrix,ds),'Color','k'); legend({'pred','actual ds'});
% plot(timepoints,smooth(neuron_data_matrix,ds),'Color','b'); legend({'pred','actual ds','actual smoo'});
title('All timepoints including ITI');
disp(['R^2 score for all timepoints: ' num2str(calculateR2(neuron_data_matrix, pre))]);

f=find(behEvents(1,:)==1);
inTrialInds=f;
for i=1:50
    inTrialInds=[inTrialInds f+i];
end
inTrialInds=sort(unique(inTrialInds));

figure();
plot(timepoints(inTrialInds),pre(inTrialInds),'Color','r'); hold on;
plot(timepoints(inTrialInds),neuron_data_matrix(inTrialInds),'Color','k'); legend({'pred','actual ds'});
title('Within-trial timepoints');
disp(['R^2 score for within-trial timepoints: ' num2str(calculateR2(neuron_data_matrix(inTrialInds), pre(inTrialInds)))]);

end

function R2 = calculateR2(y_true, y_pred)
% calculateR2 computes the coefficient of determination (R^2) between actual and predicted data.
%
% Syntax:
%   R2 = calculateR2(y_true, y_pred)
%
% Inputs:
%   y_true - Vector of actual observed values.
%   y_pred - Vector of predicted values from a model.
%
% Outputs:
%   R2 - Coefficient of determination.
%
% Example:
%   y_actual = [3; -0.5; 2; 7];
%   y_predicted = [2.5; 0.0; 2; 8];
%   R2 = calculateR2(y_actual, y_predicted)
%   % R2 = 0.9486

    % Input Validation
    if nargin ~= 2
        error('Function requires two input arguments: y_true and y_pred.');
    end
    
    % Ensure inputs are column vectors
    y_true = y_true(:);
    y_pred = y_pred(:);
    
    % Check if inputs are of the same length
    if length(y_true) ~= length(y_pred)
        error('Input vectors y_true and y_pred must be of the same length.');
    end
    
    % Check if inputs are numeric
    if ~isnumeric(y_true) || ~isnumeric(y_pred)
        error('Input vectors y_true and y_pred must be numeric.');
    end
    
    % Check if y_true has variance
    SS_tot = sum( (y_true - mean(y_true)).^2 );
    if SS_tot == 0
        warning('Total sum of squares (SS_tot) is zero. R^2 is set to NaN.');
        R2 = NaN;
        return;
    end
    
    % Calculate Residual Sum of Squares
    SS_res = sum( (y_true - y_pred).^2 );
    
    % Compute R^2
    R2 = 1 - (SS_res / SS_tot);
end