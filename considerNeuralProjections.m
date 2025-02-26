function considerNeuralProjections(cued_success_Response,uncued_success_Response,cued_failure_Response,uncued_failure_Response,cue_Response,cued_succ_to_cued,cued_succ_to_uncued,uncued_succ_to_cued,uncued_succ_to_uncued)

onlyOnD1tag=false;
onlyOnA2atag=true;

if onlyOnD1tag==true
    cued_success_Response.unitbyunit_y(cued_success_Response.D1tag(cued_success_Response.excluded==0)~=1,:)=0;
    uncued_success_Response.unitbyunit_y(uncued_success_Response.D1tag(uncued_success_Response.excluded==0)~=1,:)=0;
    cued_failure_Response.unitbyunit_y(cued_failure_Response.D1tag(cued_failure_Response.excluded==0)~=1,:)=0;
    uncued_failure_Response.unitbyunit_y(uncued_failure_Response.D1tag(uncued_failure_Response.excluded==0)~=1,:)=0;
    cue_Response.unitbyunit_y(cue_Response.D1tag(cue_Response.excluded==0)~=1,:)=0;
    cued_succ_to_cued.unitbyunit_y(cued_succ_to_cued.D1tag(cued_succ_to_cued.excluded==0)~=1,:)=0;
    cued_succ_to_uncued.unitbyunit_y(cued_succ_to_uncued.D1tag(cued_succ_to_uncued.excluded==0)~=1,:)=0;
    uncued_succ_to_cued.unitbyunit_y(uncued_succ_to_cued.D1tag(uncued_succ_to_cued.excluded==0)~=1,:)=0;
    uncued_succ_to_uncued.unitbyunit_y(uncued_succ_to_uncued.D1tag(uncued_succ_to_uncued.excluded==0)~=1,:)=0;
end
if onlyOnA2atag==true
    cued_success_Response.unitbyunit_y(cued_success_Response.A2atag(cued_success_Response.excluded==0)~=1,:)=0;
    uncued_success_Response.unitbyunit_y(uncued_success_Response.A2atag(uncued_success_Response.excluded==0)~=1,:)=0;
    cued_failure_Response.unitbyunit_y(cued_failure_Response.A2atag(cued_failure_Response.excluded==0)~=1,:)=0;
    uncued_failure_Response.unitbyunit_y(uncued_failure_Response.A2atag(uncued_failure_Response.excluded==0)~=1,:)=0;
    cue_Response.unitbyunit_y(cue_Response.A2atag(cue_Response.excluded==0)~=1,:)=0;
    cued_succ_to_cued.unitbyunit_y(cued_succ_to_cued.A2atag(cued_succ_to_cued.excluded==0)~=1,:)=0;
    cued_succ_to_uncued.unitbyunit_y(cued_succ_to_uncued.A2atag(cued_succ_to_uncued.excluded==0)~=1,:)=0;
    uncued_succ_to_cued.unitbyunit_y(uncued_succ_to_cued.A2atag(uncued_succ_to_cued.excluded==0)~=1,:)=0;
    uncued_succ_to_uncued.unitbyunit_y(uncued_succ_to_uncued.A2atag(uncued_succ_to_uncued.excluded==0)~=1,:)=0;
end

% What are the dimensions that interest me?
% cued reach vs. uncued reach
% success vs. failure
% reach vs. no reach
% initial cued response
% pre-cue baseline
% uncue baseline
% and interactions of these

r{1}=cued_success_Response;
r{2}=cued_failure_Response;
r{3}=uncued_failure_Response;
r{4}=uncued_success_Response;
r{5}=cue_Response;
r{6}=cued_succ_to_cued;
r{7}=cued_succ_to_uncued;
r{8}=uncued_succ_to_cued;
r{9}=uncued_succ_to_uncued;
r=matchAllUnits(r);
cued_success_Response=r{1};
cued_failure_Response=r{2};
uncued_failure_Response=r{3};
uncued_success_Response=r{4};
cue_Response=r{5};
cued_succ_to_cued=r{6};
cued_succ_to_uncued=r{7};
uncued_succ_to_cued=r{8};
uncued_succ_to_uncued=r{9};

times_cuedsuccess=nanmean(cued_success_Response.unitbyunit_x,1)-findZeroAtAlignComp(cued_success_Response);
times_cuedfailure=nanmean(cued_failure_Response.unitbyunit_x,1)-findZeroAtAlignComp(cued_failure_Response);
times_uncuedsuccess=nanmean(uncued_success_Response.unitbyunit_x,1)-findZeroAtAlignComp(uncued_success_Response);
times_uncuedfailure=nanmean(uncued_failure_Response.unitbyunit_x,1)-findZeroAtAlignComp(uncued_failure_Response);
times_cue=nanmean(cue_Response.unitbyunit_x,1)-findZeroAtAlignComp(cue_Response);
times_cued_succ_to_cued=nanmean(cued_succ_to_cued.unitbyunit_x,1)-findZeroAtAlignComp(cued_succ_to_cued);
times_cued_succ_to_uncued=nanmean(cued_succ_to_uncued.unitbyunit_x,1)-findZeroAtAlignComp(cued_succ_to_uncued);
times_uncued_succ_to_cued=nanmean(uncued_succ_to_cued.unitbyunit_x,1)-findZeroAtAlignComp(uncued_succ_to_cued);
times_uncued_succ_to_uncued=nanmean(uncued_succ_to_uncued.unitbyunit_x,1)-findZeroAtAlignComp(uncued_succ_to_uncued);

% First, consider "reinforcement"
% Get the "max coding directions" that represent
% 1. cued vs. uncued
% 2. outcome
% 3. given current behavior, cued or uncued reach on next trial
% get orthonormal basis spanning this same subspace
% then project neural activity onto this subspace
% could be a different subspace that contains "reinforcement" after a
% cued or uncued success, so do these independently
w1=plotNeuralProjection(cued_success_Response.unitbyunit_y,uncued_success_Response.unitbyunit_y,times_cuedsuccess,times_uncuedsuccess,[1 5],[]);
w2=plotNeuralProjection(cued_success_Response.unitbyunit_y,cued_failure_Response.unitbyunit_y,times_cuedsuccess,times_cuedfailure,[1 5],[]);
w3=plotNeuralProjection(cued_succ_to_cued.unitbyunit_y,cued_succ_to_uncued.unitbyunit_y,times_cued_succ_to_cued,times_cued_succ_to_uncued,[1 5],[]);
w1(isnan(w1))=0;
w2(isnan(w2))=0;
w3(isnan(w2))=0;
% U=gramschmidt([w1 w2 w3]);
U=[w1 w2 w3]; P=U*inv((U'*U))*U';
A=cued_success_Response.unitbyunit_y;
A(isnan(A))=0;
% A_proj=U*(U'*A);
A_proj=P*A;
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
figure();
plot(times_cuedsuccess,magnitudes); title('cued success projected onto cued_succ_high subspace'); hold all;
A=uncued_success_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
plot(times_uncuedsuccess,magnitudes);
A=cued_failure_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
plot(times_cuedfailure,magnitudes);
A=uncued_failure_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
plot(times_uncuedfailure,magnitudes);

w1=plotNeuralProjection(uncued_success_Response.unitbyunit_y,cued_success_Response.unitbyunit_y,times_uncuedsuccess,times_cuedsuccess,[1 5],[]);
w2=plotNeuralProjection(uncued_success_Response.unitbyunit_y,uncued_failure_Response.unitbyunit_y,times_uncuedsuccess,times_uncuedfailure,[1 5],[]);
w3=plotNeuralProjection(uncued_succ_to_uncued.unitbyunit_y,uncued_succ_to_cued.unitbyunit_y,times_uncued_succ_to_uncued,times_uncued_succ_to_cued,[1 5],[]);
w1(isnan(w1))=0;
w2(isnan(w2))=0;
w3(isnan(w2))=0;
U=gramschmidt([w1 w2 w3]);
A=uncued_success_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
figure();
plot(times_uncuedsuccess,magnitudes); title('uncued success projected onto uncued_succ_high subspace');

temp=cued_success_Response.unitbyunit_y;
% drop all nan rows
[R,P]=corr_timepoints(temp(~all(temp==0,2),:),25);
t=times_cuedsuccess;
f=find(~isnan(t),1,'first');
t(1:f-1)=t(f) - (t(f+1) - t(f)) * ((f-1):-1:1);
figure(); imagesc(t,t,R);

% CCA

end

function plot_imagesc(R,timerange)

figure(); imagesc(t,t,R);
% Get the current axes handle
ax = gca;
% Get current tick positions (these correspond to data indices)
xticks_current = ax.XTick;
yticks_current = ax.YTick;
% Define a new range for the x- and y-axes. For instance,
% map the x-axis from 1 to 10 into a new range [0, 100],
% and the y-axis from 1 to 10 into a new range [0, 50].
new_x_range = linspace(timerange(1), timerange(2), numel(xticks_current));
new_y_range = linspace(timerange(1), timerange(2), numel(yticks_current));
% Update the tick labels (this changes only the displayed numbers)
ax.XTickLabel = new_x_range;
ax.YTickLabel = new_y_range;

end

function [R,P] = corr_timepoints(data,ds)
% CORR_TIMEPOINTS computes the correlation between timepoints based on the
% neural population activity.
%
%   R = corr_timepoints(data)
%
%   Input:
%       data - an array of size [neurons x timepoints] where each row is the
%              time series of a neuron.
%
%   Steps:
%       1. Z-score each neuron's activity across time (i.e., along rows).
%       2. Compute the correlation coefficient matrix among timepoints by
%          treating each column as a "population vector" (i.e., activity across neurons).
%
%   Output:
%       R - a [timepoints x timepoints] correlation matrix.

    % Z-score each row (neuron) across timepoints
    data_z=data-repmat(mean(data,2,'omitnan'),1,size(data,2));
    data_z=data_z./repmat(std(data_z,0,2,'omitnan'),1,size(data_z,2));

    % Downsample
    data_z=downSampMatrix(data_z,ds);
    data_z(isnan(data_z))=0;
    
    [R,P] = corrcoef(data_z);
end


function U = gramschmidt(V)
    [n, k] = size(V);
    U = zeros(n,k);
    U(:,1) = V(:,1) / norm(V(:,1));
    for i = 2:k
        U(:,i) = V(:,i);
        for j = 1:i-1
            U(:,i) = U(:,i) - (U(:,j)'*U(:,i)) * U(:,j);
        end
        U(:,i) = U(:,i) / norm(U(:,i));
    end
end

function shiftby=findZeroAtAlignComp(response)

[~,ma]=max(mean(response.aligncomp_y,1,'omitnan'),[],2,'omitnan');
temp=mean(response.aligncomp_x,1,'omitnan');
shiftby=temp(ma);

end
