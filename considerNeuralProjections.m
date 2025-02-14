function considerNeuralProjections(cued_success_Response,uncued_success_Response,cued_failure_Response,uncued_failure_Response,cue_Response)

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
r=matchAllUnits(r);
cued_success_Response=r{1};
cued_failure_Response=r{2};
uncued_failure_Response=r{3};
uncued_success_Response=r{4};
cue_Response=r{5};

times_cuedsuccess=nanmean(cued_success_Response.unitbyunit_x,1)-2.86168;
times_cuedfailure=nanmean(cued_failure_Response.unitbyunit_x,1)-2.58234;
times_uncuedsuccess=nanmean(uncued_success_Response.unitbyunit_x,1)-6.58879;
times_uncuedfailure=nanmean(uncued_failure_Response.unitbyunit_x,1)-6.1889;

% First, consider "reinforcement"
% Get the "max coding directions" that represent
% 1. cued vs. uncued
% 2. outcome
% 3. given current behavior, cued or uncued reach on next trial
% get orthonormal basis spanning this same subspace
% then project neural activity onto this subspace
% could be a different subspace that contains "reinforcement" after a
% success or failure, so do these independently
w1=plotNeuralProjection(cued_success_Response.unitbyunit_y,uncued_success_Response.unitbyunit_y,times_cuedsuccess,times_uncuedsuccess,[1 5],[]);
w2=plotNeuralProjection(cued_success_Response.unitbyunit_y,cued_failure_Response.unitbyunit_y,times_cuedsuccess,times_cuedfailure,[1 5],[]);
w1(isnan(w1))=0;
w2(isnan(w2))=0;
U=gramschmidt([w1 w2]);
A=cued_success_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);

figure();


% CCA

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
