function considerNeuralProjections(cued_success_Response,uncued_success_Response,cued_failure_Response,uncued_failure_Response,cue_Response,cued_succ_to_cued,cued_succ_to_uncued,uncued_succ_to_cued,uncued_succ_to_uncued)

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

times_cuedsuccess=nanmean(cued_success_Response.unitbyunit_x,1)-2.86168;
times_cuedfailure=nanmean(cued_failure_Response.unitbyunit_x,1)-2.58234;
times_uncuedsuccess=nanmean(uncued_success_Response.unitbyunit_x,1)-6.58879;
times_uncuedfailure=nanmean(uncued_failure_Response.unitbyunit_x,1)-6.1889;
times_cue=nanmean(cue_Response.unitbyunit_x,1)-a;
times_cued_succ_to_cued=nanmean(cued_succ_to_cued.unitbyunit_x,1)-a;
times_cued_succ_to_uncued=nanmean(cued_succ_to_uncued.unitbyunit_x,1)-a;
times_uncued_succ_to_cued=nanmean(uncued_succ_to_cued.unitbyunit_x,1)-a;
times_uncued_succ_to_uncued=nanmean(uncued_succ_to_uncued.unitbyunit_x,1)-a;

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
U=gramschmidt([w1 w2 w3]);
A=cued_success_Response.unitbyunit_y;
A(isnan(A))=0;
A_proj=U*(U'*A);
magnitudes=sqrt(sum(A_proj.^2, 1,'omitnan'));  % vector of norms
figure();
plot(times_cuedsuccess,magnitudes); title('cued success projected onto cued_succ_high subspace');

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
