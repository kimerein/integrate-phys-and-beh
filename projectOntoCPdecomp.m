function trialWeightsOntoTypeFactors=projectOntoCPdecomp(loc,test_data_matrix,whichFactors)

% dsby=225;
% test_data_matrix must be in format neurons X times X trials

% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\allconditions_cpmodel.mat');
load(loc);
trialWeightsOntoTypeFactors=nan(length(whichFactors),size(allconditions_cpmodel.U{3},1)*size(test_data_matrix,1),size(test_data_matrix,3));
% trialWeightsOntoTypeFactors=nan(length(whichFactors),size(allconditions_cpmodel.U{3},1)*size(test_data_matrix,2),size(test_data_matrix,3));
% trialWeightsOntoTypeFactors=nan(length(whichFactors),size(allconditions_cpmodel.U{3},1)*((size(test_data_matrix,2)/dsby)*size(test_data_matrix,1)),size(test_data_matrix,3));
% trialWeightsOntoTypeFactors=nan(length(whichFactors),1*((size(test_data_matrix,2)/dsby)*size(test_data_matrix,1)),size(test_data_matrix,3));
% trialWeightsOntoTypeFactors is in format factor X (trial type X neuron) X trial 
for i=1:length(whichFactors)
    whichFactor=whichFactors(i);
    temp=allconditions_cpmodel.U{1}; fac1_vec1=temp(:,whichFactor);
    temp=allconditions_cpmodel.U{2}; fac1_vec2=temp(:,whichFactor);
    temp=allconditions_cpmodel.U{3}; fac1_vec3=temp(:,whichFactor);
    % Neuron factors do not match
    fac1_vec1=ones(size(test_data_matrix,1),1);
    fac1_ktens=ktensor({fac1_vec1,fac1_vec2,fac1_vec3});
    fac1=outerProduct(outerProduct(fac1_vec1,fac1_vec2),fac1_vec3);
    % Project each trial onto this factor, save loading onto trial type
    % factors
    for j=1:size(test_data_matrix,3)
        for k=1:size(fac1,3)
            projected=(test_data_matrix(:,:,j).*fac1(:,:,k))./norm(fac1_ktens);
%             projected=(test_data_matrix(:,:,j).*nanmean(fac1,3))./norm(fac1_ktens);
            % collapse across time but not across neurons
            curbloinds=(k-1)*size(test_data_matrix,1)+1:(k-1)*size(test_data_matrix,1)+size(test_data_matrix,1);
            trialWeightsOntoTypeFactors(i,curbloinds,j)=nansum(projected,2);
            % collapse across neurons but not across time
%             curbloinds=(k-1)*size(test_data_matrix,2)+1:(k-1)*size(test_data_matrix,2)+size(test_data_matrix,2);
%             trialWeightsOntoTypeFactors(i,curbloinds,j)=nansum(projected,1);
            % no collapse, retain neurons X time
%             curbloinds=(k-1)*((size(test_data_matrix,2)/dsby)*size(test_data_matrix,1))+1:(k-1)*((size(test_data_matrix,2)/dsby)*size(test_data_matrix,1))+((size(test_data_matrix,2)/dsby)*size(test_data_matrix,1));
%             % down samp
%             projected=downSampMatrix(projected,dsby);
%             projected=projected';
%             trialWeightsOntoTypeFactors(i,curbloinds,j)=projected(1:end);

%             ptens=tensor(projected);
%             T=hosvd(ptens,sqrt(3e-1),'rank',[1 1 1]);
%             neuron_loadings=T.U{1}./T.core(1,1,1);
%             figure(); bar(T.U{1}./T.core(1,1,1),'r');
%             figure(); subplot(1,3,1); plot(fac1_vec1,'Color','k'); hold on; plot(T.U{1}./T.core(1,1,1),'Color','r'); legend({'existing CP','new data'});
%             subplot(1,3,2); plot(fac1_vec2,'Color','k'); hold on; plot(T.U{2},'Color','r');
%             subplot(1,3,3); plot(fac1_vec3,'Color','k'); hold on; plot(T.U{3},'Color','r');
        end
    end
end

end

function out=outerProduct(vec1,vec2)

if size(vec1,2)==1 && size(vec2,2)==1
    out=vec1*transpose(vec2);
    if size(out,1)~=length(vec1) || size(out,2)~=length(vec2)
        error('problem in outerProduct in principaledCA.m');
    end
elseif length(size(vec1))==2 && size(vec2,2)==1
    out=nan(size(vec1,1),size(vec1,2),length(vec2));
    for i=1:size(vec1,1)
        for j=1:size(vec1,2)
            for k=1:length(vec2)
                out(i,j,k)=vec1(i,j)*vec2(k);
            end
        end
    end
    if size(out,1)~=size(vec1,1) || size(out,2)~=size(vec1,2) || size(out,3)~=length(vec2)
        error('problem in outerProduct in principaledCA.m');
    end
else
    error('outerProduct in principaledCA.m not implemented for inputs of these sizes');
end

end