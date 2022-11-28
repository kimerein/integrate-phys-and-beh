function principaledCA(data,dimensionNames,plotN)

if all(strcmp(dimensionNames,{'units','time','conditions'}))
    disp('Studying units X time X conditions');
else
    error('principaledCA.m currently requires data matrix structured as units X time X conditions');
end

% Normalize each unit's PSTH, don't min-subtract here bcz assume 0 is 0
data=normalizeData(data,2);

% Flatten
flatData=flatten(data,'mean',3);

% Fill in nans and infs
flatData=noNansOrInfs(flatData);

% Eigenvalues and vectors
plotEigs(flatData,plotN);

% PCA
plotPCA(flatData,plotN);

% CCA

end

function flatData=flatten(data,method,dim)

switch method
    case 'mean'
        flatData=mean(data,dim,'omitnan');
    case 'expand'
end
sz=size(flatData);
flatData=reshape(flatData,sz(find(sz>1,1,'first')),sz(find(sz>1,1,'last')));

end

function data=noNansOrInfs(data)

data(isnan(data))=mean(data(isfinite(data)),'all');
data(isinf(data))=mean(data(isfinite(data)),'all');

end

function data=normalizeData(data,whichdim)

if whichdim==1
    sz=[size(data,1) 1 1];
elseif whichdim==2
    sz=[1 size(data,2) 1];
elseif whichdim==3
    sz=[1 1 size(data,3)];
else
    error('principaledCA.m works only for 3D data matrices');
end
data=data./repmat(max(data,[],whichdim,'omitnan'),sz);

end

function plotPCA(A,plotN)

[coeff,score,latent,tsquared,explained,mu]=pca(A);
% plot
figure();
for i=1:plotN
    subplot(plotN,1,i);
    plot(coeff(:,i));
end
title('PCs');

end

function plotEigs(A,plotN)

disp('Plotting eigenvals and vectors');
% E.g., size of matrix is N neurons X T "timepoints"
% center data both dims
A=A-repmat(mean(A,1,'omitnan'),size(A,1),1);
A=A-repmat(mean(A,2,'omitnan'),1,size(A,2));
disp('Size of matrix A');
disp(size(A));
neurons_by_neurons=A*transpose(A);
times_by_times=transpose(A)*A;
% eigenvectors and eigenvalues
[eigVec_nbyn,eigVal_nbyn]=eig(neurons_by_neurons);
[sorted_eigVal_nbyn,si]=sort(diag(eigVal_nbyn),'descend');
eigVec_nbyn=eigVec_nbyn(:,si);
[eigVec_tbyt,eigVal_tbyt]=eig(times_by_times);
[sorted_eigVal_tbyt,si]=sort(diag(eigVal_tbyt),'descend');
eigVec_tbyt=eigVec_tbyt(:,si);

figure(); 
for i=1:plotN
    subplot(plotN,1,i);
    plot(eigVec_nbyn(:,i));
end
title('eigs of A * transposeA');

figure(); 
for i=1:plotN
    subplot(plotN,1,i);
    plot(eigVec_tbyt(:,i));
end
title('eigs of transposeA * A');

end
