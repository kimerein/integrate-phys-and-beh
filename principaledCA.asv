function principaledCA(data,dimensionNames,plotN,bootstrapNTimes)

if all(strcmp(dimensionNames,{'units','time','conditions'}))
    disp('Studying units X time X conditions');
else
    error('principaledCA.m currently requires data matrix structured as units X time X conditions');
end

% Normalize each unit's PSTH, don't min-subtract here bcz assume 0 is 0
data=normalizeData(data,2);
% data=ZscoreData(data,2);
% data=smoothData(data,2);

% Flatten
flatData=flatten(data,'mean',3);

% Display flattened data matrix
figure();
imagesc(flatData); title('flatData');

% Fill in nans and infs
flatData=noNansOrInfs(flatData);

% Eigenvalues and vectors
plotEigs(flatData,plotN,true);
bootstrapEigs(flatData,bootstrapNTimes,plotN);

% PCA
plotPCA(flatData,plotN);
plotPCA(flatData',plotN);

% just timepoints after outcome
X=flatten(data(:,50:end,2:5),'mean',2)'; X=noNansOrInfs(X);
[U,S,V]=svd(X);
figure(); imagesc(U);
figure(); imagesc(S(:,1:10));
figure(); imagesc(V(1:4,:)');
figure(); imagesc(V(:,1:4));



end

function [eigVec_nbyn,sorted_eigVal_nbyn,eigVec_tbyt,sorted_eigVal_tbyt]=bootstrapEigs(flatData,bootstrapNTimes,plotN)

eigVec_nbyn=[]; sorted_eigVal_nbyn=[]; eigVec_tbyt=[]; sorted_eigVal_tbyt=[];
if bootstrapNTimes==1
    return
end

takeForBoot=0.7;
k=ceil(takeForBoot*size(flatData,2));
for i=1:bootstrapNTimes
    [eigVec,val]=plotEigs(flatData(:,randsample(size(flatData,2),k)),plotN,false);
    if i==1
        eigVec_nbyn=nan(size(eigVec,1),size(eigVec,2),bootstrapNTimes);
        sorted_eigVal_nbyn=nan(length(val),bootstrapNTimes);
    end
    eigVec_nbyn(:,:,i)=eigVec;
    sorted_eigVal_nbyn(:,i)=val;
end
k=ceil(takeForBoot*size(flatData,1));
for i=1:bootstrapNTimes
    [~,~,eigVec,val]=plotEigs(flatData(randsample(size(flatData,1),k),:),plotN,false);
    if i==1
        eigVec_tbyt=nan(size(eigVec,1),size(eigVec,2),bootstrapNTimes);
        sorted_eigVal_tbyt=nan(length(val),bootstrapNTimes);
    end
    eigVec_tbyt(:,:,i)=eigVec;
    sorted_eigVal_tbyt(:,i)=val;
end

figure();
if plotN>size(eigVec_nbyn,2)
    plotTo=size(eigVec_nbyn,2);
else
    plotTo=plotN;
end
for i=1:plotTo
    subplot(plotN,1,i);
    plot(mean(eigVec_nbyn(:,i,:),3)); hold on;
    plot(mean(eigVec_nbyn(:,i,:),3)-std(eigVec_nbyn(:,i,:),0,3)); plot(mean(eigVec_nbyn(:,i,:),3)+std(eigVec_nbyn(:,i,:),0,3));
    text(1,eigVec_nbyn(1,i),num2str(mean(sorted_eigVal_nbyn(i,:))));
    if i==1
        title('eigs of A * transposeA');
    end
end
figure(); scatter(1:size(sorted_eigVal_nbyn(1:plotTo,:),1),mean(sorted_eigVal_nbyn(1:plotTo,:),2)); hold on;
for i=1:plotTo
    line([i i],[mean(sorted_eigVal_nbyn(i,:))-std(sorted_eigVal_nbyn(i,:)) mean(sorted_eigVal_nbyn(i,:))+std(sorted_eigVal_nbyn(i,:))]);
end
title('eig vals of A * transposeA');

figure();
if plotN>size(eigVec_tbyt,2)
    plotTo=size(eigVec_tbyt,2);
else
    plotTo=plotN;
end
for i=1:plotTo
    subplot(plotN,1,i);
    plot(mean(eigVec_tbyt(:,i,:),3)); hold on;
    plot(mean(eigVec_tbyt(:,i,:),3)-std(eigVec_tbyt(:,i,:),0,3)); plot(mean(eigVec_tbyt(:,i,:),3)+std(eigVec_tbyt(:,i,:),0,3));
    text(1,eigVec_tbyt(1,i),num2str(mean(sorted_eigVal_tbyt(i,:))));
    if i==1
        title('eigs of transposeA * A');
    end
end
figure(); scatter(1:size(sorted_eigVal_tbyt(1:plotTo,:),1),mean(sorted_eigVal_tbyt(1:plotTo,:),2)); hold on;
for i=1:plotTo
    line([i i],[mean(sorted_eigVal_tbyt(i,:))-std(sorted_eigVal_tbyt(i,:)) mean(sorted_eigVal_tbyt(i,:))+std(sorted_eigVal_tbyt(i,:))]);
end
title('eig vals of transposeA * A');

end

function data=ZscoreData(data,whichdim)

if whichdim==1
    sz=[size(data,1) 1 1];
elseif whichdim==2
    sz=[1 size(data,2) 1];
elseif whichdim==3
    sz=[1 1 size(data,3)];
else
    error('principaledCA.m works only for 3D data matrices');
end
data=data-repmat(mean(data,whichdim,'omitnan'),sz);
data=data./repmat(std(data,0,whichdim,'omitnan'),sz);

end

function data=smoothData(data,whichdim)

if whichdim==1
    s=[1 2 3];
elseif whichdim==2
    s=[2 1 3];
elseif whichdim==3
    s=[3 1 2];
else
    error('principaledCA.m works only for 3D data matrices');
end
origdataorder=[1 2 3];
data=permute(data,s);
origdataorder=origdataorder(s);
for i=1:size(data,2)
    for j=1:size(data,3)
        data(:,i,j)=smoothdata(data(:,i,j),'gaussian',3);
    end
end
[~,si]=sort(origdataorder);
data=permute(data,si);

end

function flatData=flatten(data,method,dim)

switch method
    case 'mean'
        flatData=mean(data,dim,'omitnan');
        sz=size(flatData);
        flatData=reshape(flatData,sz(find(sz>1,1,'first')),sz(find(sz>1,1,'last')));
    case 'max'
        flatData=max(data,[],dim,'omitnan');
        sz=size(flatData);
        flatData=reshape(flatData,sz(find(sz>1,1,'first')),sz(find(sz>1,1,'last')));
    case 'expand'
        if dim(1)==1 
            s=[2 3 1];
        elseif dim(1)==2
            s=[1 3 2];
        elseif dim(1)==3
            s=[1 2 3];
        else
            error('principaledCA.m works only for 3D data matrices');
        end
        data=permute(data,s);
        flatData=[];
        for i=1:size(data,3)
            if dim(2)==1
                flatData=[flatData; data(:,:,i)];
            elseif dim(2)==2
                flatData=[flatData data(:,:,i)];
            else
                error('failed to flatten');
            end
        end        
end

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
if plotN>size(coeff,2)
    plotTo=size(coeff,2);
else
    plotTo=plotN;
end
figure();
for i=1:plotTo
    subplot(plotN,1,i);
    plot(coeff(:,i));
    if i==1
        title('PCs');
    end
    text(1,coeff(1,i),num2str(latent(i)));
end

end

function [eigVec_nbyn,sorted_eigVal_nbyn,eigVec_tbyt,sorted_eigVal_tbyt]=plotEigs(A,plotN,doPlots)

if doPlots==true
    disp('Plotting eigenvals and vectors');
    disp('Size of matrix A');
    disp(size(A));
end
% E.g., size of matrix is N neurons X T "timepoints"
% center data both dims
A=A-repmat(mean(A,1,'omitnan'),size(A,1),1);
A=A-repmat(mean(A,2,'omitnan'),1,size(A,2));
neurons_by_neurons=A*transpose(A);
times_by_times=transpose(A)*A;
% eigenvectors and eigenvalues
[eigVec_nbyn,eigVal_nbyn]=eig(neurons_by_neurons);
[sorted_eigVal_nbyn,si]=sort(diag(eigVal_nbyn),'descend');
eigVec_nbyn=eigVec_nbyn(:,si);
[eigVec_tbyt,eigVal_tbyt]=eig(times_by_times);
[sorted_eigVal_tbyt,si]=sort(diag(eigVal_tbyt),'descend');
eigVec_tbyt=eigVec_tbyt(:,si);

if doPlots==true
    figure();
    if plotN>size(eigVec_nbyn,2)
        plotTo=size(eigVec_nbyn,2);
    else
        plotTo=plotN;
    end
    for i=1:plotTo
        subplot(plotN,1,i);
        plot(eigVec_nbyn(:,i));
        text(1,eigVec_nbyn(1,i),num2str(sorted_eigVal_nbyn(i)));
        if i==1
            title('eigs of A * transposeA');
        end
    end

    figure();
    if plotN>size(eigVec_tbyt,2)
        plotTo=size(eigVec_tbyt,2);
    else
        plotTo=plotN;
    end
    for i=1:plotTo
        subplot(plotN,1,i);
        plot(eigVec_tbyt(:,i));
        text(1,eigVec_tbyt(1,i),num2str(sorted_eigVal_tbyt(i)));
        if i==1
            title('eigs of transposeA * A');
        end
    end
end

end
