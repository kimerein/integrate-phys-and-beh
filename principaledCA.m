function principaledCA(data,dimensionNames,plotN,bootstrapNTimes)

if all(strcmp(dimensionNames,{'units','time','conditions'}))
    disp('Studying units X time X conditions');
else
    error('principaledCA.m currently requires data matrix structured as units X time X conditions');
end

% Normalize each unit's PSTH, don't min-subtract here bcz assume 0 is 0
data=normalizeData(data,[2 3],'sd'); % last argument either 'sd' or 'max'
% data=ZscoreData(data,[2 3]);
% data=smoothData(data,2,3); % last arg is gaussian smooth bin

% Flatten
flatData=flatten(data,'expand',[3 1]); % for temporal factors

% Display flattened data matrix
figure();
imagesc(flatData); title('flatData');

% Fill in nans and infs
flatData=noNansOrInfs(flatData);

% Eigenvalues and vectors
[eigVec_1,sorted_eigVal_1,eigVec_2,sorted_eigVal_2]=plotEigs(flatData,plotN,true);
bootstrapEigs(flatData,bootstrapNTimes,plotN);

% PCA
plotPCA(flatData,plotN);
plotPCA(flatData',plotN);

% Project onto eig space
% projectOntoEigSpace(flatData,2,eigVec_2(:,1:plotN),sorted_eigVal_2(1:plotN));

% TCA
% R_guess=5; % guess matrix rank
% allconditions_cpmodel=plotTCA(noNansOrInfs(data),10,R_guess);
R_guess=4; % guess matrix rank
allconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[2:5])),10,R_guess);
[allcell_PCs,allcell_archetypeCells,dimOrdering]=studyCPmodel(allconditions_cpmodel);
% get cue responses
[~,cuescores]=plotPCA(data(:,:,1),4);
% cueweights=ones(size(data,1),1);
cuePC1_archetypeCells=projectDirectionOntoCPmodel(cuescores(:,1),allconditions_cpmodel,allcell_PCs); 
% cueweights should be a vector, same length as units, with each unit's
% degree of cue responsiveness (or loading to cue factor)
cuePC2_archetypeCells=projectDirectionOntoCPmodel(cuescores(:,2),allconditions_cpmodel,allcell_PCs); 
tits={'cueSuccGrp1','cueFailGrp1','uncueSuccGrp1','uncueFailGrp1';...
      'cueSuccGrp2','cueFailGrp2','uncueSuccGrp2','uncueFailGrp2';...
      'cueSuccGrp3','cueFailGrp3','uncueSuccGrp3','uncueFailGrp3';...
      'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4'};
plotArchetypeCells(cat(4,allcell_archetypeCells,cuePC1_archetypeCells,cuePC2_archetypeCells),{'all','cuePC1','cuePC2'},tits);
% uncuedconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[4:5])),10,R_guess);

% just timepoints after outcome
% X=flatten(data(:,50:end,2:5),'mean',2)'; X=noNansOrInfs(X);
% [U,S,V]=svd(X);
% figure(); imagesc(U);
% figure(); imagesc(S(:,1:10));
% figure(); imagesc(V(1:4,:)');
% figure(); imagesc(V(:,1:4));

end

function plotArchetypeCells(acells,leg1,tits)

% expects dimensions of acells to be 
% temporal factors x trial factors x unit loadings PCs x cell groups
% dimension 4 is cell groups

figure(); set(gcf, 'Name', ['Sum of weighted tfs comparing cell group archetypes']);
cmap=colormap('jet');
colororder(cmap(1:floor((size(cmap,1)-1)/(size(acells,4))):size(cmap,1),:));
k=1;
for i=1:size(acells,3)
    for j=1:size(acells,2)
        subplot(size(acells,3),size(acells,2),k);
        plot(reshape(acells(:,j,i,:),size(acells,1),size(acells,4)));
        if k==1
            legend(leg1);
        end
        title(tits{i,j});
        k=k+1;
    end
end

end

function [TFsForUnitsAndTrialTypes,dimOrdering]=projectDirectionOntoCPmodel(cueweights,allconditions_cpmodel,allcell_PCs)

unitweights=allconditions_cpmodel.U{1}; tf=allconditions_cpmodel.U{2}; trialfactors=allconditions_cpmodel.U{3};

% first just project onto existing "archetypal neurons", i.e., onto the PCs
% from PCA of all neuron factors
inputToPCspace=repmat(cueweights,1,size(unitweights,2)).*unitweights;
newscore=(inputToPCspace-repmat(allcell_PCs.mu,size(allcell_PCs.score,1),1))*allcell_PCs.coeff;
figure(); scatter3(newscore(:,1),newscore(:,2),newscore(:,3)); set(gcf, 'Name', 'cueweights-transformed PC space');

% then do new PCA to get directions describing cue neuron responses
[pcaout.coeff,pcaout.score,pcaout.latent,pcaout.tsquared,pcaout.explained,pcaout.mu]=plotPCA(repmat(cueweights,1,size(unitweights,2)).*unitweights,size(unitweights,2));
coeff=pcaout.coeff;
figure(); set(gcf, 'Name', ['Sum of weighted temporal factors describing unit PCA -- CUE WEIGHTS']);
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(coeff,2)
    weightedTFforUnit=zeros(size(tf(:,1)));
    for j=1:size(tf,2)
        weightedTFforUnit=weightedTFforUnit+coeff(j,i)*tf(:,j)*allconditions_cpmodel.lambda(j);
    end
    plot(weightedTFforUnit); hold all;  
end

% for each trial type, get how archetypal neuron responses expected to change
[TFsForUnitsAndTrialTypes,dimOrdering]=archetypal_neuron_by_neuron(coeff,tf,trialfactors,allconditions_cpmodel);

end

function [pcaout,TFsForUnitsAndTrialTypes,dimOrdering]=studyCPmodel(allconditions_cpmodel)

% pca unit loadings and reconstruct top unit responses
unitweights=allconditions_cpmodel.U{1}; tf=allconditions_cpmodel.U{2};
[pcaout.coeff,pcaout.score,pcaout.latent,pcaout.tsquared,pcaout.explained,pcaout.mu]=plotPCA(unitweights,size(unitweights,2));
coeff=pcaout.coeff;
figure(); set(gcf, 'Name', ['Sum of weighted temporal factors describing unit PCA']);
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(coeff,2)
    weightedTFforUnit=zeros(size(tf(:,1)));
    for j=1:size(tf,2)
        weightedTFforUnit=weightedTFforUnit+coeff(j,i)*tf(:,j)*allconditions_cpmodel.lambda(j);
    end
    plot(weightedTFforUnit); hold all;  
end

% for each trial type, get average temporal factor
trialfactors=allconditions_cpmodel.U{3};
figure(); set(gcf, 'Name', ['Sum of weighted temporal factors describing each trial type averaged across neurons']);
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(trialfactors,1)
    weightedTFforTrialtype=zeros(size(tf(:,1)));
    for j=1:size(tf,2)
        weightedTFforTrialtype=weightedTFforTrialtype+trialfactors(i,j)*tf(:,j)*allconditions_cpmodel.lambda(j);
    end
    plot(weightedTFforTrialtype); hold all; 
end
 
% for each trial type, get how archetypal neuron responses expected to change
[TFsForUnitsAndTrialTypes,dimOrdering]=archetypal_neuron_by_neuron(coeff,tf,trialfactors,allconditions_cpmodel);

end

function [TFsForUnitsAndTrialTypes,dimOrdering]=archetypal_neuron_by_neuron(coeff,tf,trialfactors,allconditions_cpmodel)

TFsForUnitsAndTrialTypes=nan(length(tf(:,1)),size(trialfactors,1),size(coeff,2));
dimOrdering='temporal factors x trial factors x unit loadings PCs';
for k=1:size(trialfactors,1)
    for i=1:size(coeff,2)
        weightedTFforUnit=zeros(size(tf(:,1)));
        for j=1:size(tf,2)
            weightedTFforUnit=weightedTFforUnit+coeff(j,i)*tf(:,j)*trialfactors(k,j)*allconditions_cpmodel.lambda(j);
        end
        TFsForUnitsAndTrialTypes(:,k,i)=weightedTFforUnit; 
    end
end
figure(); set(gcf, 'Name', ['Sum of weighted temporal factors describing each trial type archetypal neuron by neuron']);
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(coeff,2)
    subplot(size(coeff,2),1,i);
    plot(reshape(TFsForUnitsAndTrialTypes(:,:,i),length(tf(:,1)),size(trialfactors,1)));
end

end

function matsum=tensorToMatrix(ktens)

if length(ktens.U)~=3
    error('current code tensorToMatrix in principaledCA.m expect 3D tensor');
end
a=ktens.U{1}; b=ktens.U{2}; c=ktens.U{3};
matsum=zeros(size(a,1),size(b,1),size(c,1));
for i=1:length(ktens.lambda)
    matsum=matsum+ktens.lambda(i)*outerProduct(outerProduct(a(:,i),b(:,i)),c(:,i));
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

function bestmodel=plotTCA(inputdata,n_fits,R)

% Fit CP Tensor Decomposition
% these commands require that you download Sandia Labs' tensor toolbox:
% http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
% data is neurons x times x trials
% convert data to a tensor object
data = tensor(inputdata);

% fit the cp decomposition from random initial guesses
% n_fits = 30;
err = zeros(n_fits,1);
% choose model with lowest error
bestmodel=[];
for n = 1:n_fits
    % fit model
    disp(n);
%     est_factors = cp_als(tensor(data),R,'printitn',0,'fixsigns',true,'maxiters',300);
    est_factors = cp_opt(tensor(data),R,'init','rand','lower',0,'printitn',0,'fixsigns',true,'maxiters',10000);
    % store error
    err(n) = norm(full(est_factors) - data)/norm(data);
    if all(err(1:n-1)>err(n))
        % best so far
        bestmodel=est_factors;
        besterr=err(n);
        % reconstruct initial matrix
        data_reconstruct=tensorToMatrix(est_factors);
    end

    % visualize fit for first several fits
    if n < 10
        % plot the estimated factors
        viz_ktensor(est_factors, ... 
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'})
        set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
        text(0.05,0.05,['lambdas ' num2str(est_factors.lambda')]);
    end
end

figure(); hold on;
plot(randn(n_fits,1), err, 'ob');
xlim([-10,10]);
ylim([0 1.0]);
set(gca,'xtick',[]);
ylabel('model error');

viz_ktensor(bestmodel, ...
    'Plottype', {'bar', 'line', 'scatter'}, ...
    'Modetitles', {'neurons', 'time', 'trials'})
set(gcf, 'Name', ['best fit w err ' num2str(besterr)])
text(0.05,0.05,['lambdas ' num2str(bestmodel.lambda')])

% temp=bestmodel.U{1};
% tsnetemp=tsne(temp,'Algorithm','exact','Distance','chebychev'); figure(); scatter(tsnetemp(:,1),tsnetemp(:,2));

end

function projectOntoEigSpace(A,whichdim,eigVecs,eigVals)

if whichdim==2
elseif whichdim==1
    A=A';
else
    error('expected 2D matrix for projectOntoEigSpace in principaledCA.m');
end
if size(A,2)~=size(eigVecs,1)
    error('array passed to projectOntoEigSpace size does not match eigVecs size in principaledCA.m');
end
% project row by row onto eig space
projections=nan(size(A,1),size(eigVecs,2));
for i=1:size(A,1)
    for j=1:size(eigVecs,2)
        projections(i,j)=dot(A(i,:),eigVecs(:,j)/norm(eigVecs(:,j)));
    end
end
scaleByEigVals=true;
if scaleByEigVals==true
    projections=projections./repmat(eigVals',size(A,1),1);
end
figure(); scatter(projections(:,1),projections(:,2)); xlabel('Project onto eig vec 1'); ylabel('Project onto eig vec 2');
figure(); scatter(projections(:,1),projections(:,3)); xlabel('Project onto eig vec 1'); ylabel('Project onto eig vec 3');
figure(); scatter(projections(:,2),projections(:,3)); xlabel('Project onto eig vec 2'); ylabel('Project onto eig vec 3');
figure(); scatter3(projections(:,1),projections(:,2),projections(:,3)); xlabel('Proj1'); ylabel('Proj2'); zlabel('Proj3');

end

function [eigVec_nbyn,sorted_eigVal_nbyn,eigVec_tbyt,sorted_eigVal_tbyt]=bootstrapEigs(flatData,bootstrapNTimes,plotN)

eigVec_nbyn=[]; sorted_eigVal_nbyn=[]; eigVec_tbyt=[]; sorted_eigVal_tbyt=[];
if bootstrapNTimes==1
    return
end

takeForBoot=0.7;
k=ceil(takeForBoot*size(flatData,2));
for i=1:bootstrapNTimes
    disp(i);
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
    disp(i);
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

if length(whichdim)>2
    data=data-repmat(mean(data,'all','omitnan'),[size(data,1) size(data,2) size(data,3)]);
    data=data./repmat(std(data,0,'all','omitnan'),[size(data,1) size(data,2) size(data,3)]);
elseif length(whichdim)>1
    me=mean(mean(data,whichdim(1),'omitnan'),whichdim(2),'omitnan');
    v=mean(var(data,0,whichdim(1),'omitnan'),whichdim(2),'omitnan');
    d=[1 2 3]; missingdim=d(~ismember(d,whichdim));
    if missingdim==1
        sz=[1 size(data,2) size(data,3)];
    elseif missingdim==2
        sz=[size(data,1) 1 size(data,3)];
    elseif missingdim==3
        sz=[size(data,1) size(data,2) 1];
    else
        error('principaledCA.m works only for 3D data matrices');
    end
    data=data-repmat(me,sz);
    data=data./repmat(sqrt(v),sz);
else
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

end

function data=smoothData(data,whichdim,gaussBin)

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
        data(:,i,j)=smoothdata(data(:,i,j),'gaussian',gaussBin);
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

function data=normalizeData(data,whichdim,bywhat)

if length(whichdim)>1
    % expand along dimension whichdim(2)
    % then max normalize along dimension whichdim(1)
    ed=flatten(data,'expand',[whichdim(2) whichdim(1)]);
    switch bywhat
        case 'max'
            ma=max(ed,[],whichdim(1),'omitnan');
        case 'sd'
            ma=std(ed,0,whichdim(1),'omitnan');
    end
    d=[1 2 3]; missingdim=d(~ismember(d,whichdim));
    if missingdim==1
        sz=[1 size(data,2) size(data,3)];
    elseif missingdim==2
        sz=[size(data,1) 1 size(data,3)];
    elseif missingdim==3
        sz=[size(data,1) size(data,2) 1];
    else
        error('principaledCA.m works only for 3D data matrices');
    end
    data=data./repmat(ma,sz);
else
    if whichdim==1
        sz=[size(data,1) 1 1];
    elseif whichdim==2
        sz=[1 size(data,2) 1];
    elseif whichdim==3
        sz=[1 1 size(data,3)];
    else
        error('principaledCA.m works only for 3D data matrices');
    end
    switch bywhat
        case 'max'
           ma=max(data,[],whichdim,'omitnan');
        case 'sd'
           ma=std(data,0,whichdim,'omitnan');
    end
    data=data./repmat(ma,sz);
end

end

function [coeff,score,latent,tsquared,explained,mu]=plotPCA(A,plotN)

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

figure(); scatter3(score(:,1),score(:,2),score(:,3));
xlabel('Proj1'); ylabel('Proj2'); zlabel('Proj3');

end

function [eigVec_nbyn,sorted_eigVal_nbyn,eigVec_tbyt,sorted_eigVal_tbyt]=plotEigs(A,plotN,doPlots)

if doPlots==true
    disp('Plotting eigenvals and vectors');
    disp('Size of matrix A');
    disp(size(A));
end
% E.g., size of matrix is N neurons X T "timepoints"
% center data both dims
% A=A-repmat(mean(A,1,'omitnan'),size(A,1),1); 
% A=A-repmat(mean(A,2,'omitnan'),1,size(A,2));
neurons_by_neurons=A*transpose(A);
times_by_times=transpose(A)*A;
% center data
neurons_by_neurons=neurons_by_neurons-mean(neurons_by_neurons,'all','omitnan');
times_by_times=times_by_times-mean(times_by_times,'all','omitnan');
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
