function [idx,out_cueresponse]=principaledCA(data,dimensionNames,plotN,bootstrapNTimes,doingHighOrLowRank)

% doingHighOrLowRank='high';

if all(strcmp(dimensionNames,{'units','time','conditions'}))
    disp('Studying units X time X conditions');
else
    error('principaledCA.m currently requires data matrix structured as units X time X conditions');
end

switch doingHighOrLowRank
    case 'high'
        backupdata=data;

%         smoodata=backupdata;
%         tensie=smoodata(:,:,2:5);
%         tensie(tensie>10)=10;
%         tensie=noNansOrInfs(tensie);
%         whichFac=1;
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_0.mat']); facvec1=ones(size(factor(:,whichFac)));
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_1.mat']); facvec2=factor(:,whichFac);
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_2.mat']); facvec3=factor(:,whichFac);
%         [T,neuron_loadings_fac1]=projectCurrentDataOntoExistingCP_justvecs(facvec1,facvec2,facvec3,tensie);
%         whichFac=2;
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_0.mat']); facvec1=ones(size(factor(:,whichFac)));
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_1.mat']); facvec2=factor(:,whichFac);
%         load(['C:\Users\sabatini\Documents\currtens train at least 20 trials\factor_2.mat']); facvec3=factor(:,whichFac);
%         [T,neuron_loadings_fac2]=projectCurrentDataOntoExistingCP_justvecs(facvec1,facvec2,facvec3,tensie);

        data=getTrialTypeDependentResidual(data);
        % Normalize each unit's PSTH, don't min-subtract here bcz assume 0 is 0
        data=normalizeData(data,[2 3],'sd'); % last argument either 'sd', 'max' or 'mean'
        % Remove outliers, only necessary if norm by the mean
        % data=rmOutliersFromNormedPSTH(data);
        % data=ZscoreData(data,[2 3]);
        % data=smoothData(data,2,10); % last arg is gaussian smooth bin

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
        R_guess=6; % guess matrix rank
        allconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[2:5])),20,R_guess);
        [allcell_PCs,allcell_archetypeCells,dimOrdering]=studyCPmodel(allconditions_cpmodel);
        % project onto existing CP model
        %test_data_matrix=noNansOrInfs(data(:,:,[2:5]));
        %loc='Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\TCA\allconditions_cpmodel.mat';
        %whichFactor=1; [T,temp]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        % meanMinusTemp=nan(length(temp),R_guess); meanMinusTemp(:,1)=temp;
        %whichFactor=2; [T,meanMinusTemp(:,2)]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        %whichFactor=3; [T,meanMinusTemp(:,3)]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        %whichFactor=4; [T,meanMinusTemp(:,4)]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        %whichFactor=5; [T,meanMinusTemp(:,5)]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        %whichFactor=6; [T,meanMinusTemp(:,6)]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor);
        %load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\sortedbyneuronloadingsPC1.mat');
        %meanMinusTemp=meanMinusTemp-repmat(nanmean(meanMinusTemp,2),1,size(meanMinusTemp,2));
        %plotLikeTrainingSet(meanMinusTemp,sipc);
    case 'low'
        backupdata=data;
        % Normalize each unit's PSTH, don't min-subtract here bcz assume 0 is 0
        data=normalizeData(data,[2 3],'sd'); % last argument either 'sd', 'max' or 'mean'
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
        R_guess=4; % guess matrix rank
        allconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[2:3])),20,R_guess);
        [allcell_PCs,allcell_archetypeCells,dimOrdering]=studyCPmodel(allconditions_cpmodel);
end

temp=allconditions_cpmodel.U{1}; meanMinusTemp=temp-repmat(nanmean(temp,2),1,size(temp,2));
[~,si]=sort(meanMinusTemp(:,1));
[~,sipc]=sort(allcell_PCs.score(:,1)); % alternate sorting
ngroups=2;
currc={[0, 0.75, 0.75], [0.4940, 0.1840, 0.5560]}; %currc={'k','r','c','g','b'};
%Y=pdist(meanMinusTemp);
%Z=linkage(Y);
%idx=cluster(Z,'maxclust',3);
idx=kmeans(meanMinusTemp,ngroups,'Replicates',50);
% replicate idx
idx_replicates=nan(length(idx),100);
for i=1:100
    temp=kmeans(meanMinusTemp,ngroups,'Replicates',5);
    if nansum(temp==idx)/length(idx)>0.5
        % same label
    else % flip labels
        newtemp=temp; newtemp(temp==1)=2; newtemp(temp==2)=1; temp=newtemp;
    end
    idx_replicates(:,i)=temp;
end
figure();
for i=1:size(meanMinusTemp,2)
subplot(1,size(meanMinusTemp,2),i);
bar(meanMinusTemp(sipc,i),'k');
thesevals=meanMinusTemp(sipc,i);
hold on;
for j=1:ngroups
scatter(find(idx(sipc)==j),thesevals(idx(sipc)==j),2,currc{j});
end
end
figure(); scatter3(allcell_PCs.score(idx==1,2),allcell_PCs.score(idx==1,3),allcell_PCs.score(idx==1,4),[],'k'); hold on;
scatter3(allcell_PCs.score(idx==2,2),allcell_PCs.score(idx==2,3),allcell_PCs.score(idx==2,4),[],'r');
xlabel('PC2'); ylabel('PC3'); zlabel('PC4'); % PC2 separates clusters

load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\for_data_matrix_D1vA2a.mat')
figure();
for i=1:size(meanMinusTemp,2)
    subplot(1,size(meanMinusTemp,2),i);
    bar(meanMinusTemp(sipc,i));
    thesevals=meanMinusTemp(sipc,i);
    hold on;
    scatter(find(D1tag(sipc)==1),thesevals(find(D1tag(sipc)==1)),2,'r');
    scatter(find(A2atag(sipc)==1),thesevals(find(A2atag(sipc)==1)),2,'b');
end

% Is there a way to regress off the trial type-independent turning off
% after the reach? 
% 1. Could choose only cells that significantly represent the success v
% failure or cued v uncued
% 2. But is there a way to get the trial type-dependent component for each
% cell? How about subtracting the averaged across all cells, averaged across all trials, 
% Z-scored rate, and then only consider deviations from this average?
score=allcell_PCs.score;
[n,x]=histcounts(score(idx==1,1),[-0.06:0.005:0.1]+0.0025); [n,x]=cityscape_hist(n,x); figure(); plot(x,n./nansum(n),'Color',currc{1});
[n,x]=histcounts(score(idx==2,1),[-0.06:0.005:0.1]+0.0025); [n,x]=cityscape_hist(n,x); hold on; plot(x,n./nansum(n),'Color',currc{2});
figure(); scatter3(score(:,1),score(:,2),score(:,3),[],[0.5 0.5 0.5]); hold on; scatter3(score(D1tag==1,1),score(D1tag==1,2),score(D1tag==1,3),[],'r','filled');
scatter3(score(A2atag==1,1),score(A2atag==1,2),score(A2atag==1,3),[],'b','filled');
xlabel('Proj1'); ylabel('Proj2'); zlabel('Proj3');
% tsnetemp=tsne(score,'Algorithm','exact','Distance','chebychev','Exaggeration',10,'NumDimensions',2,'Perplexity',60,'Standardize',false);
tsnetemp=tsne(score,'Algorithm','exact','Distance','chebychev','Exaggeration',10,'NumDimensions',2,'Perplexity',900,'Standardize',false);
load('Z:\MICROSCOPE\Kim\old batch Physiology Final Data Sets\training\CP model\data_loc_array.mat'); [~,~,uic]=unique(data_loc_array(:,2));
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\cued_success_Response.mat'); mousebymouse=uic(cued_success_Response.fromWhichSess); 
[~,~,umbym]=unique(mousebymouse);
cmap=[]; for i=1:length(unique(umbym)) cmap=[cmap; rand(1,3)]; end %cmap=colormap(parula(108)); 
figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],cmap(umbym,:),'filled'); title('Control'); 
figure(); scatter(1:length(unique(umbym)),ones(size(1:length(unique(umbym)))),[],cmap,'filled');
% cmap=colormap('jet');
% ccft=cmap((idx-1)*100+1,:);
cmap=[0, 0.75, 0.75; 0.4940, 0.1840, 0.5560];
ccft=cmap(idx,:); 
adata=nansum(idx_replicates==1,2)./size(idx_replicates,2); adata2=nansum(idx_replicates==2,2)./size(idx_replicates,2);
adata(idx==1)=adata(idx==1); adata(idx==2)=adata2(idx==2);
figure(); s=scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft,'filled'); s.AlphaData=adata; s.MarkerFaceAlpha='flat';
figure(); scatter(tsnetemp(D1tag==1,1),tsnetemp(D1tag==1,2),[],repmat(cmap(1,:),nansum(D1tag==1),1)); hold on;
scatter(tsnetemp(A2atag==1,1),tsnetemp(A2atag==1,2),[],repmat(cmap(100,:),nansum(A2atag==1),1));
tempiedata=normalizeData(backupdata,[2 3],'sd'); currcuescore=reshape(max(tempiedata(:,1:150,1),[],2,'omitnan'),size(backupdata,1),1);
currcuescore=log(currcuescore); currcuescore(currcuescore<-2)=-2; figure(); histogram(currcuescore,200); 
out_cueresponse=currcuescore;
cuecolsfortsne=currcuescore-min(currcuescore,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne./max(cuecolsfortsne,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne*255;
cuecolsfortsne=ceil(cuecolsfortsne); cuecolsfortsne(cuecolsfortsne==0)=1; dontshow(isnan(cuecolsfortsne))=true;  cuecolsfortsne(isnan(cuecolsfortsne))=1;
cmap=colormap('jet');
ccft=cmap(cuecolsfortsne,:); ccft(dontshow==true,:)=repmat([1 1 1],sum(dontshow==true),1);
figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft); title('cue responsive');

% currcuescore=cuez;
% currcuescore(currcuescore<-2)=-2; currcuescore(currcuescore>4)=4;
% [~,sicc]=sort(currcuescore,'ascend'); r=1:length(currcuescore); r(sicc)=r; currcuescore=r;
% cuecolsfortsne=currcuescore-min(currcuescore,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne./max(cuecolsfortsne,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne*255;
% cuecolsfortsne=ceil(cuecolsfortsne); cuecolsfortsne(cuecolsfortsne==0)=1; dontshow(isnan(cuecolsfortsne))=true;  cuecolsfortsne(isnan(cuecolsfortsne))=1;
% cmap=colormap('jet');
% ccft=cmap(cuecolsfortsne,:); ccft(dontshow==true,:)=repmat([1 1 1],sum(dontshow==true),1);
% figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft,'filled'); title('cue responsive');

return

response_forThisArchetype(allconditions_cpmodel,idx==1);
response_forThisArchetype(allconditions_cpmodel,idx==2);
% Subdivide along cue dimension
% response_forThisArchetype(allconditions_cpmodel,idx==2 & (meanMinusTemp(:,5)>meanMinusTemp(:,6))); % cued
% response_forThisArchetype(allconditions_cpmodel,idx==2 & (meanMinusTemp(:,6)>meanMinusTemp(:,5))); % uncued
% response_forThisArchetype(allconditions_cpmodel,idx==1 & (meanMinusTemp(:,2)>meanMinusTemp(:,3))); % cued
% response_forThisArchetype(allconditions_cpmodel,idx==1 & (meanMinusTemp(:,3)>meanMinusTemp(:,2))); % uncued


temp=allcell_PCs.score(:,1); figure(); histogram(temp,400);
response_forThisArchetype(allconditions_cpmodel,temp<-0.01532);
response_forThisArchetype(allconditions_cpmodel,temp>=-0.01532);
% get cue responses, need to dump first cue PC, because just picks up cells
% that spike a lot generally, second cue PC is better
[~,cuescores]=plotPCA(data(:,:,1),6);
cuescores=cuescores-repmat(nanmean(cuescores(:,1:6),2),1,size(cuescores,2));
% cueweights=ones(size(data,1),1);
% drop first PC because represents DC offset?
temp=allcell_PCs.score(:,1:5);
% throw out cells with very low weights on all PCs
dontshow=zeros(size(allcell_PCs.score,1),1); 
% dontshow=sum(allcell_PCs.score,2,'omitnan')<-0.055 & sum(allcell_PCs.score,2,'omitnan')>-0.065;
tsnetemp=tsne(temp,'Algorithm','exact','Distance','chebychev','Exaggeration',10,'NumDimensions',2,'Perplexity',60,'Standardize',false);
currcuescore=cuescores(:,2);
cuecolsfortsne=currcuescore-min(currcuescore,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne./max(cuecolsfortsne,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne*255;
cuecolsfortsne=ceil(cuecolsfortsne); cuecolsfortsne(cuecolsfortsne==0)=1; dontshow(isnan(cuecolsfortsne))=true;  cuecolsfortsne(isnan(cuecolsfortsne))=1;
cmap=colormap('jet');
ccft=cmap(cuecolsfortsne,:); ccft(dontshow==true,:)=repmat([1 1 1],sum(dontshow==true),1);
figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft); title('cue responsive');
currcuescore=mean(cuescores(:,[1 2 3]),2,'omitnan');
cuecolsfortsne=currcuescore-min(currcuescore,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne./max(cuecolsfortsne,[],'all','omitnan'); cuecolsfortsne=cuecolsfortsne*255;
cuecolsfortsne=ceil(cuecolsfortsne); cuecolsfortsne(cuecolsfortsne==0)=1; dontshow(isnan(cuecolsfortsne))=true;  cuecolsfortsne(isnan(cuecolsfortsne))=1;
cmap=colormap('jet');
ccft=cmap(cuecolsfortsne,:); ccft(dontshow==true,:)=repmat([1 1 1],sum(dontshow==true),1);
figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft); title('pc 1 2 3 responsive');
% ccft=ones(size(D1tag,1),3)*0.9;
% ccft(A2atag==1,:)=repmat([0 0 1],sum(A2atag==1),1);
% ccft(D1tag==1,:)=repmat([1 0 0],sum(D1tag==1),1);
% figure(); scatter(tsnetemp(:,1),tsnetemp(:,2),[],ccft);

cuePC1_archetypeCells=projectDirectionOntoCPmodel(cuescores(:,1),allconditions_cpmodel,allcell_PCs); 
% cueweights should be a vector, same length as units, with each unit's
% degree of cue responsiveness (or loading to cue factor)
cuePC2_archetypeCells=projectDirectionOntoCPmodel(cuescores(:,2),allconditions_cpmodel,allcell_PCs); 
tits={'cueSuccGrp1','cueFailGrp1','uncueSuccGrp1','uncueFailGrp1';...
      'cueSuccGrp2','cueFailGrp2','uncueSuccGrp2','uncueFailGrp2';...
      'cueSuccGrp3','cueFailGrp3','uncueSuccGrp3','uncueFailGrp3';...
      'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4'};
%       'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4';...
%       'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4';...
%       'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4';...
%       'cueSuccGrp4','cueFailGrp4','uncueSuccGrp4','uncueFailGrp4'};
plotArchetypeCells(cat(4,allcell_archetypeCells,cuePC1_archetypeCells,cuePC2_archetypeCells),{'all','cuePC1','cuePC2'},tits);
%TCAneuronloadings(allconditions_cpmodel, data);

% higher rank decomposition to study cue coding after outcome
R_guess=8; % rank 8 grabs the cue v uncue specific components
allconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[2:5])),10,R_guess);
[allcell_PCs,allcell_archetypeCells,dimOrdering]=studyCPmodel(allconditions_cpmodel);
contextscores=mean(mean(data(:,1:5,2:3),3,'omitnan'),2,'omitnan')-mean(mean(data(:,1:5,4:5),3,'omitnan'),2,'omitnan');
projectDirectionOntoCPmodel(cuescores(:,2),allconditions_cpmodel,allcell_PCs); % cue direction
projectDirectionOntoCPmodel(contextscores,allconditions_cpmodel,allcell_PCs); % context direction

% uncued trying for cell type classification
R_guess=6; % guess matrix rank
uncuedconditions_cpmodel=plotTCA(noNansOrInfs(data(:,:,[2:5])),10,R_guess);
us=uncuedconditions_cpmodel.U{1};
figure(); scatter(mean(us(:,2:3),2),mean(us(:,4:5),2));
[uncuedconditions_PCs,uncuedconditions_archetypeCells]=studyCPmodel(uncuedconditions_cpmodel);

% just timepoints after outcome
% X=flatten(data(:,50:end,2:5),'mean',2)'; X=noNansOrInfs(X);
% [U,S,V]=svd(X);
% figure(); imagesc(U);
% figure(); imagesc(S(:,1:10));
% figure(); imagesc(V(1:4,:)');
% figure(); imagesc(V(:,1:4));

end

function [T,neuron_loadings]=projectCurrentDataOntoExistingCP_justvecs(fac1_vec1,fac1_vec2,fac1_vec3,test_data_matrix)

% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\allconditions_cpmodel.mat');
% load(loc);
% temp=allconditions_cpmodel.U{1}; fac1_vec1=temp(:,whichFactor);
% temp=allconditions_cpmodel.U{2}; fac1_vec2=temp(:,whichFactor);
% temp=allconditions_cpmodel.U{3}; fac1_vec3=temp(:,whichFactor);
fac1_ktens=ktensor({fac1_vec1,fac1_vec2,fac1_vec3});
fac1=outerProduct(outerProduct(fac1_vec1,fac1_vec2),fac1_vec3);
% B=fac1(1:end);
% A=test_data_matrix(1:end);
% projected=(nansum(A.*B)/nansum(B.*B)).*fac1;
% projected=(test_data_matrix.*fac1)./norm(fac1_ktens);
projected=((test_data_matrix-repmat(mean(test_data_matrix,2),1,size(test_data_matrix,2),1)).*(fac1-repmat(mean(fac1,2),1,size(fac1,2),1)))./norm(fac1_ktens);
ptens=tensor(projected);
T=hosvd(ptens,sqrt(3e-1),'rank',[1 1 1]);
neuron_loadings=T.U{1}./T.core(1,1,1);
figure(); bar(T.U{1}./T.core(1,1,1),'r');
figure(); subplot(1,3,1); plot(fac1_vec1,'Color','k'); hold on; plot(T.U{1}./T.core(1,1,1),'Color','r'); legend({'existing CP','new data'});
subplot(1,3,2); plot(fac1_vec2,'Color','k'); hold on; plot(T.U{2},'Color','r');
subplot(1,3,3); plot(fac1_vec3,'Color','k'); hold on; plot(T.U{3},'Color','r');

end

function [T,neuron_loadings]=projectCurrentDataOntoExistingCP(loc,test_data_matrix,whichFactor)

% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\allconditions_cpmodel.mat');
load(loc);
temp=allconditions_cpmodel.U{1}; fac1_vec1=temp(:,whichFactor);
temp=allconditions_cpmodel.U{2}; fac1_vec2=temp(:,whichFactor);
temp=allconditions_cpmodel.U{3}; fac1_vec3=temp(:,whichFactor);
fac1_ktens=ktensor({fac1_vec1,fac1_vec2,fac1_vec3});
fac1=outerProduct(outerProduct(fac1_vec1,fac1_vec2),fac1_vec3);
projected=(test_data_matrix.*fac1)./norm(fac1_ktens);
ptens=tensor(projected);
T=hosvd(ptens,sqrt(3e-1),'rank',[1 1 1]);
neuron_loadings=T.U{1}./T.core(1,1,1);
figure(); bar(T.U{1}./T.core(1,1,1),'r');
figure(); subplot(1,3,1); plot(fac1_vec1,'Color','k'); hold on; plot(T.U{1}./T.core(1,1,1),'Color','r'); legend({'existing CP','new data'});
subplot(1,3,2); plot(fac1_vec2,'Color','k'); hold on; plot(T.U{2},'Color','r');
subplot(1,3,3); plot(fac1_vec3,'Color','k'); hold on; plot(T.U{3},'Color','r');

end

function plotLikeTrainingSet(meanMinusTemp,sipc)

%load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\sortedbyneuronloadingsPC1.mat');
% temp=T.U{1}./T.core(1,1,1); meanMinusTemp=temp-repmat(nanmean(temp,2),1,size(temp,2));
load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\training\CP model\idx_groupLabelsFromTCAkmeans.mat');
figure();
for i=1:size(meanMinusTemp,2)
subplot(1,size(meanMinusTemp,2),i);
bar(meanMinusTemp(sipc,i),'k');
thesevals=meanMinusTemp(sipc,i);
hold on;
ngroups=2;
currc={[0, 0.75, 0.75], [0.4940, 0.1840, 0.5560]};
for j=1:ngroups
scatter(find(idx(sipc)==j),thesevals(idx(sipc)==j),2,currc{j});
end
end
load('Z:\MICROSCOPE\Kim\20221129 lab meeting\responses unit by unit\for_data_matrix_D1vA2a.mat');
figure();
for i=1:size(meanMinusTemp,2)
    subplot(1,size(meanMinusTemp,2),i);
    bar(meanMinusTemp(sipc,i));
    thesevals=meanMinusTemp(sipc,i);
    hold on;
    scatter(find(D1tag(sipc)==1),thesevals(find(D1tag(sipc)==1)),2,'r');
    scatter(find(A2atag(sipc)==1),thesevals(find(A2atag(sipc)==1)),2,'b');
end

end

function data=getTrialTypeDependentResidual(data)

nonspecific=reshape(min(noNansOrInfs(data(:,:,2:5)),[],3,'omitnan'),size(data,1),size(data,2));
% Subtract off average across all neurons, all trial types
% nonspecific=repmat(reshape(mean(mean(noNansOrInfs(data(:,:,2:5)),3,'omitnan'),1,'omitnan'),1,size(data,2)),size(data,1),1);
% Subtract off from each neuron its average time course across trial types
% nonspecific=noNansOrInfs(reshape(mean(data(:,:,2:5),3,'omitnan'),size(data,1),size(data,2)));
for i=1:5
    data(:,:,i)=data(:,:,i)-nonspecific;
end

end

function data=rmOutliersFromNormedPSTH(data)

for i=1:size(data,3)
    temp=reshape(data(:,:,i),size(data,1),size(data,2));
    for j=1:size(temp,1)
        rmd=temp(j,:)>10*mean(temp(j,:),2,'omitnan');
        temp(j,rmd)=mean(temp(j,~rmd),2,'omitnan');
    end
    data(:,:,i)=temp;
end

end

function TCAneuronloadings(allconditions_cpmodel, data)

% study neurons x trial type factors -- do cells group?
nld = allconditions_cpmodel.U{1} ./ repmat(sum(allconditions_cpmodel.U{1},2,'omitnan'),1,size(allconditions_cpmodel.U{1},2));
est_factors = cp_opt(tensor(nld),2,'init','rand','lower',0,'printitn',0,'fixsigns',true,'maxiters',10000);
figure(); temp=est_factors.U{2}; plot(temp(:,1),'Color','r'); hold on; plot(temp(:,2),'Color','b'); legend({'fac1','fac2'})
ttf=allconditions_cpmodel.U{3}; tf=allconditions_cpmodel.U{2};
%neuron_classifications=sum(nld.*repmat(temp(:,1)',size(nld,1),1).*repmat(allconditions_cpmodel.lambda',size(nld,1),1),2)>sum(nld.*repmat(temp(:,2)',size(nld,1),1).*repmat(allconditions_cpmodel.lambda',size(nld,1),1),2);
neuron_classifications=max(nld.*repmat(temp(:,1)',size(nld,1),1).*repmat(allconditions_cpmodel.lambda',size(nld,1),1),[],2)>max(nld.*repmat(temp(:,2)',size(nld,1),1).*repmat(allconditions_cpmodel.lambda',size(nld,1),1),[],2);
figure(); 
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(data(:,:,2:5),3)
    plot(reshape(mean(data(neuron_classifications,:,i),1,'omitnan'),1,size(data,2)));
    hold all;
end
title('Classified neurons average first factor');
figure(); 
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(data(:,:,2:5),3)
    plot(reshape(mean(data(~neuron_classifications,:,i),1,'omitnan'),1,size(data,2)));
    hold all;
end
title('Classified neurons average second factor');
figure(); 
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(ttf,1)
    sumfortt=zeros(size(tf(:,1)));
    for j=1:size(ttf,2)
        sumfortt=sumfortt+ttf(i,j)*tf(:,j)*temp(j,1)*allconditions_cpmodel.lambda(j);
    end
    plot(sumfortt); hold all;
end
title('First neuron loadings factor');
figure(); 
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127]./255);
for i=1:size(ttf,1)
    sumfortt=zeros(size(tf(:,1)));
    for j=1:size(ttf,2)
        sumfortt=sumfortt+ttf(i,j)*tf(:,j)*temp(j,2)*allconditions_cpmodel.lambda(j);
    end
    plot(sumfortt); hold all;
end
title('Second neuron loadings factor');

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
%         title(tits{i,j});
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
colororder([0 204 0; 255 153 153; 0 102 0; 255 0 127; 50 120 70; 100 100 100]./255);
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
% Don't need to mean subtract because matlab pca function already does this
% unitweights=unitweights-repmat(mean(unitweights,1,'omitnan'),size(unitweights,1),1);
[pcaout.coeff,pcaout.score,pcaout.latent,pcaout.tsquared,pcaout.explained,pcaout.mu]=plotPCA(unitweights,size(unitweights,2));

% drop first PC because represents DC offset?
% temp=pcaout.score(:,2:end);
% tsnetemp=tsne(temp,'Algorithm','exact','Distance','chebychev'); figure(); scatter(tsnetemp(:,1),tsnetemp(:,2));

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

function response_forThisArchetype(allconditions_cpmodel,useTheseCells)

unitweights=allconditions_cpmodel.U{1}; tf=allconditions_cpmodel.U{2};
trialfactors=allconditions_cpmodel.U{3};
TFsForUnitsAndTrialTypes=nan(length(tf(:,1)),size(trialfactors,1),size(unitweights,1));
for k=1:size(trialfactors,1)
    for i=1:size(unitweights,1)
        weightedTFforUnit=zeros(size(tf(:,1)));
        for j=1:size(tf,2)
            weightedTFforUnit=weightedTFforUnit+unitweights(i,j)*tf(:,j)*trialfactors(k,j)*allconditions_cpmodel.lambda(j);
        end
        TFsForUnitsAndTrialTypes(:,k,i)=weightedTFforUnit; 
    end
end
TFsForUnitsAndTrialTypes=TFsForUnitsAndTrialTypes(:,:,useTheseCells==1);
figure();
cs=[0 204 0; 255 153 153; 0 102 0; 255 0 127]./255;
for i=1:size(TFsForUnitsAndTrialTypes,2)
    me=reshape(mean(TFsForUnitsAndTrialTypes(:,i,:),3,'omitnan'),length(tf(:,1)),1);
    plot(me,'Color',cs(i,:),'LineWidth',2);
    hold on;
    se=reshape(std(TFsForUnitsAndTrialTypes(:,i,:),[],3,'omitnan'),length(tf(:,1)),1)./sqrt(sum(useTheseCells==1));
    plot(me+se,'Color',cs(i,:));
    plot(me-se,'Color',cs(i,:));
end

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

userChooseDecomp=true;

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
%     est_factors = cp_apr(tensor(data),R,'init','rand','lower',0,'printitn',0,'fixsigns',true,'maxiters',10000);
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
    if n < 20
        % plot the estimated factors
        viz_ktensor(est_factors, ... 
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'})
        set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
        text(0.05,0.05,['lambdas ' num2str(est_factors.lambda')]);
    end
    if userChooseDecomp==true
%         temp=est_factors.U{1};
%         D=0; fd1=find(D1tag==1); fd2=find(A2atag==1);
%         for d1=1:length(fd1)
%             for d2=1:length(fd2)
%                 D=D+pdist([temp(fd1(d1),:); temp(fd2(d2),:)]);
%             end
%         end
%         disp(D);

%         temp=est_factors.U{1};
%         meanMinusTemp=temp-repmat(nanmean(temp,2),1,size(temp,2));
%         [~,si]=sort(meanMinusTemp(:,3));
%         ngroups=2;
%         currc={'k','r','c','g','b'};
%         %Y=pdist(meanMinusTemp);
%         %Z=linkage(Y);
%         %idx=cluster(Z,'maxclust',3);
%         idx=kmeans(meanMinusTemp,ngroups,'Replicates',100);
%         figure();
%         for i=1:size(meanMinusTemp,2)
%             subplot(1,size(meanMinusTemp,2),i);
%             bar(meanMinusTemp(si,i));
%             thesevals=meanMinusTemp(si,i);
%             hold on; 
%             for j=1:ngroups
%                 scatter(find(idx(si)==j),thesevals(idx(si)==j),2,currc{j});
%             end
%         end
% load('Z:\MICROSCOPE\Kim\20221129 lab meeting\responses unit by unit\for_data_matrix_D1vA2a.mat')
% figure();
% for i=1:size(meanMinusTemp,2)
%     subplot(1,size(meanMinusTemp,2),i);
%     bar(meanMinusTemp(si,i));
%     thesevals=meanMinusTemp(si,i);
%     hold on;
%     scatter(find(D1tag(si)==1),thesevals(find(D1tag(si)==1)),2,'r');
%     scatter(find(A2atag(si)==1),thesevals(find(A2atag(si)==1)),2,'b');
% end
        
        answer=questdlg('Use this decomp?');
        switch answer
            case 'Yes'
                bestmodel=est_factors;
                besterr=err(n);
                return
            case 'No'
            case 'Cancel'
        end
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
        case 'mean'
            ma=mean(ed,whichdim(1),'omitnan');
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
        case 'mean'
           ma=mean(data,whichdim,'omitnan');
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
