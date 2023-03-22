function studyGLMcoef(coef,times,whichCoefToUse)

smoobin=1;

if size(coef,2)>length(times)
    whichCoefs=[];
    j=1;
    for i=1:length(times):size(coef,2)
        whichCoefs=[whichCoefs j*ones([1,length(times)])];
        j=j+1;
    end
else
    whichCoefs=ones(1,length(times));
end

temp=coef(:,ismember(whichCoefs,whichCoefToUse));
% Zscore temp
temp=temp-repmat(mean(temp,2,'omitnan'),1,size(temp,2));
temp=temp-repmat(std(temp,[],2,'omitnan'),1,size(temp,2));

if smoobin~=1
    for i=1:size(temp,1)
        temp(i,:)=smooth(temp(i,:),smoobin);
    end
end
% make square matrix
% temp is currently neurons by coef times
% to study temporal factors, get times by times matrix
A=temp'*temp;
[coeff,score,latent,tsquared,explained,mu]=plotPCA(A,10);

whichPC=7;
figure();
for i=1:length(unique(whichCoefToUse))
    plot(times,coeff([1:length(times)] +(i-1)*length(times),whichPC));
    hold all;
end
% project onto PC7
projectionsOntoPC=nan(size(temp,1),1);
for i=1:size(temp,1)
    projectionsOntoPC(i)=dot(temp(i,:),coeff(:,whichPC));
end
figure();
histogram(projectionsOntoPC,1000);

load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\GLM coeffs\indexGLMcellsIntoUnitNames.mat')
load('Z:\MICROSCOPE\Kim\20230205 all SU alignments\GLM coeffs\cue_Response.mat')
labeltemp=cue_Response.D1tag(cue_Response.excluded==0); toplotr=projectionsOntoPC(labeltemp(indexGLMcellsIntoUnitNames)==1);
labeltemp=cue_Response.A2atag(cue_Response.excluded==0); toplotb=projectionsOntoPC(labeltemp(indexGLMcellsIntoUnitNames)==1);

[n,x]=histcounts(toplotr,200);
[n,x]=cityscape_hist(n,x);
figure(); plot(x,n,'Color','r');
hold on;
[n,x]=histcounts(toplotb,200);
[n,x]=cityscape_hist(n,x);
plot(x,n,'Color','b');

labeltempr=cue_Response.D1tag(cue_Response.excluded==0); figure(); plot(nanmean(coef(labeltempr(indexGLMcellsIntoUnitNames)==1,:),1),'Color','r');
labeltempb=cue_Response.A2atag(cue_Response.excluded==0); hold on; plot(nanmean(coef(labeltempb(indexGLMcellsIntoUnitNames)==1,:),1),'Color','b');

% figure(); histogram(mean(coef(:,122:142),2,'omitnan')-mean(coef(:,162:182),2,'omitnan'),1000);
% cuer=nanmax(all_glm_coef(:,10:40)-repmat(mean(all_glm_coef(1:9,:),'all','omitnan'),1,size(all_glm_coef(:,10:40),2)),[],2);


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

function bestmodel=plotTCA(inputdata,n_fits,R)

userChooseDecomp=false;

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
    if n < 10
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