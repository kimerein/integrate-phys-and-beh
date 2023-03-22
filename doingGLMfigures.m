function doingGLMfigures(all_glm_coef,metrics,mat_all_glm_coef,mat_metrics,fromWhichSess_glm)

indsdelay=44;
endminusinds=5;
smooby=21;

% some basic plots
figure(); 
scatter(metrics.postCueAmp_over1sec-metrics.preCueAmp,metrics.cXfail_sustained); 
xlabel('post minus pre cue'); ylabel('FAILURE cue interaction');
figure(); 
scatter(metrics.postCueAmp_at1sec-metrics.preCueAmp,metrics.cXsucc_sustained); 
xlabel('post minus pre cue'); ylabel('SUCCESS cue interaction');
figure(); 
scatter(metrics.allFail_sustained,metrics.allSucc_sustained); 
xlabel('failure sustained'); ylabel('success sustained'); % this works

% cluster neurons based on glm
% coefs_after_outcome=all_glm_coef(:,[71*3+1+20:71*3+1+70 71*4+1+20:71*4+1+70 71*5+1+20:71*5+1+70 71*6+1+20:71*6+1+70 71*7+1+20:71*7+1+70 71*8+1+20:71*8+1+70]);
coefs_after_outcome=[];
% coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70])+all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70]),smooby)];
% coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70])+all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70])+all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds])+all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds])+all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70-endminusinds])+all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70-endminusinds]),smooby)];
clear trialTypeIndependent
% trialTypeIndependent(:,:,1)=smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70])+all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70]),10); 
% trialTypeIndependent(:,:,2)=smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70])+all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70])+all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70]),10); 
trialTypeIndependent(:,:,1)=smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,2)=smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,3)=smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds])+all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,4)=smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds])+all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70-endminusinds])+all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70-endminusinds]),smooby); 
% trialTypeIndependent(:,:,2)=smoothMatrix(all_glm_coef(:,[71*4+1+20:71*4+1+70]),10); 
% trialTypeIndependent(:,:,3)=smoothMatrix(all_glm_coef(:,[71*5+1+20:71*5+1+70]),10); 
% trialTypeIndependent(:,:,4)=smoothMatrix(all_glm_coef(:,[71*6+1+20:71*6+1+70]),10); 
% trialTypeIndependent(:,:,5)=smoothMatrix(all_glm_coef(:,[71*7+1+20:71*7+1+70]),10); 
% trialTypeIndependent(:,:,6)=smoothMatrix(all_glm_coef(:,[71*8+1+20:71*8+1+70]),10);
ttInd=nanmin(trialTypeIndependent,[],3);
temp=smoothMatrix(coefs_after_outcome,1)-repmat(ttInd,1,floor(size(coefs_after_outcome,2)/size(ttInd,2))); temp=temp./nanstd(temp,[],2);
idx_from_glm=kmeans(temp,3,'Replicates',50);

% MAT
py_temp=temp;
coefs_after_outcome=[];
% coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70])+mat_all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70]),smooby)];
% coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70])+mat_all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70])+mat_all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds])+mat_all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds])+mat_all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70-endminusinds])+mat_all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70-endminusinds]),smooby)];
clear trialTypeIndependent
% trialTypeIndependent(:,:,1)=smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70])+mat_all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70]),10); 
% trialTypeIndependent(:,:,2)=smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70])+mat_all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70])+mat_all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70]),10); 
trialTypeIndependent(:,:,1)=smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,2)=smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,3)=smoothMatrix(mat_all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds])+mat_all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70-endminusinds]),smooby); 
trialTypeIndependent(:,:,4)=smoothMatrix(mat_all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+mat_all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds])+mat_all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70-endminusinds])+mat_all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70-endminusinds]),smooby); 
% trialTypeIndependent(:,:,2)=smoothMatrix(all_glm_coef(:,[71*4+1+20:71*4+1+70]),10); 
% trialTypeIndependent(:,:,3)=smoothMatrix(all_glm_coef(:,[71*5+1+20:71*5+1+70]),10); 
% trialTypeIndependent(:,:,4)=smoothMatrix(all_glm_coef(:,[71*6+1+20:71*6+1+70]),10); 
% trialTypeIndependent(:,:,5)=smoothMatrix(all_glm_coef(:,[71*7+1+20:71*7+1+70]),10); 
% trialTypeIndependent(:,:,6)=smoothMatrix(all_glm_coef(:,[71*8+1+20:71*8+1+70]),10);
ttInd=nanmin(trialTypeIndependent,[],3);
temp=smoothMatrix(coefs_after_outcome,1)-repmat(ttInd,1,floor(size(coefs_after_outcome,2)/size(ttInd,2))); temp=temp./nanstd(temp,[],2);
temp(~isnan(idx_from_glm),:)=py_temp(~isnan(idx_from_glm),:);

idx_from_glm=kmeans(temp,2,'Replicates',50); 

cmap=[0, 0.75, 0.75; 0.4940, 0.1840, 0.5560];

% Make tsne plot
Y=tsne(temp); 
Ynanned=all(isnan(temp),2);

% Plot output of clustering
temp(temp<-0.4)=-0.4; figure(); imagesc(temp(idx_from_glm==2,:)); figure(); imagesc(temp(idx_from_glm==1,:));
nanma=nanmax(all_glm_coef(:,214:end),[],2); nanma(nanma<0.1)=0.1;
nama_all_glm_coef=all_glm_coef./repmat(nanma,1,size(all_glm_coef,2)); nama_all_glm_coef(nama_all_glm_coef<0)=0;
figure(); imagesc(smoothMatrix(nama_all_glm_coef(idx_from_glm==2 & ~any(abs(nama_all_glm_coef)>10,2),214:end),1)); 
figure(); imagesc(smoothMatrix(nama_all_glm_coef(idx_from_glm==1 & ~any(abs(nama_all_glm_coef)>10,2),214:end),1));
figure(); scatter(metrics.allSucc_sustained(idx_from_glm==1),metrics.allFail_sustained(idx_from_glm==1),[],cmap(2,:),'filled'); hold on;
scatter(metrics.allSucc_sustained(idx_from_glm==2),metrics.allFail_sustained(idx_from_glm==2),[],cmap(1,:),'filled');

% small group belongs to success-continuing
% figure(); scatter(metrics.allSucc_sustained(idx_from_glm==2 & metrics.allFail_sustained<=metrics.allSucc_sustained),metrics.allFail_sustained(idx_from_glm==2 & metrics.allFail_sustained<=metrics.allSucc_sustained),[],'k'); hold on;
% scatter(metrics.allSucc_sustained(idx_from_glm==2 & metrics.allFail_sustained>metrics.allSucc_sustained),metrics.allFail_sustained(idx_from_glm==2 & metrics.allFail_sustained>metrics.allSucc_sustained),[],'r');
% idx_from_glm(idx_from_glm==2)=1; % if success-continuing is group 1 
% idx_from_glm(idx_from_glm==3)=2;
% figure(); scatter(metrics.allSucc_sustained(idx_from_glm==1),metrics.allFail_sustained(idx_from_glm==1),[],'k'); hold on;
% scatter(metrics.allSucc_sustained(idx_from_glm==2),metrics.allFail_sustained(idx_from_glm==2),[],'r');

% Sort order of units according to failSustained minuns succSustained 
% [n,x]=histcounts(metrics.allFail_sustained(idx_from_glm==1)-metrics.allSucc_sustained(idx_from_glm==1),[-1 -0.1-0.00125:0.0025:0.1+0.00125 1]);
% [n,x]=cityscape_hist(n,x); figure(); plot(x,n,'Color','k'); hold on; 
% [n,x]=histcounts(metrics.allFail_sustained(idx_from_glm==2)-metrics.allSucc_sustained(idx_from_glm==2),[-1 -0.1-0.00125:0.0025:0.1+0.00125 1]);
% [n,x]=cityscape_hist(n,x); plot(x,n,'Color','r'); 
[~,orderingSuccVFail]=sort(metrics.allFail_sustained-metrics.allSucc_sustained);
nanma=nanmax(all_glm_coef(:,214:end),[],2); nanma(nanma<0.1)=0.1;
nama_all_glm_coef=all_glm_coef./repmat(nanma,1,size(all_glm_coef,2));
nama_all_glm_coef=nama_all_glm_coef(orderingSuccVFail,:); 
nama_all_glm_coef(nama_all_glm_coef<0)=0;
succ_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*3+1:71*3+71),4); 
cueXsucc_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*6+1:71*6+71),4);
fail_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*4+1:71*4+71)+nama_all_glm_coef(:,71*5+1:71*5+71),4); 
cueXfail_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*7+1:71*7+71)+nama_all_glm_coef(:,71*8+1:71*8+71),4);

% Plot tsne
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\idx_from_glm.mat'); 
figure(); scatter(Y(idx_from_glm(~isnan(idx_from_glm))==1,1),Y(idx_from_glm(~isnan(idx_from_glm))==1,2),[],cmap(2,:),'filled'); 
hold on; scatter(Y(idx_from_glm(~isnan(idx_from_glm))==2,1),Y(idx_from_glm(~isnan(idx_from_glm))==2,2),[],cmap(1,:),'filled');
% Mouse by mouse tsne control
load('Z:\MICROSCOPE\Kim\old batch Physiology Final Data Sets\training\CP model\data_loc_array.mat'); [~,~,uic]=unique(data_loc_array(:,2));
mousebymouse=uic(fromWhichSess_glm); 
[~,~,umbym]=unique(mousebymouse);
umbym=umbym(~Ynanned);
cmapm=[]; for i=1:length(unique(umbym)) cmapm=[cmapm; rand(1,3)]; end %cmap=colormap(parula(108)); 
figure(); scatter(Y(:,1),Y(:,2),[],cmapm(umbym,:),'filled'); title('Control'); 
figure(); scatter(1:length(unique(umbym)),ones(size(1:length(unique(umbym)))),[],cmapm,'filled');

% Plot histogram of tsne y axis
[n,x]=histcounts(Y(idx_from_glm(~isnan(idx_from_glm))==1,2)-Y(idx_from_glm(~isnan(idx_from_glm))==1,1),-60-0.5:1:60+0.5); 
[n,x]=cityscape_hist(n,x); figure(); plot(x,n,'Color',cmap(2,:));
[n,x]=histcounts(Y(idx_from_glm(~isnan(idx_from_glm))==2,2)-Y(idx_from_glm(~isnan(idx_from_glm))==2,1),-60-0.5:1:60+0.5); 
[n,x]=cityscape_hist(n,x); hold on; plot(x,n,'Color',cmap(1,:));

% Plot coefs ordered according to failSustained minus succSustained 
idx_from_glm=idx_from_glm(orderingSuccVFail);
figure(); imagesc(smoothMatrix(nama_all_glm_coef(idx_from_glm==1,214:end),1)); 
figure(); imagesc(smoothMatrix(nama_all_glm_coef(idx_from_glm==2,214:end),1));
temptog=[succ_nama_all_glm_coef(idx_from_glm==1,:) fail_nama_all_glm_coef(idx_from_glm==1,:) cueXsucc_nama_all_glm_coef(idx_from_glm==1,:) cueXfail_nama_all_glm_coef(idx_from_glm==1,:)];
nanmatog=nanmax(temptog,[],2); 
nanmatog(nanmatog<0.1)=0.1; 
temptog=temptog./repmat(nanmatog,1,size(temptog,2));
figure(); imagesc(smoothMatrix(temptog,1)); title('grp 1 coefs');
temptog=[succ_nama_all_glm_coef(idx_from_glm==2,:) fail_nama_all_glm_coef(idx_from_glm==2,:) cueXsucc_nama_all_glm_coef(idx_from_glm==2,:) cueXfail_nama_all_glm_coef(idx_from_glm==2,:)];
nanmatog=nanmax(temptog,[],2); 
nanmatog(nanmatog<0.1)=0.1; 
temptog=temptog./repmat(nanmatog,1,size(temptog,2));
figure(); imagesc(smoothMatrix(temptog,1)); title('grp 2 coefs');

% Plot cue x reach and cue x outcome
[~,orderingSuccVFail]=sort(metrics.allFail_sustained-metrics.allSucc_sustained);
nanma=nanmax(all_glm_coef(:,214:end),[],2); nanma(nanma<0.1)=0.1;
nama_all_glm_coef=all_glm_coef./repmat(nanma,1,size(all_glm_coef,2));
nama_all_glm_coef=nama_all_glm_coef(orderingSuccVFail,:); 
nama_all_glm_coef(nama_all_glm_coef<0)=0;
succ_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*3+1:71*3+71),4); 
cueXsucc_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*6+1:71*6+71),4);
fail_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*4+1:71*4+71)+nama_all_glm_coef(:,71*5+1:71*5+71),4); 
cueXfail_nama_all_glm_coef=smoothMatrix(nama_all_glm_coef(:,71*7+1:71*7+71)+nama_all_glm_coef(:,71*8+1:71*8+71),4);
% load('Z:\MICROSCOPE\Kim\Physiology Final Data Sets\GLM\python glm training set\idx_from_glm.mat'); 
% idx_from_glm=idx_from_glm(orderingSuccVFail);
ds=10;
temptog=[smoothMatrix(succ_nama_all_glm_coef,ds) smoothMatrix(fail_nama_all_glm_coef,ds) smoothMatrix(cueXsucc_nama_all_glm_coef,ds) smoothMatrix(cueXfail_nama_all_glm_coef,ds)];
% nanmatog=nanmax(temptog,[],2); 
% temptog=temptog./repmat(nanmatog,1,size(temptog,2));
figure(); imagesc(temptog);
% cue x reach
cuexreach=smoothMatrix(temptog(:,71*3+1:71*3+71)+temptog(:,71*2+1:71*2+71),1); % failure plus success
cuexreach(cuexreach>0.85)=0.85;
figure(); imagesc([cuexreach(idx_from_glm==1,12:end-12); nan(1,size(cuexreach(idx_from_glm==1,12:end-12),2)); cuexreach(idx_from_glm==2,12:end-12)]); title('grp 1 then 2 cue x reach');
% cue x outcome
cuexoutcome=smoothMatrix(temptog(:,71*3+1:71*3+71)-temptog(:,71*2+1:71*2+71),1); % failure minus success
figure(); imagesc([cuexoutcome(idx_from_glm==1,12:end-12); nan(1,size(cuexoutcome(idx_from_glm==1,12:end-12),2)); cuexoutcome(idx_from_glm==2,12:end-12)]); title('grp 1 then 2 cue x outcome');
% scatter plots
% figure(); scatter(metrics.allSucc_sustained(idx_from_glm==1),nanmean(cuexreach(idx_from_glm==1,31:end),2),[],'k'); hold on;
% scatter(metrics.allSucc_sustained(idx_from_glm==2),nanmean(cuexreach(idx_from_glm==2,31:end),2),[],'r'); 
% xlabel('succ sustained'); ylabel('cuexreach sustained');

% Plot other coefs
nanma=nanmax(all_glm_coef(:,[71*0+1:71*0+71 71*2+1:71*2+71]),[],2); nanma(nanma<0.1)=0.1;
nama_all_glm_coef_othercoefs=all_glm_coef./repmat(nanma,1,size(all_glm_coef,2)); 
nama_all_glm_coef_othercoefs=nama_all_glm_coef_othercoefs(orderingSuccVFail,:);
nama_all_glm_coef_othercoefs(nama_all_glm_coef_othercoefs<0)=0;
figure(); imagesc(nama_all_glm_coef_othercoefs(idx_from_glm==1,[71*0+1:71*0+71 71*2+1:71*2+71])); 
figure(); imagesc(nama_all_glm_coef_othercoefs(idx_from_glm==2,[71*0+1:71*0+71 71*2+1:71*2+71]));

% Plot with reference to TCA idx
% figure(); scatter(metrics.allSucc_sustained(idx(indexGLMcellsIntoUnitNames)==2),metrics.allMiss_sustained(idx(indexGLMcellsIntoUnitNames)==2),[],'r');
% hold on; scatter(metrics.allSucc_sustained(idx(indexGLMcellsIntoUnitNames)==1),metrics.allMiss_sustained(idx(indexGLMcellsIntoUnitNames)==1),[],'k');
% xlabel('succ sus'); ylabel('miss sus');

end