function plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,isSig)

% for cue tuned plots
basesubtract=true;
basetimewindow=[5 9]; %[4 9];

% plot all SU
% although ugly, the raw raw data actually shows effects (maybe for
% supplement)
% cuezbins=prctile(cuez,0:5:100); 
% cuezbins=prctile(cuez,[0:10:90 92 94 96 97 98 99 100]); cuezbins(1)=cuezbins(1)-0.0001; cuezbins(end)=cuezbins(end)+0.0001;

% temp=prctile(cuez(groupLabelsFromTCA==1),[0:25:75 80 85 90 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==2),[0:25:75 80 85 90 95 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;
temp=prctile(cuez(groupLabelsFromTCA==1),[0:25:75 85 92 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
temp=prctile(cuez(groupLabelsFromTCA==2),[0:25:75 85 92 100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;

% temp=prctile(cuez(groupLabelsFromTCA==1),[1:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{1}=temp;
% temp=prctile(cuez(groupLabelsFromTCA==2),[1:100]); temp(1)=temp(1)-0.0001; temp(end)=temp(end)+0.0001; cuezbins{2}=temp;

plotAll=false;
Zscore=false;
minmaxnorm=false;
smoo=12; %6; %smoo=3; %smoo=42;
getResiduals=false; % but need this to get rid of mid-range
dsForCuez=6;

if isempty(isSig)
    isSig=ones(size(groupLabelsFromTCA));
end
groupLabelsFromTCA(isSig~=1)=-10; % omit these in plotting

% else plotBinRange
% cuezbins=prctile(cuez,[0 20 40 60 70 80 82 84 86 88 90 91 92 93 95 97 100]);
% plotAll=false;
% Zscore=false; % no because miss cued success activity for grp 2
% smoo=6;
% getResiduals=true; % but need this to get rid of mid-range




ds=6;
cued_success_Response.idx=groupLabelsFromTCA; exclu=cued_success_Response.excluded; gpLab=1;
sub1=subResponse(cued_success_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(cued_success_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
outCuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedSucc,'k','r'); title('Cued success');
cued_failure_Response.idx=groupLabelsFromTCA; exclu=cued_failure_Response.excluded; gpLab=1;
sub1=subResponse(cued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(cued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
outCuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedFail,'k','r'); title('Cued failure');
uncued_success_Response.idx=groupLabelsFromTCA; exclu=uncued_success_Response.excluded; gpLab=1;
sub1=subResponse(uncued_success_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(uncued_success_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
outUncuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedSucc,'k','r'); title('Uncued success');
uncued_failure_Response.idx=groupLabelsFromTCA; exclu=uncued_failure_Response.excluded; gpLab=1;
sub1=subResponse(uncued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0; gpLab=2; sub2=subResponse(uncued_failure_Response,'idx',gpLab); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=gpLab))=0;
outUncuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedFail,'k','r'); title('Uncued failure');

r.response1=outCuedSucc.response1; r.response2=outCuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Cued grp 1');
r.response1=outCuedSucc.response2; r.response2=outCuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Cued grp 2');
r.response1=outUncuedSucc.response1; r.response2=outUncuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Uncued grp 1');
r.response1=outUncuedSucc.response2; r.response2=outUncuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Uncued grp 2');

% failure_off={'unit91onCh1_A2atagged','unit97onCh28_A2atagged','unit98onCh31_A2atagged','unit99onCh28_A2atagged','unit147onCh22_A2atagged','unit159onCh27_A2atagged','unit160onCh27_A2atagged','unit162onCh27_A2atagged','unit163onCh27_A2atagged','unit208onCh23_A2atagged','unit208onCh25_A2atagged','unit209onCh21_A2atagged','unit209onCh23_A2atagged','unit215onCh27_A2atagged','unit217onCh27_A2atagged','unit227onCh30_A2atagged','unit233onCh30_A2atagged'};
% plotSU_contextAndOutcome('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_1\20210803\SU aligned to behavior',failure_off);

smoothBeforeResids=false;
if smoothBeforeResids==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=smoothResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,6);
end

if getResiduals==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=getTrialTypeSpecificResiduals(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end

if Zscore==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=ZscoreResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
elseif minmaxnorm==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=maxNorm(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end

if basesubtract==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=baseSubResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,basetimewindow);
end

% cuezbins=-2:0.5:3;
% cuezbins(1)=-2.0001; cuezbins(end)=3.0001;
% cuezbins=prctile(cuez,0:10:100);
% cuezbins=prctile(cuez,[0 10 20 30 40 50 60 70 75 80 85 87.5 90 92.5 95 97.5 100]);
% cuezbins=prctile(cuez,[0 20 40 60 70 80 82 84 86 88 90 91 92 93 95 97 100]); 
basesubtract=false;
plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,'cued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll); 
plotByCuez(cued_failure_Response,cuez,groupLabelsFromTCA,'cued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
plotByCuez(uncued_success_Response,cuez,groupLabelsFromTCA,'uncued success',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);
plotByCuez(uncued_failure_Response,cuez,groupLabelsFromTCA,'uncued failure',dsForCuez,smoo,'jet',cuezbins,basesubtract,basetimewindow,plotAll);

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=maxNorm(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

minmaxwindow=[-2 5]; % in sec wrt cue
[cuetime,t]=getCueTimeForThisResponse(cued_success_Response);
cued_success_Response.unitbyunit_y(cued_success_Response.unitbyunit_y<0.001)=0;
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y./repmat(max(cued_success_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(cued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(cued_failure_Response);
cued_failure_Response.unitbyunit_y(cued_failure_Response.unitbyunit_y<0.001)=0;
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y./repmat(max(cued_failure_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(cued_failure_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_success_Response);
uncued_success_Response.unitbyunit_y(uncued_success_Response.unitbyunit_y<0.001)=0;
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y./repmat(max(uncued_success_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(uncued_success_Response.unitbyunit_y,2));
[cuetime,t]=getCueTimeForThisResponse(uncued_failure_Response);
uncued_failure_Response.unitbyunit_y(uncued_failure_Response.unitbyunit_y<0.001)=0;
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y./repmat(max(uncued_failure_Response.unitbyunit_y(:,t>minmaxwindow(1) & t<minmaxwindow(2)),[],2,'omitnan'),1,size(uncued_failure_Response.unitbyunit_y,2));
      
end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=smoothResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,smoo)

temp=cued_success_Response.unitbyunit_y;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
end
cued_success_Response.unitbyunit_y=temp;

temp=cued_failure_Response.unitbyunit_y;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
end
cued_failure_Response.unitbyunit_y=temp;

temp=uncued_success_Response.unitbyunit_y;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
end
uncued_success_Response.unitbyunit_y=temp;

temp=uncued_failure_Response.unitbyunit_y;
for i=1:size(temp,1)
    temp(i,:)=smoothdata(temp(i,:),'gaussian',smoo);
end
uncued_failure_Response.unitbyunit_y=temp;

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=baseSubResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,basewindow)

% basewindow=[5 9]; % in sec wrt cue
t=nanmean(cued_success_Response.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(cued_success_Response.aligncomp_y,1));
talign=nanmean(cued_success_Response.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;
t2=nanmean(cued_failure_Response.unitbyunit_x,1); t2=t2-cuetime;
t3=nanmean(uncued_success_Response.unitbyunit_x,1); t3=t3-cuetime;
t4=nanmean(uncued_failure_Response.unitbyunit_x,1); t4=t4-cuetime;
base=mean([cued_success_Response.unitbyunit_y(:,t>basewindow(1) & t<basewindow(2)) cued_failure_Response.unitbyunit_y(:,t2>basewindow(1) & t2<basewindow(2)) uncued_success_Response.unitbyunit_y(:,t3>basewindow(1) & t3<basewindow(2)) uncued_failure_Response.unitbyunit_y(:,t4>basewindow(1) & t4<basewindow(2))],2,'omitnan');
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y-repmat(base,1,size(cued_success_Response.unitbyunit_y,2));
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y-repmat(base,1,size(cued_failure_Response.unitbyunit_y,2));
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y-repmat(base,1,size(uncued_success_Response.unitbyunit_y,2));
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y-repmat(base,1,size(uncued_failure_Response.unitbyunit_y,2));

end

function [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=ZscoreResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

sdwindow=[4 9]; % in sec wrt cue
t=nanmean(cued_success_Response.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(cued_success_Response.aligncomp_y,1));
talign=nanmean(cued_success_Response.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;
t2=nanmean(cued_failure_Response.unitbyunit_x,1); t2=t2-cuetime;
t3=nanmean(uncued_success_Response.unitbyunit_x,1); t3=t3-cuetime;
t4=nanmean(uncued_failure_Response.unitbyunit_x,1); t4=t4-cuetime;
sd=std([cued_success_Response.unitbyunit_y(:,t>sdwindow(1) & t<sdwindow(2)) cued_failure_Response.unitbyunit_y(:,t2>sdwindow(1) & t2<sdwindow(2)) uncued_success_Response.unitbyunit_y(:,t3>sdwindow(1) & t3<sdwindow(2)) uncued_failure_Response.unitbyunit_y(:,t4>sdwindow(1) & t4<sdwindow(2))],[],2,'omitnan');
sd(sd<0.01)=0.01;
cued_success_Response.unitbyunit_y=cued_success_Response.unitbyunit_y./repmat(sd,1,size(cued_success_Response.unitbyunit_y,2));
cued_success_Response.unitbyunit_y(cued_success_Response.unitbyunit_y>60)=60;
cued_failure_Response.unitbyunit_y=cued_failure_Response.unitbyunit_y./repmat(sd,1,size(cued_failure_Response.unitbyunit_y,2));
cued_failure_Response.unitbyunit_y(cued_failure_Response.unitbyunit_y>60)=60;
uncued_success_Response.unitbyunit_y=uncued_success_Response.unitbyunit_y./repmat(sd,1,size(uncued_success_Response.unitbyunit_y,2));
uncued_success_Response.unitbyunit_y(uncued_success_Response.unitbyunit_y>60)=60;
uncued_failure_Response.unitbyunit_y=uncued_failure_Response.unitbyunit_y./repmat(sd,1,size(uncued_failure_Response.unitbyunit_y,2));
uncued_failure_Response.unitbyunit_y(uncued_failure_Response.unitbyunit_y>60)=60;

end

function [cuetime,t]=getCueTimeForThisResponse(r1)

t=nanmean(r1.unitbyunit_x,1);
[~,cueind]=nanmax(nanmean(r1.aligncomp_y,1));
talign=nanmean(r1.aligncomp_x,1);
cuetime=talign(cueind);
t=t-cuetime;

end

function [r1,r2,r3,r4]=getTrialTypeSpecificResiduals(r1,r2,r3,r4)

smoothBeforeNonspecific=true;
smooBefore=200;

[cuetime,t]=getCueTimeForThisResponse(r1);

% assumes same binning for all
timesteps=[mode(diff(nanmean(r1.unitbyunit_x,1))) mode(diff(nanmean(r2.unitbyunit_x,1))) mode(diff(nanmean(r3.unitbyunit_x,1))) mode(diff(nanmean(r4.unitbyunit_x,1)))];
if ~all(diff(timesteps)<0.001)
    error('Binning must match in all responses');
end

resp1=r1.unitbyunit_y(:,t>=-3 & t<=9); 
[cuetime,t]=getCueTimeForThisResponse(r2); resp2=r2.unitbyunit_y(:,t>=-3 & t<=9.01);
[cuetime,t]=getCueTimeForThisResponse(r3); resp3=r3.unitbyunit_y(:,t>=-3 & t<=9.01);
[cuetime,t]=getCueTimeForThisResponse(r4); resp4=r4.unitbyunit_y(:,t>=-3 & t<=9.01);
resp2=resp2(:,1:size(resp1,2));
resp3=resp3(:,1:size(resp1,2));
resp4=resp4(:,1:size(resp1,2));
allResp=nan(size(resp1,1),size(resp1,2),4);
allResp(:,:,1)=resp1; allResp(:,:,2)=resp2; allResp(:,:,3)=resp3; allResp(:,:,4)=resp4;

forNonspec=allResp;
if smoothBeforeNonspecific==true
    for i=1:size(forNonspec,3)
        for j=1:size(forNonspec,1)
            forNonspec(j,:,i)=smoothdata(reshape(forNonspec(j,:,i),1,size(forNonspec,2)),'gaussian',smooBefore);
        end
    end
end
nonspecific=reshape(min(noNansOrInfs(forNonspec),[],3,'omitnan'),size(forNonspec,1),size(forNonspec,2));
for i=1:4
    allResp(:,:,i)=allResp(:,:,i)-nonspecific;
end

% put back
[cuetime,t]=getCueTimeForThisResponse(r1); r1.unitbyunit_y(:,t>=-3 & t<=9)=reshape(allResp(:,:,1),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r2); f=find(t>=-3,1,'first'); r2.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,2),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r3); f=find(t>=-3,1,'first'); r3.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,3),size(allResp,1),size(allResp,2));
[cuetime,t]=getCueTimeForThisResponse(r4); f=find(t>=-3,1,'first'); r4.unitbyunit_y(:,f:f+size(allResp,2)-1)=reshape(allResp(:,:,4),size(allResp,1),size(allResp,2));

end

function data=noNansOrInfs(data)

data(isnan(data))=mean(data(isfinite(data)),'all');
data(isinf(data))=mean(data(isfinite(data)),'all');

end

function plotOverlayedResponses(out,c1,c2)

figure();
plot(out.response1.t,out.response1.me,'Color',c1); hold on; 
plot(out.response1.t,out.response1.plusSe,'Color',c1); plot(out.response1.t,out.response1.minusSe,'Color',c1);
plot(out.response2.t,out.response2.me,'Color',c2); hold on; 
plot(out.response2.t,out.response2.plusSe,'Color',c2); plot(out.response2.t,out.response2.minusSe,'Color',c2);

end

function cmap=getCmapWithRed(cuezbins)

cmap=colormap(jet(255));
cmapstep=floor(255/(length(cuezbins)-1));
cmap=cmap(end:-cmapstep:1,:);
cmap=cmap(1:length(cuezbins)-1,:);
cmap=flipud(cmap);

end

function plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,addToTit,ds,smoo,cmapname,cuezbins,basesubtract,basetimewindow,plotAll)

alph=0.8;
plotBackwards=false;
doMaxInstead=false;
maxTimeBin=0.1;
noGaussSmooth=true;

backup_cuezbins=cuezbins;
if iscell(backup_cuezbins)
    cuezbins=backup_cuezbins{1};
end
temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
time=cell(length(cuezbins)-1,1);
gpLab=1;
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0;
    sub1.incuezrange=incuezrange(temp.idx==gpLab); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    if all(isnan(sub1.unitbyunit_y),'all')
        continue
    end
    if doMaxInstead==true
        out=plotVariousSUResponsesAlignedToBeh('maxAcrossUnits',sub1,ds,maxTimeBin,true);
    else
        out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    end
    if ~isempty(smoo)
        if noGaussSmooth==true
            out.me=smooth(out.me,smoo); out.plusSe=smooth(out.plusSe,smoo); out.minusSe=smooth(out.minusSe,smoo); 
        else
            out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
        end
        if plotAll==true
            for j=1:size(out.unitbyunit,1)
                if noGaussSmooth==true
                    out.unitbyunit(j,:)=smooth(out.unitbyunit(j,:),smoo);
                else
                    out.unitbyunit(j,:)=smoothdata(out.unitbyunit(j,:),'gaussian',smoo);
                end
            end
        end
    end
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    if plotAll==true
        bycuez{i}=out.unitbyunit;
        time{i}=out.unittimes;
    else
        bycuez{i}=out.me;
        time{i}=out.t;
    end
end
figure();
if strcmp(cmapname,'jet')
%     cmap=colormap(jet(length(cuezbins)-1));
    cmap=getCmapWithRed(cuezbins);
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
if plotBackwards==true
    for i=length(cuezbins)-1:-1:1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
else
    for i=1:length(cuezbins)-1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
end
title(['groupLabel is 1 ' addToTit]);

if iscell(backup_cuezbins)
    cuezbins=backup_cuezbins{2};
end
temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
gpLab=2;
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',gpLab); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=gpLab))=0;
    sub1.incuezrange=incuezrange(temp.idx==gpLab); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    if all(isnan(sub1.unitbyunit_y),'all')
        continue
    end
    if doMaxInstead==true
        out=plotVariousSUResponsesAlignedToBeh('maxAcrossUnits',sub1,ds,maxTimeBin,true);
    else
        out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    end
    if ~isempty(smoo)
        if noGaussSmooth==true
            out.me=smooth(out.me,smoo); out.plusSe=smooth(out.plusSe,smoo); out.minusSe=smooth(out.minusSe,smoo); 
        else
            out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
        end
        if plotAll==true
            for j=1:size(out.unitbyunit,1)
                if noGaussSmooth==true
                    out.unitbyunit(j,:)=smooth(out.unitbyunit(j,:),smoo);
                else
                    out.unitbyunit(j,:)=smoothdata(out.unitbyunit(j,:),'gaussian',smoo);
                end
            end
        end
    end
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    if plotAll==true
        bycuez{i}=out.unitbyunit;
        time{i}=out.unittimes;
    else
        bycuez{i}=out.me;
        time{i}=out.t;
    end
end
figure();
if strcmp(cmapname,'jet')
%     cmap=colormap(jet(length(cuezbins)-1));
    cmap=getCmapWithRed(cuezbins);
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
if plotBackwards==true
    for i=length(cuezbins)-1:-1:1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
    end
else
    for i=1:length(cuezbins)-1
        if isempty(bycuez{i})
            continue
        end
        h=plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
        if plotAll==true
            % h = findobj(gca,'Type','line');
            for j=1:length(h)
                h(j).Color=[cmap(i,:),alph];
            end
        end
        hold on;
    end
end
title(['groupLabel is 2 ' addToTit]);

end