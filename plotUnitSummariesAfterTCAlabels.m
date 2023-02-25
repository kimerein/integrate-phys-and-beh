function plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

% for cue tuned plots
basesubtract=false;
basetimewindow=[-3 -2];
Zscore=true;
getResiduals=true;
smoothBeforeResids=true;

ds=6;
cued_success_Response.idx=groupLabelsFromTCA; exclu=cued_success_Response.excluded;
sub1=subResponse(cued_success_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(cued_success_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
outCuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedSucc,'k','r'); title('Cued success');
cued_failure_Response.idx=groupLabelsFromTCA; exclu=cued_failure_Response.excluded;
sub1=subResponse(cued_failure_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(cued_failure_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
outCuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outCuedFail,'k','r'); title('Cued failure');
uncued_success_Response.idx=groupLabelsFromTCA; exclu=uncued_success_Response.excluded;
sub1=subResponse(uncued_success_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(uncued_success_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
outUncuedSucc=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedSucc,'k','r'); title('Uncued success');
uncued_failure_Response.idx=groupLabelsFromTCA; exclu=uncued_failure_Response.excluded;
sub1=subResponse(uncued_failure_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(uncued_failure_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
outUncuedFail=plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',ds,true); plotOverlayedResponses(outUncuedFail,'k','r'); title('Uncued failure');

r.response1=outCuedSucc.response1; r.response2=outCuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Cued grp 1');
r.response1=outCuedSucc.response2; r.response2=outCuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Cued grp 2');
r.response1=outUncuedSucc.response1; r.response2=outUncuedFail.response1; plotOverlayedResponses(r,'g','r'); title('Uncued grp 1');
r.response1=outUncuedSucc.response2; r.response2=outUncuedFail.response2; plotOverlayedResponses(r,'g','r'); title('Uncued grp 2');

% failure_off={'unit91onCh1_A2atagged','unit97onCh28_A2atagged','unit98onCh31_A2atagged','unit99onCh28_A2atagged','unit147onCh22_A2atagged','unit159onCh27_A2atagged','unit160onCh27_A2atagged','unit162onCh27_A2atagged','unit163onCh27_A2atagged','unit208onCh23_A2atagged','unit208onCh25_A2atagged','unit209onCh21_A2atagged','unit209onCh23_A2atagged','unit215onCh27_A2atagged','unit217onCh27_A2atagged','unit227onCh30_A2atagged','unit233onCh30_A2atagged'};
% plotSU_contextAndOutcome('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_1\20210803\SU aligned to behavior',failure_off);

if smoothBeforeResids==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=smoothResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,10);
end

if getResiduals==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=getTrialTypeSpecificResiduals(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end

if Zscore==true
    [cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response]=ZscoreResponses(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response);
end


smoo=3;
% cuezbins=-2:0.5:3;
% cuezbins(1)=-2.0001; cuezbins(end)=3.0001;
% cuezbins=prctile(cuez,0:10:100);
% cuezbins=prctile(cuez,[0 10 20 30 40 50 60 70 75 80 85 87.5 90 92.5 95 97.5 100]);
cuezbins=prctile(cuez,[0 20 40 60 70 80 90 95 100]);
cuezbins(1)=cuezbins(1)-0.0001; cuezbins(end)=cuezbins(end)+0.0001; 
plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,'cued success',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow); 
plotByCuez(cued_failure_Response,cuez,groupLabelsFromTCA,'cued failure',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);
plotByCuez(uncued_success_Response,cuez,groupLabelsFromTCA,'uncued success',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);
plotByCuez(uncued_failure_Response,cuez,groupLabelsFromTCA,'uncued failure',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);

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

nonspecific=reshape(min(noNansOrInfs(allResp),[],3,'omitnan'),size(allResp,1),size(allResp,2));
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

function plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,addToTit,ds,smoo,cmapname,cuezbins,basesubtract,basetimewindow)

temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
time=cell(length(cuezbins)-1,1);
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0;
    sub1.incuezrange=incuezrange(temp.idx==1); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    if ~isempty(smoo)
        out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
    end
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    bycuez{i}=out.me;
    time{i}=out.t;
end
figure();
if strcmp(cmapname,'jet')
    cmap=colormap(jet(length(cuezbins)-1));
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
for i=1:length(cuezbins)-1
    plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
end
title(['groupLabel is 1 ' addToTit]);

temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',2); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=2))=0;
    sub1.incuezrange=incuezrange(temp.idx==2); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,ds,true);
    if ~isempty(smoo)
        out.me=smoothdata(out.me,'gaussian',smoo); out.plusSe=smoothdata(out.plusSe,'gaussian',smoo); out.minusSe=smoothdata(out.minusSe,'gaussian',smoo); 
    end
    if basesubtract==true
        out.me=out.me-nanmean(out.me(out.t>=basetimewindow(1) & out.t<=basetimewindow(2)));
    end
    bycuez{i}=out.me;
    time{i}=out.t;
end
figure();
if strcmp(cmapname,'jet')
    cmap=colormap(jet(length(cuezbins)-1));
else
    cmap=colormap(brewermap((length(cuezbins)-1)*2,cmap));
    cmap=cmap(length(cuezbins):end,:);
end
for i=1:length(cuezbins)-1
    plot(time{i},bycuez{i},'Color',cmap(i,:)); hold on;
end
title(['groupLabel is 2 ' addToTit]);

end