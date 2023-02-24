function plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

% for cue tuned plots
basesubtract=false;
basetimewindow=[-3 -2];

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

% pause;
ds=1;
smoo=42;
% cuezbins=-2:0.5:3;
% cuezbins(1)=-2.0001; cuezbins(end)=3.0001;
cuezbins=prctile(cuez,0:5:100);
cuezbins(1)=cuezbins(1)-0.0001; cuezbins(end)=cuezbins(end)+0.0001; 
plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA,'cued success',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow); 
plotByCuez(cued_failure_Response,cuez,groupLabelsFromTCA,'cued failure',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);
plotByCuez(uncued_success_Response,cuez,groupLabelsFromTCA,'uncued success',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);
plotByCuez(uncued_failure_Response,cuez,groupLabelsFromTCA,'uncued failure',ds,smoo,'jet',cuezbins,basesubtract,basetimewindow);

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