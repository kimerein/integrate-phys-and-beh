function plotUnitSummariesAfterTCAlabels(groupLabelsFromTCA,cuez,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response)

cued_success_Response.idx=groupLabelsFromTCA; exclu=cued_success_Response.excluded;
sub1=subResponse(cued_success_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(cued_success_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',1); title('Cued success');
cued_failure_Response.idx=groupLabelsFromTCA; exclu=cued_failure_Response.excluded;
sub1=subResponse(cued_failure_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(cued_failure_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',1); title('Cued failure');
uncued_success_Response.idx=groupLabelsFromTCA; exclu=uncued_success_Response.excluded;
sub1=subResponse(uncued_success_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(uncued_success_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',1); title('Uncued success');
uncued_failure_Response.idx=groupLabelsFromTCA; exclu=uncued_failure_Response.excluded;
sub1=subResponse(uncued_failure_Response,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0; sub2=subResponse(uncued_failure_Response,'idx',2); f=find(exclu~=0); sub2.excluded(f(groupLabelsFromTCA~=1))=0;
plotVariousSUResponsesAlignedToBeh('scatterResponseVsResponse',sub1,sub2,'meanAcrossUnits',1); title('Uncued failure');
% failure_off={'unit91onCh1_A2atagged','unit97onCh28_A2atagged','unit98onCh31_A2atagged','unit99onCh28_A2atagged','unit147onCh22_A2atagged','unit159onCh27_A2atagged','unit160onCh27_A2atagged','unit162onCh27_A2atagged','unit163onCh27_A2atagged','unit208onCh23_A2atagged','unit208onCh25_A2atagged','unit209onCh21_A2atagged','unit209onCh23_A2atagged','unit215onCh27_A2atagged','unit217onCh27_A2atagged','unit227onCh30_A2atagged','unit233onCh30_A2atagged'};
% plotSU_contextAndOutcome('Z:\MICROSCOPE\Kim\WHISPER recs\Mar_1\20210803\SU aligned to behavior',failure_off);

pause;

plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA);
plotByCuez(cued_failure_Response,cuez,groupLabelsFromTCA);
plotByCuez(uncued_success_Response,cuez,groupLabelsFromTCA);
plotByCuez(uncued_failure_Response,cuez,groupLabelsFromTCA);

end

function plotByCuez(cued_success_Response,cuez,groupLabelsFromTCA)

cuezbins=-2:1:3;
cuezbins(1)=-2.0001; cuezbins(end)=3.0001;
temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',1); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=1))=0;
    sub1.incuezrange=incuezrange(temp.idx==1); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,1,true);
    bycuez{i}=out.me;
end
figure();
cmap=colormap(brewermap((length(cuezbins)-1)*2,"Purples"));
for i=1:length(cuezbins)-1
    plot(bycuez{i},'Color',cmap(i+(length(cuezbins)-1),:)); hold on;
end
title('groupLabel is 1');

temp=cued_success_Response;
bycuez=cell(length(cuezbins)-1,1);
for i=1:length(cuezbins)-1
    temp.idx=groupLabelsFromTCA; exclu=temp.excluded;
    incuezrange=cuez>cuezbins(i) & cuez<=cuezbins(i+1);
    sub1=subResponse(temp,'idx',2); f=find(exclu~=0); sub1.excluded(f(groupLabelsFromTCA~=2))=0;
    sub1.incuezrange=incuezrange(temp.idx==2); sub1=subResponse(sub1,'incuezrange',1); f=find(exclu~=0); sub1.excluded(f(incuezrange~=1))=0;
    out=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',sub1,1,true);
    bycuez{i}=out.me;
end
figure();
cmap=colormap(brewermap((length(cuezbins)-1)*2,"Purples"));
for i=1:length(cuezbins)-1
    plot(bycuez{i},'Color',cmap(i+(length(cuezbins)-1),:)); hold on;
end
title('groupLabel is 2');

end