function plotSU_contextAndOutcome(loc,un)

a=load([loc '\cued_success\' un '_cuedSuccess.mat']);
cuedSuccess.dataout=a.dataout;
cuedSuccess.alignComp=a.alignComp;
a=load([loc '\cued_drop\' un '_cuedDrop.mat']);
cuedFailure.dataout=a.dataout;
cuedFailure.alignComp=a.alignComp;
a=load([loc '\uncued_success\' un '_uncuedSuccess.mat']);
uncuedSuccess.dataout=a.dataout;
uncuedSuccess.alignComp=a.alignComp;
a=load([loc '\uncued_drop\' un '_uncuedDrop.mat']);
uncuedFailure.dataout=a.dataout;
uncuedFailure.alignComp=a.alignComp;

figure(); 
subplot(2,2,1);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g'); 
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r');
ylim([0, max(maxy1,maxy2)]);
text(3,0,'cued');

subplot(2,2,2);
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g'); 
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r');
ylim([0, max(maxy1,maxy2)]);
text(3,0,'uncued');

subplot(2,2,3);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'c'); 
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy2=plotit(dataout,alignComp,'k');
ylim([0, max(maxy1,maxy2)]);
text(3,0,'success');

subplot(2,2,4);
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy1=plotit(dataout,alignComp,'c'); 
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'k');
ylim([0, max(maxy1,maxy2)]);
text(3,0,'failure');

end

function maxy=plotit(dataout,alignComp,c)

[~,ma]=nanmax(nanmean(alignComp.y,1));
alignmentTime=alignComp.x(ma);
plot(dataout.x-alignmentTime,smoothdata(nanmean(dataout.y,1),'gaussian',40),'Color',c);
maxy=max(ylim);
hold on; 
line([alignComp.x(ma)-alignmentTime alignComp.x(ma)-alignmentTime],[0 1],'Color','b');
xlim([-3 +10]);

end
