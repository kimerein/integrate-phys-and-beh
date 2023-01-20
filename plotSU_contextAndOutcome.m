function plotSU_contextAndOutcome(loc,un)

plotDerivativeInstead=false;
gaussmooth=80;

if iscell(un)
    units=un;
    un=un{1};
else
    units=[];
end

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
h1=subplot(2,2,1);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); 
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth);
h1max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h1max]);
end
text(3,0,'cued');

h2=subplot(2,2,2);
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); 
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth);
h2max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h2max]);
end
text(3,0,'uncued');

h3=subplot(2,2,3);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); 
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy2=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth);
h3max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h3max]);
end
text(3,0,'success');

h4=subplot(2,2,4);
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy1=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); 
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth);
h4max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h4max]);
end
text(3,0,'failure');

if ~isempty(units)
    for i=1:length(units)
        un=units{i};
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

        axes(h1);
        dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth);
        dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth);
        h1max=max(h1max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h1max]);
        end

        axes(h2);
        dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth);
        dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth);
        h2max=max(h2max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h2max]);
        end

        axes(h3);
        dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; maxy1=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth);
        dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; maxy2=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth);
        h3max=max(h3max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h3max]);
        end

        axes(h4);
        dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; maxy1=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth);
        dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; maxy2=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth);
        h4max=max(h4max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h4max]);
        end
    end
end

end

function maxy=plotit(dataout,alignComp,c,plotDerivativeInstead,gaussmooth)

[~,ma]=nanmax(nanmean(alignComp.y,1));
alignmentTime=alignComp.x(ma);
temp=nanmean(dataout.y,1);
if plotDerivativeInstead==true
    temp=diff(temp);
    temp=[temp temp(end)];
end
temp=smoothdata(temp,'gaussian',gaussmooth);
plot(dataout.x-alignmentTime,temp,'Color',c);
maxy=max(temp(dataout.x-alignmentTime>-3 & dataout.x-alignmentTime<+10),[],'omitnan');
% maxy=max(temp(dataout.x-alignmentTime>-1 & dataout.x-alignmentTime<+4),[],'omitnan');
hold on; 
line([alignComp.x(ma)-alignmentTime alignComp.x(ma)-alignmentTime],[0 1],'Color','b');
xlim([-3 +10]);
% xlim([-1 +4]);

end
