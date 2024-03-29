function plotSU_contextAndOutcome(loc,un)

plotDerivativeInstead=false;
gaussmooth=80;
dropZeros=true;
takeWhichUnits='all'; % can be 'cueTuned','noncueTuned' or 'all'
cueStartCutoff=-0.25;

if iscell(un)
    units=un;
    un=un{1};
else
    units=[];
end

a=load([loc '\cue\' un '_cueAligned.mat']);
cue.dataout=a.dataout;
cue.alignComp=a.alignComp;
a=load([loc '\uncued_reach\' un '_uncuedReach.mat']);
uncuedReach.dataout=a.dataout;
uncuedReach.alignComp=a.alignComp;
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

resp_cue=[];
resp_uncued_reach=[];
resp_h1_green=[];
resp_h1_red=[];
resp_h2_green=[];
resp_h2_red=[];
resp_h3_cyan=[];
resp_h3_black=[];
resp_h4_cyan=[];
resp_h4_black=[];

% Get cue tuning
f=figure(); 
dataout=cue.dataout; alignComp=cue.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_cue=[resp_cue; r];
dataout=uncuedReach.dataout; alignComp=uncuedReach.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_uncued_reach=[resp_uncued_reach; r];
t=linspace(-3,10,size(resp_cue,2)); beforecue=max(resp_cue(:,t>-0.9 & t<cueStartCutoff),[],'omitnan'); duringcue=max(resp_cue(:,t>cueStartCutoff & t<0.8),[],'omitnan'); cuer=duringcue-beforecue;
beforereach=max(resp_uncued_reach(:,t<-1),[],'omitnan'); duringreach=max(resp_uncued_reach(:,t>-1 & t<1),[],'omitnan'); reachr=duringreach-beforereach;
h1max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h1max]);
end
switch takeWhichUnits
    case 'cueTuned'
        if (cuer-reachr)/(cuer+reachr)>0.5
        else
            close(f); figure();
            resp_cue=[]; resp_uncued_reach=[];
        end
    case 'noncueTuned'
        if (cuer-reachr)/(cuer+reachr)>0.5
            close(f); figure();
            resp_cue=[]; resp_uncued_reach=[];
        else
        end
    case 'all'
end
cuetune_hax=gca;
text(3,0,'cue aligned v uncued reach');

% Plot other conditions
figure(); 
h1=subplot(2,2,1);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); resp_h1_green=[resp_h1_green; r];
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth); resp_h1_red=[resp_h1_red; r];
h1max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h1max]);
end
text(3,0,'cued');

h2=subplot(2,2,2);
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); resp_h2_green=[resp_h2_green; r];
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth); resp_h2_red=[resp_h2_red; r];
h2max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h2max]);
end
text(3,0,'uncued');

h3=subplot(2,2,3);
dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_h3_cyan=[resp_h3_cyan; r];
dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_h3_black=[resp_h3_black; r];
h3max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h3max]);
end
text(3,0,'success');

h4=subplot(2,2,4);
dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_h4_cyan=[resp_h4_cyan; r];
dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_h4_black=[resp_h4_black; r];
h4max=max(maxy1,maxy2); 
if plotDerivativeInstead==false 
    ylim([0, h4max]);
end
text(3,0,'failure');

if ~isempty(units)
    for i=1:length(units)
        un=units{i};
        a=load([loc '\cue\' un '_cueAligned.mat']);
        cue.dataout=a.dataout;
        cue.alignComp=a.alignComp;
        a=load([loc '\uncued_reach\' un '_uncuedReach.mat']);
        uncuedReach.dataout=a.dataout;
        uncuedReach.alignComp=a.alignComp;
        dataout=cue.dataout; alignComp=cue.alignComp; [~,rcue]=calcit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); 
        dataout=uncuedReach.dataout; alignComp=uncuedReach.alignComp; [~,runcuedreach]=calcit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); 
        t=linspace(-3,10,size(rcue,2)); beforecue=max(rcue(:,t>-0.9 & t<cueStartCutoff),[],'omitnan'); duringcue=max(rcue(:,t>cueStartCutoff & t<0.8),[],'omitnan'); cuer=duringcue-beforecue;
        beforereach=max(runcuedreach(:,t<-1),[],'omitnan'); duringreach=max(runcuedreach(:,t>-1 & t<1),[],'omitnan'); reachr=duringreach-beforereach;
        switch takeWhichUnits
            case 'cueTuned'
                if (cuer-reachr)/(cuer+reachr)>0.5
                else
                    continue % skip this unit
                end
            case 'noncueTuned'
                if (cuer-reachr)/(cuer+reachr)>0.5
                    continue % skip this unit
                else
                end
            case 'all'
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

        axes(cuetune_hax);
        dataout=cue.dataout; alignComp=cue.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_cue=[resp_cue; r];
        dataout=uncuedReach.dataout; alignComp=uncuedReach.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_uncued_reach=[resp_uncued_reach; r];
        h1max=max(h1max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h1max]);
        end

        axes(h1);
        dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); resp_h1_green=[resp_h1_green; r];
        dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth); resp_h1_red=[resp_h1_red; r];
        h1max=max(h1max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h1max]);
        end

        axes(h2);
        dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'g',plotDerivativeInstead,gaussmooth); resp_h2_green=[resp_h2_green; r];
        dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'r',plotDerivativeInstead,gaussmooth); resp_h2_red=[resp_h2_red; r];
        h2max=max(h2max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h2max]);
        end

        axes(h3);
        dataout=cuedSuccess.dataout; alignComp=cuedSuccess.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_h3_cyan=[resp_h3_cyan; r];
        dataout=uncuedSuccess.dataout; alignComp=uncuedSuccess.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_h3_black=[resp_h3_black; r];
        h3max=max(h3max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h3max]);
        end

        axes(h4);
        dataout=cuedFailure.dataout; alignComp=cuedFailure.alignComp; [maxy1,r]=plotit(dataout,alignComp,'c',plotDerivativeInstead,gaussmooth); resp_h4_cyan=[resp_h4_cyan; r];
        dataout=uncuedFailure.dataout; alignComp=uncuedFailure.alignComp; [maxy2,r]=plotit(dataout,alignComp,'k',plotDerivativeInstead,gaussmooth); resp_h4_black=[resp_h4_black; r];
        h4max=max(h4max,max(maxy1,maxy2)); 
        if plotDerivativeInstead==false 
            ylim([0, h4max]);
        end
    end
end

dodprime=true;
if dodprime==true
%     figure();
%     plot(linspace(-3,10,size(resp_cue,2)),dprime(resp_cue,resp_uncued_reach,dropZeros)); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);
%     title('cue aligned v uncued reach');

    figure();
    subplot(2,2,1);
    plot(linspace(-3,10,size(resp_h1_green,2)),dprime(resp_h1_green,resp_h1_red,dropZeros)); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,2);
    plot(linspace(-3,10,size(resp_h2_green,2)),dprime(resp_h2_green,resp_h2_red,dropZeros)); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,3);
    plot(linspace(-3,10,size(resp_h3_cyan,2)),dprime(resp_h3_cyan,resp_h3_black,dropZeros)); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,4);
    plot(linspace(-3,10,size(resp_h4_cyan,2)),dprime(resp_h4_cyan,resp_h4_black,dropZeros)); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    figure();
    subplot(2,2,1);
    if dropZeros==true
        include1=~(abs(resp_h1_green)<0.0001);
        include2=~(abs(resp_h1_red)<0.0001);
        resp_h1_green(include1~=1)=nan;
        resp_h1_red(include2~=1)=nan;
    else
        include1=ones(size(resp_h1_green)); include2=ones(size(resp_h1_red));
    end
    plot(linspace(-3,10,size(resp_h1_green,2)),mean(resp_h1_green,1,'omitnan'),'Color','g'); hold on; 
    plot(linspace(-3,10,size(resp_h1_red,2)),mean(resp_h1_red,1,'omitnan'),'Color','r'); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,2);
    if dropZeros==true
        include1=~(abs(resp_h2_green)<0.0001);
        include2=~(abs(resp_h2_red)<0.0001);
        resp_h2_green(include1~=1)=nan;
        resp_h2_red(include2~=1)=nan;
    else
        include1=ones(size(resp_h2_green)); include2=ones(size(resp_h2_red));
    end
    plot(linspace(-3,10,size(resp_h2_green,2)),mean(resp_h2_green,1,'omitnan'),'Color','g'); hold on; 
    plot(linspace(-3,10,size(resp_h2_red,2)),mean(resp_h2_red,1,'omitnan'),'Color','r'); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,3);
    if dropZeros==true
        include1=~(abs(resp_h3_cyan)<0.0001);
        include2=~(abs(resp_h3_black)<0.0001);
        resp_h3_cyan(include1~=1)=nan;
        resp_h3_black(include2~=1)=nan;
    else
        include1=ones(size(resp_h3_cyan)); include2=ones(size(resp_h3_black));
    end
    plot(linspace(-3,10,size(resp_h3_cyan,2)),mean(resp_h3_cyan,1,'omitnan'),'Color','c'); hold on; 
    plot(linspace(-3,10,size(resp_h3_black,2)),mean(resp_h3_black,1,'omitnan'),'Color','k'); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);

    subplot(2,2,4);
    if dropZeros==true
        include1=~(abs(resp_h4_cyan)<0.0001);
        include2=~(abs(resp_h4_black)<0.0001);
        resp_h4_cyan(include1~=1)=nan;
        resp_h4_black(include2~=1)=nan;
    else
        include1=ones(size(resp_h4_cyan)); include2=ones(size(resp_h4_black));
    end
    plot(linspace(-3,10,size(resp_h4_cyan,2)),mean(resp_h4_cyan,1,'omitnan'),'Color','c'); hold on; 
    plot(linspace(-3,10,size(resp_h4_black,2)),mean(resp_h4_black,1,'omitnan'),'Color','k'); line([0 0],[0 1],'Color','b'); xlim([-3 +10]);
end

end

function dprimes=dprime(resp1,resp2,dropZeros)

accuracies=nan(1,size(resp1,2));
dprimes=nan(1,size(resp1,2));
for tp=1:length(accuracies)
    % for each timepoint
    trythresh=0:0.01:max(max(resp1(:,tp),[],1,'omitnan'),max(resp2(:,tp),[],1,'omitnan'));
    if length(trythresh)==0
        continue
    end
    a1=nan(size(trythresh));
    a2=nan(size(trythresh));
    for i=1:length(trythresh)
        grp1=resp1(:,tp); grp2=resp2(:,tp);
        grp1_guess=grp1<trythresh(i); grp2_guess=grp2>=trythresh(i);
        if dropZeros==true
            isZero_gp1=abs(grp1)<0.0001;
            isZero_gp2=abs(grp2)<0.0001;
        else
            isZero_gp1=zeros(size(grp1));
            isZero_gp2=zeros(size(grp2));
        end
        a1(i)=(sum(grp1_guess(isZero_gp1~=1)==1)+sum(grp2_guess(isZero_gp2~=1)==1))/(length(grp1(isZero_gp1~=1))+length(grp2(isZero_gp2~=1)));
        grp1_guess=grp1>=trythresh(i); grp2_guess=grp2<trythresh(i);
        a2(i)=(sum(grp1_guess(isZero_gp1~=1)==1)+sum(grp2_guess(isZero_gp2~=1)==1))/(length(grp1(isZero_gp1~=1))+length(grp2(isZero_gp2~=1)));
    end
    accuracies(tp)=max([a1; a2],[],'all','omitnan');
    dprimes(tp)=2*norminv(accuracies(tp));
end

end

function [maxy,yout]=plotit(dataout,alignComp,c,plotDerivativeInstead,gaussmooth)

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
yout=temp(dataout.x-alignmentTime>-3 & dataout.x-alignmentTime<+10);

end

function [maxy,yout]=calcit(dataout,alignComp,c,plotDerivativeInstead,gaussmooth)

[~,ma]=nanmax(nanmean(alignComp.y,1));
alignmentTime=alignComp.x(ma);
temp=nanmean(dataout.y,1);
if plotDerivativeInstead==true
    temp=diff(temp);
    temp=[temp temp(end)];
end
temp=smoothdata(temp,'gaussian',gaussmooth);
maxy=max(temp(dataout.x-alignmentTime>-3 & dataout.x-alignmentTime<+10),[],'omitnan');
% maxy=max(temp(dataout.x-alignmentTime>-1 & dataout.x-alignmentTime<+4),[],'omitnan');
yout=temp(dataout.x-alignmentTime>-3 & dataout.x-alignmentTime<+10);

end
