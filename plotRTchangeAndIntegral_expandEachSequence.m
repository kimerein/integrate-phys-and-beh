function plotRTchangeAndIntegral_expandEachSequence(outputOfCompareCombos,metadata,sessid,binRT,sequenceLength)

plotFirstFig=0;

rt_pairs1_contingent=outputOfCompareCombos.rt_pairs1_contingent;
rt_pairs2_contingent=outputOfCompareCombos.rt_pairs2_contingent;
rt1=outputOfCompareCombos.all_rt1;
sequenceMatchStarts1=outputOfCompareCombos.sequenceMatchStarts1;
rt2=outputOfCompareCombos.all_rt2;
sequenceMatchStarts2=outputOfCompareCombos.sequenceMatchStarts2;

% expand each sequence using sequenceMatchStarts1 and sequenceLength
takeInds=nan(1,sum(sequenceMatchStarts1)*sequenceLength);
j=1;
for i=1:length(sequenceMatchStarts1)
    if sequenceMatchStarts1(i)==1
        takeInds(j:j+sequenceLength-1)=i:i+sequenceLength-1;
        j=j+sequenceLength;
    end
end
if any(takeInds>length(rt1))
    takeInds(takeInds>length(rt1))=nan;
    takeInds=takeInds(~isnan(takeInds));
end

% expand reaction times by sequences
sessids=metadata.sessid;
backup_rt1=rt1;
rt1=rt1(takeInds);
sessids=sessids(takeInds);
sequenceMatchStarts_backup=sequenceMatchStarts1;
sequenceMatchStarts1=ones(1,length(sessids));

if plotFirstFig==1
    figure();
    subplot(2,1,1);
    plot(diff(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1)),'Color','k');
    subplot(2,1,2);
    inds=downSampAv(1:length(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1)),binRT);
    plot(inds,downSampAv(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1),binRT),'Color','k');
    hold on;
    temp=downSampAv(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1),binRT);
    for i=1:length(temp)
        if temp(i)<1
            scatter(inds(i),temp(i),[],'k');
        end
        if mod(i,sequenceLength)==1
            scatter(inds(i),temp(i),[],'r');
        end
    end
    temp=rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1);
    tempdiff=[0 diff(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1))];
    trySteps=0:0.01:10;
    sumdiffs=100000*ones(1,length(trySteps));
    for i=1:length(trySteps)
        sumdiffs(i)=nansum(abs(temp-(tempdiff-trySteps(i))));
    end
    [~,bestStep]=nanmin(sumdiffs);
    plot(2:length(temp),trySteps(bestStep)+diff(rt1(ismember(sessids,sessid) & sequenceMatchStarts1'==1)),'Color','b');
    ylim([-2 4]);
end
    
% average RT across sequence
sequenceMatchStarts1=sequenceMatchStarts_backup==1 & ismember(metadata.sessid',sessid);

% Make sequence averaged RT
acrossSequences=nan(sum(sequenceMatchStarts1),sequenceLength);
f=find(sequenceMatchStarts1==1);
for i=1:size(acrossSequences,1)
    if f(i)+sequenceLength-1>length(backup_rt1)
        break
    end
    acrossSequences(i,:)=backup_rt1(f(i):f(i)+sequenceLength-1);
end

if binRT>1
    figure();
    acrossSequences=downSampMatrix(acrossSequences,binRT);
    temp=nanmean(acrossSequences,1);
    plot(downSampAv(1:sequenceLength,binRT),temp);
    hold on;
    plot(downSampAv(1:sequenceLength,binRT),temp+nanstd(acrossSequences,[],1)./sqrt(size(acrossSequences,1)),'Color','k');
    plot(downSampAv(1:sequenceLength,binRT),temp-nanstd(acrossSequences,[],1)./sqrt(size(acrossSequences,1)),'Color','k');
    binSteps=downSampAv(1:sequenceLength,binRT);
    for i=1:size(temp,2)
        if temp(i)<1
            scatter(binSteps(i),temp(i),[],'k');
        end
    end
    scatter(1,temp(1),[],'r');
else
    figure();
    temp=nanmean(acrossSequences,1);
    plot(1:sequenceLength,temp);
    hold on;
    plot(1:sequenceLength,temp+nanstd(acrossSequences,[],1)./sqrt(size(acrossSequences,1)),'Color','k');
    plot(1:sequenceLength,temp-nanstd(acrossSequences,[],1)./sqrt(size(acrossSequences,1)),'Color','k');
    for i=1:size(temp,2)
        if temp(i)<1
            scatter(i,temp(i),[],'k');
        end
    end
    scatter(1,temp(1),[],'r');
end
