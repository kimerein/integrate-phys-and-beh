function plotRTchangeAndIntegral_expandEachSequence(outputOfCompareCombos,metadata,sessid,binRT,sequenceLength)

rt_pairs1_contingent=outputOfCompareCombos.rt_pairs1_contingent;
rt_pairs2_contingent=outputOfCompareCombos.rt_pairs2_contingent;
rt1=outputOfCompareCombos.all_rt1;
sequenceMatchStarts1=outputOfCompareCombos.sequenceMatchStarts1;
rt2=outputOfCompareCombos.all_rt2;
sequenceMatchStarts2=outputOfCompareCombos.sequenceMatchStarts2;

% get which points to plot from sequenceMatchStarts
if ~isempty(sequenceLength)
    useTrials=zeros(size(sequenceMatchStarts1));
    for i=1:length(sequenceMatchStarts1)
        if sequenceMatchStarts1(i)==1
            useTrials(i:i+sequenceLength-1)=1;
        end
    end
end

figure();
subplot(2,1,1);
plot(diff(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1)),'Color','k');
subplot(2,1,2);
inds=downSampAv(1:length(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1)),binRT);
plot(inds,downSampAv(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1),binRT),'Color','k');
hold on;
temp=downSampAv(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1),binRT);
for i=1:length(temp)
    if temp(i)<1
        scatter(inds(i),temp(i),[],'k');
    end
end
temp=rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1);
tempdiff=[0 diff(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1))];
trySteps=0:0.01:10;
sumdiffs=100000*ones(1,length(trySteps));
for i=1:length(trySteps)
    sumdiffs(i)=nansum(abs(temp-(tempdiff-trySteps(i))));
end
[~,bestStep]=nanmin(sumdiffs);
plot(2:length(temp),trySteps(bestStep)+diff(rt1(metadata.sessid==sessid & sequenceMatchStarts1'==1)),'Color','b');
ylim([-2 4]);
