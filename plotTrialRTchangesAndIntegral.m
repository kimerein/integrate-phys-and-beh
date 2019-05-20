function plotTrialRTchangesAndIntegral(outputOfCompareCombos,nApart,outputForIntegral,integral_nApart)

plotAverage=1;

rt_pairs1_contingent=outputOfCompareCombos.rt_pairs1_contingent;
rt_pairs2_contingent=outputOfCompareCombos.rt_pairs2_contingent;
rt1=outputOfCompareCombos.all_rt1;
sequenceMatchStarts1=outputOfCompareCombos.sequenceMatchStarts1;
rt2=outputOfCompareCombos.all_rt2;
sequenceMatchStarts2=outputOfCompareCombos.sequenceMatchStarts2;

rt_pairs1_integral=outputForIntegral.rt_pairs1_contingent;
rt_pairs2_integral=outputForIntegral.rt_pairs2_contingent;

if plotAverage==1
    f=find(sequenceMatchStarts1==1);
    addUpRTchange1=nan(length(f),integral_nApart);
    for i=1:length(f)
        if f(i)+integral_nApart-1>length(rt1)
            break
        end
        addUpRTchange1(i,:)=rt1(f(i):f(i)+integral_nApart-1);
    end
    figure();
    plot(nanmean(addUpRTchange1,1),'Color','k');
    hold on;
    plot(nanmean(addUpRTchange1,1)+nanstd(addUpRTchange1,[],1)./sqrt(length(f)),'Color','k');
    plot(nanmean(addUpRTchange1,1)-nanstd(addUpRTchange1,[],1)./sqrt(length(f)),'Color','k');
    
    f=find(sequenceMatchStarts2==1);
    addUpRTchange2=nan(length(f),integral_nApart);
    for i=1:length(f)
        if f(i)+integral_nApart-1>length(rt2)
            break
        end
        addUpRTchange2(i,:)=rt2(f(i):f(i)+integral_nApart-1);
    end
    hold on;
    plot(nanmean(addUpRTchange2,1),'Color','r');
    hold on;
    plot(nanmean(addUpRTchange2,1)+nanstd(addUpRTchange2,[],1)./sqrt(length(f)),'Color','r');
    plot(nanmean(addUpRTchange2,1)-nanstd(addUpRTchange2,[],1)./sqrt(length(f)),'Color','r');
    
else
    trials=1:nApart:nApart*length(rt_pairs1_contingent);
    figure();
    subplot(2,1,1);
    plot(rt_pairs1_contingent,'Color','k');
    axis tight;
    ylim([-7 7]);
    % xlim([0 max([nApart*length(rt_pairs1_contingent) integral_nApart*length(rt_pairs1_integral)])]);
    subplot(2,1,2);
    trials=1:integral_nApart:integral_nApart*length(rt_pairs1_integral);
    plot(rt1(sequenceMatchStarts1==1),'Color','k');
    axis tight;
%     ylim([-5 5]);
    % xlim([0 max([nApart*length(rt_pairs1_contingent) integral_nApart*length(rt_pairs1_integral)])]);
    
    trials=1:nApart:nApart*length(rt_pairs2_contingent);
    figure();
    subplot(2,1,1);
    plot(rt_pairs2_contingent,'Color','r');
    axis tight;
    ylim([-7 7]);
    % xlim([0 max([nApart*length(rt_pairs1_contingent) integral_nApart*length(rt_pairs1_integral)])]);
    subplot(2,1,2);
    trials=1:integral_nApart:integral_nApart*length(rt_pairs2_integral);
    plot(rt2(sequenceMatchStarts2==1),'Color','r');
    axis tight;
%     ylim([-5 5]);
    % xlim([0 max([nApart*length(rt_pairs1_contingent) integral_nApart*length(rt_pairs1_integral)])]);
end

