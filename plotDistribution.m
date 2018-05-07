function plotDistribution(data,nbins,condition,useTheseTrials)

data=data(useTheseTrials==1);
condition=condition(useTheseTrials==1);

u=unique(condition);
u=u(~isnan(u));
n=cell(1,length(u));
trialsPerU=nan(1,length(u));
for i=1:length(u)
    [n{i},x]=histcounts(data(condition==u(i)),nbins);
    trialsPerU(i)=sum(condition==u(i));
end

figure();
leg=cell(1,length(u));
for i=1:length(u)
    plot(nanmean([x(1:end-1); x(2:end)],1),n{i});
    hold all;
    leg{i}=['condition: ' num2str(u(i))];
end
xlabel('data');
ylabel('counts');
title('histogram');
legend(leg);

figure();
leg=cell(1,length(u));
for i=1:length(u)
    plot(nanmean([x(1:end-1); x(2:end)],1),n{i}./trialsPerU(i));
    hold all;
    leg{i}=['condition: ' num2str(u(i))];
end
xlabel('data');
ylabel('counts');
title('histogram normalized by # trials per condition');
legend(leg);