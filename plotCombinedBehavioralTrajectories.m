function plotCombinedBehavioralTrajectories(cued_reach_rates,noncued_reach_rates)

% cued_reach_rates is a cell array where each element is a vector of cued
% reach rates for one mouse, organized in time from the first to the last
% session
% noncued_reach_rates is a cell array where each element is a vector of non-cued
% reach rates for one mouse, organized in time from the first to the last
% session

alignToFirstSess=true; % true if want to align all mice to begin at the same point

maxLength=0;
for i=1:length(cued_reach_rates)
    if length(cued_reach_rates{i})>maxLength
        maxLength=length(cued_reach_rates{i});
    end
end
tog_cued=nan(length(cued_reach_rates),maxLength);
tog_uncued=nan(length(cued_reach_rates),maxLength);
for i=1:length(cued_reach_rates)
    if alignToFirstSess==true
        temp=cued_reach_rates{i};
        tog_cued(i,1:length(cued_reach_rates{i}))=cued_reach_rates{i}-temp(1);
        temp=noncued_reach_rates{i};
        tog_uncued(i,1:length(cued_reach_rates{i}))=noncued_reach_rates{i}-temp(1);
    else
        tog_cued(i,1:length(cued_reach_rates{i}))=cued_reach_rates{i};
        tog_uncued(i,1:length(cued_reach_rates{i}))=noncued_reach_rates{i};
    end
end
mean_cued=nanmean(tog_cued,1);
mean_uncued=nanmean(tog_uncued,1);
se_cued=nanstd(tog_cued,[],1)./sqrt(size(tog_cued,1));
se_uncued=nanstd(tog_uncued,[],1)./sqrt(size(tog_uncued,1));

figure();
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(mean_cued));
scatter(mean_uncued(1),mean_cued(1),[],'k');
hold on;
for i=1:length(mean_uncued)
    line([mean_uncued(i)-se_uncued(i) mean_uncued(i)+se_uncued(i)],[mean_cued(i) mean_cued(i)],'Color',cmap(k,:));
    line([mean_uncued(i) mean_uncued(i)],[mean_cued(i)-se_cued(i) mean_cued(i)+se_cued(i)],'Color',cmap(k,:));
    hold on;
    if i>length(mean_uncued)-1
        break
    end
    quiver(mean_uncued(i),mean_cued(i),mean_uncued(i+1)-mean_uncued(i),mean_cued(i+1)-mean_cued(i),'Color',cmap(k,:));
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
xlabel('Non-cued reaching rate (Hz)');
ylabel('Cued reaching rate (Hz)');

end