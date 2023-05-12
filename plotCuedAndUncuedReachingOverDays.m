function plotCuedAndUncuedReachingOverDays(all_rr_cued,all_rr_uncued,days,whichdays)

colorsUpTo=20;
% colorsUpTo=40;

cmap=colormap(cool(colorsUpTo));

figure();
for i=1:size(all_rr_cued,2)
    allrrcued_forthisday=all_rr_cued(:,i);
    allrruncued_forthisday=all_rr_uncued(:,i);
    if ~ismember(days(i),whichdays)
        continue
    end
    me_cued=nanmean(allrrcued_forthisday);
    me_uncued=nanmean(allrruncued_forthisday);
    se_cued=nanstd(allrrcued_forthisday)./sqrt(nansum(~isnan(allrrcued_forthisday)));
    se_uncued=nanstd(allrruncued_forthisday)./sqrt(nansum(~isnan(allrruncued_forthisday)));
    if days(i)<1
        currc=cmap(1,:);
    elseif days(i)>colorsUpTo
        currc=cmap(colorsUpTo,:);
    else
        currc=cmap(days(i),:);
    end
%     scatter(me_uncued,me_cued,75,currc,'filled'); hold on;
    line([me_uncued-se_uncued me_uncued+se_uncued],[me_cued me_cued],'Color',currc); hold on;
    line([me_uncued me_uncued],[me_cued-se_cued me_cued+se_cued],'Color',currc);
end
plot(nanmean(all_rr_uncued(:,ismember(days,whichdays)),1),nanmean(all_rr_cued(:,ismember(days,whichdays)),1),'Color','k');

% daspect([1 1 1]);
xlabel('Uncued reaching (reaches per sec)');
ylabel('Cued reaching (reaches per sec)');

figure();
colormap('cool');
imagesc(1:colorsUpTo);
title(['Colormap days 1 to ' num2str(colorsUpTo)]);

end