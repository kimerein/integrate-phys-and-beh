function plotOutcomeDependentShift(firstTrial_rrs,secondTrial_rrs,c)

doQuiver=1;

figure();
temp_uncued=nanmean(firstTrial_rrs.approach2_alltrials_uncued,2);
temp_cued=nanmean(firstTrial_rrs.approach2_alltrials_cued,2);
line([nanmean(temp_uncued) nanmean(temp_uncued)],[nanmean(temp_cued)-nanstd(temp_cued)./sqrt(length(temp_cued)) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(length(temp_cued))],'Color','k');
hold on;
line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(length(temp_uncued)) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(length(temp_uncued))],[nanmean(temp_cued) nanmean(temp_cued)],'Color','k');
if any(isnan([nanmean(temp_uncued) nanmean(temp_cued)])) || any(isinf([nanmean(temp_uncued) nanmean(temp_cued)]))
    doQuiver=0;
end

temp_uncued=nanmean(secondTrial_rrs.approach2_alltrials_uncued,2);
temp_cued=nanmean(secondTrial_rrs.approach2_alltrials_cued,2);
line([nanmean(temp_uncued) nanmean(temp_uncued)],[nanmean(temp_cued)-nanstd(temp_cued)./sqrt(length(temp_cued)) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(length(temp_cued))],'Color',c);
hold on;
line([nanmean(temp_uncued)-nanstd(temp_uncued)./sqrt(length(temp_uncued)) nanmean(temp_uncued)+nanstd(temp_uncued)./sqrt(length(temp_uncued))],[nanmean(temp_cued) nanmean(temp_cued)],'Color',c);
if any(isnan([nanmean(temp_uncued) nanmean(temp_cued)])) || any(isinf([nanmean(temp_uncued) nanmean(temp_cued)]))
    doQuiver=0;
end

if doQuiver==1
    quiver(nanmean(nanmean(firstTrial_rrs.approach2_alltrials_uncued,2)),nanmean(nanmean(firstTrial_rrs.approach2_alltrials_cued,2)),...
          nanmean(temp_uncued)-nanmean(nanmean(firstTrial_rrs.approach2_alltrials_uncued,2)),nanmean(temp_cued)-nanmean(nanmean(firstTrial_rrs.approach2_alltrials_cued,2)),'Color','k');
end

end



