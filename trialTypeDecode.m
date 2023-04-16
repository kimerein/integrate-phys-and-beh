function trialTypeDecode(tensor,allLabels,timepoints_for_tensor)

disp(['Rows of tensor must be ']);
disp('Cued grp 1 cells');
disp('Cued grp 2 cells');
disp('Uncued grp 1 cells');
disp('Uncued grp 2 cells');

timebin1=[0 8*0.06*5]; % in secs
timebin2=[8*0.06*5 18*0.06*5]; % in secs
timebin3=[8*0.06*5 18*0.06*5]; % in secs

% Fill nans w zeros
tensor(isnan(tensor))=0;

% Get inds in tensor corresponding to these time bins
indbin1=find(timepoints_for_tensor>=timebin1(1),1,'first'):find(timepoints_for_tensor<=timebin1(2),1,'last');
indbin2=find(timepoints_for_tensor>=timebin2(1),1,'first'):find(timepoints_for_tensor<=timebin2(2),1,'last');
indbin3=find(timepoints_for_tensor>=timebin3(1),1,'first'):find(timepoints_for_tensor<=timebin3(2),1,'last');

% Get cue axis values
cue_axis_vals_part1=mean(mean(tensor(1:2,indbin1,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor(3:4,indbin1,:),1,'omitnan'),2,'omitnan');
cue_axis_vals_part2=mean(mean(tensor(1:2,indbin2,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor(3:4,indbin2,:),1,'omitnan'),2,'omitnan');
cue_axis_vals_part2=-cue_axis_vals_part2;
cue_axis_vals=cue_axis_vals_part1-cue_axis_vals_part2;

% Get outcome axis values
out_axis_vals=mean(mean(tensor([1 3],indbin3,:),1,'omitnan'),2,'omitnan')-mean(mean(tensor([2 4],indbin3,:),1,'omitnan'),2,'omitnan');

% Plot
c{1}='b'; c{2}='g'; c{3}='r'; c{4}='m';
figure();
unique_allLabels=unique(allLabels);
meanOfAll=[nanmean(squeeze(out_axis_vals(:,:,:))) nanmean(squeeze(cue_axis_vals(:,:,:)))];
for i=1:length(unique_allLabels)
    scatter(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i))),squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i))),[],c{i}); hold on;
    scatter(nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i)))),nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i)))),[],c{i},'filled');
    line([meanOfAll(1) nanmean(squeeze(out_axis_vals(:,:,allLabels==unique_allLabels(i))))],[meanOfAll(2) nanmean(squeeze(cue_axis_vals(:,:,allLabels==unique_allLabels(i))))],'LineWidth',2,'Color',c{i});
end
xlabel('Outcome axis'); ylabel('Cue axis');
legend({'cued succ','','','uncued succ','','','cued fail','','','uncued fail','',''});

end