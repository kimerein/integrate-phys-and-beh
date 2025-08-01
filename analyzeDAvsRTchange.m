function analyzeDAvsRTchange(alltbt,trialTypes,metadata,DAbaselinewindow,DApostwindowstart,wtt,meanWindow)

if isempty(wtt) % which to take
    wtt='ones(size(alltbt.RTs))==1';
end

timestep=mode(diff(nanmean(alltbt.times,1)));
cueOffsetInds=0;
% cueOffset if didn't realign to start of cue
% cueOffset=-0.16;
% cueOffsetInds=ceil(abs(cueOffset)/timestep);

% calculate reaction times
usewhichreachfield='all_reachBatch';
[~,startCol]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
startCol=startCol-cueOffsetInds;
[nRows, ~] = size(alltbt.(usewhichreachfield));
% Preallocate result
firstIdx = nan(nRows, 1);
% Loop over each row
for i = 1:nRows
    % Look for the first column (≥startCol) where value > 0.05
    temp=alltbt.(usewhichreachfield);
    rel = find(temp(i, startCol:end) > 0.05, 1, 'first');
    if ~isempty(rel)
        firstIdx(i) = rel + startCol - 1;
    end
end
alltbt.RTs=(firstIdx-startCol).*timestep;

usewhichreachfield='reachStarts';
[~,startCol]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
startCol=startCol-cueOffsetInds;
[nRows, ~] = size(alltbt.(usewhichreachfield));
% Preallocate result
firstIdx = nan(nRows, 1);
% Loop over each row
for i = 1:nRows
    % Look for the first column (≥startCol) where value > 0.05
    temp=alltbt.(usewhichreachfield);
    rel = find(temp(i, startCol:end) > 0.05, 1, 'first');
    if ~isempty(rel)
        firstIdx(i) = rel + startCol - 1;
    end
end
alltbt.RTs_fromReachStarts=(firstIdx-startCol).*timestep;
% dump trials where reach is too extended in time
alltbt.RTs(abs(alltbt.RTs-alltbt.RTs_fromReachStarts)>1)=nan;

alltbt.RTs(alltbt.RTs<0.05)=nan;

% make consistent with taking cued reaches only
% alltbt.cued_successChunks(alltbt.RTs>3,:)=nan;
% alltbt.cued_failureChunks(alltbt.RTs>3,:)=nan;

% dump artifacts of infrequent LabJack problem
alltbt.cued_successChunks(alltbt.cued_successChunks<-4)=nan;
alltbt.cued_failureChunks(alltbt.cued_failureChunks<-4)=nan;

% get max and min for each
alltbt.cued_successChunks_ma=nanmax(alltbt.cued_successChunks,[],2);
alltbt.cued_successChunks_mi=nanmin(alltbt.cued_successChunks,[],2);
alltbt.cued_failureChunks_ma=nanmax(alltbt.cued_failureChunks,[],2);
alltbt.cued_failureChunks_mi=nanmin(alltbt.cued_failureChunks,[],2);

% DA signals for successes and failures combined
alltbt.cued_allChunks=alltbt.cued_successChunks;
alltbt.cued_allChunks(~isnan(alltbt.cued_failureChunks_ma),:)=alltbt.cued_failureChunks(~isnan(alltbt.cued_failureChunks_ma),:);

% DA derivatives for successes and failures combined
% alltbt.cued_derivChunks=[diff(alltbt.cued_allChunks,1,2) nan(size(alltbt.cued_allChunks,1),1)];

% plot DA responses 
nTime = size(alltbt.cued_successChunks, 2);
x     = (0:nTime-1) * timestep;
meanSuc = mean(alltbt.cued_successChunks, 1, 'omitnan');
semSuc  = std(alltbt.cued_successChunks, 0, 1, 'omitnan') ./ ...
          sqrt(sum(~isnan(alltbt.cued_successChunks), 1));
meanFail = mean(alltbt.cued_failureChunks, 1, 'omitnan');
semFail  = std(alltbt.cued_failureChunks, 0, 1, 'omitnan') ./ ...
           sqrt(sum(~isnan(alltbt.cued_failureChunks), 1));
figure; hold on;
% success (green)
plot(x, meanSuc, 'g-', 'LineWidth', 1.5);
patch([x, fliplr(x)], ...
      [meanSuc+semSuc, fliplr(meanSuc-semSuc)], ...
      'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% failure (red)
plot(x, meanFail, 'r-', 'LineWidth', 1.5);
patch([x, fliplr(x)], ...
      [meanFail+semFail, fliplr(meanFail-semFail)], ...
      'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('dLight signal');
legend({'Success mean','Success ± SEM','Failure mean','Failure ± SEM'}, ...
       'Location','Best');
title('Cued Success vs. Failure (±SEM)');
hold off;

% separate misses (paw never touches pellet) and drops (paw initially grabs
% pellet, then pellet falls before mouse eats pellet)
isdrop=any(alltbt.reachBatch_drop_reachStarts>0.05,2);
ismiss=any(alltbt.reachBatch_miss_reachStarts>0.05,2) & ~any(alltbt.reachBatch_drop_reachStarts>0.05,2);

% plot DA 3 outcomes
figure; hold on;
% success (green)
plot(x, meanSuc, 'g-', 'LineWidth', 1.5);
patch([x, fliplr(x)], [meanSuc+semSuc, fliplr(meanSuc-semSuc)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% drop (red)
meanFail = mean(alltbt.cued_failureChunks(isdrop,:), 1, 'omitnan');
semFail  = std(alltbt.cued_failureChunks(isdrop,:), 0, 1, 'omitnan') ./ ...
           sqrt(sum(~isnan(alltbt.cued_failureChunks(isdrop,:)), 1));
plot(x, meanFail, 'r-', 'LineWidth', 1.5);
patch([x, fliplr(x)], [meanFail+semFail, fliplr(meanFail-semFail)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% miss (cyan)
meanFail = mean(alltbt.cued_failureChunks(ismiss,:), 1, 'omitnan');
semFail  = std(alltbt.cued_failureChunks(ismiss,:), 0, 1, 'omitnan') ./ ...
           sqrt(sum(~isnan(alltbt.cued_failureChunks(ismiss,:)), 1));
plot(x, meanFail, 'c-', 'LineWidth', 1.5);
patch([x, fliplr(x)], [meanFail+semFail, fliplr(meanFail-semFail)], 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('dLight signal');
legend({'Success mean','Success ± SEM','Drop mean','Drop ± SEM','Miss mean','Miss ± SEM'}, ...
       'Location','Best');
title('Cued Success vs. Drop vs. Miss (±SEM)');
hold off;

% discard drops, may confuse things
alltbt.cued_failureChunks(isdrop==1,:)=nan;
alltbt.cued_allChunks(isdrop==1,:)=nan;

% plot DA derivs 
% nTime = size(alltbt.cued_derivChunks, 2);
% x     = (0:nTime-1) * timestep;
% meanSuc = mean(alltbt.cued_derivChunks, 1, 'omitnan');
% semSuc  = std(alltbt.cued_derivChunks, 0, 1, 'omitnan') ./ ...
%           sqrt(sum(~isnan(alltbt.cued_derivChunks), 1));
% figure; hold on;
% plot(x, meanSuc, 'k-', 'LineWidth', 1.5);
% patch([x, fliplr(x)], ...
%       [meanSuc+semSuc, fliplr(meanSuc-semSuc)], ...
%       'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% title('Derivative of DA (±SEM)');
% hold off;

DAbaseToInd=find(x<DAbaselinewindow,1,'last');
DAbaseAtEnd=find(x>DApostwindowstart,1,'first');
% DAmeanInds=x>2.4 & x<=4.41;
DAmeanInds=x>meanWindow(1) & x<=meanWindow(2);

% plot DA 3 outcomes
meanSuc = mean(alltbt.cued_successChunks-repmat(mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan'),1,size(alltbt.cued_successChunks,2)), 1, 'omitnan');
semSuc  = std(alltbt.cued_successChunks-repmat(mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan'),1,size(alltbt.cued_successChunks,2)), 0, 1, 'omitnan') ./ ...
          sqrt(sum(~isnan(alltbt.cued_successChunks), 1));
figure; hold on;
% success (green)
plot(x, meanSuc, 'g-', 'LineWidth', 1.5);
patch([x, fliplr(x)], [meanSuc+semSuc, fliplr(meanSuc-semSuc)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% miss (cyan)
meanFail = mean(alltbt.cued_failureChunks(ismiss,:)-repmat(mean(alltbt.cued_failureChunks(ismiss,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan'),1,size(alltbt.cued_failureChunks,2)), 1, 'omitnan');
semFail  = std(alltbt.cued_failureChunks(ismiss,:)-repmat(mean(alltbt.cued_failureChunks(ismiss,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan'),1,size(alltbt.cued_failureChunks,2)), 0, 1, 'omitnan') ./ ...
           sqrt(sum(~isnan(alltbt.cued_failureChunks(ismiss,:)), 1));
plot(x, meanFail, 'c-', 'LineWidth', 1.5);
patch([x, fliplr(x)], [meanFail+semFail, fliplr(meanFail-semFail)], 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('dLight signal');
legend({'Success mean','Success ± SEM','Drop mean','Drop ± SEM','Miss mean','Miss ± SEM'}, ...
       'Location','Best');
title('Cued Success vs. Miss (±SEM) MINUS BASELINE');
hold off;

% get peak DA minus baseline
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,x>DAbaselinewindow & x<DApostwindowstart),[],2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks,[],2,'omitnan')-mean(alltbt.cued_successChunks,2,'omitnan');
alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,x>DAbaselinewindow & x<DApostwindowstart),[],2,'omitnan')-mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,:),[],2,'omitnan');
% alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,:),[],2,'omitnan');
% cued failure peaks
% alltbt.cued_failureChunks_peak=max(alltbt.cued_failureChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_failureChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_failureChunks_peak=mean(alltbt.cued_failureChunks(:,DAmeanInds),2,'omitnan');

% alltbt.cued_failureChunks_peak=max(alltbt.cued_failureChunks(:,:),[],2,'omitnan');
% get mean DA
% alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,DAmeanInds),2,'omitnan');
alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,DAmeanInds),2,'omitnan');
% alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,DAmeanInds),2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,DAmeanInds),2,'omitnan')-mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% get derivative max
% alltbt.cued_derivChunks_peak=max(alltbt.cued_derivChunks(:,:),[],2,'omitnan');

% plot scatter of peaks minus baseline for successes v failures
figure(); scatter(ones(size(alltbt.cued_successChunks_peak))+rand(size(alltbt.cued_successChunks_peak)),alltbt.cued_successChunks_peak,[],[107 76 154]./255);
title('cued success peaks minus baseline');
figure(); scatter(ones(size(alltbt.cued_failureChunks_peak))+rand(size(alltbt.cued_failureChunks_peak)),alltbt.cued_failureChunks_peak,[],'r');
title('cued failure peaks minus baseline');

% get change in RT from trial to trial = deltaRTs
% put RT_n+1 - RT_n at the nth position
alltbt.deltaRTs=alltbt.RTs(2:end)-alltbt.RTs(1:end-1);
alltbt.deltaRTs=[alltbt.deltaRTs; nan];

% then exclude any deltaRTs across session boundaries
for i=2:size(alltbt.deltaRTs,1)
    if metadata.sessid(i)~=metadata.sessid(i-1)
        alltbt.deltaRTs(i)=nan;
        alltbt.deltaRTs(i-1)=nan;
    end
end

% plot DA signal vs. deltaRTs
whichToTake=~isnan(alltbt.cued_successChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
figure(); scatter(alltbt.cued_successChunks_peak(),alltbt.deltaRTs);
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued successes');
figure(); scatter(alltbt.cued_successChunks_mean,alltbt.deltaRTs);
xlabel('DA mean minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA means vs. change in RT for cued successes');

figure(); scatter(alltbt.cued_successChunks_peak,alltbt.deltaRTs,[],[107 76 154]./255); hold on;
scatter(alltbt.cued_failureChunks_peak,alltbt.deltaRTs,[],'r');
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued successes and failures');
alltbt.cued_allChunks_peak=alltbt.cued_failureChunks_peak;
alltbt.cued_allChunks_peak(~isnan(alltbt.cued_successChunks_peak))=alltbt.cued_successChunks_peak(~isnan(alltbt.cued_successChunks_peak));

whichToTake=~isnan(alltbt.cued_failureChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
figure(); scatter(alltbt.cued_failureChunks_peak,alltbt.deltaRTs);
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued failures');

figure(); scatter(alltbt.cued_allChunks_peak,alltbt.deltaRTs);
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued all reaches');
figure(); scatter(alltbt.cued_allChunks_mean,alltbt.deltaRTs);
xlabel('DA mean minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA means vs. change in RT for cued all reaches');

% figure(); scatter(alltbt.cued_derivChunks_peak,alltbt.deltaRTs);
% xlabel('DA deriv peak -- Z-scored fluorescence units');
% ylabel('Change in RT (trial n+1 minus trial n)');
% title('DA derivative vs. change in RT for cued all reaches');

alltbt.isLongITI=trialTypes.isLongITI;
% shuffle_tbt=shuffleTbtTrialOrder(alltbt,metadata,trialTypes);

% test correlations
whichToTake=~isnan(alltbt.cued_successChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for cued successes: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for cued successes: ']); disp(p(1,2));
whichToTake=~isnan(alltbt.cued_successChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA means and deltaRTs for cued successes: ']); disp(r(1,2));
disp(['pval correlation between DA means and deltaRTs for cued successes: ']); disp(p(1,2));

% whichToTake=~isnan(shuffle_tbt.cued_successChunks_peak) & ~isnan(shuffle_tbt.deltaRTs) & shuffle_tbt.isLongITI==1; % this won't take wtt
% [r,p]=corrcoef(shuffle_tbt.cued_successChunks_peak(whichToTake),shuffle_tbt.deltaRTs(whichToTake));
% disp(['SHUFFLE R correlation between DA peaks and deltaRTs for cued successes: ']); disp(r(1,2));
% disp(['SHUFFLE pval correlation between DA peaks and deltaRTs for cued successes: ']); disp(p(1,2));

whichToTake=~isnan(alltbt.cued_failureChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_failureChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for cued failures: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for cued failures: ']); disp(p(1,2));

whichToTake=~isnan(alltbt.cued_allChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_allChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for all reaches: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for all reaches: ']); disp(p(1,2));

% session by session cued_successChunks_peak
whichToTake=~isnan(alltbt.cued_successChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
u=unique(metadata.sessid);
rs=nan(length(u),1); ps=nan(length(u),1);
for i=1:length(u)
[r,p]=corrcoef(alltbt.cued_successChunks_peak(metadata.sessid==u(i) & whichToTake==1),alltbt.deltaRTs(metadata.sessid==u(i) & whichToTake==1));
if size(r,1)<2
continue
end
rs(i)=r(1,2); ps(i)=p(1,2);
end

figure(); 
histogram(rs,100);

whichToTake=~isnan(alltbt.cued_allChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_allChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA means and deltaRTs for all reaches: ']); disp(r(1,2));
disp(['pval correlation between DA means and deltaRTs for all reaches: ']); disp(p(1,2));

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_allChunks_mean');
figure(); scatter(alltbt.cued_allChunks_mean_normed,alltbt.deltaRTs); title('normed all reaches DA mean');
whichToTake=~isnan(alltbt.cued_allChunks_mean_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_allChunks_mean_normed(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between normed DA means and deltaRTs for all reaches: ']); disp(r(1,2));
disp(['pval correlation between normed DA means and deltaRTs for all reaches: ']); disp(p(1,2));

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_successChunks_peak');
figure(); scatter(alltbt.cued_successChunks_peak_normed,alltbt.deltaRTs); title('normed success reaches DA peak');
whichToTake=~isnan(alltbt.cued_successChunks_peak_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_peak_normed(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between normed DA peaks and deltaRTs for success reaches: ']); disp(r(1,2));
disp(['pval correlation between normed DA peaks and deltaRTs for success reaches: ']); disp(p(1,2));
fitDeltaRTtoDA(alltbt,whichToTake,'cued_successChunks_peak_normed','deltaRTs');
alltbt.fidgetSum=sum(alltbt.fidgetData,2,'omitnan');
[DA_coeff,DA_pValue,mvmt_coeff,mvmt_pValue]=fitDeltaRTtoDA_wMovement(alltbt,whichToTake,'cued_successChunks_peak_normed','deltaRTs','fidgetSum');
[rs,ps]=getCoefSessbySess(alltbt,trialTypes,metadata,'cued_successChunks_peak',wtt);

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_failureChunks_peak');
figure(); scatter(alltbt.cued_failureChunks_peak_normed,alltbt.deltaRTs); title('normed failure reaches DA peak');
whichToTake=~isnan(alltbt.cued_failureChunks_peak_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_failureChunks_peak_normed(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between normed DA peaks and deltaRTs for failure reaches: ']); disp(r(1,2));
disp(['pval correlation between normed DA peaks and deltaRTs for failure reaches: ']); disp(p(1,2));

% alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_allChunks_peak');
% figure(); scatter(alltbt.cued_allChunks_peak_normed,alltbt.deltaRTs); title('normed all reaches DA peak');
% whichToTake=~isnan(alltbt.cued_allChunks_peak_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_allChunks_peak_normed(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between normed DA peaks and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between normed DA peaks and deltaRTs for all reaches: ']); disp(p(1,2));

% whichToTake=~isnan(alltbt.cued_derivChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_derivChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA derivs and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between DA derivs and deltaRTs for all reaches: ']); disp(p(1,2));

% new_baseline_window_sec = 300; 
% alltbt.Zscored_DA=alltbt.cued_successChunks;
% alltbt=recalculateZScore(alltbt, metadata, new_baseline_window_sec);
% trial_to_plot = find(all(~isnan(alltbt.Zscored_DA),2));
% trial_to_plot=trial_to_plot(randperm(length(trial_to_plot)));
% trial_to_plot=trial_to_plot(1);
% figure; hold on;
% plot(alltbt.Zscored_DA(trial_to_plot, :), 'b-', 'DisplayName', 'Original Z-Score');
% plot(alltbt.Zscored_DA_recalculated(trial_to_plot, :), 'r--', 'DisplayName', 'Recalculated Z-Score');
% title(['Comparison for Trial ' num2str(trial_to_plot)]);
% xlabel('Timepoints in Trial');
% ylabel('Z-Scored Dopamine');
% legend;
% hold off;
% alltbt.Zscored_DA_mean=mean(alltbt.Zscored_DA_recalculated(:,DAmeanInds),2,'omitnan');
% whichToTake=~isnan(alltbt.cued_successChunks_mean) & ~isnan(alltbt.Zscored_DA_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.Zscored_DA_mean(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA means and deltaRTs for cued successes RE-Z-SCORED: ']); disp(r(1,2));
% disp(['pval correlation between DA means and deltaRTs for cued successes RE-Z-SCORED: ']); disp(p(1,2));


end

function [rs,ps]=getCoefSessbySess(alltbt,trialTypes,metadata,whichDAField,wtt)

da=alltbt.(whichDAField);
whichToTake=~isnan(da) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
u=unique(metadata.sessid);
rs=nan(length(u),1); ps=nan(length(u),1);
for i=1:length(u)
    [r,p]=fitDeltaRTtoDA(alltbt,metadata.sessid==u(i) & whichToTake==1,whichDAField,'deltaRTs');
    rs(i)=r; ps(i)=p;
end

end

function [DA_coeff,DA_pValue,mvmt_coeff,mvmt_pValue]=fitDeltaRTtoDA_wMovement(alltbt,whichToTake,fieldDA,field2,fieldMvmt)

tempDA=alltbt.(fieldDA);
deltaRT=alltbt.(field2);
tempMvmt=alltbt.(fieldMvmt);
DA_data = tempDA(whichToTake,:);
Mvmt_data = tempMvmt(whichToTake,:);
RT_data = deltaRT(whichToTake,:);
tbl = table(DA_data, Mvmt_data, RT_data, 'VariableNames', {'DA', 'Mvmt', 'deltaRT'});
mdl = fitlm(tbl, 'deltaRT ~ DA + Mvmt');
DA_coeff = mdl.Coefficients.Estimate('DA');
DA_pValue = mdl.Coefficients.pValue('DA');
mvmt_coeff = mdl.Coefficients.Estimate('Mvmt');
mvmt_pValue = mdl.Coefficients.pValue('Mvmt');

end

function [slope,pValue]=fitDeltaRTtoDA(alltbt,whichToTake,fieldDA,field2)

tempDA=alltbt.(fieldDA);
temp2=alltbt.(field2);
mdl=fitlm(tempDA(whichToTake,:),temp2(whichToTake,:));
% Get the table of coefficients
coeff_table = mdl.Coefficients;
% Extract the slope (the 2nd row, 'Estimate' column)
slope = coeff_table.Estimate(2);
% Extract the p-value for the slope (the 2nd row, 'pValue' column)
pValue = coeff_table.pValue(2);
% Display the results
% fprintf('The linear regression coefficient (slope) is: %.4f\n', slope);
% fprintf('The p-value for the slope is: %e\n', pValue);

end

function shuffle_tbt=shuffleTbtTrialOrder(alltbt,metadata,trialTypes)

u=unique(metadata.sessid);
shuffle_tbt=alltbt; 
thissize=size(alltbt.cueZone_onVoff,1);
fie=fieldnames(alltbt);
for i=1:length(u)
    f=find(metadata.sessid==u(i));
    rapper=randperm(length(f));
    for j=1:length(fie)
        if size(alltbt.(fie{j}),1)==thissize
            temp=alltbt.(fie{j});
            tempshuf=shuffle_tbt.(fie{j});
            tempshuf(f,:)=temp(f(rapper),:);
            shuffle_tbt.(fie{j})=tempshuf;
        end
    end
end
shuffle_tbt.deltaRTs=shuffle_tbt.RTs(2:end)-shuffle_tbt.RTs(1:end-1);
shuffle_tbt.deltaRTs=[shuffle_tbt.deltaRTs; nan];

end

function alltbt=normalizeWithinEachSession(metadata,alltbt,whichField)

% Assume `metadata.sessid` and `alltbt.cued_successChunks_ma` already exist,
% and both have the same number of rows (one row per trial).

% Extract session IDs and dopamine values:
sessids = metadata.sessid;                           % [numTrials × 1]
daValues = alltbt.(whichField);             % [numTrials × 1]

% Preallocate an array to hold the median for each trial’s session:
maxPerTrial = nan(size(daValues));
minPerTrial = nan(size(daValues));

% Loop over each unique session:
uniqueSess = unique(sessids);
for k = 1:numel(uniqueSess)
    s = uniqueSess(k);
    idx = (sessids == s);                   % logical indices of trials in session s
    maxVal = nanmax( daValues(idx) );     % median DA for that session
    maxPerTrial(idx) = maxVal;         % assign to all rows of that session
    minVal = nanmin( daValues(idx) );     % median DA for that session
    minPerTrial(idx) = minVal;         % assign to all rows of that session
end

% Save as a new field in alltbt:
alltbt.([whichField '_normed'])=alltbt.(whichField)-minPerTrial;
maxPerTrial=maxPerTrial-minPerTrial;
alltbt.([whichField '_normed'])=alltbt.([whichField '_normed'])./maxPerTrial;

end

function alltbt = recalculateZScore(alltbt, metadata, newBaselineDuration_sec)
% NOTE THAT THIS FUNCTION DOESN'T REALLY DO WHAT IT SAYS
% THIS IS A PURELY APPROXIMATE RE-SCORING, TO AVOID HAVING TO REPROCESS ALL
% DATA
%recalculateZScore Re-Z-scores dopamine data based on an approximate number of initial trials.
%
%   alltbt = recalculateZScore(alltbt, metadata, newBaselineDuration_sec)
%
%   This function recalculates the Z-score for dopamine data on a per-session
%   basis. It approximates a baseline duration by assuming each trial is 12
%   seconds long, and uses the first N trials of each session to compute a
%   new mean and standard deviation for Z-scoring.
%
%   INPUTS:
%   - alltbt: A structure containing trial-based data. Must include the field:
%       - 'Zscored_DA': A matrix of Z-scored dopamine data where rows are
%         trials and columns are timepoints.
%   - metadata: A structure containing metadata. Must include the field:
%       - 'sessid': A vector where each element corresponds to a trial (row)
%         in 'alltbt' and gives the unique integer ID for that session.
%   - newBaselineDuration_sec: A scalar specifying the desired duration of the
%     new baseline window in seconds (e.g., 300 for the first 5 minutes).
%
%   OUTPUT:
%   - alltbt: The input structure with a new field added:
%       - 'Zscored_DA_recalculated': The re-Z-scored dopamine data.

% --- Input Validation ---
if ~isfield(alltbt, 'Zscored_DA')
    error('Input structure ''alltbt'' must contain the field ''Zscored_DA''.');
end
if ~isfield(metadata, 'sessid')
    error('Input structure ''metadata'' must contain the field ''sessid''.');
end
fprintf('Starting re-Z-scoring process based on initial trials...\n');
% --- Calculate Approximate Number of Baseline Trials ---
% Assume each trial is 12 seconds long
TRIAL_DURATION_SEC = 12;
num_baseline_trials = round(newBaselineDuration_sec / TRIAL_DURATION_SEC);
% fprintf('Using the first ~%d trials of each session as the baseline.\n', num_baseline_trials);
% Create a new field to store the recalculated data, preserving the original
alltbt.Zscored_DA_recalculated = nan(size(alltbt.Zscored_DA));
% Get the list of unique session IDs
unique_sessions = unique(metadata.sessid);
num_sessions = length(unique_sessions);
% --- Loop Through Each Session ---
for i = 1:num_sessions
    current_sessid = unique_sessions(i);
%     fprintf('Processing Session ID: %d (%d of %d)\n', current_sessid, i, num_sessions);

    % Find the indices of all trials belonging to the current session
    session_trial_indices = find(metadata.sessid == current_sessid);
    
    % Extract all dopamine data for this session
    session_DA_data = alltbt.Zscored_DA(session_trial_indices, :);
    
    % Determine how many trials are actually available in this session
    num_trials_in_session = size(session_DA_data, 1);
    
    % Use the smaller of the two values to avoid errors in short sessions
    actual_baseline_trials = min(num_baseline_trials, num_trials_in_session);
    
    if actual_baseline_trials < 1
        warning('Session %d has no trials to use for baseline. Skipping.', current_sessid);
        continue;
    end

    % --- Calculate New Baseline Mean and Std ---
    % Extract the first N trials for the baseline calculation
    baseline_data_matrix = session_DA_data(1:actual_baseline_trials, :);
    
    % Reshape the baseline data into a single vector to get the overall mean and std
    baseline_data_vector = baseline_data_matrix(:);
    
    new_mean = mean(baseline_data_vector, 'omitnan');
    new_std = std(baseline_data_vector, 0, 'omitnan');

    % Avoid division by zero if baseline is flat
    if new_std == 0
        warning('Standard deviation of baseline for session %d is zero. Z-scores will be NaN.', current_sessid);
        new_std = NaN; % This will result in NaN Z-scores
    end

    % --- Re-Z-score All Data for the Current Session ---
    recalculated_session_DA = (session_DA_data - new_mean) / new_std;

    % Place the re-calculated Z-scores back into the new field
    alltbt.Zscored_DA_recalculated(session_trial_indices, :) = recalculated_session_DA;
end

fprintf('Re-Z-scoring complete.\n');

end