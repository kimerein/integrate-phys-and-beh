function analyzeDAvsRTchange(alltbt,trialTypes,metadata,DAbaselinewindow,DApostwindowstart,wtt)


timestep=mode(diff(nanmean(alltbt.times,1)));
cueOffsetInds=0;
% cueOffset if didn't realign to start of cue
% cueOffset=-0.16;
% cueOffsetInds=ceil(abs(cueOffset)/timestep);

% calculate reaction times
[~,startCol]=nanmax(nanmean(alltbt.cueZone_onVoff,1));
startCol=startCol-cueOffsetInds;
[nRows, ~] = size(alltbt.reachStarts);
% Preallocate result
firstIdx = nan(nRows, 1);
% Loop over each row
for i = 1:nRows
    % Look for the first column (≥startCol) where value > 0.05
    rel = find(alltbt.reachStarts(i, startCol:end) > 0.05, 1, 'first');
    if ~isempty(rel)
        firstIdx(i) = rel + startCol - 1;
    end
end
alltbt.RTs=(firstIdx-startCol).*timestep;

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
DAmeanInds=x>1.5 & x<=2.4;

% get peak DA minus baseline
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,x>DAbaselinewindow & x<DApostwindowstart),[],2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks,[],2,'omitnan')-mean(alltbt.cued_successChunks,2,'omitnan');
alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,:),[],2,'omitnan');
% alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,:),[],2,'omitnan');
% cued failure peaks
alltbt.cued_failureChunks_peak=max(alltbt.cued_failureChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_failureChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_failureChunks_peak=max(alltbt.cued_failureChunks(:,:),[],2,'omitnan');
% get mean DA
% alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,DAmeanInds),2,'omitnan');
alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,DAmeanInds),2,'omitnan');
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
figure(); scatter(alltbt.cued_successChunks_peak,alltbt.deltaRTs);
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

figure(); scatter(alltbt.cued_allChunks_peak,alltbt.deltaRTs);
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued all reaches');
figure(); scatter(alltbt.cued_allChunks_mean,alltbt.deltaRTs);
xlabel('DA mean minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA means vs. change in RT for cued all reaches');

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_allChunks_peak');
figure(); scatter(alltbt.cued_allChunks_peak_normed,alltbt.deltaRTs);

% figure(); scatter(alltbt.cued_derivChunks_peak,alltbt.deltaRTs);
% xlabel('DA deriv peak -- Z-scored fluorescence units');
% ylabel('Change in RT (trial n+1 minus trial n)');
% title('DA derivative vs. change in RT for cued all reaches');

% test correlations
whichToTake=~isnan(alltbt.cued_successChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for cued successes: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for cued successes: ']); disp(p(1,2));
whichToTake=~isnan(alltbt.cued_successChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA means and deltaRTs for cued successes: ']); disp(r(1,2));
disp(['pval correlation between DA means and deltaRTs for cued successes: ']); disp(p(1,2));

whichToTake=~isnan(alltbt.cued_allChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_allChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for all reaches: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for all reaches: ']); disp(p(1,2));

% session by session cued_allChunks_peak
u=unique(metadata.sessid);
rs=nan(length(u),1); ps=nan(length(u),1);
for i=1:length(u)
[r,p]=corrcoef(alltbt.cued_allChunks_peak(metadata.sessid==u(i) & whichToTake==1),alltbt.deltaRTs(metadata.sessid==u(i) & whichToTake==1));
if size(r,1)<2
continue
end
rs(i)=r(1,2); ps(i)=p(1,2);
end

whichToTake=~isnan(alltbt.cued_allChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_allChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA means and deltaRTs for all reaches: ']); disp(r(1,2));
disp(['pval correlation between DA means and deltaRTs for all reaches: ']); disp(p(1,2));

% whichToTake=~isnan(alltbt.cued_derivChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_derivChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA derivs and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between DA derivs and deltaRTs for all reaches: ']); disp(p(1,2));


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