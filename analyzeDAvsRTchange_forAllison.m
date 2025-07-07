function analyzeDAvsRTchange_forAllison(alltbt,trialTypes,metadata,DAbaselinewindow,DApostwindowstart,wtt,meanWindow)

% Kim is using 
% analyzeDAvsRTchange(alltbt,trialTypes,metadata,1.5,12,[],[2.45 4.2]);
% This sets the baseline for calculating the DA peak as the time window up
% until the arm is outstretched (at 1.5 secs) 
% 5th argument set to 12 secs indicates NO post-reach baseline time window
% According to these arguments, the time window for calculating the DA mean 
% is 2.45 to 4.2 secs, wrt the arm outstretched at 1.5 secs

% if want to further subsample
% otherwise leave wtt empty
if isempty(wtt) % which to take
    wtt='ones(size(alltbt.RTs))==1';
end

timestep=mode(diff(nanmean(alltbt.times,1)));
cueOffsetInds=0;
% only need cueOffset if didn't realign to start of cue
% cueOffset=-0.16;
% cueOffsetInds=ceil(abs(cueOffset)/timestep);

% calculate reaction times
alltbt=getReactionTimes(alltbt,cueOffsetInds,timestep);

alltbt.RTs(alltbt.RTs<0.05)=nan;

% dump artifacts of infrequent LabJack problem
alltbt.cued_successChunks(alltbt.cued_successChunks<-4)=nan;
alltbt.cued_failureChunks(alltbt.cued_failureChunks<-4)=nan;

% get max and min for each
alltbt.cued_successChunks_ma=nanmax(alltbt.cued_successChunks,[],2);
alltbt.cued_successChunks_mi=nanmin(alltbt.cued_successChunks,[],2);
alltbt.cued_failureChunks_ma=nanmax(alltbt.cued_failureChunks,[],2);
alltbt.cued_failureChunks_mi=nanmin(alltbt.cued_failureChunks,[],2);

% DA signals for successes and failures combined
% alltbt.cued_allChunks=alltbt.cued_successChunks;
alltbt.cued_allChunks(~isnan(alltbt.cued_failureChunks_ma),:)=alltbt.cued_failureChunks(~isnan(alltbt.cued_failureChunks_ma),:);

% separate misses (paw never touches pellet) and drops (paw initially grabs
% pellet, then pellet falls before mouse eats pellet)
isdrop=any(alltbt.reachBatch_drop_reachStarts>0.05,2);
ismiss=any(alltbt.reachBatch_miss_reachStarts>0.05,2) & ~any(alltbt.reachBatch_drop_reachStarts>0.05,2);

% plot DA responses 
nTime = size(alltbt.cued_successChunks, 2);
x     = (0:nTime-1) * timestep;
meanSuc = mean(alltbt.cued_successChunks, 1, 'omitnan');
semSuc  = std(alltbt.cued_successChunks, 0, 1, 'omitnan') ./ ...
          sqrt(sum(~isnan(alltbt.cued_successChunks), 1));
meanFail = mean(alltbt.cued_failureChunks, 1, 'omitnan');
semFail  = std(alltbt.cued_failureChunks, 0, 1, 'omitnan') ./ ...
           sqrt(sum(~isnan(alltbt.cued_failureChunks), 1));
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
disp('OMITTING DROPS!')

DAbaseToInd=find(x<DAbaselinewindow,1,'last');
DAbaseAtEnd=find(x>DApostwindowstart,1,'first');
DAmeanInds=x>meanWindow(1) & x<=meanWindow(2);

% plot DA 3 outcomes with baseline subtracted
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
alltbt.cued_successChunks_peak=max(alltbt.cued_successChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_successChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% alltbt.cued_allChunks_peak=max(alltbt.cued_allChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_allChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
alltbt.cued_failureChunks_peak=max(alltbt.cued_failureChunks(:,:),[],2,'omitnan')-mean(alltbt.cued_failureChunks(:,[1:DAbaseToInd DAbaseAtEnd:end]),2,'omitnan');
% get DA means over initial response window
alltbt.cued_failureChunks_mean=mean(alltbt.cued_failureChunks(:,DAmeanInds),2,'omitnan');
alltbt.cued_successChunks_mean=mean(alltbt.cued_successChunks(:,DAmeanInds),2,'omitnan');
% alltbt.cued_allChunks_mean=mean(alltbt.cued_allChunks(:,DAmeanInds),2,'omitnan');

% plot scatter of peaks minus baseline for successes v failures
figure(); scatter(ones(size(alltbt.cued_successChunks_peak))+rand(size(alltbt.cued_successChunks_peak)),alltbt.cued_successChunks_peak,[],[107 76 154]./255);
hold on; scatter(2*ones(size(alltbt.cued_failureChunks_peak))+rand(size(alltbt.cued_failureChunks_peak)),alltbt.cued_failureChunks_peak,[],'r');
title('DA peaks minus baseline, purple is success, red is failure');

% plot scatter of DA means for successes v failures
figure(); scatter(ones(size(alltbt.cued_successChunks_mean))+rand(size(alltbt.cued_successChunks_mean)),alltbt.cued_successChunks_mean,[],[107 76 154]./255);
hold on; scatter(2*ones(size(alltbt.cued_failureChunks_mean))+rand(size(alltbt.cued_failureChunks_mean)),alltbt.cued_failureChunks_mean,[],'r');
title('DA pmeans, purple is success, red is failure');

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
xlabel('DA mean -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA means vs. change in RT for cued successes');

figure(); scatter(alltbt.cued_successChunks_peak,alltbt.deltaRTs,[],[107 76 154]./255); hold on;
scatter(alltbt.cued_failureChunks_peak,alltbt.deltaRTs,[],'r');
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued successes and failures');

% alltbt.cued_allChunks_peak=alltbt.cued_failureChunks_peak;
% alltbt.cued_allChunks_peak(~isnan(alltbt.cued_successChunks_peak))=alltbt.cued_successChunks_peak(~isnan(alltbt.cued_successChunks_peak));

whichToTake=~isnan(alltbt.cued_failureChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
figure(); scatter(alltbt.cued_failureChunks_peak,alltbt.deltaRTs);
xlabel('DA peak minus baseline -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA peaks vs. change in RT for cued failures');
figure(); scatter(alltbt.cued_failureChunks_mean,alltbt.deltaRTs);
xlabel('DA mean -- Z-scored fluorescence units');
ylabel('Change in RT (trial n+1 minus trial n)');
title('DA means vs. change in RT for cued failures');

alltbt.isLongITI=trialTypes.isLongITI;
% shuffle_tbt=shuffleTbtTrialOrder(alltbt,metadata,trialTypes);

% test correlations
whichToTake=~isnan(alltbt.cued_successChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA peaks and deltaRTs for cued successes: ']); disp(r(1,2));
disp(['pval correlation between DA peaks and deltaRTs for cued successes: ']); disp(p(1,2));
% whichToTake=~isnan(alltbt.cued_successChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_successChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA means and deltaRTs for cued successes: ']); disp(r(1,2));
% disp(['pval correlation between DA means and deltaRTs for cued successes: ']); disp(p(1,2));

% whichToTake=~isnan(shuffle_tbt.cued_successChunks_peak) & ~isnan(shuffle_tbt.deltaRTs) & shuffle_tbt.isLongITI==1; % this won't take wtt
% [r,p]=corrcoef(shuffle_tbt.cued_successChunks_peak(whichToTake),shuffle_tbt.deltaRTs(whichToTake));
% disp(['SHUFFLE R correlation between DA peaks and deltaRTs for cued successes: ']); disp(r(1,2));
% disp(['SHUFFLE pval correlation between DA peaks and deltaRTs for cued successes: ']); disp(p(1,2));

whichToTake=~isnan(alltbt.cued_failureChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_failureChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between DA means and deltaRTs for cued failures: ']); disp(r(1,2));
disp(['pval correlation between DA means and deltaRTs for cued failures: ']); disp(p(1,2));

% whichToTake=~isnan(alltbt.cued_allChunks_peak) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_allChunks_peak(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA peaks and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between DA peaks and deltaRTs for all reaches: ']); disp(p(1,2));

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
figure(); bin_edges = getCenteredBinEdges(rs,1);
histogram(rs,bin_edges); xlabel('Pearson corrcoef for successes session by session'); ylabel('Count');

% session by session cued_failureChunks_mean
whichToTake=~isnan(alltbt.cued_failureChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
u=unique(metadata.sessid);
rs=nan(length(u),1); ps=nan(length(u),1);
for i=1:length(u)
    [r,p]=corrcoef(alltbt.cued_failureChunks_mean(metadata.sessid==u(i) & whichToTake==1),alltbt.deltaRTs(metadata.sessid==u(i) & whichToTake==1));
    if size(r,1)<2
        continue
    end
    rs(i)=r(1,2); ps(i)=p(1,2);
end
figure(); bin_edges = getCenteredBinEdges(rs,1);
histogram(rs,bin_edges); xlabel('Pearson corrcoef for failures session by session'); ylabel('Count');

% whichToTake=~isnan(alltbt.cued_allChunks_mean) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_allChunks_mean(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between DA means and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between DA means and deltaRTs for all reaches: ']); disp(p(1,2));

% alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_allChunks_mean');
% figure(); scatter(alltbt.cued_allChunks_mean_normed,alltbt.deltaRTs); title('normed all reaches DA mean');
% whichToTake=~isnan(alltbt.cued_allChunks_mean_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
% [r,p]=corrcoef(alltbt.cued_allChunks_mean_normed(whichToTake),alltbt.deltaRTs(whichToTake));
% disp(['R correlation between normed DA means and deltaRTs for all reaches: ']); disp(r(1,2));
% disp(['pval correlation between normed DA means and deltaRTs for all reaches: ']); disp(p(1,2));

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_successChunks_peak');
figure(); scatter(alltbt.cued_successChunks_peak_normed,alltbt.deltaRTs); title('normed successes DA peaks');
xlabel('DA peak minus baseline -- Normed within each session');
ylabel('Change in RT (trial n+1 minus trial n)');
whichToTake=~isnan(alltbt.cued_successChunks_peak_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_successChunks_peak_normed(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between normed DA peaks and deltaRTs for success reaches: ']); disp(r(1,2));
disp(['pval correlation between normed DA peaks and deltaRTs for success reaches: ']); disp(p(1,2));
[slope,pValue]=fitDeltaRTtoDA(alltbt,whichToTake,'cued_successChunks_peak_normed','deltaRTs');
disp(['Regression coefficient for normed DA peaks and success reaches: ']); disp(slope);
disp(['pval of regression coefficient for normed DA peaks and success reaches: ']); disp(pValue);
alltbt.fidgetSum=sum(alltbt.fidgetData,2,'omitnan');
% keep movement as regressor
[DA_coeff,DA_pValue,mvmt_coeff,mvmt_pValue]=fitDeltaRTtoDA_wMovement(alltbt,whichToTake,'cued_successChunks_peak_normed','deltaRTs','fidgetSum');
% session by session regression
[rs,ps,std_DA,std_deltaRT,ns]=getCoefSessbySess(alltbt,trialTypes,metadata,'cued_successChunks_peak',wtt);
figure(); bin_edges = getCenteredBinEdges(rs,3);
histogram(rs,bin_edges); xlabel('Regression coefficients for successes session by session'); ylabel('Count');
figure(); bin_edges = getCenteredBinEdges(rs(ns>=20),1);
histogram(rs(ns>=20),bin_edges); xlabel('Regression coefficients for successes session by session EXCLUDING SESSIONS W n<20'); ylabel('Count');

alltbt=normalizeWithinEachSession(metadata,alltbt,'cued_failureChunks_mean');
figure(); scatter(alltbt.cued_failureChunks_mean_normed,alltbt.deltaRTs); title('normed failures DA means');
xlabel('DA means -- Normed within each session');
ylabel('Change in RT (trial n+1 minus trial n)');
whichToTake=~isnan(alltbt.cued_failureChunks_mean_normed) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
[r,p]=corrcoef(alltbt.cued_failureChunks_mean_normed(whichToTake),alltbt.deltaRTs(whichToTake));
disp(['R correlation between normed DA means and deltaRTs for failure reaches: ']); disp(r(1,2));
disp(['pval correlation between normed DA means and deltaRTs for failure reaches: ']); disp(p(1,2));

[slope,pValue]=fitDeltaRTtoDA(alltbt,whichToTake,'cued_failureChunks_mean_normed','deltaRTs');
disp(['Regression coefficient for normed DA means and failure reaches: ']); disp(slope);
disp(['pval of regression coefficient for normed DA means and failure reaches: ']); disp(pValue);
alltbt.fidgetSum=sum(alltbt.fidgetData,2,'omitnan');
% keep movement as regressor
[DA_coeff,DA_pValue,mvmt_coeff,mvmt_pValue]=fitDeltaRTtoDA_wMovement(alltbt,whichToTake,'cued_failureChunks_mean_normed','deltaRTs','fidgetSum');
% session by session regression
[rs,ps,std_DA,std_deltaRT,ns]=getCoefSessbySess(alltbt,trialTypes,metadata,'cued_failureChunks_mean',wtt);
figure(); bin_edges = getCenteredBinEdges(rs,1);
histogram(rs,bin_edges); xlabel('Regression coefficients for failures session by session'); ylabel('Count');
figure(); bin_edges = getCenteredBinEdges(rs(ns>=20),1);
histogram(rs(ns>=20),bin_edges); xlabel('Regression coefficients for failures session by session EXCLUDING SESSIONS W n<20'); ylabel('Count');

end

function [rs,ps,std_DA,std_deltaRT,ns]=getCoefSessbySess(alltbt,trialTypes,metadata,whichDAField,wtt)

da=alltbt.(whichDAField);
whichToTake=~isnan(da) & ~isnan(alltbt.deltaRTs) & trialTypes.isLongITI==1 & eval(wtt);
u=unique(metadata.sessid);
rs=nan(length(u),1); ps=nan(length(u),1);
std_DA=nan(length(u),1); std_deltaRT=nan(length(u),1);
ns=nan(length(u),1); 
for i=1:length(u)
    % don't keep movement as regressor
%     [r,p]=fitDeltaRTtoDA(alltbt,metadata.sessid==u(i) & whichToTake==1,whichDAField,'deltaRTs');
    % keep movement as regressor
    if nansum(metadata.sessid==u(i) & whichToTake==1)==0
        rs(i)=nan; ps(i)=nan;
    else
        [r,p,mvmt_coeff,mvmt_pValue]=fitDeltaRTtoDA_wMovement(alltbt,metadata.sessid==u(i) & whichToTake==1,whichDAField,'deltaRTs','fidgetSum');
        rs(i)=r; ps(i)=p;
        temp=alltbt.(whichDAField);
        std_DA(i)=std(temp(metadata.sessid==u(i) & whichToTake==1,:),[],1,'omitnan');
        temp=alltbt.deltaRTs;
        std_deltaRT(i)=std(temp(metadata.sessid==u(i) & whichToTake==1,:),[],1,'omitnan');
        ns(i)=nansum(metadata.sessid==u(i) & whichToTake==1);
    end
end

end

function alltbt=getReactionTimes(alltbt,cueOffsetInds,timestep)

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

function bin_edges = getCenteredBinEdges(data,multiplyBins)
%getCenteredBinEdges Calculates optimal histogram bin edges with one bin
%   centered on zero.
%
%   bin_edges = getCenteredBinEdges(data)
%
%   INPUTS:
%   - data: A numeric vector of data points.
%
%   OUTPUT:
%   - bin_edges: A vector containing the calculated bin edges.

% --- 1. Determine a "Reasonable" Bin Width using Freedman-Diaconis Rule ---

n = length(data);
iqr_val = iqr(data);

% Handle edge case where IQR is zero (e.g., all data points are the same)
if iqr_val == 0
    % Fallback to a simpler rule, like Scott's rule, if IQR is 0
    bin_width = 3.5 * std(data) / (n^(1/3));
else
    % Freedman-Diaconis rule for bin width
    bin_width = 2 * iqr_val / (n^(1/3));
end

% If bin_width is still zero or NaN (e.g., for constant data), set a default.
if bin_width <= 0 || isnan(bin_width)
    bin_width = 1; % A sensible default
end

bin_width=bin_width/multiplyBins;


% --- 2. Construct Bin Edges Centered Around Zero ---

% Find the limits of the data
min_val = min(data);
max_val = max(data);

% Create the positive-side bin edges, starting from the center
pos_edges = (bin_width/2) : bin_width : (max_val + bin_width);

% Create the negative-side bin edges, starting from the center
neg_edges = (-bin_width/2) : -bin_width : (min_val - bin_width);

% Combine and sort the edges to form the complete set of bins.
bin_edges = unique([neg_edges, pos_edges]);

end