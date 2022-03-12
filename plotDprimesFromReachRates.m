function plotDprimesFromReachRates(reachrates)

settings.stopPlottingTrialsAfterN=175;
settings.binTrialsForAvAcrossSess=true;
settings.binThisManyTrials=6;
settings.stopPlottingBinsAfterN=[];
settings.furtherBinBins=false;
settings.binThisManyBins=1;

dprimes=getdprimes(reachrates,settings.binThisManyTrials,settings);
makePlot(settings,dprimes);
xlabel('Trial #');
ylabel('d-prime');

end

function makePlot(settings,approach_alltrials_dprime)

figure();
cmap=colormap('cool');
hold on;
k=1;
if ~isempty(settings.stopPlottingBinsAfterN)
    kstep=ceil(size(cmap,1)/settings.stopPlottingBinsAfterN);
else
    kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(approach_alltrials_dprime,1))));
end
currbincued=zeros(size(approach_alltrials_dprime,1),1);
currbintrialnum=0;
currbincounter=0;
trial_nums=1:settings.binThisManyTrials:settings.binThisManyTrials*(size(approach_alltrials_dprime,2)+1);
for i=1:size(approach_alltrials_dprime,2) % across trials
    % rows are different sessions, columns are different trials in each
    % session
    % so ACROSS ALL SESSIONS, take each trial in session
    temp_cued=approach_alltrials_dprime(:,i); % reach rate in trial n+i (last trial) of sequence
    if i==1 
        scatter(trial_nums(i),nanmean(temp_cued),[],'k'); % first trial in SESSION, last trial in sequence
    end
    if ~isempty(settings.stopPlottingBinsAfterN)
        if i<=settings.stopPlottingBinsAfterN
            goAhead=true;
        else
            goAhead=false;
        end
    else
        goAhead=true;
    end
    if goAhead
        if settings.furtherBinBins==true
            currbincued=currbincued+temp_cued;
            currbintrialnum=currbintrialnum+trial_nums(i);
            currbincounter=currbincounter+1;
            if currbincounter==settings.binThisManyBins
                % plot and reset
                currbincued=currbincued/settings.binThisManyBins;
                currbintrialnum=currbintrialnum/settings.binThisManyBins;
                currbincounter=0;
                backup_temp_cued=temp_cued;
                temp_cued=currbincued;
                temp_currbintrialnum=currbintrialnum;
                currbincued=zeros(size(approach_alltrials_dprime,1),1);
                currbintrialnum=0;
                scatter(temp_currbintrialnum,nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
                hold on;
                % plot mean and s.e. across session
                line([temp_currbintrialnum temp_currbintrialnum],...
                     [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                temp_cued=backup_temp_cued;
            end
        else
            scatter(trial_nums(i),nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
            hold on;
            line([trial_nums(i) trial_nums(i)],...
                [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
        end
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

end

function dprimes_over_sess=getdprimes(reachrates,trialBinSize,settings)

% get dprimes per average trial in session
if isempty(reachrates)
    return
end
dprimes_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
currbincounter=0;
currtrial=0;
goAhead=true;
for i=1:size(reachrates.alltrials_uncued,2)
    if goAhead==false
        break
    end
    currbincounter=currbincounter+1;
    takeInds=i:i+trialBinSize-1;
    if any(takeInds>size(reachrates.alltrials_uncued,2))
        takeInds=i:size(reachrates.alltrials_uncued,2);
    end
    currtrial=currtrial+length(takeInds);
    if ~isempty(settings.stopPlottingTrialsAfterN)
        if currtrial>settings.stopPlottingTrialsAfterN
            % stop after this bin
            goAhead=false;
        end
    end            
    currbincued=reachrates.alltrials_cued(:,takeInds); % reach rate in trial n+i (last trial) of sequence
    currbinuncued=reachrates.alltrials_uncued(:,takeInds); % reach rate in trial n+i (last trial) of sequence
    if all(isnan(currbincued(1:end))) && all(isnan(currbinuncued(1:end)))
        currbincounter=currbincounter-1;
        break
    end
    dprimes_lasttrial=calc_dprimes(currbinuncued,currbincued);
    dprimes_over_sess(:,currbincounter)=dprimes_lasttrial;
end
dprimes_over_sess=dprimes_over_sess(:,1:currbincounter);
fi=find(all(isnan(dprimes_over_sess),1),2,'first');
if ~isempty(fi)
    dprimes_over_sess=dprimes_over_sess(:,1:fi-1);
end

end

function dprimes=calc_dprimes(uncued_events,cued_events)

hit_rates=nansum(cued_events>0,2)./nansum(~isnan(cued_events),2);
fa_rates=nansum(uncued_events>0,2)./nansum(~isnan(uncued_events),2);
% closest we can get to 1 or zero is defined by number of trials
ns=nansum(~isnan(cued_events),2);
hit_rates(ns<3)=nan;
fa_rates(ns<3)=nan;
hit_rates(hit_rates==1)=1-(1./ns(hit_rates==1));
hit_rates(hit_rates==0)=0+(1./ns(hit_rates==0));
fa_rates(fa_rates==1)=1-(1./ns(fa_rates==1));
fa_rates(fa_rates==0)=0+(1./ns(fa_rates==0));
dprimes=dprime(hit_rates,fa_rates);

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end
