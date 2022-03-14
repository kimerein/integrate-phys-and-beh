function [dprimes,fracs_over_sess]=plotDprimesFromReachRates(reachrates,suppressPlots,plotVersusFrac)

% settings.stopPlottingTrialsAfterN=286;
settings.stopPlottingTrialsAfterN=10000; %2500;
settings.binTrialsForAvAcrossSess=true;
settings.binThisManyTrials=200; %50; %4; %10; % somehow this makes dprime bigger, SO sensitive to this
settings.stopPlottingBinsAfterN=60; %55;
settings.furtherBinBins=true;
settings.binThisManyBins=5;
settings.plotVersusFrac=plotVersusFrac; % if is true, will plot dprime versus fraction through session instead of trial count
settings.plotChangeInDprimes=true;

[dprimes,fracs_over_sess]=getdprimes(reachrates,settings.binThisManyTrials,settings);
if settings.plotChangeInDprimes==true
    dprimes=dprimes-repmat(dprimes(:,1),1,size(dprimes,2));
end
if suppressPlots==false
    makePlot(settings,dprimes,fracs_over_sess);
    if plotVersusFrac==true
        xlabel('Fraction through session');
    else
        xlabel('Trial #');
    end
    ylabel('d-prime');
end

end

function makePlot(settings,approach_alltrials_dprime,fracs_over_sess)

plotVersusFrac=settings.plotVersusFrac;
figure();
cmap=colormap('cool');
hold on;
k=1;
if ~isempty(settings.stopPlottingBinsAfterN)
    kstep=ceil(size(cmap,1)/settings.stopPlottingBinsAfterN);
else
    kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(approach_alltrials_dprime,1))));
end
currbincued=[];
currbinfracs=[];
currbintrialnum=0;
currbincounter=0;
trial_nums=1:settings.binThisManyTrials:settings.binThisManyTrials*(size(approach_alltrials_dprime,2)+1);
for i=1:size(approach_alltrials_dprime,2) % across trials
    % rows are different sessions, columns are different trials in each
    % session
    % so ACROSS ALL SESSIONS, take each trial in session
    temp_cued=approach_alltrials_dprime(:,i); % reach rate in trial n+i (last trial) of sequence
    temp_fracs=fracs_over_sess(:,i);
    if i==1 
        if plotVersusFrac==true
            xToPlot=nanmean(temp_fracs);
        else
            xToPlot=trial_nums(i);
        end
        scatter(xToPlot,nanmean(temp_cued),[],'k'); % first trial in SESSION, last trial in sequence
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
            currbincued=[currbincued temp_cued];
            currbinfracs=[currbinfracs temp_fracs];
            currbintrialnum=currbintrialnum+trial_nums(i);
            currbincounter=currbincounter+1;
            if currbincounter==settings.binThisManyBins
                % plot and reset
%                 currbincued=nanmean(currbincued,2);
%                 currbinfracs=nanmean(currbinfracs,2);
                currbincued=currbincued(1:end);
                currbinfracs=currbinfracs(1:end);
                currbintrialnum=currbintrialnum/settings.binThisManyBins;
                currbincounter=0;
                backup_temp_cued=temp_cued;
                backup_temp_fracs=temp_fracs;
                temp_cued=currbincued;
                temp_fracs=currbinfracs;
                temp_currbintrialnum=currbintrialnum;
                currbincued=[];
                currbinfracs=[];
                currbintrialnum=0;
                if plotVersusFrac==true
                    xToPlot=nanmean(temp_fracs);
                else
                    xToPlot=temp_currbintrialnum;
                end
                scatter(xToPlot,nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
                hold on;
                % plot mean and s.e. across session
                line([xToPlot xToPlot],...
                     [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
                temp_cued=backup_temp_cued;
                temp_fracs=backup_temp_fracs;
            end
        else
            if plotVersusFrac==true
                xToPlot=nanmean(temp_fracs);
            else
                xToPlot=trial_nums(i);
            end
            scatter(xToPlot,nanmean(temp_cued),[],cmap(k,:),'filled'); % later trials in SESSION, last trial in sequence
            hold on;
            line([xToPlot xToPlot],...
                [nanmean(temp_cued)-nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued))) nanmean(temp_cued)+nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)))],'Color',cmap(k,:),'LineWidth',1);
        end
    end
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

end

function [dprimes_over_sess,fracs_over_sess]=getdprimes(reachrates,trialBinSize,settings)

% get dprimes per average trial in session
if isempty(reachrates)
    return
end
dprimes_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
fracs_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
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
    curr_bin_fracs=reachrates.fracsThroughSess(:,takeInds);
    if all(isnan(currbincued(1:end))) && all(isnan(currbinuncued(1:end)))
        currbincounter=currbincounter-1;
        break
    end
    fracs_av=nanmean(curr_bin_fracs,2);
    dprimes_lasttrial=calc_dprimes(currbinuncued,currbincued);
    dprimes_over_sess(:,currbincounter)=dprimes_lasttrial;
    fracs_over_sess(:,currbincounter)=fracs_av;
end
dprimes_over_sess=dprimes_over_sess(:,1:currbincounter);
fracs_over_sess=fracs_over_sess(:,1:currbincounter);
fi=find(all(isnan(dprimes_over_sess),1),2,'first');
if ~isempty(fi)
    dprimes_over_sess=dprimes_over_sess(:,1:fi(1)-1);
    fracs_over_sess=fracs_over_sess(:,1:fi(1)-1);
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