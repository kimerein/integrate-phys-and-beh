function [dprimes,fracs_over_sess,firstBinDprimes]=plotDprimesFromReachRates(varargin)

if length(varargin)==3
    reachrates=varargin{1};
    suppressPlots=varargin{2};
    plotVersusFrac=varargin{3};
    initialDprimes=[];
elseif length(varargin)==4
    reachrates=varargin{1};
    suppressPlots=varargin{2};
    plotVersusFrac=varargin{3};
    initialDprimes=varargin{4};
else
    error('Wrong number of arguments to plotDprimesFromReachRates');
end

% settings.stopPlottingTrialsAfterN=286;
settings.stopPlottingTrialsAfterN=10000; %2500;
settings.binTrialsForAvAcrossSess=true;
settings.binThisManyTrials=10; %50; %200; %50; %4; %10; % somehow this makes dprime bigger, SO sensitive to this
settings.stopPlottingBinsAfterN=200; %60; %55;
settings.furtherBinBins=false; %true; %false; %true;
settings.binThisManyBins=5;
settings.plotVersusFrac=plotVersusFrac; % if is true, will plot dprime versus fraction through session instead of trial count
settings.plotChangeInDprimes=false;

[dprimes,fracs_over_sess]=getdprimes(reachrates,settings.binThisManyTrials,settings);
firstBinDprimes=dprimes(:,1);
if settings.plotChangeInDprimes==true
    if ~isempty(initialDprimes)
        dprimes=dprimes-repmat(initialDprimes,1,size(dprimes,2));
    else
        dprimes=dprimes-repmat(dprimes(:,1),1,size(dprimes,2));
    end
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
    if settings.stopPlottingBinsAfterN>nansum(~isnan(nanmean(approach_alltrials_dprime,1)))
        settings.stopPlottingBinsAfterN=nansum(~isnan(nanmean(approach_alltrials_dprime,1)));
    end
    kstep=ceil(size(cmap,1)/settings.stopPlottingBinsAfterN);
else
    kstep=ceil(size(cmap,1)/nansum(~isnan(nanmean(approach_alltrials_dprime,1))));
end
currbincued=[];
currbinfracs=[];
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
            currbintrialnum=trial_nums(i);
            currbincounter=currbincounter+1;
            if currbincounter==settings.binThisManyBins
                % plot and reset
%                 currbincued=nanmean(currbincued,2);
%                 currbinfracs=nanmean(currbinfracs,2);
                currbincued=currbincued(1:end);
                currbinfracs=currbinfracs(1:end);
                currbincounter=0;
                backup_temp_cued=temp_cued;
                backup_temp_fracs=temp_fracs;
                temp_cued=currbincued;
                temp_fracs=currbinfracs;
                temp_currbintrialnum=currbintrialnum;
                currbincued=[];
                currbinfracs=[];
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
hit_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
fa_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
fracs_over_sess=nan(size(reachrates.alltrials_uncued,1),ceil(size(reachrates.alltrials_uncued,2)/trialBinSize));
currbincounter=0;
currtrial=0;
goAhead=true;
for i=1:ceil(size(reachrates.alltrials_uncued,2)/trialBinSize)
    if goAhead==false
        break
    end
    currbincounter=currbincounter+1;
    takeInds=(i-1)*trialBinSize+1:(i-1)*trialBinSize+trialBinSize;
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
    [dprimes_lasttrial,hit_lasttrial,fa_lasttrial]=calc_dprimes(currbinuncued,currbincued);
    dprimes_over_sess(:,currbincounter)=dprimes_lasttrial;
    hit_over_sess(:,currbincounter)=hit_lasttrial;
    fa_over_sess(:,currbincounter)=fa_lasttrial;
    fracs_over_sess(:,currbincounter)=fracs_av;
end
dprimes_over_sess=dprimes_over_sess(:,1:currbincounter);
hit_over_sess=hit_over_sess(:,1:currbincounter);
fa_over_sess=fa_over_sess(:,1:currbincounter);
fracs_over_sess=fracs_over_sess(:,1:currbincounter);
fi=find(all(isnan(dprimes_over_sess),1),2,'first');
if ~isempty(fi)
    dprimes_over_sess=dprimes_over_sess(:,1:fi(1)-1);
    fracs_over_sess=fracs_over_sess(:,1:fi(1)-1);
end

end

function [dprimes,hit_rates,fa_rates]=calc_dprimes(uncued_events,cued_events)

useBayes=false; % Bayes estimator helps to ameliorate SOME of the shift in d-prime that results from simply having too few trials
flipContingency=false; % if want to get probability that a reach preceded or followed by cue, rather than probability that cue preceded or followed by reach

if flipContingency==true
    disp('DOING FLIP CONTINGENCY!');
    % https://www.researchgate.net/publication/251102295_Corrections_for_extreme_proportions_and_their_biasing_effects_on_estimated_values_of_d_'
    % Stanislaw, Harold, and Natasha Todorov. 1999. "Calculation of Signal Detection Theory Measures." Behavior Research Methods, Instruments, & Computers 31 (1): 137–49. http://link.springer.com/article/10.3758/BF03207704.
    % Could also get
    % Probability that reach preceded by cue
    hit_rates=nansum(cued_events>0,2)./nansum((cued_events>0) + (uncued_events>0),2); % 3 reaches same as 1 reach
%     hit_rates=nansum(cued_events,2)./nansum(nansum([cued_events; uncued_events],1),2);
    % probability that reach followed by cue (this is from the uncued time bin,
    % when precedes cue)
    fa_rates=nansum(uncued_events>0,2)./nansum((cued_events>0) + (uncued_events>0),2); % 3 reaches same as 1 reach
%     fa_rates=nansum(uncued_events,2)./nansum(nansum([cued_events; uncued_events],1),2);
    % Correct for extreme values, i.e., 0 or 1
    n=nansum((cued_events>0) + (uncued_events>0),2);
    if any(hit_rates==0)
        hit_rates(hit_rates==0)=1./(2.*n(hit_rates==0));
    end
    if any(hit_rates==1)
        hit_rates(hit_rates==1)=1-(1./(2.*n(hit_rates==1)));
    end
    if any(fa_rates==0)
        fa_rates(fa_rates==0)=1./(2.*n(fa_rates==0));
    end
    if any(fa_rates==1)
        fa_rates(fa_rates==1)=1-(1./(2.*n(fa_rates==1)));
    end
elseif useBayes==true
    disp('USING BAYES!');
    % This is calculating probability that cue is followed by reach
    hit_rates=(nansum(cued_events>0,2)+1)./(nansum(~isnan(cued_events),2)+2);
    % probability that uncue is followed by reach
    fa_rates=(nansum(uncued_events>0,2)+1)./(nansum(~isnan(uncued_events),2)+2);
else
    disp('Neither Bayes nor flip contingency!');
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
end
dprimes=dprime(hit_rates,fa_rates);

end

function out=dprime(hit_rates,FA_rates)

out=norminv(hit_rates)-norminv(FA_rates);

end
