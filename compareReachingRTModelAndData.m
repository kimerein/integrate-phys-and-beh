function compareReachingRTModelAndData()

% Model
% Rate-based term is fit from data
% Depends on behavioral event and ITI
% RPE-based term is built from principle
% Need to fit alpha, the learning rate



%n_steps_to_update=5000;
% sample reaction times from current RT PDF
curr_rts=randsample(bin_centers(bins),n_steps_to_update,true,rt_pdf(rts,bins));
suppressOutput=true;
[rt_pdf_out,bins]=update_RPE_term(rt_pdf(rts,bins),bins,curr_rts(1),alpha,false);
% for i=2:n_steps_to_update-1
%     if mod(i,200)==0
%         disp(i);
%     end
%     [rt_pdf_out,bins]=update_RPE_term(rt_pdf_out,bins,curr_rts(i),alpha,suppressOutput);
% end
% [rt_pdf_out,bins]=update_RPE_term(rt_pdf_out,bins,curr_rts(8),alpha,false);










sum2exp = @(x,tau1,tau2,a,b) a*exp(-x/tau1)+b*exp(-x/tau2);
getTimeIndShift = @(x_step,timeShift) floor(timeShift/x_step);
applyDriveBasedTerm = @(reach_pdf,reach_suppress,timeIndShift) subplus(reach_pdf(1:end-timeIndShift)-reach_suppress(timeIndShift+1:end));


rt_pdf = @(rts,bins) histcounts(rts,bins)./nansum(histcounts(rts,bins));
bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);
nTrials=1000;
fake_rts1=nan(1,nTrials);
fake_rts2=nan(1,nTrials);
fake_rts3=nan(1,nTrials);
fake_rts4=nan(1,nTrials);
fake_rts5=nan(1,nTrials);
for i=1:nTrials
    % Get expected number of reaches per trial from pdf?
    sampleFromFunc1=randsample(bin_centers(bins),1,true,rt_pdf(rts,bins));
    sampleFromFunc2=randsample(bin_centers(bins),5,true,rt_pdf(rts,bins));
    sampleFromFunc3=randsample(bin_centers(bins),10,true,rt_pdf(rts,bins));
    sampleFromFunc4=randsample(bin_centers(bins),20,true,rt_pdf(rts,bins));
    sampleFromFunc5=randsample(bin_centers(bins),40,true,rt_pdf(rts,bins));
    % Get "reaction time" of first reach
    fake_rts1(i)=min(sampleFromFunc1);
    fake_rts2(i)=min(sampleFromFunc2);
    fake_rts3(i)=min(sampleFromFunc3);
    fake_rts4(i)=min(sampleFromFunc4);
    fake_rts5(i)=min(sampleFromFunc5);
end
histo_nbins=plotHist(fake_rts1,fake_rts2,2000,'Sample and get RTs','RT','c','m');
plotHist(fake_rts3,fake_rts4,histo_nbins,'Sample and get RTs','RT','k','r');
plotHist(fake_rts5,fake_rts4,histo_nbins,'Sample and get RTs','RT','b','r');



% Question 1: Can I fit RTs as Poisson?
% As sum of 2 Poissons plus a baseline rate, yes
% x=0:0.07:40;
% poisspdf(x,17)/30+poisspdf(x,6)/4+0.00135
% Question 1: Can I fit RTs as gammas?
% As sum of 3 gammas, yes

% Plot reaction times distribution
% plot_rt_pdf=true;
% histo_nbins=200; % number of bins for reaction time histogram
% backup_histo_nbins=histo_nbins;
% if plot_rt_pdf==true
%     for i=1:length(dataset.allTrialsSequence_RT_trial1InSeq)
%         [histo_nbins,counts1,counts2]=plotHist(dataset.allTrialsSequence_RT_trial1InSeq{i},dataset.allTrialsSequence_RT_trialiInSeq{i},histo_nbins,['Reaction Times all trials reference: trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)','k','r');
%     end
%     % Attempt to fit gamma 
%     temp=dataset.allTrialsSequence_RT_trial1InSeq{i};
%     %lambdahat=poissfit(temp(~isnan(temp))/0.07);
%     figure(); 
%     [n,x]=cityscape_hist(counts1,histo_nbins);
%     plot(x,n./nansum(n),'Color','k');
%     %plot(x/0.07,n./nansum(n),'Color','k');
%     %hold on;
%     %plot(x/0.07,poisspdf(x/0.07,lambdahat),'Color','b');
%     y=gampdf(x,3,1/6)/1;
%     y2=gampdf(x,3,5)/2;
%     y3=gampdf(x,3,1)/4;
%     hold on;
%     plot(x,(y)./nansum(y),'Color','m');
%     plot(x,(y2)./nansum(y2),'Color','c');
%     plot(x,(y3)./nansum(y3),'Color','g');
%     plot(x,(y+y2+y3)./nansum(y+y2+y3),'Color','b');
%     
%     temp=dataset.event_name;
%     temp(regexp(temp,'_'))=' ';
%     for i=1:length(dataset.event_RT_trial1InSeq)
%         histo_nbins=plotHist(dataset.event_RT_trial1InSeq{i},dataset.event_RT_trialiInSeq{i},histo_nbins,['Reaction Times fx of ' temp ': trial 1 (black) vs ' num2str(dataset.nInSequence(i)-1) ' later (magenta)'],'RT (sec)','k','r');
%     end
% end












% Drive-based term
% Acts primarily on reach rate
% Which will impact RT measurements



















% Model change in RT as a shift in rate of a Poisson process

% From relative suppression of reaching after mouse drops pellet wrt
% missing pellet, which can be fit as sum of 2 exponentials
sum2exp = @(x,tau1,tau2,a,b) a*exp(-x/tau1)+b*exp(-x/tau2);
getTimeIndShift = @(x_step,timeShift) floor(timeShift/x_step);
applyDriveBasedTerm = @(reach_pdf,reach_suppress,timeIndShift) subplus(reach_pdf(1:end-timeIndShift)-reach_suppress(timeIndShift+1:end));

% Transform to change in reaction time, accounting for any decrease in rate



x=0:0.001:15; y=0.1*lognpdf(x,0,0.5); figure(); plot(x,y);
tau=0.5;
tau2=5;
a=0.7;
b=0.3;
x_exp=0:0.001:15; y_exp=a*exp(-x_exp/tau)+b*exp(-x_exp/tau2); figure(); plot(x_exp,y_exp);
figure(); plot(x,y);
addBaseline=0.005;
% hold on; i=7000; func1=subplus(y(1:end-i)-y_exp(i+1:end))+ones(size(y_exp(i+1:end))).*addBaseline; plot(x(1:end-i),func1,'Color','c');
% hold on; i=9000; func2=subplus(y(1:end-i)-y_exp(i+1:end))+ones(size(y_exp(i+1:end))).*addBaseline; plot(x(1:end-i),func2,'Color','m');
hold on; i=7000; func1=subplus(y(1:end-i)-y_exp(i+1:end)); plot(x(1:end-i),func1,'Color','c');
hold on; i=9000; func2=subplus(y(1:end-i)-y_exp(i+1:end)); plot(x(1:end-i),func2,'Color','m');

sample_n=10000;
sampleFromFunc1=randsample(x(1:length(func1)),sample_n,true,func1);
sampleFromFunc2=randsample(x(1:length(func2)),sample_n,true,func2);
plotHist(sampleFromFunc1,sampleFromFunc2,50,'Sampling','Reaches','c','m');

nTrials=1000;
rts1=nan(1,nTrials);
rts2=nan(1,nTrials);
for i=1:nTrials
    sample_n=50;
    % Get expected number of reaches per trial from pdf?
    sampleFromFunc1=randsample(x(1:length(func1)),sample_n,true,func1);
    sampleFromFunc2=randsample(x(1:length(func2)),sample_n,true,func2);
    % Get "reaction time" of first reach
    rts1(i)=min(sampleFromFunc1);
    rts2(i)=min(sampleFromFunc2);
end
[histo_nbins,counts1,counts2]=plotHist(rts1,rts2,50,'Sample and get RTs','RT','c','m');

figure();
plot(x(1:length(func1)),accumulateDistribution(func1),'Color','c');
hold on;
plot(x(1:length(func2)),accumulateDistribution(func2),'Color','m');

return