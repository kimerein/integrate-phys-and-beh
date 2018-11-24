function RTdataAcrossDays(alltbt,out,metadata)

unique_sess=unique(metadata.sessid);

f=fieldnames(alltbt);
fout=fieldnames(out);

allout.allreactiontimes.ledFalse=cell(1,length(unique_sess));
allout.allreactiontimes.ledTrue=cell(1,length(unique_sess));
all_mu=nan(1,length(unique_sess));
all_sigma=nan(1,length(unique_sess));
all_mu_red=nan(1,length(unique_sess));
all_sigma_red=nan(1,length(unique_sess));
allmouseids=nan(1,length(unique_sess));
for i=1:length(unique_sess)
    currsessid=unique_sess(i);
    for j=1:length(f)
        temp=alltbt.(f{j});
        currtbt.(f{j})=temp(metadata.sessid==currsessid,:);
    end
    for j=1:length(fout)
        temp=out.(fout{j});
        currout.(fout{j})=temp(metadata.sessid==currsessid);
    end
    output=getRTanalysis(currtbt,currout);
    fracThrough=metadata.fractionThroughSess(metadata.sessid==currsessid);
%     allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse;
    allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse(fracThrough>=0.1 & fracThrough<=0.5);
    hax=axes(); 
    parmhat=fitAndPlot(hax,allout.allreactiontimes.ledFalse{i},'lognfit',0-0.03:0.06:15-0.03,'k');
    all_mu(i)=parmhat(1);
    all_sigma(i)=parmhat(2);
    nth_session=nanmean(metadata.nth_session(metadata.sessid==currsessid));
    whichmouse=nanmean(metadata.mouseid(metadata.sessid==currsessid));
    allmouseids(i)=whichmouse;
%     disp(nth_session);
    close all;
    allout.allreactiontimes.ledTrue{i}=output.reactionTimes.ledFalse(fracThrough>=0.5 & fracThrough<=0.9);
    hax=axes();
    parmhat=fitAndPlot(hax,allout.allreactiontimes.ledTrue{i},'lognfit',0-0.03:0.06:15-0.03,'k');
    all_mu_red(i)=parmhat(1);
    all_sigma_red(i)=parmhat(2);
    % allout.allreactiontimes.ledFalse{i}=output.reactionTimes.ledFalse;
    % allout.allreactiontimes.ledTrue{i}=output.reactionTimes.ledTrue;
end

allRT_ledFalse=allout.allreactiontimes.ledFalse;
allRT_ledTrue=allout.allreactiontimes.ledTrue;

for i=1:length(unique_sess)
    all_med_ledFalse(i)=nanmedian(allRT_ledFalse{i});
    all_med_ledTrue(i)=nanmedian(allRT_ledTrue{i});
end

useMetric=all_mu_red;

% useMetric=exp(all_mu_red-all_sigma_red.^2);
% all_mu=exp(all_mu-all_sigma.^2);
% all_mu_red=exp(all_mu_red-all_sigma_red.^2);

% useMetric=log2(all_sigma_red.*exp(all_mu_red+0.5).*sqrt(2*pi));
% all_mu=log2(all_sigma.*exp(all_mu+0.5).*sqrt(2*pi));
% all_mu_red=log2(all_sigma_red.*exp(all_mu_red+0.5).*sqrt(2*pi));

% useMetric=all_mu_red+0.5*all_sigma_red.^2;
% all_mu=all_mu+0.5*all_sigma.^2;
% all_mu_red=all_mu_red+0.5*all_sigma_red.^2;

me=nanmean(useMetric);
sd=nanstd(useMetric,[],2);
useMetric(useMetric>me+10*sd | useMetric<me-10*sd)=nan;

figure(); 
all_nth_sessions=nan(1,length(unique_sess));
for i=1:length(unique_sess)
    currsessid=unique_sess(i);
    nth_session=nanmean(metadata.nth_session(metadata.sessid==currsessid));
    scatter(nth_session,all_mu(i),[],'r');
%     scatter(nth_session,useMetric(i),[],'r');
    all_nth_sessions(i)=nth_session;
    hold on;
    scatter(nth_session+0.5,all_mu_red(i),[],'c');
    line([nth_session nth_session+0.5],[all_mu(i) all_mu_red(i)],'Color','k');
end
disp('p-value comparing beginning vs end of session');
disp('mean beginning');
disp(nanmean(all_mu));
disp('mean end');
disp(nanmean(all_mu_red));
p=signrank(all_mu,all_mu_red);
disp(p);

hold all;
unique_mouse=unique(metadata.mouseid);
for i=1:length(unique_mouse)
    plot(all_nth_sessions(allmouseids==unique_mouse(i))+0.5,useMetric(allmouseids==unique_mouse(i)));       
end

useNthSess=zeros(size(all_nth_sessions));
useNthSess=ones(size(all_nth_sessions));
[r,p]=corrcoef(all_nth_sessions(useNthSess & ~isnan(all_nth_sessions) & ~isnan(useMetric)),useMetric(useNthSess & ~isnan(all_nth_sessions) & ~isnan(useMetric)));
disp('corr coef');
disp(r);
disp('p-value');
disp(p);

end

function output=getRTanalysis(alltbt,out)

nbins=[];

% Get reaction times for all trials where mouse reached after cue onset
[reactionTimes,alltbt]=plotOnlyFirstReach(alltbt,1,'reachStarts_noPawOnWheel','cueZone_onVoff',out,'led',0);
% output.reactionTimes.ledFalse=reactionTimes(out.led_1back==0 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0);
% output.reactionTimes.ledTrue=reactionTimes(out.led_1back==1 & out.cued_reach_1back==1 & out.touched_pellet_1back==1 & out.chewing_at_trial_start==0);
temp=reactionTimes;
temp(~(out.chewing_at_trial_start==0 & out.led==0))=nan;
output.reactionTimes.ledFalse=temp;
temp=reactionTimes;
temp(~(out.chewing_at_trial_start==0 & out.led==0))=nan;
output.reactionTimes.ledTrue=temp;
% temp(out.led_1back==1)=nan;
% output.reactionTimes.ledFalse=temp;
% temp=reactionTimes;
% temp(out.led_1back==0)=nan;
% output.reactionTimes.ledTrue=temp;

end