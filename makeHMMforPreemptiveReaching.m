function [rt_pdf_outs,all_pdfs]=makeHMMforPreemptiveReaching(preemptive_process,cued_process,rt_bins)

% % If in cued reaching mode, probability of switching to preemptive reaching mode
% p_cp=0.05;
% % If in preemptive reaching mode, probability of switching to cued reaching mode
% p_pc=0.02;
% % If in cued reaching mode, probability of staying in cued reaching mode
% p_cc=0.95;
% % If in preemptive reaching mode, probability of staying in preemptive
% % reaching mode
% p_pp=0.98;

% % If in cued reaching mode, probability of switching to preemptive reaching mode
% p_cp=0.05;
% % If in preemptive reaching mode, probability of switching to cued reaching mode
% p_pc=0.01;
% % If in cued reaching mode, probability of staying in cued reaching mode
% p_cc=0.95;
% % If in preemptive reaching mode, probability of staying in preemptive
% % reaching mode
% p_pp=0.99;

% % If in cued reaching mode, probability of switching to preemptive reaching mode
% p_cp=0.3;
% % If in preemptive reaching mode, probability of switching to cued reaching mode
% p_pc=0.005;
% % If in cued reaching mode, probability of staying in cued reaching mode
% p_cc=0.7;
% % If in preemptive reaching mode, probability of staying in preemptive
% % reaching mode
% p_pp=0.995;

% If in cued reaching mode, probability of switching to preemptive reaching mode
p_cp=0.30;
% If in preemptive reaching mode, probability of switching to cued reaching mode
p_pc=0.005;
% If in cued reaching mode, probability of staying in cued reaching mode
p_cc=0.70;
% If in preemptive reaching mode, probability of staying in preemptive
% reaching mode
p_pp=0.995;


% Transition probabilities from each state must sum to 1
temp1=p_pp/(p_pp+p_pc);
temp2=p_pc/(p_pp+p_pc);
p_pp=temp1;
p_pc=temp2;

temp1=p_cc/(p_cc+p_cp);
temp2=p_cp/(p_cc+p_cp);
p_cc=temp1;
p_cp=temp2;

disp('sum of transition probabilities from preemptive reaching');
disp(num2str(p_pp+p_pc));
disp('sum of transition probabilities from cued reaching');
disp(num2str(p_cc+p_cp));

bin_centers = @(bins) nanmean([bins(1:end-1); bins(2:end)],1);

% PDF of reaction times for preemptive reaching
% preempt_pdf=gampdf(bin_centers(rt_bins),preemptive_process.shape,preemptive_process.rate);
% preempt_pdf=gampdf(bin_centers(rt_bins),1.15,0.25);
% preempt_pdf=gampdf(bin_centers(rt_bins),1.5,0.25);
% preempt_pdf=gampdf(bin_centers(rt_bins),1.7,0.25);
if isfield(preemptive_process,'preempt_pdf')
    preempt_pdf=preemptive_process.preempt_pdf;
else
    preempt_pdf=gampdf(bin_centers(rt_bins),1.7,0.25);
    preempt_pdf=[preempt_pdf(4:end) zeros(1,3)];
    preempt_pdf=preempt_pdf./nansum(preempt_pdf);
    preempt_pdf(preempt_pdf<=0)=0.00000001;
end

% PDF of reaction times for cued reaching
cued_pdf=cued_process.rt_pdf;
cued_pdf=cued_pdf./nansum(cued_pdf);
cued_pdf(cued_pdf<=0)=0.00000001;

% Simulate
nSamples=30000;
% Choose initial state randomly
state=rand<0.5; % if state is 0, cued reaching; if state is 1, preemptive reaching
all_rts=nan(1,nSamples);
all_pdfs=zeros(length(bin_centers(rt_bins)),length(bin_centers(rt_bins)));
for i=1:nSamples
    if state==0
        % cued reaching
        % draw curr_rt from cued rt pdf
        all_rts(i)=randsample(bin_centers(rt_bins),1,true,cued_pdf);
        all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)=all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)-cued_pdf;
        % test for transition
        % also transition in proportion to how certain I am that in this
        % state
        if rand<p_cc.*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))/(p_cc.*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))+p_cp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first')))
            % stay in this state
        else
            state=1;
        end
    elseif state==1
        % preemptive reaching
        % draw curr_rt from preemptive rt pdf
        all_rts(i)=randsample(bin_centers(rt_bins),1,true,preempt_pdf);
        all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)=all_pdfs(find(bin_centers(rt_bins)>=all_rts(i),1,'first'),:)+preempt_pdf;
        % test for transition
        % also transition in proportion to how certain I am that in this
        % state
        if rand<p_pp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))/(p_pp.*preempt_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first'))+p_pc*cued_pdf(find(bin_centers(rt_bins)>=all_rts(i),1,'first')))
            % stay in this state
        else
            state=0;
        end
    end
end

% Make plot of output rt pdfs as a function of previous trial rt
rt_pdf_outs=nan(length(bin_centers(rt_bins)),length(bin_centers(rt_bins)));
for i=1:length(bin_centers(rt_bins))
    curr_rt_bin=[rt_bins(i) rt_bins(i+1)];
    f=find(all_rts>=curr_rt_bin(1) & all_rts<curr_rt_bin(2));
    f=f+1;
    f=f(f<=length(all_rts));
    next_rts=all_rts(f);
    [n,x]=histcounts(next_rts,rt_bins);
    rt_pdf_outs(i,:)=n;
end
figure()
imagesc(bin_centers(rt_bins),bin_centers(rt_bins),rt_pdf_outs);

figure();
imagesc(bin_centers(rt_bins),bin_centers(rt_bins),all_pdfs);

% figure();
% plot(bin_centers(rt_bins),rt_pdf_outs(1,:)./nansum(rt_pdf_outs(1,:)),'Color','k');
% hold on; 
% plot(bin_centers(rt_bins),preempt_pdf./nansum(preempt_pdf),'Color','r');

    
    
    