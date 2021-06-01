function [new_x,new_n]=plotPairedChangeMinusSatiety(reachRates)

trial1_reachRates_uncued=reachRates.trial1_alltrials_uncued(1:end);
trial1_reachRates_cued=reachRates.trial1_alltrials_cued(1:end);
trialn_reachRates_uncued=reachRates.alltrials_uncued(1:end);
trialn_reachRates_cued=reachRates.alltrials_cued(1:end);
fracsThroughSess=reachRates.fracsThroughSess(1:end);

vec1=[trialn_reachRates_uncued-trial1_reachRates_uncued; trialn_reachRates_cued-trial1_reachRates_cued];
vec1_trial1=[trial1_reachRates_uncued; trial1_reachRates_cued];
vec1_trialn=[trialn_reachRates_uncued; trialn_reachRates_cued];
curr_m=(vec1_trial1(2,:)./vec1_trial1(1,:));
actual_y_change=vec1(2,:);
expected_y_change=curr_m.*vec1(1,:);
actual_x_change=vec1(1,:);
expected_x_change=actual_y_change./curr_m;
expected_x_change(isnan(curr_m))=nan;
expected_x_change(isinf(curr_m))=0;
expected_x_change(curr_m==0)=nan;
y_comps_change=actual_y_change-expected_y_change;
x_comps_change=actual_x_change-expected_x_change;
% Get projections
projected_to_orth=nan(1,length(actual_y_change));
projected_to_par=nan(1,length(actual_y_change));
for i=1:length(actual_y_change)
    V_par=[1 1];
    V_par=V_par./norm(V_par);
    % Get orthogonal to this
    V_orth=[-V_par(2)/V_par(1) 1];
    V_orth=V_orth./norm(V_orth);
    % Project
    projected_to_orth(i)=x_comps_change(i)*V_orth(1)+y_comps_change(i)*V_orth(2);
    projected_to_par(i)=x_comps_change(i)*V_par(1)+y_comps_change(i)*V_par(2);
end
    
% Deviation from expected
xlab='UNCUED actual rate change minus expected rate change (1/sec)';
ylab='CUED actual rate change minus expected rate change (1/sec)';
plotVersus(fracsThroughSess,x_comps_change,y_comps_change,xlab,ylab,0.01);
daspect([1 1 1]);

% Expected
% xlab='UNCUED expected rate change (1/sec)';
% ylab='CUED expected rate change (1/sec)';
% plotVersus(fracsThroughSess,expected_x_change,expected_y_change,xlab,ylab,0.1);

% X real vs. Y deviation from expected
% xlab='UNCUED actual rate change (1/sec)';
% ylab='CUED actual rate change minus expected rate change (1/sec)';
% plotVersus(fracsThroughSess,actual_x_change,y_comps_change,xlab,ylab,0.1);

% Y real vs. X deviation from expected
% ylab='CUED actual rate change (1/sec)';
% xlab='UNCUED actual rate change minus expected rate change (1/sec)';
% plotVersus(fracsThroughSess,x_comps_change,actual_y_change,xlab,ylab,0.1);

% Projections
xlab='Deviation projected onto 45 deg line increase cue, increase uncue (1/sec)';
ylab='Deviation projected onto 45 deg line increase cue, decrease uncue (1/sec)';
plotVersus(fracsThroughSess,projected_to_par,projected_to_orth,xlab,ylab,0.01);
daspect([1 1 1]);
[n,edges]=histcounts(projected_to_orth,-10-0.125:0.25:10+0.125);
[new_n,new_x]=cityscape_hist(n,edges);
figure();
plot(new_x,new_n,'Color','k');
xlabel('Deviation projected onto 45 deg line increase cue, increase uncue (1/sec)');
ylabel('Counts');

end

function plotVersus(fracsThroughSess,x_comps_change,y_comps_change,xlab,ylab,randscale)

figure();
cmap=colormap('cool');
temp=round(fracsThroughSess.*size(cmap,1));
temp(temp>size(cmap,1))=size(cmap,1);
temp(temp<1)=1;
r1=(-1 + (1+1) .* rand(size(fracsThroughSess)))*randscale;
r2=(-1 + (1+1) .* rand(size(fracsThroughSess)))*randscale;
scatter(x_comps_change+r1,y_comps_change+r2,30,cmap(temp,:),'LineWidth',1);
keeps=~isnan(x_comps_change) & ~isnan(y_comps_change) & ~isinf(x_comps_change) & ~isinf(y_comps_change);
x_comps_change=x_comps_change(keeps);
y_comps_change=y_comps_change(keeps);
howmany=length(x_comps_change);
line([nanmean(x_comps_change)-nanstd(x_comps_change,[],2)./sqrt(howmany) nanmean(x_comps_change)+nanstd(x_comps_change,[],2)./sqrt(howmany)],...
     [nanmean(y_comps_change) nanmean(y_comps_change)],'Color','k','LineWidth',2);
line([nanmean(x_comps_change) nanmean(x_comps_change)],...
     [nanmean(y_comps_change)-nanstd(y_comps_change,[],2)./sqrt(howmany) nanmean(y_comps_change)+nanstd(y_comps_change,[],2)./sqrt(howmany)],'Color','k','LineWidth',2);
xlim([min(x_comps_change,[],2,'omitnan') max(x_comps_change,[],2,'omitnan')]);
ylim([min(y_comps_change,[],2,'omitnan') max(y_comps_change,[],2,'omitnan')]);
line([0 0],[min(y_comps_change,[],2,'omitnan') max(y_comps_change,[],2,'omitnan')],'Color','k');
line([min(x_comps_change,[],2,'omitnan') max(x_comps_change,[],2,'omitnan')],[0 0],'Color','k');
xlabel(xlab);
ylabel(ylab);

end