function plotReachRateChangeMinusSatiety()

% Changes from paired trials (trial to trial)
whichSess=ones(size(reachRates.trial1_alltrials_cued)).*repmat([1:size(reachRates.trial1_alltrials_cued,1)]',1,size(reachRates.trial1_alltrials_cued,2));
whichSess=whichSess(1:end);
fracsThroughSess=reachRates.fracsThroughSess(1:end);
trial1_reachRates_uncued=reachRates.trial1_alltrials_uncued(1:end);
trial1_reachRates_cued=reachRates.trial1_alltrials_cued(1:end);
trialn_reachRates_uncued=reachRates.alltrials_uncued(1:end);
trialn_reachRates_cued=reachRates.alltrials_cued(1:end);
m_binned_inds1to100=1:(length(m_binned)/102):length(m_binned);
orth_comps_change=nan(1,length(fracsThroughSess));
par_comps_change=nan(1,length(fracsThroughSess));
orth_comps_trial1=nan(1,length(fracsThroughSess));
par_comps_trial1=nan(1,length(fracsThroughSess));
orth_comps_trialn=nan(1,length(fracsThroughSess));
par_comps_trialn=nan(1,length(fracsThroughSess));
all_angles=nan(1,length(fracsThroughSess));
all_mags=nan(1,length(fracsThroughSess));
all_angles_minusSatiety=nan(1,length(fracsThroughSess));
all_mags_minusSatiety=nan(1,length(fracsThroughSess));
figure();
for i=1:length(fracsThroughSess) % for each trial to trial change from each session from each mouse
    if isnan(trial1_reachRates_uncued(i)) || isnan(trialn_reachRates_uncued(i))
        continue
    end
    % Get closest m based on the fraction through this session
    % Note that have already averaged across mice and sessions to calculate
    % these m's
    temp=round(fracsThroughSess(i)*100);
    if temp>length(m_binned_inds1to100)
        temp=length(m_binned_inds1to100);
    elseif temp<1
        temp=1;
    end
    temp=ceil(m_binned_inds1to100(temp));
    if temp>length(m_binned)
        temp=length(m_binned);
    end
    curr_m=m_binned(temp);
    
    vec3=[1 curr_m];
    V_par=vec3./norm(vec3);
    % Get m orthogonal to this
    vec2=[-V_par(2)/V_par(1) 1];
    V_orth=vec2./norm(vec2);
    
    
    % Current trial-to-trial change
    vec1=[trialn_reachRates_uncued(i)-trial1_reachRates_uncued(i) trialn_reachRates_cued(i)-trial1_reachRates_cued(i)];
    vec1_trial1=[trial1_reachRates_uncued(i) trial1_reachRates_cued(i)];
    vec1_trialn=[trialn_reachRates_uncued(i) trialn_reachRates_cued(i)];
    
    % Get cued change minus expected cued change for this uncued change
    currusem=(vec1_trial1(2)/vec1_trial1(1));
    curr_m=currusem;
    actual_y_change=vec1(2);
%     expected_y_change=(vec1_trial1(2)/vec1_trial1(1))*vec1(1);
    expected_y_change=curr_m*vec1(1);
    actual_x_change=vec1(1);
    
    if isnan(curr_m) | isinf(curr_m) | curr_m==0
        expected_x_change=0;
    else
        expected_x_change=actual_y_change/(curr_m);
    end
    
    line([0 vec1(1)+rand(1)/50],[0 vec1(2)],'Color','k');
    hold all;
    
    all_angles(i)=angle(complex(vec1(1),vec1(2)));
    all_mags(i)=abs(complex(vec1(1),vec1(2)));
%     orth_comps_change(i)=vec1(1)*V_orth(1)+vec1(2)*V_orth(2);
%     par_comps_change(i)=vec1(1)*V_par(1)+vec1(2)*V_par(2);
    orth_comps_change(i)=actual_y_change-expected_y_change;
    par_comps_change(i)=actual_x_change-expected_x_change;

    veep=[-1 1];
    Veepers_orth=veep./norm(veep);
    Veepers_par=[-Veepers_orth(2)/Veepers_orth(1) 1];
    Veepers_par=Veepers_par./norm(Veepers_par);
    tempers1=par_comps_change(i)*Veepers_orth(1)+orth_comps_change(i)*Veepers_orth(2);
    tempers2=par_comps_change(i)*Veepers_par(1)+orth_comps_change(i)*Veepers_par(2);
    orth_comps_change(i)=tempers1;
    par_comps_change(i)=tempers2;
    
    orth_comps_trial1(i)=vec1_trial1(1)*V_orth(1)+vec1_trial1(2)*V_orth(2);
    par_comps_trial1(i)=vec1_trial1(1)*V_par(1)+vec1_trial1(2)*V_par(2);
    orth_comps_trialn(i)=vec1_trialn(1)*V_orth(1)+vec1_trialn(2)*V_orth(2);
    par_comps_trialn(i)=vec1_trialn(1)*V_par(1)+vec1_trialn(2)*V_par(2);
    all_angles_minusSatiety(i)=angle(complex(par_comps_change(i),orth_comps_change(i)));
    all_mags_minusSatiety(i)=abs(complex(par_comps_change(i),orth_comps_change(i)));
end

figure();
for i=1:length(par_comps_change)
if isnan(fracsThroughSess(i)) | abs(orth_comps_change(i))>100
continue
end
temp=ceil(fracsThroughSess(i)*100);
if temp>size(cmap,1)
temp=size(cmap,1);
end
scatter(par_comps_change(i)+(-1 + (1+1) .* rand(1,1))*0.3,orth_comps_change(i)+(-1 + (1+1) .* rand(1,1))*0.3,[],'k');
hold on;
end
line([0 0],[-5 5],'Color','k');
line([-5 5],[0 0],'Color','k');
throwout=isinf(par_comps_change) | isinf(orth_comps_change);
par_comps_change=par_comps_change(throwout==false);
orth_comps_change=orth_comps_change(throwout==false);
outies=(par_comps_change>nanmean(par_comps_change)+3*mad(par_comps_change) | par_comps_change<nanmean(par_comps_change)-3*mad(par_comps_change)) | ...
       (orth_comps_change>nanmean(orth_comps_change)+3*mad(orth_comps_change) | orth_comps_change<nanmean(orth_comps_change)-3*mad(orth_comps_change));
par_comps_change=par_comps_change(outies==false);
orth_comps_change=orth_comps_change(outies==false);
line([0 nanmean(par_comps_change)],[0 nanmean(orth_comps_change)],'Color','m');


 figure();
for i=1:length(par_comps_change)
if isnan(fracsThroughSess(i))
continue
end
temp=ceil(fracsThroughSess(i)*100);
if temp>size(cmap,1)
temp=size(cmap,1);
end
vec1=[trialn_reachRates_uncued(i)-trial1_reachRates_uncued(i) trialn_reachRates_cued(i)-trial1_reachRates_cued(i)];
scatter(vec1(1)+(-1 + (1+1) .* rand(1,1))*1,vec1(2)+(-1 + (1+1) .* rand(1,1))*1,[],'k');
hold on;
end
line([0 0],[-5 5],'Color','k');
line([-5 5],[0 0],'Color','k');

angs=all_angles(~isnan(all_angles))';
mags=all_mags(~isnan(all_angles))';
angs(mags>nanmean(mags)+3*mad(mags) | mags<nanmean(mags)-3*mad(mags))=nan;
mags=mags(~isnan(angs));
angs=angs(~isnan(angs));
% Bin angles
% edges=0:10:360;
% edges=deg2rad(edges);
% countsperbin=nan(1,length(edges)-1);
% for i=1:length(edges)-1
%     countsperbin(i)=nansum(angs>=edges(i) & angs<edges(i+1))*(nansum(mags(angs>=edges(i) & angs<edges(i+1)))/0.33);
% end
% [mu,ul,ll]=circ_mean(((edges(1:end-1)+edges(2:end))/2)',countsperbin');
% edges=0:10:360;
% n=histcounts(angs,deg2rad(edges));
[mu,ul,ll]=circ_mean(angs,(mags/0.33)*10,1);
% [mu,ul,ll]=circ_mean(deg2rad((edges(1:end-1)+edges(2:end))/2)',n');
disp(['mean angle: ' num2str(rad2deg(mu))]);
disp(['mean angle lower confidence bound: ' num2str(rad2deg(ul))]);
disp(['mean angle upper confidence bound: ' num2str(rad2deg(ll))]);
disp(['mean magnitude: ' num2str(nanmean(mags))]);
disp(['se magnitude: ' num2str(nanstd(mags,[],1)./sqrt(length(mags)))]);

angs=all_angles_minusSatiety(~isnan(all_angles_minusSatiety))';
mags=all_mags_minusSatiety(~isnan(all_angles_minusSatiety))';
angs(mags>nanmean(mags)+3*mad(mags) | mags<nanmean(mags)-3*mad(mags))=nan;
dontuse=isnan(angs) | isinf(mags);
mags=mags(~dontuse);
angs=angs(~dontuse);
[mu,ul,ll]=circ_mean(angs,(mags/0.33)*10,1); % mean mag is about 0.3
% edges=0:10:360;
% edges=deg2rad(edges);
% countsperbin=nan(1,length(edges)-1);
% for i=1:length(edges)-1
%     countsperbin(i)=nansum(angs>=edges(i) & angs<edges(i+1))*(nansum(mags(angs>=edges(i) & angs<edges(i+1)))/0.33);
% end
% [mu,ul,ll]=circ_mean(((edges(1:end-1)+edges(2:end))/2)',countsperbin');
% n=histcounts(angs,deg2rad(edges));
% [mu,ul,ll]=circ_mean(deg2rad((edges(1:end-1)+edges(2:end))/2)',n');
disp(['mean angle MINUS SATIETY: ' num2str(rad2deg(mu))]);
disp(['mean angle lower confidence bound MINUS SATIETY: ' num2str(rad2deg(ul))]);
disp(['mean angle upper confidence bound MINUS SATIETY: ' num2str(rad2deg(ll))]);
disp(['mean magnitude: ' num2str(nanmean(mags))]);

% Average across each session, then plot the average shift session by
% session
u=unique(whichSess);
orth_trial1=nan(1,length(u));
par_trial1=nan(1,length(u));
orth_trialn=nan(1,length(u));
par_trialn=nan(1,length(u));
for i=1:length(u)
    currSess=u(i);
    orth_trial1(i)=nanmean(orth_comps_trial1(whichSess==currSess));
    par_trial1(i)=nanmean(par_comps_trial1(whichSess==currSess));
    orth_trialn(i)=nanmean(orth_comps_trialn(whichSess==currSess));
    par_trialn(i)=nanmean(par_comps_trialn(whichSess==currSess));
end
figure(); 
for i=1:length(orth_trial1)
    if any(isnan([orth_trial1(i) par_trial1(i) orth_trialn(i) par_trialn(i)]))
    else
%         quiver(par_trial1(i),orth_trial1(i),par_trialn(i)-par_trial1(i),orth_trialn(i)-orth_trial1(i),'Color','k');
        quiver(0,0,par_trialn(i)-par_trial1(i),orth_trialn(i)-orth_trial1(i),'Color','k');
        hold on;
    end
end
changesAcrossSess_x=par_trialn-par_trial1;
changesAcrossSess_y=orth_trialn-orth_trial1;
outlie=par_trial1>nanmean(par_trial1)+3*mad(par_trial1) | par_trial1<nanmean(par_trial1)-3*mad(par_trial1) | ...
       par_trialn>nanmean(par_trialn)+3*mad(par_trialn) | par_trialn<nanmean(par_trialn)-3*mad(par_trialn) | ...
       orth_trial1>nanmean(orth_trial1)+3*mad(orth_trial1) | orth_trial1<nanmean(orth_trial1)-3*mad(orth_trial1) | ...
       orth_trialn>nanmean(orth_trialn)+3*mad(orth_trialn) | orth_trialn<nanmean(orth_trialn)-3*mad(orth_trialn);
changesAcrossSess_x(outlie==1)=nan;
changesAcrossSess_y(outlie==1)=nan;
figure(); 
for i=1:length(changesAcrossSess_x)
    if any(isnan([changesAcrossSess_x(i) changesAcrossSess_y(i)]))
    else
        quiver(0,0,changesAcrossSess_x(i),changesAcrossSess_y(i),'Color','k');
        hold on;
    end
end
figure();
if isnan(nanmean(changesAcrossSess_x)) | isnan(nanmean(changesAcrossSess_y))
else
    quiver(0,0,nanmean(changesAcrossSess_x),nanmean(changesAcrossSess_y));
    disp(['Magnitude of vector: ' num2str(sqrt(nanmean(changesAcrossSess_x)^2+nanmean(changesAcrossSess_y)^2))]);
end

return

end