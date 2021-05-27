function plotTrialToTrialChangeMinusSatiety(reachRates,m_binned,cued_binned,uncued_binned)

% THIS FUNCTION ASSUMES THAT M_BINNED ARE EQUALLY SPACED ACROSS EACH
% SESSION!

m_binned=m_binned(~isnan(m_binned));

% Plot the m's you are removing
cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(m_binned));
f=figure(); 
for i=1:length(m_binned)
    if i==length(m_binned)
    else
        length_line_segment=(cued_binned(i+1)-cued_binned(i))/m_binned(i);
    end
    line([uncued_binned(i)-length_line_segment/2 uncued_binned(i)+length_line_segment/2],[cued_binned(i)-(length_line_segment*m_binned(i))/2 cued_binned(i)+(length_line_segment*m_binned(i))/2],'Color',cmap(k,:),'LineWidth',1);
    hold on;
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

% kstep=size(cmap,1)/101;
% fracThroughColormap=1:kstep:size(cmap,1);
% temp=ceil(fracsThroughSess(i)*100);
%     temp=round(fracThroughColormap(temp));
%     if temp<1
%         temp=1;
%     elseif temp>size(cmap,1)
%         temp=size(cmap,1);
%     end
%     indToCmap=temp;

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
    all_angles(i)=angle(complex(vec1(1),vec1(2)));
    all_mags(i)=abs(complex(vec1(1),vec1(2)));
    hold on;
    orth_comps_change(i)=vec1(1)*V_orth(1)+vec1(2)*V_orth(2);
    par_comps_change(i)=vec1(1)*V_par(1)+vec1(2)*V_par(2);
    orth_comps_trial1(i)=vec1_trial1(1)*V_orth(1)+vec1_trial1(2)*V_orth(2);
    par_comps_trial1(i)=vec1_trial1(1)*V_par(1)+vec1_trial1(2)*V_par(2);
    orth_comps_trialn(i)=vec1_trialn(1)*V_orth(1)+vec1_trialn(2)*V_orth(2);
    par_comps_trialn(i)=vec1_trialn(1)*V_par(1)+vec1_trialn(2)*V_par(2);
    all_angles_minusSatiety(i)=angle(complex(par_comps_change(i),orth_comps_change(i)));
    all_mags_minusSatiety(i)=abs(complex(par_comps_change(i),orth_comps_change(i)));
end
angs=all_angles(~isnan(all_angles))';
mags=all_mags(~isnan(all_angles))';
angs(mags>nanmean(mags)+3*mad(mags) | mags<nanmean(mags)-3*mad(mags))=nan;
mags=mags(~isnan(angs));
angs=angs(~isnan(angs));
[mu,ul,ll]=circ_mean(angs,mags*100);
disp(['mean angle: ' num2str(rad2deg(mu))]);
disp(['mean angle lower confidence bound: ' num2str(rad2deg(ul))]);
disp(['mean angle upper confidence bound: ' num2str(rad2deg(ll))]);

angs=all_angles_minusSatiety(~isnan(all_angles_minusSatiety))';
mags=all_mags_minusSatiety(~isnan(all_angles_minusSatiety))';
angs(mags>nanmean(mags)+3*mad(mags) | mags<nanmean(mags)-3*mad(mags))=nan;
mags=mags(~isnan(angs));
angs=angs(~isnan(angs));
[mu,ul,ll]=circ_mean(angs,mags*100,1);
disp(['mean angle MINUS SATIETY: ' num2str(rad2deg(mu))]);
disp(['mean angle lower confidence bound MINUS SATIETY: ' num2str(rad2deg(ul))]);
disp(['mean angle upper confidence bound MINUS SATIETY: ' num2str(rad2deg(ll))]);

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
quiver(0,0,nanmean(changesAcrossSess_x),nanmean(changesAcrossSess_y));

return
% Take components orthogonal to m's capturing satiety
orth_comps=nan(1,length(reachRates.uncued));
par_comps=nan(1,length(reachRates.uncued));
orth_se=nan(1,length(reachRates.uncued));
par_se=nan(1,length(reachRates.uncued));
for i=1:length(reachRates.uncued)
    % Get closest m
    curr_m=m_binned(refsToMbinned(i));
    vec3=[1 curr_m];
    V_par=vec3./norm(vec3);
    % Get m orthogonal to this
    vec2=[-V_par(2)/V_par(1) 1];
    V_orth=vec2./norm(vec2);
    % Project current point onto this orthogonal_m
    vec1=[reachRates.uncued(i) reachRates.cued(i)];
    orth_comps(i)=vec1(1)*V_orth(1)+vec1(2)*V_orth(2);
    par_comps(i)=vec1(1)*V_par(1)+vec1(2)*V_par(2);
    % Project se
    temp_uncued=reachRates.alltrials_uncued(:,i);
    uncued_se=nanstd(temp_uncued)./sqrt(nansum(~isnan(temp_uncued)));
    temp_cued=reachRates.alltrials_cued(:,i);
    cued_se=nanstd(temp_cued)./sqrt(nansum(~isnan(temp_cued)));
    se_vec=[uncued_se cued_se];
    orth_se(i)=se_vec(1)*V_orth(1)+se_vec(2)*V_orth(2);
    par_se(i)=se_vec(1)*V_par(1)+se_vec(2)*V_par(2);
end
% Plot how components evolve over trials
figure();
k=1;
kstep=ceil(size(cmap,1)/length(reachRates.uncued));
for i=1:length(orth_comps)
    line([par_comps(i)-par_se(i) par_comps(i)+par_se(i)],[orth_comps(i) orth_comps(i)],'Color',cmap(k,:));
    hold on; 
    line([par_comps(i) par_comps(i)],[orth_comps(i)-orth_se(i) orth_comps(i)+orth_se(i)],'Color',cmap(k,:));
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
downSamp_par_comps=downSampAv(par_comps,ceil(length(reachRates.m)/nbins));
downSamp_orth_comps=downSampAv(orth_comps,ceil(length(reachRates.m)/nbins));
k=1;
kstep=ceil(size(cmap,1)/length(downSamp_orth_comps));
for i=1:length(downSamp_par_comps)
    scatter(downSamp_par_comps(i),downSamp_orth_comps(i),[],cmap(k,:),'filled');
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
downSamp_par_comps(downSamp_par_comps>nanmedian(downSamp_par_comps)+3*mad(downSamp_par_comps) | downSamp_par_comps<nanmedian(downSamp_par_comps)-3*mad(downSamp_par_comps))=nan;
downSamp_orth_comps(downSamp_orth_comps>nanmedian(downSamp_orth_comps)+3*mad(downSamp_orth_comps) | downSamp_orth_comps<nanmedian(downSamp_orth_comps)-3*mad(downSamp_orth_comps))=nan;
firstHalf=[nanmean(downSamp_par_comps(1:floor(length(downSamp_par_comps)/2))) nanmean(downSamp_orth_comps(1:floor(length(downSamp_orth_comps)/2)))];
secondHalf=[nanmean(downSamp_par_comps(floor(length(downSamp_par_comps)/2)+1:end)) nanmean(downSamp_orth_comps(floor(length(downSamp_orth_comps)/2)+1:end))];
if any(isnan(firstHalf)) | any(isnan(secondHalf)) | any(isinf(firstHalf)) | any(isinf(secondHalf))
else
    quiver(firstHalf(1),firstHalf(2),secondHalf(1)-firstHalf(1),secondHalf(2)-firstHalf(2),'Color','k');
end

end