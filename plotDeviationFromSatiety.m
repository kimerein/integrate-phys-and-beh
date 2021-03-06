function [m_binned,refsToMbinned,cued_binned,uncued_binned]=plotDeviationFromSatiety(reachRates,nbins)

% dim 1 of alltrials_cued is across sessions, dim 2 of alltrials_cued is
% across trials
% nbins is how many bins for grouping trials
m_binned=downSampAv(reachRates.m,ceil(length(reachRates.m)/nbins));
cued_binned=downSampAv(reachRates.cued,ceil(length(reachRates.m)/nbins));
uncued_binned=downSampAv(reachRates.uncued,ceil(length(reachRates.m)/nbins));
refToOrig=1:length(m_binned);
j=1;
refsToMbinned=nan(length(reachRates.m),1);
for i=1:length(refToOrig)
    if j+nbins-1>length(reachRates.m)
        refsToMbinned(j:length(reachRates.m))=ones(length(j:length(reachRates.m)),1)*refToOrig(i);
    else
        refsToMbinned(j:j+ceil(length(reachRates.m)/nbins)-1)=ones(length(j:j+ceil(length(reachRates.m)/nbins)-1),1)*refToOrig(i);
    end
    j=j+ceil(length(reachRates.m)/nbins);
end

cmap=colormap('cool');
k=1;
kstep=ceil(size(cmap,1)/length(m_binned));
figure(); 
for i=1:length(m_binned)
    if i==length(m_binned)
    else
        length_line_segment=(cued_binned(i+1)-cued_binned(i))/m_binned(i);
    end
    line([uncued_binned(i)-length_line_segment/2 uncued_binned(i)+length_line_segment/2],[cued_binned(i)-(length_line_segment*m_binned(i))/2 cued_binned(i)+(length_line_segment*m_binned(i))/2],'Color','k','LineWidth',1);
    hold on;
end
for i=1:length(uncued_binned)
    scatter(uncued_binned(i),cued_binned(i),[],cmap(k,:),'filled');
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

k=1;
kstep=ceil(size(cmap,1)/length(reachRates.uncued));
for i=1:length(reachRates.uncued)
    scatter(reachRates.uncued(i),reachRates.cued(i),[],cmap(k,:));
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end

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