function [lc,udays]=learningCurves(alltbt,trialTypes,metadata,dayField,day1is,dayNis)

fillInToEnd=true;

u=unique(metadata.mouseid);
dayField=metadata.(dayField);
learnCurves=nan(length(u),200); 
learnCurves_distract=nan(length(u),200); 
rr_cued=nan(length(u),200); 
rr_uncued=nan(length(u),200); 
udays=-50:1:149;
for i=1:length(u)
    currmouseid=u(i);
    subday=dayField(metadata.mouseid==currmouseid);
    subdprimes=metadata.dprimes(metadata.mouseid==currmouseid);
    subdprimesdistract=metadata.distract_dprimes(metadata.mouseid==currmouseid);
    sub_rr_cued=metadata.reachrate_cued(metadata.mouseid==currmouseid);
    sub_rr_uncued=metadata.reachrate_uncued(metadata.mouseid==currmouseid);
    [udays_for_mouse,ui]=unique(subday);
    dp=subdprimes(ui);
    distract_dp=subdprimesdistract(ui);
    rcue=sub_rr_cued(ui);
    runcue=sub_rr_uncued(ui);
    for j=1:length(udays_for_mouse)
        f=find(udays==udays_for_mouse(j));
        if isempty(f)
            error(['Array missing day ' num2str(udays_for_mouse(j))]);
        end
        learnCurves(i,f)=dp(j);
        learnCurves_distract(i,f)=distract_dp(j);
        rr_cued(i,f)=rcue(j);
        rr_uncued(i,f)=runcue(j);
    end
end
if fillInToEnd==true
    for i=1:size(learnCurves,1)
        temp=learnCurves(i,:);
        f=find(~isnan(temp),1,'last');
        learnCurves(i,end)=temp(f);

        temp=rr_cued(i,:);
        f=find(~isnan(temp),1,'last');
        rr_cued(i,end)=temp(f);

        temp=rr_uncued(i,:);
        f=find(~isnan(temp),1,'last');
        rr_uncued(i,end)=temp(f);
    end
end
lc=learnCurves;
lcminusdistract=learnCurves-learnCurves_distract;
rr_cued_interp=rr_cued;
rr_uncued_interp=rr_uncued;
% linearly interpolate
for i=1:size(lc,1)
    lc(i,:)=interpMissing(lc(i,:));
    lcminusdistract(i,:)=interpMissing(lcminusdistract(i,:));
    rr_cued_interp(i,:)=interpMissing(rr_cued_interp(i,:));
    rr_uncued_interp(i,:)=interpMissing(rr_uncued_interp(i,:));
end

figure();
plot(udays,lc'); xlabel('days'); ylabel('dprime'); hold on;
plot(udays,mean(lc,1,'omitnan'),'Color','k','LineWidth',2);

figure();
plot(udays,mean(lc,1,'omitnan'),'Color','k','LineWidth',2); hold on; xlabel('days'); ylabel('dprime');
plot(udays,mean(lc,1,'omitnan')+std(lc,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);
plot(udays,mean(lc,1,'omitnan')-std(lc,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);

figure();
plot(udays,mean(learnCurves_distract,1,'omitnan'),'Color','k','LineWidth',2); hold on; xlabel('days'); ylabel('distract dprime');
plot(udays,mean(learnCurves_distract,1,'omitnan')+std(learnCurves_distract,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);
plot(udays,mean(learnCurves_distract,1,'omitnan')-std(learnCurves_distract,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1); 

figure(); 
plot(udays,lcminusdistract'); xlabel('days'); ylabel('dprime minus distract dprime'); hold on;
plot(udays,mean(lcminusdistract,1,'omitnan'),'Color','k','LineWidth',2);

% figure(); 
% plot(udays,max(cat(3,lc,lcminusdistract),[],3,'omitnan')'); xlabel('days'); ylabel('max(dprime,dprime minus distract)'); hold on;
% plot(udays,mean(max(cat(3,lc,lcminusdistract),[],3,'omitnan'),1,'omitnan'),'Color','k','LineWidth',2);

[lc_day1,lc_dayN]=getFirstAndLastRR(lc,lc,day1is,dayNis,udays);   
figure(); histogram(lc_dayN-lc_day1,10); xlabel('Change in dprime'); ylabel('Count');

[day1_rr_cued,dayN_rr_cued,day1_rr_uncued,dayN_rr_uncued]=getFirstAndLastRR(rr_cued_interp,rr_uncued_interp,day1is,dayNis,udays);
quiverPlot(day1_rr_cued,dayN_rr_cued,day1_rr_uncued,dayN_rr_uncued);

doSmooth=true;
if doSmooth
    smoothbin=5;
    for i=1:size(rr_cued_interp,1)
        rr_cued_interp(i,ismember(udays,2:20))=smooth(rr_cued_interp(i,ismember(udays,2:20)),smoothbin);
        rr_uncued_interp(i,ismember(udays,2:20))=smooth(rr_uncued_interp(i,ismember(udays,2:20)),smoothbin);
    end
end
figure();
plot(udays,mean(rr_cued_interp,1,'omitnan'),'Color','k','LineWidth',2); 
hold on; plot(udays,mean(rr_uncued_interp,1,'omitnan'),'Color','r','LineWidth',2);
plot(udays,mean(rr_cued_interp,1,'omitnan')+std(rr_cued_interp,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);
plot(udays,mean(rr_cued_interp,1,'omitnan')-std(rr_cued_interp,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);
plot(udays,mean(rr_uncued_interp,1,'omitnan')+std(rr_uncued_interp,[],1,'omitnan')./sqrt(size(lc,1)),'Color','r','LineWidth',1);
plot(udays,mean(rr_uncued_interp,1,'omitnan')-std(rr_uncued_interp,[],1,'omitnan')./sqrt(size(lc,1)),'Color','r','LineWidth',1);

end

function quiverPlot(day1_rr_cued,dayN_rr_cued,day1_rr_uncued,dayN_rr_uncued)

figure();
for i=1:length(day1_rr_cued)
    quiver(0,0,dayN_rr_uncued(i)-day1_rr_uncued(i),dayN_rr_cued(i)-day1_rr_cued(i),'Color','k'); hold on;
end
xlabel('Change in uncued reach rate');
ylabel('Change in cued reach rate');

end

function [day1_rr_cued,dayN_rr_cued,day1_rr_uncued,dayN_rr_uncued]=getFirstAndLastRR(rr_cued,rr_uncued,day1is,dayNis,udays)

day1_rr_cued=nan(size(rr_cued,1),1);
dayN_rr_cued=nan(size(rr_cued,1),1);
day1_rr_uncued=nan(size(rr_cued,1),1);
dayN_rr_uncued=nan(size(rr_cued,1),1);
for i=1:size(rr_cued,1)
    if ischar(day1is)
        switch day1is
            case 'first'
                temp=rr_cued(i,:);
                f=find(~isnan(temp),1,'first');
                day1_rr_cued(i)=temp(f);
                temp=rr_uncued(i,:);
                f=find(~isnan(temp),1,'first');
                day1_rr_uncued(i)=temp(f);
            case 'last'
                temp=rr_cued(i,:);
                f=find(~isnan(temp),1,'last');
                day1_rr_cued(i)=temp(f);
                temp=rr_uncued(i,:);
                f=find(~isnan(temp),1,'last');
                day1_rr_uncued(i)=temp(f);
        end
    else
        temp=rr_cued(i,:);
        day1_rr_cued(i)=nanmean(temp(ismember(udays,day1is)));
        temp=rr_uncued(i,:);
        day1_rr_uncued(i)=nanmean(temp(ismember(udays,day1is)));
    end

    if ischar(dayNis)
        switch dayNis
            case 'first'
                temp=rr_cued(i,:);
                f=find(~isnan(temp),1,'first');
                dayN_rr_cued(i)=temp(f);
                temp=rr_uncued(i,:);
                f=find(~isnan(temp),1,'first');
                dayN_rr_uncued(i)=temp(f);
            case 'last'
                temp=rr_cued(i,:);
                f=find(~isnan(temp),1,'last');
                dayN_rr_cued(i)=temp(f);
                temp=rr_uncued(i,:);
                f=find(~isnan(temp),1,'last');
                dayN_rr_uncued(i)=temp(f);
        end
    else
        temp=rr_cued(i,:);
        dayN_rr_cued(i)=nanmean(temp(ismember(udays,dayNis)));
        temp=rr_uncued(i,:);
        dayN_rr_uncued(i)=nanmean(temp(ismember(udays,dayNis)));
    end
end

end

function data=interpMissing(data)

% Find the indices of the nans
nan_indices = find(isnan(data));

% Find the indices of the non-nans
non_nan_indices = find(~isnan(data));

% Use interp1 to interpolate the missing values
interpolated_data = interp1(non_nan_indices, data(non_nan_indices), nan_indices, 'linear');

% Replace the nans with the interpolated values
data(nan_indices) = interpolated_data;

end