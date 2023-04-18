function learningCurves(alltbt,trialTypes,metadata,dayField)

fillInToEnd=true;

u=unique(metadata.mouseid);
dayField=metadata.(dayField);
learnCurves=nan(length(u),200); 
learnCurves_distract=nan(length(u),200); 
udays=-50:1:149;
for i=1:length(u)
    currmouseid=u(i);
    subday=dayField(metadata.mouseid==currmouseid);
    subdprimes=metadata.dprimes(metadata.mouseid==currmouseid);
    subdprimesdistract=metadata.distract_dprimes(metadata.mouseid==currmouseid);
    [udays_for_mouse,ui]=unique(subday);
    dp=subdprimes(ui);
    distract_dp=subdprimesdistract(ui);
    for j=1:length(udays_for_mouse)
        f=find(udays==udays_for_mouse(j));
        if isempty(f)
            error(['Array missing day ' num2str(udays_for_mouse(j))]);
        end
        learnCurves(i,f)=dp(j);
        learnCurves_distract(i,f)=distract_dp(j);
    end
end
if fillInToEnd==true
    for i=1:size(learnCurves,1)
        temp=learnCurves(i,:);
        f=find(~isnan(temp),1,'last');
        learnCurves(i,end)=temp(f);
    end
end
lc=learnCurves;
lcminusdistract=learnCurves-learnCurves_distract;
% linearly interpolate
for i=1:size(lc,1)
    lc(i,:)=interpMissing(lc(i,:));
    lcminusdistract(i,:)=interpMissing(lcminusdistract(i,:));
end

figure();
plot(udays,lc'); xlabel('days'); ylabel('dprime'); hold on;
plot(udays,mean(lc,1,'omitnan'),'Color','k','LineWidth',2);

figure();
plot(udays,mean(lc,1,'omitnan'),'Color','k','LineWidth',2); hold on;
plot(udays,mean(lc,1,'omitnan')+std(lc,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);
plot(udays,mean(lc,1,'omitnan')-std(lc,[],1,'omitnan')./sqrt(size(lc,1)),'Color','k','LineWidth',1);

figure(); 
plot(udays,lcminusdistract'); xlabel('days'); ylabel('dprime minus distract dprime'); hold on;
plot(udays,mean(lcminusdistract,1,'omitnan'),'Color','k','LineWidth',2);

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