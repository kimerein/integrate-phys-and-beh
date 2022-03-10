function plotChangeInReachCDF_throughSession(alltbt,trialTypes,metadata)

if ~isfield(trialTypes,'optoGroup')
    trialTypes.optoGroup=nan(size(trialTypes.led));
end

doTheseMice=unique(metadata.mouseid);

%fractionThroughSess_steps={[0-0.01 0.1+0.01],[0.1-0.01 0.2+0.01],[0.2-0.01 0.3+0.01],[0.3-0.01 0.4+0.01],[0.4-0.01 0.5+0.01],[0.5-0.01 0.6+0.01],[0.6-0.01 0.7+0.01],[0.7-0.01 0.8+0.01],[0.8-0.01 0.9+0.01],[0.9-0.01 1+0.01]};
% fractionThroughSess_steps={[0-0.01 0.2+0.01],[0.2-0.01 0.4+0.01],[0.4-0.01 0.6+0.01],[0.6-0.01 0.8+0.01],[0.8-0.01 1+0.01]};
% fractionThroughSess_steps={[0-0.01 0.2+0.01],[0.2-0.01 0.4+0.01],[0.4-0.01 0.6+0.01],[0.6-0.01 1+0.01]};
fractionThroughSess_steps={[0-0.01 0.2+0.01],[0.2-0.01 0.4+0.01],[0.3-0.01 0.5+0.01],[0.5-0.01 0.75+0.01]};
%fractionThroughSess_steps={[0-0.01 0.3+0.01],[0.3-0.01 0.6+0.01],[0.6-0.01 0.9+0.01]};
%fractionThroughSess_steps={[0-0.01 0.5+0.01],[0.5-0.01 1+0.01]};

backupbackup.alltbt=alltbt; backupbackup.trialTypes=trialTypes; backupbackup.metadata=metadata;

cuedComponentAcrossMice=nan(length(fractionThroughSess_steps),length(doTheseMice));
uncuedComponentAcrossMice=nan(length(fractionThroughSess_steps),length(doTheseMice));
% for j=1:length(doTheseMice)
for j=1:1
%     currMouse=doTheseMice(j);
%     
%     alltbt=backupbackup.alltbt; trialTypes=backupbackup.trialTypes; metadata=backupbackup.metadata;
%     %%%%%%%%%%%%%%%%%%%%% perform any filtering on alltbt
%     temp=datestr(datetime('now'));
%     temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
%     temp=temp(~isspace(temp));
%     saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
%     % filter settings
%     tbt_filter.sortField='mouseid';
%     tbt_filter.range_values=[currMouse-0.5 currMouse+0.5];
%     tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
%     temp=tbt_filter.name;
%     temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
%     temp=temp(~isspace(temp));
%     tbt_filter.name=temp;
%     tbt_filter.clock_progress=true;
%     % filter alltbt
%     [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
    
    backup.alltbt=alltbt; backup.trialTypes=trialTypes; backup.metadata=metadata;    
    fracsForPlot=nan(1,length(fractionThroughSess_steps));
    for i=1:length(fractionThroughSess_steps)
        curr_fracthroughsess=fractionThroughSess_steps{i};
        % exclude bad fits
        %%%% str sil bad fits
%         if j==1 && (i==length(fractionThroughSess_steps))
%             continue
%         end
%         if j==3 && (i==length(fractionThroughSess_steps))
%             continue
%         end
        %if j==4 && (i==length(fractionThroughSess_steps))
        %    continue
        %end
        %%%% control bad fits
%         if j==1 && (i==length(fractionThroughSess_steps)) %|| i==length(fractionThroughSess_steps)-1)
%             continue
%         end
%         if j==3 && (i==2)
%             continue
%         end
%         if j==2
%             continue
%         end
%         if j==8 && (i==length(fractionThroughSess_steps)) %|| i==length(fractionThroughSess_steps)-1)
%             continue
%         end
%         if j==9 
%             continue 
%         end
        fracsForPlot(i)=nanmean(curr_fracthroughsess);
        
        alltbt=backup.alltbt; trialTypes=backup.trialTypes; metadata=backup.metadata;
        
        %%%%%%%%%%%%%%%%%%%%% perform any filtering on alltbt
        temp=datestr(datetime('now'));
        temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
        temp=temp(~isspace(temp));
        saveDir=['/Volumes/Neurobio/MICROSCOPE/Kim/RT pairs data sets/' temp]; % where to save details of alltbt filtering and RT pairs data set
        
        % filter settings
        tbt_filter.sortField='fractionThroughSess';
        tbt_filter.range_values=curr_fracthroughsess;
        tbt_filter.name=[tbt_filter.sortField num2str(tbt_filter.range_values(1)) 'to' num2str(tbt_filter.range_values(2))];
        temp=tbt_filter.name;
        temp(~ismember(temp,['A':'Z' 'a':'z' '0':'9']))='';
        temp=temp(~isspace(temp));
        tbt_filter.name=temp;
        tbt_filter.clock_progress=true;
        
        % filter alltbt
        [alltbt,trialTypes,metadata]=filtTbt(alltbt,trialTypes,tbt_filter.sortField,tbt_filter.range_values,metadata,tbt_filter.clock_progress);
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%% build relevant data sets
        
        % settings for paired RT data set
        test.nInSequence=[2]; % defines trial pairs, e.g., 2 means will compare each trial with its subsequent trial, 3 means will compare each trial with the trial after next, etc.
        % requirement for first trial in pair
        %trial1='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
        trial1='trialTypes.optoGroup~=1';
        trialTypes.isLongITI_1forward=[trialTypes.isLongITI(2:end); 0];
        trialTypes.optoGroup_1forward=[trialTypes.optoGroup(2:end); 0];
        test.trial1=trial1;
        test.templateSequence2_cond=eval(trial1);
        %trial2='trialTypes.chewing_at_trial_start==0 | trialTypes.chewing_at_trial_start==1';
        trial2='trialTypes.optoGroup~=1 & trialTypes.led==1';
        %trial2='trialTypes.led==0 & trialTypes.optoGroup~=1 & trialTypes.optoGroup~=3';
        test.trial2=trial2;
        test.templateSequence2_end=eval(trial2);
        test.fillInBetweenWithAnything=true; % if true, will allow middle trials to be anything; otherwise, middle trials must match cond1
        test.event_name=['alltrials' tbt_filter.name 'inBetweenAnything' num2str(test.fillInBetweenWithAnything)];
        saveDir2=[];
        % build paired RT data set
        fakeCueInd=50; % in indices, this is not relevant if not using PCA-based RT model
        skipCorrected=true;
        % this function builds the dataset using the trial type sequences specified above
        [dataset]=buildReachingRTModel(alltbt,trialTypes,metadata,fakeCueInd,saveDir2,test,skipCorrected);
        
        rawreachout=getRawReaches(alltbt,dataset.realDistributions);
        
        % plot trial to trial change in reach CDF
        cdfout=plotChangeInReachCDF(dataset.realDistributions,alltbt);
        close all;
        
        if i==1
            reach_pdfs=nan(length(fractionThroughSess_steps),length(rawreachout.data1_x));
            reach_pdfs_x=nan(length(fractionThroughSess_steps),length(rawreachout.data1_x));
            reach_cdfs=nan(length(fractionThroughSess_steps),length(cdfout.meanCDF)+1);
            reach_uncuedests=nan(length(fractionThroughSess_steps),length(cdfout.y_fitToUncued));
            reach_timebins=nan(length(fractionThroughSess_steps),length(cdfout.y_fitToUncued));
        end
        reach_pdfs(i,:)=rawreachout.data2_n;
        reach_pdfs_x(i,:)=rawreachout.data2_x;
        reach_cdfs(i,:)=[0 cdfout.meanCDF];
        reach_uncuedests(i,:)=cdfout.y_fitToUncued;
        reach_timebins(i,:)=cdfout.timebins;
    end
    
    % Plot raw reaching
    figure();
    cmap=colormap('cool');
    k=1;
    kstep=ceil(size(cmap,1)/size(reach_pdfs_x,1));
    line([cdfout.cueTime cdfout.cueTime],[0 1],'Color','b');
    hold on;
    for i=1:size(reach_pdfs_x,1)
        plot(reach_pdfs_x(i,:),reach_pdfs(i,:),'Color',cmap(k,:),'LineWidth',1);
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    
    % Plot CDFs
    figure();
    cmap=colormap('cool');
    k=1;
    kstep=ceil(size(cmap,1)/size(reach_cdfs,1));
    cued_integrals=nan(1,size(reach_cdfs,1));
    uncued_integrals=nan(1,size(reach_cdfs,1));
    for i=1:size(reach_cdfs,1)
        plot(reach_timebins(i,:),reach_cdfs(i,:),'Color',cmap(k,:),'LineWidth',1);
        hold on;
        plot(reach_timebins(i,:),reach_uncuedests(i,:),'--','Color',cmap(k,:),'LineWidth',1);
        temp_fit=reach_uncuedests(i,:);
        temp_fit(temp_fit>1)=1;
        temp_fit(temp_fit<0)=0;
        tempcdf=reach_cdfs(i,:);
        uncued_integrals(i)=nansum(temp_fit(reach_timebins(i,:)>cdfout.cueTime & reach_timebins(i,:)<cdfout.cueTime+3))./nansum(tempcdf(reach_timebins(i,:)>cdfout.cueTime & reach_timebins(i,:)<cdfout.cueTime+3));
        cued_integrals(i)=nansum(tempcdf(reach_timebins(i,:)>cdfout.cueTime & reach_timebins(i,:)<cdfout.cueTime+3)-temp_fit(reach_timebins(i,:)>cdfout.cueTime & reach_timebins(i,:)<cdfout.cueTime+3))./nansum(tempcdf(reach_timebins(i,:)>cdfout.cueTime & reach_timebins(i,:)<cdfout.cueTime+3));
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    line([cdfout.cueTime cdfout.cueTime],[0 1],'Color','b');
    ylim([0 1]);
    xlim([0 9.5]);
    
    figure();
    k=1;
    for i=1:size(reach_cdfs,1)
        scatter(fracsForPlot(i),cued_integrals(i),[],cmap(k,:),'filled');
        hold on;
        scatter(fracsForPlot(i),uncued_integrals(i),[],'k','filled');
        k=k+kstep;
        if k>size(cmap,1)
            k=size(cmap,1);
        end
    end
    plot(fracsForPlot,cued_integrals,'-','Color','k');
    plot(fracsForPlot,uncued_integrals,'--','Color','k');
    ylim([0 1]);
    
    cuedComponentAcrossMice(:,j)=cued_integrals;
    uncuedComponentAcrossMice(:,j)=uncued_integrals;
    close all;
end

figure();
k=1;
plot(fracsForPlot,cuedComponentAcrossMice,'Color',[0.5 0.5 0.5]);
hold on;
for i=1:size(cuedComponentAcrossMice,1)
    scatter(fracsForPlot(i),cuedComponentAcrossMice(i,:),50,cmap(k,:),'LineWidth',1);
    hold on;
    line([fracsForPlot(i)-0.05 fracsForPlot(i)+0.05],[nanmean(cuedComponentAcrossMice(i,:)) nanmean(cuedComponentAcrossMice(i,:))],'Color',cmap(k,:),'LineWidth',2);
    line([fracsForPlot(i) fracsForPlot(i)],[nanmean(cuedComponentAcrossMice(i,:))-nanstd(cuedComponentAcrossMice(i,:))./sqrt(length(cuedComponentAcrossMice(i,:))) nanmean(cuedComponentAcrossMice(i,:))+nanstd(cuedComponentAcrossMice(i,:))./sqrt(length(cuedComponentAcrossMice(i,:)))],'Color',cmap(k,:),'LineWidth',2);
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
title('cued component');

figure();
k=1;
plot(fracsForPlot,uncuedComponentAcrossMice,'Color',[0.5 0.5 0.5]);
hold on;
for i=1:size(cuedComponentAcrossMice,1)
    scatter(fracsForPlot(i),uncuedComponentAcrossMice(i,:),50,cmap(k,:),'LineWidth',1);
    hold on;
    line([fracsForPlot(i)-0.05 fracsForPlot(i)+0.05],[nanmean(uncuedComponentAcrossMice(i,:)) nanmean(uncuedComponentAcrossMice(i,:))],'Color',cmap(k,:),'LineWidth',2);
    line([fracsForPlot(i) fracsForPlot(i)],[nanmean(uncuedComponentAcrossMice(i,:))-nanstd(uncuedComponentAcrossMice(i,:))./sqrt(length(uncuedComponentAcrossMice(i,:))) nanmean(uncuedComponentAcrossMice(i,:))+nanstd(uncuedComponentAcrossMice(i,:))./sqrt(length(uncuedComponentAcrossMice(i,:)))],'Color',cmap(k,:),'LineWidth',2);
    k=k+kstep;
    if k>size(cmap,1)
        k=size(cmap,1);
    end
end
title('uncued component');

end

function out=getRawReaches(alltbt,dataset)

timeStep=mode(diff(nanmean(alltbt.times,1)));
timeBinsForReaching=0:timeStep:(size(dataset.rawReaching_allTrialsSequence_trial1InSeq{1},2)-1)*timeStep;
for i=1:1
    out=plotTimeseries(dataset.rawReaching_event_trial1InSeq{i},dataset.se_rawReaching_event_trial1InSeq{i},'k',dataset.rawReaching_event_trialiInSeq{i},dataset.se_rawReaching_event_trialiInSeq{i},'m',timeBinsForReaching);
end
    
end

function out=whetherToSuppressPlots()

out=false;

end

function out=plotTimeseries(data1_mean,data1_se,color1,data2_mean,data2_se,color2,timeBins)

suppPlots=whetherToSuppressPlots();

plotAsCityscape=true;
correctForBinSize=true;

if correctForBinSize==true
    data1_mean=data1_mean./mode(diff(timeBins));
    data2_mean=data2_mean./mode(diff(timeBins));
    data1_se=data1_se./mode(diff(timeBins));
    data2_se=data2_se./mode(diff(timeBins));
end

data1_mean=nanmean(data1_mean,1);
data2_mean=nanmean(data2_mean,1);
data1_se=sqrt(nansum(data1_se.^2,1));
data2_se=sqrt(nansum(data2_se.^2,1));

if suppPlots==false
    figure();
end
%fill([timeBins fliplr(timeBins)],[data1_mean+data1_se fliplr(data1_mean-data1_se)],[0.5 0.5 0.5]);
%hold on;
if plotAsCityscape==true
    if suppPlots==false
        [n,x]=cityscape_hist(data1_mean,timeBins);
        out.data1_x=x; 
        out.data1_n=n;
        plot(x,n,'Color',color1); hold on;
        [n,x]=cityscape_hist(data1_mean+data1_se,timeBins);
        plot(x,n,'Color',color1);
        [n,x]=cityscape_hist(data1_mean-data1_se,timeBins);
        plot(x,n,'Color',color1);
    end
else
    if suppPlots==false    
        plot(timeBins,data1_mean,'Color',color1); hold on;
        plot(timeBins,data1_mean+data1_se,'Color',color1);
        plot(timeBins,data1_mean-data1_se,'Color',color1);
    end
end

%fill([timeBins fliplr(timeBins)],[data2_mean+data2_se fliplr(data2_mean-data2_se)],[0.1 0.7 0.5]);
%hold on;
if plotAsCityscape==true
    if suppPlots==false   
        [n,x]=cityscape_hist(data2_mean,timeBins);
        out.data2_x=x; 
        out.data2_n=n;
        plot(x,n,'Color',color2); hold on;
        [n,x]=cityscape_hist(data2_mean+data2_se,timeBins);
        plot(x,n,'Color',color2);
        [n,x]=cityscape_hist(data2_mean-data2_se,timeBins);
        plot(x,n,'Color',color2);
    end
else
    if suppPlots==false    
        plot(timeBins,data2_mean,'Color',color2); hold on;
        plot(timeBins,data2_mean+data2_se,'Color',color2);
        plot(timeBins,data2_mean-data2_se,'Color',color2);
    end
end

end