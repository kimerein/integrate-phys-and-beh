function plotGLMcoef(glm_coef,glm_intercept,feature_names,timestep,nShiftsBefore)

smoobin=10;
event_types = {'cue', 'opto', 'distract', 'reachzone', 'fidget', 'success', 'drop', 'missing', 'failure', 'chew'};

l=1;
coef=[];
allcoef=[];
figure();
cm=jet(length(event_types));
for i=1:length(event_types)
    currev=event_types{i};
    k=1;
    for j=1:size(feature_names,1)
        if ~isempty(regexp(feature_names(j,:),currev))
            coef(k)=glm_coef(l);
            k=k+1;
            l=l+1;
        end
    end
    if isempty(allcoef)
        allcoef=nan(length(event_types),length(coef));
    end
    allcoef(i,:)=smooth(coef,smoobin);
end
for i=1:length(event_types)
    subplot(1,length(event_types),i);
    ts=0:timestep:(size(allcoef,2)-1)*timestep;
    ts=ts-nShiftsBefore*timestep;
    plot(ts,allcoef(i,:),'Color',cm(i,:));
    xlim([ts(1) ts(end)]);
    hold on;
    line([0 0],[-0.1 0.5],'Color',[0.3 0.3 0.3]);
    line([ts(1) ts(end)],[0 0],'Color',[0.3 0.3 0.3]);
    xlabel(event_types{i});
    ylim([nanmin(allcoef(1:end)) nanmax(allcoef(1:end))]);
end
set(gcf,'Position',[20 20 1300 400]);


