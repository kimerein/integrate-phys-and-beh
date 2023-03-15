function [ts,allcoef]=plotGLMcoef(glm_coef,glm_intercept,feature_names,timestep,nShiftsBefore,putTog)

% GLM python
% plotGLMcoef(glm_coef,glm_intercept,feature_names,10*0.01,9);
% GLM matlab
% plotGLMcoef(coef,[],fnames,10*0.01,9);

% putTog can be 'mean','median','me+se','me-se'

smoobin=1;

if iscell(feature_names)
    event_types=getFeatureNames(feature_names);
else
    event_types=getFeatureNames(feature_names);
end
if isempty(glm_intercept)
    glm_intercept=0;
end

l=1;
coef=[];
allcoef=[];
figure();
cm=jet(length(event_types));
for i=1:length(event_types)
    currev=event_types{i};
    k=1;
    for j=1:size(feature_names,1)
        if iscell(feature_names)
            if ~isempty(regexp(feature_names{j},currev))
                if size(glm_coef,1)>1 && size(glm_coef,2)>1
                    coef(k)=putTogether(glm_coef(:,l),putTog);
                else
                    coef(k)=glm_coef(l);
                end
                k=k+1;
                l=l+1;
            end
        else
            if ~isempty(regexp(feature_names(j,:),currev))
                if size(glm_coef,1)>1 && size(glm_coef,2)>1
                    coef(k)=putTogether(glm_coef(:,l),putTog);
                else
                    coef(k)=glm_coef(l);
                end
                k=k+1;
                l=l+1;
            end
        end
    end
    if isempty(allcoef)
        allcoef=nan(length(event_types),length(coef));
    end
    %allcoef(i,:)=smooth(coef,smoobin);
    if smoobin~=1
        allcoef(i,:)=smoothdata(coef,'gaussian',smoobin);
    else
        allcoef(i,:)=coef;
    end
end
slopes=nan(1,length(event_types));
for i=1:length(event_types)
    subplot(1,length(event_types),i);
    ts=0:timestep:(size(allcoef,2)-1)*timestep;
    ts=ts-nShiftsBefore*timestep;
    plot(ts,allcoef(i,:),'Color',cm(i,:));
    if ~isempty(regexp(event_types{i},'cue','once'))
        params.tapers=[1 3];
        params.Fs=1/mode(diff(ts));
        params.trialave=1;
        filter_range=[3 6]; 
        addpath(genpath('C:\Users\sabatini\Documents\GitHub\chronux_2_11'));
        [amp,S,t,f]=getAmpWithChronux(allcoef(i,:),[0.5 0.1],params,filter_range);
        rmpath(genpath('C:\Users\sabatini\Documents\GitHub\chronux_2_11'));
        hold on; plot(t-(nShiftsBefore*timestep),amp*10,'Color',[0 0 0.5]);
    end
    P=polyfit(ts,allcoef(i,:),1);
    yfit=P(1)*ts+P(2);
    hold on; plot(ts,yfit,'Color','m');
    slopes(i)=P(1);
    xlim([ts(1) ts(end)]);
    hold on;
    line([0 0],[-0.1 0.5],'Color',[0.3 0.3 0.3]);
    line([ts(1) ts(end)],[0 0],'Color',[0.3 0.3 0.3]);
    xlabel(event_types{i});
    ylim([nanmin(allcoef(1:end)) nanmax(allcoef(1:end))]);
end
set(gcf,'Position',[20 20 1300 400]);

end

function [amp,S,t,f]=getAmpWithChronux(data,movingwin,params,filter_range)

[S,t,f]=mtspecgramc(data,movingwin,params);

S=sqrt(S);

amp=nanmean(S(:,f>=filter_range(1) & f<=filter_range(2)),2)';

end

function out=putTogether(co,method)

% method can be 'mean','median','me+se','me-se'
switch method
    case 'mean'
        out=mean(co,'all','omitnan');
    case 'median'
        out=median(co,'all','omitnan');
    case 'me+se'
        out=mean(co,'all','omitnan')+std(co,[],'all','omitnan')./sqrt(length(co));
    case 'me-se'
        out=mean(co,'all','omitnan')-std(co,[],'all','omitnan')./sqrt(length(co));
    otherwise
        disp('Averaging as default');
        out=mean(co,'all','omitnan');
end

end

function event_types=getFeatureNames(feature_names)

potential_event_types = {'cue', 'opto', 'distract', 'reach', 'fidget', 'success', 'drop', 'miss', 'failure', 'chew', ...
                          'cXsuc', 'cXdro', 'cXmis','uXsuc','uXdro','uXmis'};
evcount=1; 
for i=1:length(potential_event_types)
    for j=1:length(feature_names)
        if iscell(feature_names)
            if ~isempty(regexp(feature_names{j},potential_event_types{i}))
                event_types{evcount}=potential_event_types{i};
                evcount=evcount+1;
                break
            end
        else
            if ~isempty(regexp(feature_names(j,:),potential_event_types{i}))
                event_types{evcount}=potential_event_types{i};
                evcount=evcount+1;
                break
            end
        end
    end
end

end
