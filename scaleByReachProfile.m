function transformtbt=scaleByReachProfile(alltbt,trialTypes,metadata,ledIsOn,useLongITIonly,reachName,durationOfWheelTurn,template)

warning('off');

cueName='cueZone_onVoff';
useReaches=alltbt.(reachName);

if ~isempty(ledIsOn)
    trialTypes.ledIsOn=ledIsOn;
else
    trialTypes.ledIsOn=trialTypes.led;
end
if useLongITIonly==true
    trialTypes.ITItoUse=trialTypes.isLongITI==1;
else
    trialTypes.ITItoUse=ones(size(trialTypes.led));
end

% cue turns on here
temp=nanmean(alltbt.(cueName),1);
[~,cueInd]=nanmax(temp);

% inds for wheel to turn
timeStep=mode(mode(nanmean(alltbt.times(:,2:end),1)-nanmean(alltbt.times(:,1:end-1),1)));
indsForWheelTurn=floor(durationOfWheelTurn/timeStep); % both in sec

% For each day for each mouse, get the response epochs
unique_mice=unique(metadata.mouseid);
unique_sess=unique(metadata.nth_session);

plotThis_i=randi(length(unique_mice));
plotThis_j=randi(length(unique_sess));
newReaches=nan(size(useReaches));
for i=1:length(unique_mice)
    for j=1:length(unique_sess)
        if mod(j,10)==0
            disp(j);
        end
        temp=useReaches(ismember(metadata.nth_session,unique_sess(j)) & ismember(metadata.mouseid,unique_mice(i)),:);
        avReaches=nanmean(temp,1);
        epochs=findEpochs(avReaches,cueInd,indsForWheelTurn,timeStep,alltbt);
        if i==plotThis_i && j==plotThis_j
            figure(); 
            title('Example session');
            plot(avReaches,'Color','k');
            hold on;
            plotVertLineAt(epochs.baseline1(1),[0 nanmax(avReaches)],'c');
            plotVertLineAt(epochs.preempt(1),[0 nanmax(avReaches)],'r');
            plotVertLineAt(epochs.cued(1),[0 nanmax(avReaches)],'b');
            plotVertLineAt(epochs.prereach(1),[0 nanmax(avReaches)],'m');
            plotVertLineAt(epochs.postreach(1),[0 nanmax(avReaches)],'g');
            plotVertLineAt(epochs.baseline2(1),[0 nanmax(avReaches)],[0.2 0.5 1]);
            plotVertLineAt(epochs.baseline2(2),[0 nanmax(avReaches)],'y');
        end
        
        if j==1 && i==1
            % set up template
            if isempty(template)
                template.baseline1=epochs.baseline1;
                template.preempt=epochs.preempt;
                template.cued=epochs.cued;
                template.prereach=[epochs.cued(2)+1 epochs.cued(2)+floor(1.5/timeStep)];
                template.postreach=[template.prereach(2)+1 template.prereach(2)+floor(3/timeStep)];
                template.baseline2=[template.postreach(2)+1 epochs.baseline2(2)];
            end
        end
        
        % for all trials in this session, scale TIME such that matches
        % template epochs
        % note that baselines, preempt and cued are fixed
        % template prereach extends for 1.5 s
        % template postreach extends for 3 s
        outdata=stretchToTemplate(epochs,template,temp);
        newReaches(ismember(metadata.nth_session,unique_sess(j)) & ismember(metadata.mouseid,unique_mice(i)),:)=outdata;
        
        if i==plotThis_i && j==plotThis_j
            figure(); 
            title('Example session AFTER TIME RESCALE');
            navReaches=nanmean(newReaches,1);
            plot(navReaches,'Color','k');
            hold on;
            plotVertLineAt(template.baseline1(1),[0 nanmax(navReaches)],'c');
            plotVertLineAt(template.preempt(1),[0 nanmax(navReaches)],'r');
            plotVertLineAt(template.cued(1),[0 nanmax(navReaches)],'b');
            plotVertLineAt(template.prereach(1),[0 nanmax(navReaches)],'m');
            plotVertLineAt(template.postreach(1),[0 nanmax(navReaches)],'g');
            plotVertLineAt(template.baseline2(1),[0 nanmax(navReaches)],[0.2 0.5 1]);
            plotVertLineAt(template.baseline2(2),[0 nanmax(navReaches)],'y');
        end
    end
end

transformtbt=alltbt;
transformtbt.(reachName)=newReaches;

end

function outdata=stretchToTemplate(epochs,template_epochs,data)

% up to prereach
uptoprereach_data=data(:,1:epochs.prereach(1)-1);

% prereach
prereach_data=data(:,epochs.prereach(1):epochs.prereach(2));
if duration(epochs.prereach)~=duration(template_epochs.prereach)
    % stretch or shrink
    % get scale factor
    scale_fac=duration(template_epochs.prereach)/duration(epochs.prereach);
    if scale_fac>20
        scale_fac=20;
    elseif size(prereach_data,2)*scale_fac<2
        scale_fac=1;
    end
    templateTimes=linspace(0,duration(template_epochs.prereach),size(prereach_data,2)*scale_fac);
    currTimes=linspace(0,duration(template_epochs.prereach),size(prereach_data,2));
    % resample data
    new_prereach_data=resampleData(currTimes,prereach_data,templateTimes);
else
    new_prereach_data=prereach_data;
end

% postreach
postreach_data=data(:,epochs.postreach(1):epochs.postreach(2));
if duration(epochs.postreach)~=duration(template_epochs.postreach)
    % stretch or shrink
    % get scale factor
    scale_fac=duration(template_epochs.postreach)/duration(epochs.postreach);
    if scale_fac>20
        scale_fac=20;
    elseif size(postreach_data,2)*scale_fac<2
        scale_fac=1;
    end
    templateTimes=linspace(0,duration(template_epochs.postreach),size(postreach_data,2)*scale_fac);
    currTimes=linspace(0,duration(template_epochs.postreach),size(postreach_data,2));
    % resample data
    new_postreach_data=resampleData(currTimes,postreach_data,templateTimes);
else
    new_postreach_data=postreach_data;
end

% after reach
afterreach_data=data(:,epochs.postreach(2)+1:end);

% concat
newdata=[uptoprereach_data new_prereach_data new_postreach_data afterreach_data];
if size(newdata,2)<size(data,2)
    % expand
    outdata=[newdata nan(size(data,1),size(data,2)-size(newdata,2))];
else
    % or cut
    outdata=newdata(:,1:size(data,2));
end

end

function newData=resampleData(currTimes,data,newTimes)

newData=[];

for i=1:size(data,1)
    temp=data(i,:);
    curr=timeseries(temp,currTimes);
    newTemp=resample(curr,newTimes);
    newData(i,:)=newTemp.data;
end

end

function a=duration(d)

a=d(2)-d(1);

end

function plotVertLineAt(x,y_range,c)

line([x x],[y_range(1) y_range(2)],'Color',c);

end

function epochs=findEpochs(avReaches,cueInd,indsForWheelTurn,timeStep,alltbt)

% return epochs in inds

% baseline 1 before wheel turns
% begin when more than 95% of trials have no nan
beg=find(nanmean(~isnan(alltbt.times),1)>0.95,1,'first');
epochs.baseline1=[beg cueInd-indsForWheelTurn];

% preempt is during wheel turn up to and including cueInd
epochs.preempt=[cueInd-indsForWheelTurn+1 cueInd];

% cued is 50 ms beginning at the time of the cue
epochs.cued=[cueInd+1 cueInd+ceil(0.05/timeStep)];

% pre-reach is from cueInd to the peak of reaching, given a peak at least 2
% times the standard deviation of the noise during the baseline, else
% cued window is up to 5 sec after cue
s=nanstd(avReaches(epochs.baseline1(1):epochs.baseline1(2)));
inds5sec=floor(5/timeStep);
if any(avReaches(cueInd+1:cueInd+inds5sec)>nanmean(avReaches(epochs.baseline1(1):epochs.baseline1(2)))+2*s)
    [~,ma]=nanmax(avReaches(cueInd+1:cueInd+inds5sec));
    peakInd=cueInd+ma;
else
    peakInd=cueInd+inds5sec;
end
epochs.prereach=[cueInd+ceil(0.05/timeStep) peakInd];
if epochs.prereach(2)==epochs.prereach(1)
    epochs.prereach(2)=epochs.prereach(1)+1;
end

% post-reach is from peak of reaching to baseline 2
f=find(avReaches(epochs.prereach(2)+1:end)<nanmean(avReaches(epochs.baseline1(1):epochs.baseline1(2))),1,'first');
if isempty(f)
    f=find(avReaches(epochs.prereach(2)+1:end)<nanmean(avReaches(epochs.baseline1(1):epochs.baseline1(2)))+2*s,1,'first');
    if isempty(f)
        f=100;
    end
end
epochs.postreach=[epochs.prereach(2)+1 epochs.prereach(2)+f];
if epochs.postreach(2)==epochs.postreach(1)
    epochs.postreach(2)=epochs.postreach(1)+1;
end

% baseline 2 is until <95% of trials have no nan
en=find(nanmean(~isnan(alltbt.times),1)>0.95,1,'last');
epochs.baseline2=[epochs.prereach(2)+f+1 en];

end