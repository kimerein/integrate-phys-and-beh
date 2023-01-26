function [data1,data2]=useDistractorAlignment(data1,whichTime1,whichField1,data2,whichTime2,whichField2,whichToShift,downSampData2,alignmentAnchor_data1,alignmentAnchor_data2)

% whichToShift is either 'data1' or 'data2' -- indicates which dataset to
% resample for alignment; other dataset will not be changed
% 
% whichTime1 is time matched to distractor
% whichField1 is distractor
% whichTime2 is time matched to distractor
% whichField2 is distractor
% for photometry
% settings.try_delay1=-150;
% settings.try_delay2=150;
% settings.tryinc=0.0005; 
% settings.try_scale1=0.8;
% settings.try_scale2=1.2;  
% for physiology
% INITIAL ALIGNMENT
settings.try_delay1=-10; %50; %-10;
settings.delaysteps=1;
settings.try_delay2=10; %70; %10;
settings.tryinc=0.01;
% ROW BY ROW
% this is red before black
settings.forSearchMinus=-10; %100; % inds around optimal for search for each row
% this is black before red
settings.forSearchPlus=10; %200; % inds around optimal for search for each row
% settings.try_scale1=0.6;
% settings.try_scale2=1;  
alignInd=11;
% downSampData2=true;
ds=1; %1000;

if downSampData2==true
    a=questdlg('Beware will downsample data2 and return downsampled. Do you want to continue? ');
    if ~strcmp(a,'Yes')
        return
    end
end

if downSampData2==true
    f=fieldnames(data2);
    for i=1:length(f)
        temp=data2.(f{i});
        data2.(f{i})=downSampMatrix(temp,ds);
    end
end

% Note that slow video DVR drops frames
% Thus, passage of time in behavior_tbt is not uniform
% Need to fix this before aligning to photometry or physiology, which do
% not drop times
% if strcmp(whichField1,'movie_distractor')
%     % fix data1
%     data1=makeTimeUniform(data1,'timesfromarduino',whichTime1,'cueZone_onVoff');
% %     timestep=mode(diff(nanmean(data1.(whichTime1),1)));
% %     data1.(whichTime1)=repmat(0:timestep:(size(data1.(whichTime1),2)-1)*timestep,size(data1.(whichTime1),1),1);
%     data1.(whichTime1)=data1.uniformtime-repmat(min(data1.uniformtime,[],2,'omitnan'),1,size(data1.uniformtime,2));
% elseif strcmp(whichField2,'movie_distractor')
%     % fix data2
%     data2=makeTimeUniform(data2,'timesfromarduino',whichTime2,'cueZone_onVoff');
% %     timestep=mode(diff(nanmean(data2.(whichTime2),1)));
% %     data2.(whichTime2)=repmat(0:timestep:(size(data2.(whichTime2),2)-1)*timestep,size(data2.(whichTime2),1),1);
%     data2.(whichTime2)=data2.uniformtime-repmat(min(data2.uniformtime,[],2,'omitnan'),1,size(data2.uniformtime,2));
% end

switch whichToShift
    case 'data1'
        % resample data 1
        disp('Black is first argument');
    case 'data2'
        disp('Black is fourth argument');
        % switch data 1 and data 2
        temp_data1=data1;
        temp_whichTime1=whichTime1;
        temp_whichField1=whichField1;
        temp_alignmentAnchor_data1=alignmentAnchor_data1;
        data1=data2;
        whichTime1=whichTime2;
        whichField1=whichField2;
        alignmentAnchor_data1=alignmentAnchor_data2;
        data2=temp_data1;
        whichTime2=temp_whichTime1;
        whichField2=temp_whichField1;
        alignmentAnchor_data2=temp_alignmentAnchor_data1;
    otherwise
        error('Unrecognized value passed in as whichToShift');
end

% Plot starting alignment
dis1=data1.(whichField1);
t1=data1.(whichTime1);
dis2=data2.(whichField2);
t2=data2.(whichTime2);

if ~isempty(alignmentAnchor_data1)
    anchordata1=mean(data1.(alignmentAnchor_data1),1,'omitnan');
    [~,anchor1]=find(anchordata1>0.5,1,'first');
    anchordata2=mean(data2.(alignmentAnchor_data2),1,'omitnan');
    [~,anchor2]=find(anchordata2>0.5,1,'first');
    for i=1:size(t1,1)
        % anchor 1 and anchor 2 are same time
        temp1=t1(i,anchor1);
        temp2=t2(i,anchor2);
        t1(i,:)=t1(i,:)-(temp1-temp2);
    end
    data1.(whichTime1)=t1;
end

figure(); plot(t2(alignInd,:),dis2(alignInd,:),'Color','r'); 
figure(); plot(t1(alignInd,:),dis1(alignInd,:),'Color','k');
scaleMid=input('What is scaling of red divided by black time duration? ');
settings.try_scale1=scaleMid-0.2;
settings.try_scale2=scaleMid+0.2;  

figure();
offset=1;
dis1=dis1./max(dis1(1:end),[],2,'omitnan');
dis2=dis2./max(dis2(1:end),[],2,'omitnan');
for i=1:size(data1.(whichTime1),1)
    plot(t1(i,:),offset+0.9*dis1(i,:),'Color','k');
    hold on;
    plot(t2(i,:),offset+0.9*dis2(i,:),'Color','r');
%     offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
    offset=offset+1;
end

% Use distractor to further align
% Resample and shift data1 times
disp('Now refining alignment ...');
tryinc=settings.tryinc;
guess_best_delay=0;
trydelays=guess_best_delay+settings.try_delay1:settings.delaysteps:guess_best_delay+settings.try_delay2;
tryscales=settings.try_scale1:tryinc:settings.try_scale2;
sumdiffs=nan(length(tryscales),length(trydelays));

data1_LED=dis1(alignInd,:);
data1_LED(isnan(data1_LED))=0;
data1_t=t1(alignInd,:);
data2_LED=dis2(alignInd,:);
data2_t=t2(alignInd,:);
data1_t_backup=data1_t;
% if ~isempty(alignmentAnchor_data1)
%     trydelays=0;
% end
for j=1:length(tryscales)
    if mod(j,100)==0
        disp('Processing ...');
    end
    currscale=tryscales(j);
    data1_t=data1_t_backup*currscale;
    data1_t=data1_t-nanmin(data1_t);
    data1_t=data1_t+nanmin(data1_t_backup);
    mid_backup_data1_t=data1_t;
    for k=1:length(trydelays)
        data1_t=mid_backup_data1_t;
        currdelay=trydelays(k);
        temp=data1_t(2:end)-data1_t(1:end-1);
        timestep=mode(temp(~isnan(temp)));
        if ~isempty(alignmentAnchor_data1)
            % anchor 1 and anchor 2 are same time
            temp1=data1_t(anchor1);
            temp2=data2_t(anchor2);
            data1_t=data1_t-(temp1-temp2);
        end
        if currdelay>0
            data1_t=[fliplr(nanmin(data1_t)-timestep:-timestep:nanmin(data1_t)-currdelay*timestep) data1_t];
            data1_t=data1_t(1:length(data1_LED));
        elseif currdelay<0
            data1_t=[data1_t(abs(currdelay):end) nanmax(data1_t)+timestep:timestep:nanmax(data1_t)+abs(currdelay)*timestep];
            data1_t=data1_t(1:length(data1_LED));
        end
        takeoutnan=isnan(data2_t);
        if isempty(data2_t(~takeoutnan)) || isempty(data2_LED(~takeoutnan)) || isempty(data1_t)
            error('There is a problem with this alignInd. Change it.');
        end
        y2i=interp1(data2_t(~takeoutnan),data2_LED(~takeoutnan),data1_t);
        dontuse=isnan(y2i);
%         % zero out
%         tempie_data1LED=data1_LED;
%         tempie_y2i=y2i;
%         tempie_data1LED(dontuse)=0;
%         tempie_y2i(dontuse)=0;
%         dontuse=zeros(size(dontuse));
%         [c,lags]=xcorr(data1_LED(~dontuse),y2i(~dontuse),0,'normalized');
        R=corrcoef(data1_LED(~dontuse),y2i(~dontuse));
        c=R(1,2);
        sumdiffs(j,k)=c;
    end
end
sumdiffs=1-sumdiffs;
sumdiffs(isinf(sumdiffs))=nan;
ma=max(sumdiffs(:));
sumdiffs(isnan(sumdiffs))=3*ma;
[~,mi]=min(sumdiffs(:));
[mi_row,mi_col]=ind2sub(size(sumdiffs),mi);

figure(); 
imagesc(sumdiffs);
title('Finding best alignment');
xlabel('Trying different delays');
ylabel('Trying different scales');
disp('Best delay column');
disp(mi_col);
disp('Best scale row');
disp(mi_row);

% [data1_LED,data1_t,timestep]=shiftRow(dis1(alignInd,:),t1(alignInd,:),tryscales,mi_row,trydelays,mi_col,true,[],[]);
guess_best_delay=trydelays(mi_col); 
trydelays=guess_best_delay+settings.forSearchMinus:1:guess_best_delay+settings.forSearchPlus;
mi_col=floor(length(trydelays)/2);
[data1_LED,data1_t,timestep]=shiftRow(dis1(alignInd,:),t1(alignInd,:),tryscales,mi_row,trydelays,mi_col,true,dis2(alignInd,:),t2(alignInd,:),anchor1,anchor2);
figure(); plot(data1_t,data1_LED,'Color','k'); hold on; plot(t2(alignInd,:),dis2(alignInd,:),'Color','r'); title('Alignment of test row');
pause;

% Scale/shift all rows of time structure similarly
for i=1:size(t1,1)
    if all(isnan(t1(i,:)))
        pause;
    end
    [~,data1_t]=shiftRow(dis1(i,:),t1(i,:),tryscales,mi_row,trydelays,mi_col,true,dis2(i,:),t2(i,:),anchor1,anchor2);
    t1(i,:)=data1_t;
end
data1.(whichTime1)=t1;

% Note that I have scaled and shifted the time by a set amount, so will
% apply to all fields with this same timestamp

% othertimestamp='green_time';
% Apply same scaling and shift to this othertimestamp
% delayInTime=trydelays(mi_col)*timestep; % convert indices to time
% othertimes=data1.(othertimestamp);
% dt_other=othertimes(alignInd,:);
% temp=dt_other(2:end)-dt_other(1:end-1);
% timestep2=mode(temp(~isnan(temp)));
% delayInOtherInds=floor(delayInTime/timestep2); % convert time back to indices
% for i=1:size(othertimes,1)
%     [~,data1_t]=shiftRow(othertimes(i,:),othertimes(i,:),tryscales,mi_row,delayInOtherInds,1);
%     othertimes(i,:)=data1_t;
% end
% data1.(othertimestamp)=othertimes;

% Plot final alignments
dis1=data1.(whichField1);
t1=data1.(whichTime1);
dis2=data2.(whichField2);
t2=data2.(whichTime2);
figure();
offset=1;
dis1=dis1./max(dis1(1:end),[],2,'omitnan');
dis2=dis2./max(dis2(1:end),[],2,'omitnan');
for i=1:size(data1.(whichTime1),1)
    plot(t1(i,:),offset+0.9*dis1(i,:),'Color','k');
    hold on;
    plot(t2(i,:),offset+0.9*dis2(i,:),'Color','r');
%     offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
    offset=offset+1;
end
title('Final alignment of distractor');
% dis1=data1.(othertimestamp);
% dis2=data2.(whichTime1);
% figure();
% offset=0;
% for i=1:size(data1.(othertimestamp),1)
%     plot(1:length(dis1(i,:)),offset+dis1(i,:),'Color','k');
%     hold on;
%     plot(1:length(dis2(i,:)),offset+dis2(i,:),'Color','r');
%     offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
% end
% title('Times in data1 vs. othertimes after alignment');

% Flip back if needed so output is consistent with input
switch whichToShift
    case 'data1'
        % resampled data 1
    case 'data2'
        % switched data 1 and data 2
        temp_data1=data1;
        temp_whichTime1=whichTime1;
        temp_whichField1=whichField1;
        data1=data2;
        whichTime1=whichTime2;
        whichField1=whichField2;
        data2=temp_data1;
        whichTime2=temp_whichTime1;
        whichField2=temp_whichField1;
    otherwise
        error('Unrecognized value passed in as whichToShift');
end

end

function data=makeTimeUniform(data,timefield,reset_timefield,anchorfield)

anchorthresh=0.4;
if ~isempty(anchorfield)
    anchordata=mean(data.(anchorfield),1,'omitnan');
    anchorind=find(anchordata>anchorthresh,1,'first');
    % location of this anchorind doesn't move
else
    anchorind=[];
end

timesFromArduino=data.(timefield);
data.uniformtime=nan(size(timesFromArduino));
for i=1:size(timesFromArduino,1)
    currtime=timesFromArduino(i,:);
    uniformtime=linspace(min(currtime,[],'all','omitnan'),max(currtime,[],'all','omitnan'),sum(~isnan(currtime)));
    if ~isempty(anchorind)
        curranchortime=currtime(anchorind);
        fillintime=nan(size(currtime));
        [~,mi]=nanmin(abs(uniformtime-curranchortime));
        fillintime(mi-(length(uniformtime(1:mi))-1):mi+(length(uniformtime(mi:end))-1))=uniformtime;
    else
        fillintime(~isnan(currtime))=uniformtime;
    end
    temp=data.uniformtime;
    temp(i,:)=fillintime;
    data.uniformtime=temp;
end
% resample other fields in data to match uniformtime
newdata_inds=resampleToMatchTimes(timesFromArduino,data.uniformtime);
f=fieldnames(data);
for i=1:length(f)
    temp=data.(f{i});
    if all(size(temp)==size(newdata_inds))
        newtemp=nan(1,length(temp(i,:)));
        for j=1:length(temp(i,:))
            newtemp(j)=temp(i,newdata_inds(i,j));
        end
        temp(i,:)=newtemp;
        data.(f{i})=temp;
    end
end
data.(reset_timefield)=data.uniformtime;

end

function newdata_inds=resampleToMatchTimes(datatimes,newdatatimes)

newdata_inds=nan(size(datatimes));
for i=1:size(datatimes,1)
    curr_datatimes=datatimes(i,:);
    curr_newdatatimes=newdatatimes(i,:);
    for j=1:length(curr_newdatatimes)
        [~,mi]=nanmin(abs(curr_newdatatimes(j)-curr_datatimes));
        newdata_inds(i,j)=mi;
    end
end

end

function [data1_LED,data1_t,timestep]=shiftRow(data1_LED,data1_t,tryscales,mi_row,trydelays,mi_col,getNewDelay,data2_LED,data2_t,anchor1,anchor2)

data1_t_backup=data1_t;

data1_LED(isnan(data1_LED))=0;

currscale=tryscales(mi_row);
data1_t=data1_t_backup*currscale;
data1_t=data1_t-nanmin(data1_t);
data1_t=data1_t+nanmin(data1_t_backup);
mid_backup_data1_t=data1_t;
data1_t=mid_backup_data1_t;
currdelay=trydelays(mi_col);
temp=data1_t(2:end)-data1_t(1:end-1);
timestep=mode(temp(~isnan(temp)));

if ~isempty(anchor1)
    % anchor 1 and anchor 2 are same time
    temp1=data1_t(anchor1);
    temp2=data2_t(anchor2);
    data1_t=data1_t-(temp1-temp2);
    % getNewDelay=false;
    currdelay=0;
end
if getNewDelay==true
    sumdiffs=nan(1,length(trydelays));
    for k=1:length(trydelays)
        data1_t=mid_backup_data1_t;
        currdelay=trydelays(k);
        if currdelay>0
            data1_t=[fliplr(nanmin(data1_t)-timestep:-timestep:nanmin(data1_t)-currdelay*timestep) data1_t];
            data1_t=data1_t(1:length(data1_LED));
        elseif currdelay<0
            data1_t=[data1_t(abs(currdelay):end) nanmax(data1_t)+timestep:timestep:nanmax(data1_t)+abs(currdelay)*timestep];
            if length(data1_t)<length(data1_LED)
                continue
            end
            data1_t=data1_t(1:length(data1_LED));
        end
        takeoutnan=isnan(data2_t);
        if isempty(data2_t(~takeoutnan)) || isempty(data2_LED(~takeoutnan)) || isempty(data1_t) || all(isnan(data1_t))
            continue
        end
        if length(unique(data2_t(~takeoutnan)))==1
            continue
        end
        y2i=interp1(data2_t(~takeoutnan),data2_LED(~takeoutnan),data1_t);
        dontuse=isnan(y2i);
%         [c,lags]=xcorr(data1_LED(~dontuse),y2i(~dontuse),0,'normalized');
        R=corrcoef(data1_LED(~dontuse),y2i(~dontuse));
        if isnan(R)
            c=0;
        else
            c=R(1,2);
        end
        sumdiffs(k)=c;
    end
%     sumdiffs=1-sumdiffs;
%     sumdiffs(isinf(sumdiffs))=nan;
%     ma=max(sumdiffs(:));
%     sumdiffs(isnan(sumdiffs))=3*ma;
    [~,mi]=nanmax(sumdiffs(:));
    currdelay=trydelays(mi);
    data1_t=mid_backup_data1_t;
    if currdelay>0
        data1_t=[fliplr(nanmin(data1_t)-timestep:-timestep:nanmin(data1_t)-currdelay*timestep) data1_t];
        data1_t=data1_t(1:length(data1_LED));
    elseif currdelay<0
        data1_t=[data1_t(abs(currdelay):end) nanmax(data1_t)+timestep:timestep:nanmax(data1_t)+abs(currdelay)*timestep];
        if length(data1_t)<length(data1_LED)
            data1_t=[data1_t nan(1,length(data1_LED)-length(data1_t))];
        end
        data1_t=data1_t(1:length(data1_LED));
    end
else
    if currdelay>0
        data1_t=[fliplr(nanmin(data1_t)-timestep:-timestep:nanmin(data1_t)-currdelay*timestep) data1_t];
        data1_t=data1_t(1:length(data1_LED));
    elseif currdelay<0
        data1_t=[data1_t(abs(currdelay):end) nanmax(data1_t)+timestep:timestep:nanmax(data1_t)+abs(currdelay)*timestep];
        data1_t=data1_t(1:length(data1_LED));
    end
end

end