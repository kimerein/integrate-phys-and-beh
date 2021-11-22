function [data1,data2]=useDistractorAlignment(data1,whichTime1,whichField1,data2,whichTime2,whichField2,whichToShift)

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
settings.try_delay1=-50;
settings.delaysteps=1;
settings.try_delay2=50;
settings.tryinc=0.0005; 
settings.try_scale1=0.8;
settings.try_scale2=1.2;  

switch whichToShift
    case 'data1'
        % resample data 1
    case 'data2'
        % switch data 1 and data 2
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

% Plot starting alignment
dis1=data1.(whichField1);
t1=data1.(whichTime1);
dis2=data2.(whichField2);
t2=data2.(whichTime2);
figure();
offset=0;
for i=1:size(data1.(whichTime1),1)
    plot(t1(i,:),offset+dis1(i,:),'Color','k');
    hold on;
    plot(t2(i,:),offset+dis2(i,:),'Color','r');
    offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
end

% Use distractor to further align
% Resample and shift data1 times
disp('Now refining alignment ...');
tryinc=settings.tryinc;
guess_best_delay=0;
trydelays=guess_best_delay+settings.try_delay1:settings.delaysteps:guess_best_delay+settings.try_delay2;
tryscales=settings.try_scale1:tryinc:settings.try_scale2;
sumdiffs=nan(length(tryscales),length(trydelays));

alignInd=1;
data1_LED=dis1(alignInd,:);
data1_LED(isnan(data1_LED))=0;
data1_t=t1(alignInd,:);
data2_LED=dis2(alignInd,:);
data2_t=t2(alignInd,:);
data1_t_backup=data1_t;
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
        if currdelay>0
            data1_t=[fliplr(nanmin(data1_t)-timestep:-timestep:nanmin(data1_t)-currdelay*timestep) data1_t];
            data1_t=data1_t(1:length(data1_LED));
        elseif currdelay<0
            data1_t=[data1_t(abs(currdelay):end) nanmax(data1_t)+timestep:timestep:nanmax(data1_t)+abs(currdelay)*timestep];
            data1_t=data1_t(1:length(data1_LED));
        end
        takeoutnan=isnan(data2_t);
        y2i=interp1(data2_t(~takeoutnan),data2_LED(~takeoutnan),data1_t);
        dontuse=isnan(y2i);
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

[data1_LED,data1_t,timestep]=shiftRow(dis1(alignInd,:),t1(alignInd,:),tryscales,mi_row,trydelays,mi_col,false,[],[]);
figure(); plot(data1_t,data1_LED,'Color','k'); hold on; plot(t2(alignInd,:),dis2(alignInd,:),'Color','r'); title('Alignment of test row');

% Scale/shift all rows of time structure similarly
for i=1:size(t1,1)
    [~,data1_t]=shiftRow(dis1(i,:),t1(i,:),tryscales,mi_row,trydelays,mi_col,true,dis2(i,:),t2(i,:));
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
offset=0;
for i=1:size(data1.(whichTime1),1)
    plot(t1(i,:),offset+dis1(i,:),'Color','k');
    hold on;
    plot(t2(i,:),offset+dis2(i,:),'Color','r');
    offset=offset+nanmax([dis1(i,:) dis2(i,:)])+0.1;
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

function [data1_LED,data1_t,timestep]=shiftRow(data1_LED,data1_t,tryscales,mi_row,trydelays,mi_col,getNewDelay,data2_LED,data2_t)

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
            data1_t=data1_t(1:length(data1_LED));
        end
        takeoutnan=isnan(data2_t);
        y2i=interp1(data2_t(~takeoutnan),data2_LED(~takeoutnan),data1_t);
        dontuse=isnan(y2i);
%         [c,lags]=xcorr(data1_LED(~dontuse),y2i(~dontuse),0,'normalized');
        R=corrcoef(data1_LED(~dontuse),y2i(~dontuse));
        c=R(1,2);
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