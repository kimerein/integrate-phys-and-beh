function trialtypes=classifyTrialTypes(tbt,led)

% User-defined settings in trialTypeSettings.m
settings=trialTypeSettings();

% Find and convert reach batches



for i=1:length(settings.bool_test)
    curr_test=settings.bool_test(i);
    switch curr_test.testwhat 
        case 'single reach'
        case 'reach batch'
        otherwise
    end
end
    
    bool_test(1).testwhat='single reach';
bool_test(1).fieldname='pawOnWheel';
bool_test(1).test='any';
bool_test(1).thresh=0.5;
bool_test(1).comparator='>';
bool_test(1).window='wheel_turning';
    


% Exclude trials where paw was on wheel while wheel turning
if excludePawOnWheelTrials==1
    % Find trials where paw was on wheel while wheel turning
    plot_cues=[];
    for i=1:size(tbt.(nameOfCue),1)
        presentInd=find(tbt.pelletPresented(i,:)>0.5,1,'first');
        temp=tbt.(nameOfCue);
        cueInd=find(temp(i,:)>0.5,1,'first');
        pawWasOnWheel=0;
        if any(tbt.pawOnWheel(i,presentInd:cueInd)>0.5)
            pawWasOnWheel=1;
        else
            plot_cues=[plot_cues i];
        end
    end
else
    plot_cues=1:size(tbt.(nameOfCue),1);
end
if settings.excludeFirstTrial==1
    plot_cues=plot_cues(~ismember(plot_cues,1));
end




end

function tbt=findReachBatches(tbt)

% Finds reach batches according to definition in trialTypeSettings.m
settings=trialTypeSettings();
reach_batch=settings.reach_batch;

% Get times of various reach types
secondtype_f=cell(1,length(reach_batch.secondreach_type));
for i=1:length(reach_batch.secondreach_type)
    temp=tbt.(reach_batch.secondreach_type{i});
    tbtsecondtype=temp(1:end);
    f=find(tbtsecondtype>0.5);
    secondtype_f{i}=f;
end

for i=1:length(reach_batch.firstreach_type)
    currtype=reach_batch.firstreach_type{i};
    % Find all reaches of this type
    % For each reach of this type, check whether a second reach of an appropriate type occurs
    % within window seconds
    temp=tbt.(currtype);
    tbtcurrtype=temp(1:end);
    f=find(tbtcurrtype>0.5);
    for j=1:length(f)-1
        [row,col]=ind2sub(size(tbt.(currtype)),f(j));
        candidate_secondreaches.inds=[];
        candidate_secondreaches.types=[];
        for k=1:length(reach_batch.secondreach_type)
            f2=secondtype_f{k};
            ne=find(f2>f(j),1,'first');
            % is ne in the same trial and within window?
            [row2,col2]=ind2sub(size(tbt.(currtype)),ne);
            if row==row2 % reaches occur in same trial
                % do reaches occur within time window defining reach batch?
                timediff=tbt.times(row2,col2)-tbt.times(row,col);
                if timediff<reach_batch.window
                    % candidate second reach
                    candidate_secondreaches.inds=[candidate_secondreaches.inds ne];
                    candidate_secondreaches.types=[candidate_secondreaches.types k];
                end
            end
        end
        % find reach immediately following first reach
        if ~isempty(candidate_secondreaches.inds)
            % found a reach qualifying for reach batch
            [secondreach.ind,secondreach.type]=min(candidate_secondreaches.inds);
        else
            % no reaches in reach batch
        end
        
            
            
            
            
        [row_next,col_next]=ind2sub(size(tbt.(currtype)),f(j+1));
        % To be in same batch, reaches 


reach_batch.window=0.3; % in seconds, reaches must occur within this many seconds of each other to be in same batch
reach_batch.firstreach_type={'miss_reachStarts'}; % first reach must be one of these types
reach_batch.secondreach_type={'success_reachStarts_pawOnWheel','drop_reachStarts_pawOnWheel','miss_reachStarts_pawOnWheel'}; % second reach must be one of these types
reach_batch.take_first_or_second_type=2; % number indicates whether to convert batch to the type of the first or second reach


end