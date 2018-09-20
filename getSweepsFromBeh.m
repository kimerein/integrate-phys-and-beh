function out=getSweepsFromBeh(tbt)

% Get user-defined settings from this file
settings=RTanalysis_settings();
lowThresh=settings.lowThresh;
ttsettings=trialTypeSettings();

% Fix hold lacking expts
if ~isfield(tbt,'isHold')
    tbt.isHold=zeros(size(tbt.(ttsettings.nameOfCue)));
end

% Get opto
out=getOptoTrials(tbt,ttsettings.nameOfCue,lowThresh);

% Classify trial types
tt=classifyTrialTypes(tbt);
for i=1:length(tt.trialtype)
    out.(tt.trialtype(i).name)=tt.trialtype(i).isThisType;
end

% Get states of previous trials (up to 4 back)
f=fieldnames(out);
for i=1:length(f)
    if ~(isnumeric(out.(f{i})) || islogical(out.(f{i})))
        continue
    end
    temp=out.(f{i});
    % 1 back
    newfieldname=[f{i} '_1back'];
    out.(newfieldname)=[nan; temp(1:end-1)];
    % 2 back
    newfieldname=[f{i} '_2back'];
    out.(newfieldname)=[nan; nan; temp(1:end-2)];
    % 3 back
    newfieldname=[f{i} '_3back'];
    out.(newfieldname)=[nan; nan; nan; temp(1:end-3)];
    % 4 back
    newfieldname=[f{i} '_4back'];
    out.(newfieldname)=[nan; nan; nan; nan; temp(1:end-4)];
    
    % 1 forward
    newfieldname=[f{i} '_1forward'];
    out.(newfieldname)=[temp(2:end); nan];
    % 2 forward
    newfieldname=[f{i} '_2forward'];
    out.(newfieldname)=[temp(3:end); nan; nan];
    % 3 forward
    newfieldname=[f{i} '_3forward'];
    out.(newfieldname)=[temp(4:end); nan; nan; nan];
    % 4 forward
    newfieldname=[f{i} '_4forward'];
    out.(newfieldname)=[temp(5:end); nan; nan; nan; nan];
end

end

function out=getOptoTrials(tbt,cueName,lowThresh)

% Get trial durations from tbt
figure();
plot(nanmean(tbt.(cueName),1));
xlabel('indices');
ylabel('av');
title('average of movie cue across trials');
ind=input('Enter index showing minimum trial length -- note that any opto within a trial should occur before this index. ');
if ~isnumeric(ind)
    error('Please enter a number.');
end

% Get trials with opto
led=any(tbt.optoOn(:,1:ind)>lowThresh,2);
led=single(led);

% Save to structure for output
out.led=led;

end
