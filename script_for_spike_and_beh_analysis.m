function [spikes,tbt,beh_sweeps,cuetimes,optotimes]=script_for_spike_and_beh_analysis(filename,spikes,tbt)
% for matching spikes to behavior

matchLEDs=1; % if 1, will take consensus led from behavior and spikes

% load spikes
% load tbt

% filename='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Behavior on Electrophys Rig\str_black\20180411\str.pl2';

% convert continuous spike structure to trial-by-trial
[spikes,cuetimes,optotimes]=process_spikes_to_plot(filename,spikes);

% classify behavior sweeps
beh_sweeps=getSweepsFromBeh(tbt,'cueZone_onVoff');

% align behavior sweeps to spikes sweeps
D=alignBehToSweeps(beh_sweeps,spikes.sweeps);
% if D is positive, then spikes.sweeps is delayed with respect to beh_sweeps
% align fields of beh_sweeps to match spikes.sweeps
if D>0
    f=fieldnames(beh_sweeps);
    for i=1:length(f)
        temp=beh_sweeps.(f{i});
        temp=padVector(temp,D);
        beh_sweeps.(f{i})=temp;
    end
    % also pad tbt to match spikes
    f=fieldnames(tbt);
    for i=1:length(f)
        tbt.(f{i})=[nan(D,size(tbt.(f{i}),2)); tbt.(f{i})];
    end
elseif D<0
    D=-D;
    % then need to delay spikes.sweeps to match delayed beh_sweeps
    % pad spikes to match beh_sweeps
    f=fieldnames(spikes.sweeps);
    for i=1:length(f)
        temp=spikes.sweeps.(f{i});
        temp=padVector(temp,D);
        spikes.sweeps.(f{i})=temp;
    end
    % now adjust trials in spikes
    spikes.sweeps.trials=spikes.sweeps.trials+D;
    spikes.sweeps.trials(1:D)=1:D;
    spikes.trials=spikes.trials+D;
end
% pad at ends if different lengths
max_spikes=length(spikes.sweeps.trials);
max_beh=length(beh_sweeps.led);
if max_spikes<max_beh
    % pad end of spikes
    f=fieldnames(spikes.sweeps);
    for i=1:length(f)
        temp=spikes.sweeps.(f{i});
        temp=endPadVector(temp,max_beh-max_spikes);
        spikes.sweeps.(f{i})=temp;
    end
    spikes.sweeps.trials(end-(max_beh-max_spikes)+1:end)=length(spikes.sweeps.trials)-(max_beh-max_spikes)+1:length(spikes.sweeps.trials);
elseif max_beh<max_spikes
    % pad end of beh
    f=fieldnames(beh_sweeps);
    for i=1:length(f)
        temp=beh_sweeps.(f{i});
        temp=endPadVector(temp,max_spikes-max_beh);
        beh_sweeps.(f{i})=temp;
    end
    % also pad tbt to match spikes
    f=fieldnames(tbt);
    for i=1:length(f)
        tbt.(f{i})=[tbt.(f{i}); nan(max_spikes-max_beh,size(tbt.(f{i}),2))];
    end
end

% Find consensus led for behavior and spikes
if matchLEDs==1
    % spikes is more reliable
    best_led=spikes.sweeps.led;
    beh_sweeps.led=best_led;
    best_led=best_led';
    beh_sweeps.led_previousTrial=[nan; best_led(1:end-1)];
    beh_sweeps.led_2back=[nan; nan; best_led(1:end-2)];
    beh_sweeps.led_3back=[nan; nan; nan; best_led(1:end-3)];
    beh_sweeps.led_4back=[nan; nan; nan; nan; best_led(1:end-4)];
end

end

function temp=padVector(temp,D)

if size(temp,1)>1
    % column vector
    temp=[nan(D,1); temp];
elseif size(temp,2)>1
    % row vector
    temp=[nan(1,D) temp];
else
    error('should be vector, not array');
end

end

function temp=endPadVector(temp,D)

if size(temp,1)>1
    % column vector
    temp=[temp; nan(D,1)];
elseif size(temp,2)>1
    % row vector
    temp=[temp nan(1,D)];
else
    error('should be vector, not array');
end

end
