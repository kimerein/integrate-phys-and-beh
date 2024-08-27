function fixOptoAligned(physDir,cueDir,opto_aligned_dir,optoDur)

% assumes opto came on at start of cue and was always aligned to cue
% opto starts 5 ms before cue

a=load([physDir sep 'physiology_tbt.mat']);
physiology_tbt=a.physiology_tbt;
%figure(); plot(nanmean(physiology_tbt.cuetimes_wrt_trial_start,1),nanmean(physiology_tbt.cue,1));
% find cue start
x=nanmean(physiology_tbt.cuetimes_wrt_trial_start,1);
y=nanmean(physiology_tbt.cue,1);
f=find(y>0.001,1,'first');
optotime=x(f)-0.005; % in seconds

% Get a list of all .mat files in the directory that begin with "phys_tbt_for_spikes"
mat_files = dir(fullfile(opto_aligned_dir, 'phys_tbt_for_spikes*.mat'));

% Sort the files alphabetically as per Windows sorting order
sortedStr = sort_nat({mat_files.name});

% Loop through each file
for k = 1:length(sortedStr)
    % Load the .mat file
    mat_file_path = fullfile(opto_aligned_dir, sortedStr{k});
    data = load(mat_file_path);
    data=data.optoAligned_phys_tbt;

    % Get the field names of the loaded structure
    field_names = fieldnames(data);

    % Loop through each field name
    unique_units=[];
    for i = 1:length(field_names)
        % Check if the field name starts with "unit" followed by an integer
        field_name = field_names{i};
        if startsWith(field_name, 'unit')
            % Extract the number immediately following "unit"
            num_str = regexp(field_name, '^unit(\d+)', 'tokens', 'once');
            if ~isempty(num_str)
                % Convert the extracted string to a number
                unit_num = str2double(num_str{1});

                % Add the number to the unique_units array if not already present
                if ~ismember(unit_num, unique_units)
                    unique_units = [unique_units, unit_num];
                end
            end
        end
    end

    % Sort the unique unit numbers
    unique_units = sort(unique_units);

    % Display the unique unit numbers
    disp('Unique unit numbers found:');
    disp(unique_units);

    % Fix opto alignment
    rdet=regexp(mat_file_path,'opto_aligned');
    for uni=1:length(unique_units)
        d=[mat_file_path(1:rdet-1) 'unit_details'];
        s=['unit' num2str(unique_units(uni)) 'on'];
        filesun = dir(fullfile(d, [s '*.mat']));
        if length(filesun)>1
            % ask user to disambiguate
            for sname=1:length(filesun)
                disp(filesun(sname).name);
            end
            ui=input([s ' on Ch ?']);
            filesun = dir(fullfile(d, [s 'Ch' num2str(ui) '*.mat']));
        end
        nurm=filesun(1).name;
        rarg=regexp(nurm,'_');
        subname=nurm(1:rarg-1);

        a=load(mat_file_path);
        optoAligned_phys_tbt=a.optoAligned_phys_tbt;

        a=load([cueDir sep subname '__cueAligned.mat']);
        dataout=a.dataout;
        % Replace unit data
        optoAligned_phys_tbt.(['unit' num2str(unique_units(uni))])=dataout.y;
        optoAligned_phys_tbt.unitTimes_wrt_trial_start=dataout.x;
        % opto starts at optotime
        %figure(); plot(nanmean(physiology_tbt.cuetimes_wrt_trial_start,1),nanmean(physiology_tbt.cue,1));
        f=find(dataout.x>optotime,1,'first');
        optoDur_inds=floor(optoDur/mode(diff(nanmean(dataout.x,1))));
        optoAligned_phys_tbt.optoOnInUnitTimes=zeros(1,size(dataout.y,2));
        optoAligned_phys_tbt.optoOnInUnitTimes(f:f+optoDur_inds)=1;
        optoAligned_phys_tbt.(['unit' num2str(unique_units(uni)) '_avAlignedToOpto'])=nanmean(dataout.y(optoAligned_phys_tbt.hasOpto==1,:),1);
        save(mat_file_path,'optoAligned_phys_tbt');
    end
end

end
    

function sortedStr = sort_nat(cellStr)
    % SORT_NAT Natural order sort of strings.
    % This function sorts strings containing numbers in a way that 
    % aligns with how a human would naturally sort them.

    [isPresent,si] = ismember('phys_tbt_for_spikes_sorted.mat', cellStr);
    currs=1:length(cellStr);
    cellStr=cellStr(currs(~ismember(currs,si)));

    % Split the strings into components that are numeric and non-numeric
    splitStr = regexp(cellStr, '(\d+)', 'split');
    numStr = regexp(cellStr, '\d+', 'match');
    allNums=[];
    for i=1:length(numStr)
        allNums=[allNums; str2double(numStr{i})];
    end
    [~,si]=sort(allNums,'ascend');

    % Put this one first if exists
    if isPresent
        sortedStr{1}='phys_tbt_for_spikes_sorted.mat';
    else
        sortedStr=[];
    end
    cellStr=cellStr(si);
    for i=1:length(cellStr)
        sortedStr{length(sortedStr)+1}=cellStr{i};
    end

end

