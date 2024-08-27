function [unique_units,allcon,allopt,tims] = getOptoSuppressionWHISPER_wrapper(opto_aligned_dir,pointsBeforeOptoStart,pointsAfterOptoStart)

    % Get a list of all .mat files in the directory that begin with "phys_tbt_for_spikes"
    mat_files = dir(fullfile(opto_aligned_dir, 'phys_tbt_for_spikes*.mat'));
    
    % Sort the files alphabetically as per Windows sorting order
    sortedStr = sort_nat({mat_files.name});

    % Loop through each file
    allcon=[];
    allopt=[];
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

        frchange=[];
        unitdetails=[];
        rdet=regexp(mat_file_path,'opto_aligned');
        % For each unit, get effects of opto
        for uni=1:length(unique_units)
            [confr,optfr,supp,unitName,con,opt,tim,optoOn]=getOptoSuppressionWHISPER(data,['unit' num2str(unique_units(uni))],1);
            f=find(optoOn>0.5,1,'first');
            allcon=[allcon; con(f-pointsBeforeOptoStart:f+pointsAfterOptoStart)];
            allopt=[allopt; opt(f-pointsBeforeOptoStart:f+pointsAfterOptoStart)];
            tims=tim(f-pointsBeforeOptoStart:f+pointsAfterOptoStart);
            pause
            close all;
            frchange=[frchange [confr;optfr;supp]];
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
            aunit=load(fullfile(filesun(1).folder,filesun(1).name));
            unitdetails=[unitdetails [aunit.unitdets.isFS; aunit.unitdets.isTAN; aunit.unitdets.isSPN; aunit.unitdets.isLowFRThin]];
        end
        displayArrayWithSpaces([frchange; unitdetails]);
        pause
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

function displayArrayWithSpaces(array)
    % Check if the input is a matrix or vector
    if ~ismatrix(array)
        error('Input must be a matrix or vector.');
    end

    % Get the size of the array
    [rows, cols] = size(array);

    % Loop through each row of the array
    for i = 1:rows
        % Create a string of the row's elements separated by spaces
        rowStr = sprintf('%g ', array(i, :));
        % Remove the trailing space
        rowStr = strtrim(rowStr);
        % Display the row string
        disp(rowStr);
    end
end

