function [physiology_tbt,beh2_tbt,behavior_tbt,photometry_tbt]=loadReachingExptPhysData(data_loc_array)

% a row of data_loc_array structured thusly
% Date, mouse, raw data physiology file name and directory, probe depth,
% recording location, location of trial by trial data, spikes location,
% location of SU aligned to behavior, raw data binary name

physiology_tbt=[];
beh2_tbt=[];
behavior_tbt=[];
photometry_tbt=[];

if exist([data_loc_array{6} sep 'physiology_tbt.mat'],'file')
    a=load([data_loc_array{6} sep 'physiology_tbt.mat']);
    physiology_tbt=a.physiology_tbt;
end
if exist([data_loc_array{6} sep 'beh2_tbt.mat'],'file')
    a=load([data_loc_array{6} sep 'beh2_tbt.mat']);
    beh2_tbt=a.beh2_tbt;
end
if exist([data_loc_array{6} sep 'behavior_tbt.mat'],'file')
    a=load([data_loc_array{6} sep 'behavior_tbt.mat']);
    behavior_tbt=a.behavior_tbt;
end
if exist([data_loc_array{6} sep 'photometry_tbt.mat'],'file')
    a=load([data_loc_array{6} sep 'photometry_tbt.mat']);
    photometry_tbt=a.photometry_tbt;
end

end

function separ=sep()

if ismac
    separ='/';
else
    separ='\';
end

end