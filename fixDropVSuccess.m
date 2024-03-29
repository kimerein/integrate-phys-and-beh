function fixDropVSuccess(varargin)

% script for fixing drop vs success classifications

if isempty(varargin)
    currentVid='Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20230313\Nov_test\O2 output\2013-10-09 22-17-31-C_processed_data';
    datestr='20230313';
    mousename='Nov_test';
    f_pr=regexp(currentVid,'_processed_data');
    fslash=regexp(currentVid,'\');
    aviName=currentVid(fslash(end)+1:f_pr-1);
    % placeForO2data=['Z:\MICROSCOPE\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];
    placeForO2data=['Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];
    usePreviousChewingThresh=true;
    alignment=[];
else
    if length(varargin)==3 || length(varargin)==4
        currentVid=varargin{1};
        placeForO2data=varargin{2};
        usePreviousChewingThresh=true;
        alignment=varargin{3};
    end 
    % check for whether have already run this code, in which case use
    % existing chewing threshold
    if exist([currentVid '\fixed_tbt_success_v_drop.txt'],'file')
        disp('using previously set chewing threshold');
    else
        usePreviousChewingThresh=false;
    end

    if length(varargin)==4
        usePreviousChewingThresh=varargin{4};
    end 
end

addpath(genpath('C:\Users\sabatini\Documents\GitHub\chronux_2_11'));

load([currentVid '\tbt.mat'])
backuptbt=tbt;
if isempty(alignment)
    load([currentVid '\final_aligned_data.mat'])
end
load([placeForO2data '_zoneVals.mat'])
load([placeForO2data '_eat.mat'])
load([placeForO2data '_savehandles.mat'])
newtbt=postAnalysis_checkForChewedPellet(tbt,alignment,savehandles,zoneVals,eat,usePreviousChewingThresh);
newtbt=addReachBatchesToSingleTbt(newtbt,'cueZone_onVoff',0.25,0,[]);
tbt=newtbt;

save([currentVid '\backuptbt.mat'],'backuptbt');
save([currentVid '\tbt.mat'],'tbt');

fid=fopen([currentVid '\fixed_tbt_success_v_drop.txt'],'wt');
fclose(fid);

rmpath(genpath('C:\Users\sabatini\Documents\GitHub\chronux_2_11'));