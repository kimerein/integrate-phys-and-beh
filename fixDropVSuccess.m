function fixDropVSuccess(varargin)

% script for fixing drop vs success classifications

if isempty(varargin)
    currentVid='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210312\O2 output\2021-03-12 16-44-31-C_processed_data';
    datestr='20210312';
    mousename='Oct_B';
    f_pr=regexp(currentVid,'_processed_data');
    fslash=regexp(currentVid,'\');
    aviName=currentVid(fslash(end)+1:f_pr-1);
    placeForO2data=['Z:\MICROSCOPE\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];
    usePreviousChewingThresh=false;
else
    if length(varargin)==2
        currentVid=varargin{1};
        placeForO2data=varargin{2};
        usePreviousChewingThresh=true;
    end 
    % check for whether have already run this code, in which case use
    % existing chewing threshold
    if exist([currentVid '\fixed_tbt_success_v_drop.txt'],'file')
        disp('using previously set chewing threshold');
    else
        usePreviousChewingThresh=false;
    end
end

load([currentVid '\tbt.mat'])
backuptbt=tbt;
load([currentVid '\final_aligned_data.mat'])
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