function fixDropVSuccess

% script for fixing drop vs success classifications

currentVid='Z:\MICROSCOPE\Kim\WHISPER recs\Oct_B\20210312\O2 output\2021-03-12 16-44-31-C_processed_data';
datestr='20210312';
mousename='Oct_B';

f_pr=regexp(currentVid,'_processed_data');
fslash=regexp(currentVid,'\');
aviName=currentVid(fslash(end)+1:f_pr-1);
placeForO2data=['Z:\MICROSCOPE\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];

load([currentVid '\tbt.mat'])
backuptbt=tbt;
load([currentVid '\final_aligned_data.mat'])
load([placeForO2data '_zoneVals.mat'])
load([placeForO2data '_eat.mat'])
load([placeForO2data '_savehandles.mat'])
newtbt=postAnalysis_checkForChewedPellet(tbt,alignment,savehandles,zoneVals,eat);
newtbt=addReachBatchesToSingleTbt(newtbt,'cueZone_onVoff',0.25,0,[]);
tbt=newtbt;

save([currentVid '\backuptbt.mat'],'backuptbt');
save([currentVid '\tbt.mat'],'tbt');

fid=fopen([currentVid '\fixed_tbt_success_v_drop.txt'],'wt');
fclose(fid);