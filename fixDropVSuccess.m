function fixDropVSuccess

% script for fixing drop vs success classifications

currentVid='Z:\Kim\WHISPER recs\Jan_1\20201202\O2 output\2020-12-02 16-03-28-C_processed_data';
datestr='20201202';
mousename='Jan_1';

f_pr=regexp(currentVid,'_processed_data');
fslash=regexp(currentVid,'\');
aviName=currentVid(fslash(end)+1:f_pr-1);
placeForO2data=['Z:\Kim\WHISPER recs\' mousename '\' datestr '\O2 output\' aviName];

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