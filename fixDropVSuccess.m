function fixDropVSuccess

% script for fixing drop vs success classifications

currentVid='Z:\Kim\Behavior Final Data Sets\Learning curves w interleaved silencing\March_A\2011-05-24 03-05-15-C_processed_data';
datestr='20190729';
mousename='March_A';

f_pr=regexp(currentVid,'_processed_data');
fslash=regexp(currentVid,'\');
aviName=currentVid(fslash(end)+1:f_pr-1);
placeForO2data=['Z:\Kim\KER Behavior\By date\Low speed\' datestr '\' mousename '\O2 output\' aviName];

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