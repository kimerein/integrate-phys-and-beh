function plotRecSiteLocs()

% rows are different recording sites
% column 4 is which mouse
% column 5 is distance (A-P) in microns from dye-marked recording site tip
% -- negative is more posterior
% column 6 is distance (D-V) in microns from dye-marked recording site tip
% -- positive is deeper
% column 7 is distance (M-L) in microns from dye-marked recording site tip
% -- positive is more lateral
% column 8 is A-P relative to bregma (negative is more posterior) of dye
% track tip
% column 9 is D-V relative to bregma (positive is deeper) of dye
% track tip
% column 10 is M-L relative to bregma (positive is more lateral) of dye
% track tip
% column 11 is flag (true) if track missed striatum "not_in_striatum"
% column 12 is probe_tip passed to process_units
% column 13 is top_of_structure passed to process_units
% column 14 is bottom_of_structure passed to process_units
% e.g., data_loc_array{i,11}=2600;
% make a colored point on grid for any channel in structure

% ultimately, ingest this from .csv file
dataTable='C:\Users\sabatini\Downloads\rec_site_locs.csv';
data_loc_array=table2cell(readtable(dataTable,'Format','%s%s%s%s%u%u%u%u%u%u%s%u%u%u'));

[a,b,c]=unique(data_loc_array{:,4});
whichMouse=c;

% First, calculate probe track angle for Kim's Floor 1 WHISPER rig striatum physiology data
errorAllowance=0; % in microns, expand structure to allow for error in anatomy
data_loc_array(:,13)=data_loc_array(:,13)-errorAllowance;
data_loc_array(:,14)=data_loc_array(:,14)+errorAllowance;
chSpacing=20; % in microns
entry_X=-2277;
final_X=-3857;
entry_Y=7902;
final_Y=7902;
entry_Z=-330;
final_Z=2300;
theta=17.75; % in degrees, A-P angle off midline pointing straight posterior
alpha=90-theta; % in degrees
phi=asind(abs(entry_X-final_X)/abs(entry_Z-final_Z)); % in degrees, angle of probe moving at diagonal off straight down (pointed ventral)

% site 32 is most dorsall, 1 is most ventral
% Calculate all site locations
chOffsets=([1:32]-1)*chSpacing; % ventral tip minus this
% account for angle of probe (not straight down)
chOffsets=chOffsets*cosd(phi);
chOffsetsAP=chOffsets*sind(phi)*cosd(theta);
chOffsetsML=chOffsets*sind(phi)*sind(theta);
% Get A-P positions of all sites
% already accounted for theta in google sheet
siteAPpositions=data_loc_array{:,8}-data_loc_array{:,5};
siteMLpositions=data_loc_array{:,10}+data_loc_array{:,7};
siteDVpositions=data_loc_array{:,9}+data_loc_array{:,6};
didnotuse=data_loc_array{:,6}==1; 
allSites=[siteAPpositions siteMLpositions siteDVpositions whichMouse didnotuse data_loc_array{:,12} data_loc_array{:,13} data_loc_array{:,14}]; % fifth column is whether used this recording

allChs=[];
for j=1:size(allSites,1)
    for i=1:length(chOffsets)
        currZoffset=chOffsets(i);
        currMLoffset=chOffsetsML(i);
        currAPoffset=chOffsetsAP(i);
        temp=allSites(j,:);
        temp(3)=temp(3)-currZoffset;
        temp(2)=temp(2)-currMLoffset;
        temp(1)=temp(1)-currAPoffset;
        threwout=temp(5)==1;
        temp(5)=threwout || temp(3)<temp(7) || temp(3)>temp(8);
        allChs=[allChs; temp];
    end
end

% Only plot sites in str
stepsForFigs=1:-0.12:-3; % A-P wrt bregma, in mm
stepsForFigs=stepsForFigs*1000; % convert to microns
for i=1:length(stepsForFigs)
    currAPstep=stepsForFigs(i);
    % plot all channels that were not thrown out based on anatomy that are
    % in this slice
    whichToPlot=allChs(:,)
    figure();
    grid on
    
end




end