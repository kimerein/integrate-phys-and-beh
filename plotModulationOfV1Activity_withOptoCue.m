function [allexpts_allmod_rawmod,allexpts_allchs,allexpts_allpsth_N,allbincenters]=plotModulationOfV1Activity_withOptoCue(directoryPaths,whichToRead)

if ~iscell(directoryPaths)
    temp=directoryPaths;
    clear directoryPaths
    directoryPaths{1}=temp;
end

allexpts_allpsth_N=[];
allexpts_allmod_baserate=[];
allexpts_allmod_afterrate=[];
allexpts_allmod_rawmod=[];
allexpts_allmod_fracmod=[];
allexpts_allchs=[];
for di=1:length(directoryPaths)
    disp(di);
    directoryPath=directoryPaths{di};
    % Get a list of all .mat files in the directory that contain whichToRead in their name
    files = dir(fullfile(directoryPath, ['*', whichToRead, '*.mat']));

    % Get n's
    a=load(fullfile(directoryPath,'cue_n.mat'));
    cue_n=a.n;
    a=load(fullfile(directoryPath,'cuePlusOpto_n.mat'));
    cuePlusOptoData_n=a.n;
    a=load(fullfile(directoryPath,'cueWithoutOpto_n.mat'));
    cueWithoutOptoData_n=a.n;
    a=load(fullfile(directoryPath,'distractor_n.mat'));
    distractor_n=a.n;
    a=load(fullfile(directoryPath,'opto_n.mat'));
    opto_n=a.n;

    switch whichToRead
        case 'cueAligned'
            n=cue_n;
        case 'distractorAligned'
            n=distractor_n;
        case 'cuePlusOptoData'
            n=cuePlusOptoData_n;
        case 'cueWithoutOptoData'
            n=cueWithoutOptoData_n;
        case 'optoAligned'
            n=opto_n;
    end

    % Load each file
    % Save modulation, save PSTH
    for i = 1:length(files)
        a = load(fullfile(directoryPath, files(i).name));
        if i==1
            allpsth_N=nan(length(files),length(a.alignData.N));
            allbincenters=a.alignData.bincenters;
            allmod_baserate=nan(length(files),1);
            allmod_afterrate=nan(length(files),1);
            allmod_rawmod=nan(length(files),1);
            allmod_fracmod=nan(length(files),1);
            allchs=nan(length(files),1);
        end
        allpsth_N(i,:)=a.alignData.N;
        allmod_baserate(i)=a.alignData.modulation.baserate;
        allmod_afterrate(i)=a.alignData.modulation.afterrate;
        allmod_rawmod(i)=a.alignData.modulation.rawmod;
        allmod_fracmod(i)=a.alignData.modulation.fracmod;
        allchs(i)=str2double(regexp(files(i).name, '(?<=Ch)\d+', 'match', 'once'));
    end

    % Adjust channels to match whole V1 mapping
    a=load(fullfile(directoryPath,'chAlign.mat'));
    addToChs=chOffset(a.alignTo,a.chAlign);
    allchs=allchs+addToChs;

    % Switch to firing rates
    timestep=mode(diff(allbincenters));
    allpsth_N=(allpsth_N./n)./timestep;
    allmod_baserate=(allmod_baserate./n)./timestep;
    allmod_afterrate=(allmod_afterrate./n)./timestep;
    allmod_rawmod=(allmod_rawmod./n)./timestep;
    allmod_fracmod=(allmod_fracmod./n)./timestep;

    % Concatenate to all experiments
    allexpts_allpsth_N=[allexpts_allpsth_N; allpsth_N];
    allexpts_allmod_baserate=[allexpts_allmod_baserate; allmod_baserate];
    allexpts_allmod_afterrate=[allexpts_allmod_afterrate; allmod_afterrate];
    allexpts_allmod_rawmod=[allexpts_allmod_rawmod; allmod_rawmod];
    allexpts_allmod_fracmod=[allexpts_allmod_fracmod; allmod_fracmod];
    allexpts_allchs=[allexpts_allchs; allchs];
end

% Put any very superficial channels at top
allexpts_allchs(allexpts_allchs>40)=40;
% Add random offset to channels just to spread out dots
randforplot=0.5*rand(size(allexpts_allchs))-0.5;

% Fix infinities and too big
allexpts_allmod_fracmod(isinf(allexpts_allmod_fracmod))=15;
allexpts_allmod_fracmod(allexpts_allmod_fracmod>15)=15;

% Plot average PSTH and se
figure();
plot(allbincenters,mean(allexpts_allpsth_N,1,'omitnan'),'Color','k');
hold on;
plot(allbincenters,mean(allexpts_allpsth_N,1,'omitnan')-std(allexpts_allpsth_N,[],1,'omitnan')./sqrt(length(allexpts_allmod_baserate)),'Color','k');
plot(allbincenters,mean(allexpts_allpsth_N,1,'omitnan')+std(allexpts_allpsth_N,[],1,'omitnan')./sqrt(length(allexpts_allmod_baserate)),'Color','k');

% Plot modulation as a function of unit depth
figure();
subplot(1,4,1);
scatter(allexpts_allmod_rawmod,randforplot+allexpts_allchs); xlabel('Raw modulation'); ylabel('Channel, 32 dorsal, 1 ventral'); 
temp=nan(40,1);
for i=1:40
    temp(i)=mean(allexpts_allmod_rawmod(allexpts_allchs==i),'all','omitnan');
end
hold on; plot(smooth(temp,10),1:40,'Color','k');

subplot(1,4,2);
scatter(allexpts_allmod_fracmod,randforplot+allexpts_allchs); xlabel('Frac modulation'); ylabel('Channel, 32 dorsal, 1 ventral'); 
temp=nan(40,1);
for i=1:40
    temp(i)=mean(allexpts_allmod_fracmod(allexpts_allchs==i),'all','omitnan');
end
hold on; plot(smooth(temp,10),1:40,'Color','k');

subplot(1,4,3);
scatter(allexpts_allmod_baserate,randforplot+allexpts_allchs); xlabel('Baseline firing rate'); ylabel('Channel, 32 dorsal, 1 ventral'); 
temp=nan(40,1);
for i=1:40
    temp(i)=mean(allexpts_allmod_baserate(allexpts_allchs==i),'all','omitnan');
end
hold on; plot(smooth(temp,10),1:40,'Color','k');

subplot(1,4,4);
scatter(allexpts_allmod_afterrate,randforplot+allexpts_allchs); xlabel('Evoked firing rate'); ylabel('Channel, 32 dorsal, 1 ventral'); 
temp=nan(40,1);
for i=1:40
    temp(i)=mean(allexpts_allmod_afterrate(allexpts_allchs==i),'all','omitnan');
end
hold on; plot(smooth(temp,10),1:40,'Color','k');

end

function differ=chOffset(alignTo,chAlign)

% Map all to a scheme where #40 is the most superficial channel in cortex
% and #1 is the deepest channel in cortex, given that probe
% has 20 micron spacing for all recordings

% If alignTo='dorsal'
% Then chAlign should be 40
% So get difference between chAlign and 40, add this to all chs

% If alignTo='ventral'
% Then chAlign should be 1
% So get chAlign-1, and add this to all chs

switch alignTo
    case 'dorsal'
        differ=40-chAlign;
    case 'ventral'
        differ=chAlign-1;
end

end
