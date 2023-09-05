function plotModulationOfV1Activity_withOptoCue(directoryPath,whichToRead)

% Get a list of all .mat files in the directory that contain whichToRead in their name
files = dir(fullfile(directoryPath, ['*', whichToRead, '*.mat']));

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

% Plot average PSTH and se
figure();
plot(allbincenters,mean(allpsth_N,1,'omitnan'),'Color','k');
hold on;
plot(allbincenters,mean(allpsth_N,1,'omitnan')-std(allpsth_N,[],1,'omitnan')./sqrt(length(files)),'Color','k');
plot(allbincenters,mean(allpsth_N,1,'omitnan')+std(allpsth_N,[],1,'omitnan')./sqrt(length(files)),'Color','k');

% Plot modulation as a function of unit depth
figure();
scatter(allmod_rawmod,allchs); xlabel('Modulation'); ylabel('Channel, 32 dorsal, 1 ventral'); title('Raw modulation');

figure();
scatter(allmod_fracmod,allchs); xlabel('Modulation'); ylabel('Channel, 32 dorsal, 1 ventral'); title('Frac modulation');

end
