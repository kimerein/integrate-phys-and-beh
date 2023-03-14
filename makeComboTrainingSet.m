function makeComboTrainingSet(whichToCombine,whichToCombineFileNames,saveTo,saveToFileName,currDir)

ls=dir([currDir sep whichToCombine{1}]);
for j=3:length(ls)
    if contains(ls(j).name,'trainingSet')
        % Get this training set
        a=load([currDir sep whichToCombine{1} sep ls(j).name]);
        trainingSetTrials=a.trainingSetTrials;
        unitbasename=ls(j).name(1:regexp(ls(j).name,['_' whichToCombineFileNames{1} '_trainingSet'],'once')-1);
        disp(['grabbing ' currDir sep whichToCombine{1} sep unitbasename]);
        % Find this unit in other whichToCombine folders
        for i=2:length(whichToCombine)
            a=load([currDir sep whichToCombine{i} sep unitbasename '_' whichToCombineFileNames{i} '_trainingSet.mat']);
            trainingSetTrials=[trainingSetTrials; a.trainingSetTrials];
        end
        save([currDir sep saveTo sep unitbasename '_' saveToFileName '_trainingSet.mat'],'trainingSetTrials');
    end
end
