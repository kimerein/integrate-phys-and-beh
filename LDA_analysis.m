function LDA_analysis(whichSess,downSampBy,takeNPointsAfterEvent,response_to_plot1,response_to_plot2,dd)
% wrapper for code from B

% get all units recorded simultaneously from one session
% whichSess=3;
% takeNPointsAfterEvent=100;

whichUnitsToGrab='_'; plotUnitCriteria=[-100 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
setForUn=settingsForStriatumUnitPlots;
if setForUn.keepAllSingleTrials~=true
    error('need trial by trial data for LDA analysis');
end
response_to_plot=response_to_plot1; %'cued_success'; 
ResponseCued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
ResponseCued.unitbyunit_x=downSampMatrix(ResponseCued.unitbyunit_x,downSampBy);
ResponseCued.unitbyunit_y=downSampMatrix(ResponseCued.unitbyunit_y,downSampBy);
response_to_plot=response_to_plot2; %'uncued_success'; 
ResponseUncued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
ResponseUncued.unitbyunit_x=downSampMatrix(ResponseUncued.unitbyunit_x,downSampBy);
ResponseUncued.unitbyunit_y=downSampMatrix(ResponseUncued.unitbyunit_y,downSampBy);

temp=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',ResponseCued,1); timesAlignedToEvent=temp.t; [~,zeroind]=nanmin(abs(timesAlignedToEvent));
us=unique(ResponseCued.fromWhichUnit); allTri=unique(ResponseCued.fromWhichTrial);
mergedData=nan(length(allTri),length(us)*takeNPointsAfterEvent);
for i=1:length(allTri)
    currTri=allTri(i);
    for j=1:length(us)
        currU=us(j);
        if sum(ResponseCued.fromWhichTrial==currTri & ResponseCued.fromWhichUnit==currU,'all','omitnan')==0
            % unit was not active on this trial, fill w nans
            mergedData(i,(j-1)*takeNPointsAfterEvent+1:j*takeNPointsAfterEvent)=nan;
        else
            mergedData(i,(j-1)*takeNPointsAfterEvent+1:j*takeNPointsAfterEvent)=ResponseCued.unitbyunit_y(ResponseCued.fromWhichTrial==currTri & ResponseCued.fromWhichUnit==currU,zeroind:zeroind+takeNPointsAfterEvent-1);
        end
    end
end

mergedAllData=mergedData;
labels=zeros(size(mergedData,1),1);

temp=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',ResponseUncued,1); timesAlignedToEvent=temp.t; [~,zeroind]=nanmin(abs(timesAlignedToEvent));
us=unique(ResponseUncued.fromWhichUnit); allTri=unique(ResponseUncued.fromWhichTrial);
mergedData=nan(length(allTri),length(us)*takeNPointsAfterEvent);
for i=1:length(allTri)
    currTri=allTri(i);
    for j=1:length(us)
        currU=us(j);
        if sum(ResponseUncued.fromWhichTrial==currTri & ResponseUncued.fromWhichUnit==currU,'all','omitnan')==0
            % unit was not active on this trial, fill w nans
            mergedData(i,(j-1)*takeNPointsAfterEvent+1:j*takeNPointsAfterEvent)=nan;
        else
            mergedData(i,(j-1)*takeNPointsAfterEvent+1:j*takeNPointsAfterEvent)=ResponseUncued.unitbyunit_y(ResponseUncued.fromWhichTrial==currTri & ResponseUncued.fromWhichUnit==currU,zeroind:zeroind+takeNPointsAfterEvent-1);
        end
    end
end

mergedAllData=[mergedAllData; mergedData];
labels=[labels; ones(size(mergedData,1),1)];

% throw out rows where everything is nan (no spikes on these trials in this time window)


LDA_analysis_sub(mergedAllData,labels+1); % in B code, group labels start at 1, not 0

end

function LDA_analysis_sub(mergedAllData,mergedDataLabels)

% each row of mergedAllData is a different observation (e.g., trial)
% and each column is a predictor variable, e.g., firing rate at a different
% time point aligned to some event
% labels are the group labels for each row

testFraction=0.3;
trainFraction=1-testFraction;
nModels=10;

% takes, for example, physiology data points from post-outcome period
% and tries to classify according to labels

% Set up the statsTable to hold results
for sCounter=1:1
    statsTable.(['moment1_train_cond' num2str(sCounter) '_avg'])=zeros(0,1); % the average and SD of the first 4 moments of all the random train data sets
    statsTable.(['moment1_train_cond' num2str(sCounter) '_sd'])=zeros(0,1);
    statsTable.(['moment2_train_cond' num2str(sCounter) '_avg'])=zeros(0,1);
    statsTable.(['moment2_train_cond' num2str(sCounter) '_sd'])=zeros(0,1);
    statsTable.(['moment3_train_cond' num2str(sCounter) '_avg'])=zeros(0,1);
    statsTable.(['moment3_train_cond' num2str(sCounter) '_sd'])=zeros(0,1);
    statsTable.(['moment4_train_cond' num2str(sCounter) '_avg'])=zeros(0,1);
    statsTable.(['moment4_train_cond' num2str(sCounter) '_sd'])=zeros(0,1);
    statsTable.(['MSE_train_cond' num2str(sCounter) '_avg'])=zeros(0,1);
    statsTable.(['MSE_train_cond' num2str(sCounter) '_sd'])=zeros(0,1);
    statsTable.(['MSE_test_cond' num2str(sCounter) '_avg'])=zeros(0,1); % the average and SD of the MSE of the test data compared to the average of the train data
    statsTable.(['MSE_test_cond' num2str(sCounter) '_sd'])=zeros(0,1);

    dummyFill={0 0 0 0 0 0 0 0 0 0 0 0 };
end

% loop through all the conditions to test
tableCounter=0;

conditionSets=1;
howManyTrials=sum(mergedDataLabels==1,'all','omitnan');
nConditions=length(unique(mergedDataLabels));
for conditionCounter=1:length(conditionSets)

    if numel(mergedAllData)>0
        % LDA analysis
        minTrainSize=min(floor(trainFraction * howManyTrials));
        minTestSize=min(floor(testFraction * howManyTrials));

        trainSets=zeros(1, minTrainSize);
        testSets=zeros(1, minTestSize);

        trainPred=zeros(nModels, 1);
        testPred=zeros(nModels, 1);
        testPredHyp=zeros(nModels, 1);

        testMSE_Train=zeros(1, nModels);
        testMSE_Test=zeros(1, nModels);

        trainMoments=zeros(1, 4, nModels);

        for randCounter=1:2
            drawnow

            if randCounter==1
                randomShuffle=false;
            else
                randomShuffle=true;
            end
            tableCounter=tableCounter+1;

%             statsTable(tableCounter,:)=dummyFill;
%             statsTable.condition(tableCounter)=conditionCounter;
%             statsTable.mouseID(tableCounter)=mouseID;
%             statsTable.analysisMode(tableCounter)=analysisMode;
            statsTable.randomShuffle(tableCounter)=randomShuffle;
            statsTable.testFraction(tableCounter)=testFraction;
            statsTable.trainFraction(tableCounter)=trainFraction;
%             statsTable.channel(tableCounter)=channel;
            statsTable.LDA_nModels(tableCounter)=nModels;
%             statsTable.channelSwapped(tableCounter)=swappedChannels;

%             disp(' ')
%             disp(['Running models Channel ' num2str(channel) ' rand ' num2str(randCounter-1)])

            if randomShuffle
                disp('   with random shuffling');
            end

            for mCounter=1:nModels
                %          close all
                if nModels>100 && mCounter/(nModels/10)==floor(mCounter/(nModels/10))
                    disp(mCounter)
                end

                % set up indices of each test/train set
                for sCounter=1:nConditions
                    sIndices=find(mergedDataLabels==sCounter);

                    nIndices=length(sIndices);
                    testII=randperm(nIndices, minTestSize);
                    testI=sort(sIndices(testII));
                    testSets(sCounter,:)=testI;

                    nonTestIndices=setdiff(sIndices, testI);
                    trainII=randperm(length(nonTestIndices), minTrainSize);
                    trainSets(sCounter,:)=sort(nonTestIndices(trainII));
                end

                testIndices=reshape(testSets', numel(testSets), 1);
                trainIndices=reshape(trainSets', numel(trainSets), 1);

                % pull data - not really necessary but makes it easier
                mergedTestData=mergedAllData(testIndices,:);
                mergedTrainData=mergedAllData(trainIndices,:);
                mergedTestLabels=mergedDataLabels(testIndices);
                mergedTrainLabels=mergedDataLabels(trainIndices);

                % randomize labels and data on this balanced set
                originalTestLabels=mergedTestLabels;
                if randomShuffle % shuffle the labels maintaining balance
                    mergedTestLabels=mergedTestLabels(randperm(length(mergedTestLabels), length(mergedTestLabels)));
                    mergedTrainLabels=mergedTrainLabels(randperm(length(mergedTrainLabels), length(mergedTrainLabels)));
                end

                % run LDA model
                % with default - may overfit and reduce fit quality on test data
                try
                    modelType='quadratic'; % linear may be better for the means (vs. all)
                    Mdl = fitcdiscr(mergedTrainData, mergedTrainLabels, 'DiscrimType', modelType);
                catch
                    modelType='pseudoQuadratic'; % linear may be better for the means (vs. all)
                    Mdl = fitcdiscr(mergedTrainData, mergedTrainLabels, 'DiscrimType', modelType);
                end

                yy=predict(Mdl, mergedTrainData);
                trainPred(mCounter)=mean((mergedTrainLabels==yy));
                yy=predict(Mdl, mergedTestData);
                testPred(mCounter)=mean((mergedTestLabels==yy));

                % comparisons to means and calc moments
                for sCounter=1:nConditions
                    cTrainData=mergedTrainData(mergedTrainLabels==sCounter, :);
                    cTestData=mergedTestData(originalTestLabels==sCounter, :); %mergedTestLabels==sCounter,:);

                    testMSE_Train(sCounter, mCounter)=...
                        mean(mean((cTestData-mean(cTrainData, 1)).^2, 1).^0.5);

                    testMSE_Test(sCounter, mCounter)=...
                        mean(mean((cTestData-mean(cTestData, 1)).^2, 1).^0.5);

                    % moments
                    cTrainData=reshape(cTrainData, 1, numel(cTrainData));
                    trainMoments(sCounter, 1, mCounter)=mean(cTrainData);
                    trainMoments(sCounter, 2, mCounter)=var(cTrainData);
                    trainMoments(sCounter, 3, mCounter)=skewness(cTrainData);
                    trainMoments(sCounter, 4, mCounter)=kurtosis(cTrainData);
                end
            end

            %   with hyper parameter sweep and cross-validation - better fits
            Mdl = fitcdiscr(mergedTrainData,mergedTrainLabels, ...
                'OptimizeHyperparameters', 'auto',...
                'HyperparameterOptimizationOptions', ...
                struct('ShowPlots', false, 'Verbose', 0, ...
                'Repartition', true, ...
                'AcquisitionFunctionName','expected-improvement-plus'));
            %    close all
            yy=predict(Mdl, mergedTestData);
            testPredHyp=mean((mergedTestLabels==yy));
            yy=predict(Mdl, mergedTrainData);
            trainPredHyp=mean((mergedTrainLabels==yy));

            statsTable.modelType(tableCounter)=modelType;
            statsTable.LDA_trainPred_avg(tableCounter)=mean(trainPred);
            statsTable.LDA_trainPred_sd(tableCounter)=std(trainPred);
            disp(['OVERALL train set accuracy: ' ...
                num2str(statsTable.LDA_trainPred_avg(tableCounter))...
                ' +/- ' num2str(statsTable.LDA_trainPred_sd(tableCounter)) ])

            statsTable.LDA_testPred_avg(tableCounter)=mean(testPred);
            statsTable.LDA_testPred_sd(tableCounter)=std(testPred);
            disp(['OVERALL test set accuracy: '  ...
                num2str(statsTable.LDA_testPred_avg(tableCounter))...
                ' +/- ' num2str(statsTable.LDA_testPred_sd(tableCounter)) ])

            disp(['OVERALL hyp test set accuracy: '  num2str(testPredHyp) ' train ' num2str(trainPredHyp)])
            statsTable.LDA_hyp_trainPred(tableCounter)=trainPredHyp;
            statsTable.LDA_hyp_testPred(tableCounter)=testPredHyp;

            for sCounter=1:nConditions
                statsTable.(['cond' num2str(sCounter)])(tableCounter)=conditionStrings{sCounter};

                statsTable.(['MSE_test_cond' num2str(sCounter) '_avg'])(tableCounter)=mean(testMSE_Test(sCounter, :));
                statsTable.(['MSE_test_cond' num2str(sCounter) '_sd'])(tableCounter)=std(testMSE_Test(sCounter,:));

                statsTable.(['MSE_train_cond' num2str(sCounter) '_avg'])(tableCounter)=mean(testMSE_Train(sCounter, :));
                statsTable.(['MSE_train_cond' num2str(sCounter) '_sd'])(tableCounter)=std(testMSE_Train(sCounter,:));

                disp(['OVERALL MSE Cond ' num2str(sCounter) ' : test vs test '  ...
                    num2str(mean(testMSE_Test(sCounter,:)))  ' +/- ' num2str(std(testMSE_Test(sCounter,:))) ])
                disp(['OVERALL MSE Cond ' num2str(sCounter) ' : test vs train '  ...
                    num2str(mean(testMSE_Train(sCounter,:)))  ' +/- ' num2str(std(testMSE_Train(sCounter,:))) ])

                for momentCounter=1:4
                    statsTable.(['moment' num2str(momentCounter) '_train_cond' num2str(sCounter) '_avg'])(tableCounter)=...
                        mean(squeeze(trainMoments(sCounter, momentCounter, :)));
                    statsTable.(['moment' num2str(momentCounter) '_train_cond' num2str(sCounter) '_sd'])(tableCounter)=...
                        std(squeeze(trainMoments(sCounter, momentCounter, :)));
                end
            end
        end
    end
end

% calcalate dPrimes
dPrimeTable=table;
dPrimeTable.mouseID=statsTable.mouseID(1:2:end);
dPrimeTable.condition=statsTable.condition(1:2:end);
dPrimeTable.channel=statsTable.channel(1:2:end);

colNames=statsTable.Properties.VariableNames;
avgColumns=find(contains(colNames, '_avg'));

for cCounter=1:nConditions
    dPrimeTable.(['cond' num2str(cCounter)])=statsTable.(['cond' num2str(cCounter)])(1:2:end);
end

dPrimeArray=zeros(size(dPrimeTable, 1), 0);
dPrimeNames={};

for cCounter=avgColumns
    cName=statsTable.Properties.VariableNames{cCounter};
    cNameDash=strfind(cName, '_avg');
    cNameStub=cName(1:cNameDash(end)-1);
    cNameDP=[cNameStub '_dp'];
    a1=statsTable{1:2:end, cCounter};
    a2=statsTable{2:2:end, cCounter};
    s1=statsTable{1:2:end, cCounter+1};
    s2=statsTable{2:2:end, cCounter+1};
    dPrimeTable.(cNameDP)=abs(a1-a2)./((s1.^2+s2.^2).^0.5);
    dPrimeArray(:,end+1)=dPrimeTable.(cNameDP);
    dPrimeNames{end+1}=removeDash(cNameStub);
end

end