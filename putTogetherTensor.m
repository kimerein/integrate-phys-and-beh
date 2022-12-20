function putTogetherTensor(whichSess,downSampBy,takeNPointsAfterEvent,nPointsBeforeEvent,response_to_plot1,response_to_plot2,dd,saveDir)
% wrapper for code from B

% get all units recorded simultaneously from one session
% whichSess=3;
% takeNPointsAfterEvent=100;

if length(whichSess)>1
    dd=dd(whichSess);
end

whichUnitsToGrab='_'; plotUnitCriteria=[-100 0 0 1 0]; getCriteriaForUnitsToPlot(plotUnitCriteria);
setForUn=settingsForStriatumUnitPlots;
if setForUn.keepAllSingleTrials~=true
    error('need trial by trial data for LDA analysis');
end
response_to_plot=response_to_plot1; %'cued_success'; 
if length(whichSess)>1
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
    end
    ResponseCued=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
else
    ResponseCued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
end
ResponseCued.unitbyunit_x=downSampMatrix(ResponseCued.unitbyunit_x,downSampBy);
ResponseCued.unitbyunit_y=downSampMatrix(ResponseCued.unitbyunit_y,downSampBy);
ResponseCued=makeUnitsUnique(ResponseCued);
response_to_plot=response_to_plot2; %'uncued_success'; 
if length(whichSess)>1
    dd_more=cell(1,length(dd));
    for i=1:length(dd)
        dd_more{i}=[dd{i} sep response_to_plot];
    end
    ResponseUncued=getAndSaveResponse(dd_more,whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
else
    ResponseUncued=getAndSaveResponse([dd{whichSess} sep response_to_plot],whichUnitsToGrab,settingsForStriatumUnitPlots,[]);
end
ResponseUncued.unitbyunit_x=downSampMatrix(ResponseUncued.unitbyunit_x,downSampBy);
ResponseUncued.unitbyunit_y=downSampMatrix(ResponseUncued.unitbyunit_y,downSampBy);
ResponseUncued=makeUnitsUnique(ResponseUncued);

temp=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',ResponseCued,1); timesAlignedToEvent=temp.t; [~,zeroind]=nanmin(abs(timesAlignedToEvent));
us=unique(ResponseCued.fromWhichUnit); allTri=unique(ResponseCued.fromWhichTrial);
mergedData=nan(length(allTri),length(us)*(takeNPointsAfterEvent+nPointsBeforeEvent));
didEventOccur=nan(length(allTri),1);
for i=1:length(allTri)
    currTri=allTri(i);
    didEventOccur(i)=mode(ResponseCued.isEventInThisTrial(ResponseCued.fromWhichTrial==currTri));
    for j=1:length(us)
        currU=us(j);
        if sum(ResponseCued.fromWhichTrial==currTri & ResponseCued.fromWhichUnit==currU,'all','omitnan')==0
            % unit was not active on this trial, fill w nans
            mergedData(i,(j-1)*(takeNPointsAfterEvent+nPointsBeforeEvent)+1:j*(takeNPointsAfterEvent+nPointsBeforeEvent))=nan;
        else
            temp=ResponseCued.unitbyunit_y(ResponseCued.fromWhichTrial==currTri & ResponseCued.fromWhichUnit==currU,zeroind-nPointsBeforeEvent:zeroind+takeNPointsAfterEvent-1);
%             temp=smoothdata(temp,'gaussian',7);
            mergedData(i,(j-1)*(takeNPointsAfterEvent+nPointsBeforeEvent)+1:j*(takeNPointsAfterEvent+nPointsBeforeEvent))=temp;
        end
    end
end
% take only trials where event occurred
mergedData=mergedData(didEventOccur==1,:);
% fill in nans with zeros, because no spiking
% mergedData(isnan(mergedData))=0;

mergedAllData=mergedData;
labels=zeros(size(mergedData,1),1);

temp=plotVariousSUResponsesAlignedToBeh('meanAcrossUnits',ResponseUncued,1); timesAlignedToEvent=temp.t; [~,zeroind]=nanmin(abs(timesAlignedToEvent));
us=unique(ResponseUncued.fromWhichUnit); allTri=unique(ResponseUncued.fromWhichTrial);
mergedData=nan(length(allTri),length(us)*(takeNPointsAfterEvent+nPointsBeforeEvent));
didEventOccur=nan(length(allTri),1);
for i=1:length(allTri)
    currTri=allTri(i);
    didEventOccur(i)=mode(ResponseUncued.isEventInThisTrial(ResponseUncued.fromWhichTrial==currTri));
    for j=1:length(us)
        currU=us(j);
        if sum(ResponseUncued.fromWhichTrial==currTri & ResponseUncued.fromWhichUnit==currU,'all','omitnan')==0
            % unit was not active on this trial, fill w nans
            mergedData(i,(j-1)*(takeNPointsAfterEvent+nPointsBeforeEvent)+1:j*(takeNPointsAfterEvent+nPointsBeforeEvent))=nan;
        else
            temp=ResponseUncued.unitbyunit_y(ResponseUncued.fromWhichTrial==currTri & ResponseUncued.fromWhichUnit==currU,zeroind-nPointsBeforeEvent:zeroind+takeNPointsAfterEvent-1);
%             temp=smoothdata(temp,'gaussian',7);
            mergedData(i,(j-1)*(takeNPointsAfterEvent+nPointsBeforeEvent)+1:j*(takeNPointsAfterEvent+nPointsBeforeEvent))=temp;
        end
    end
end
disp(['format of mergedAllData is ' num2str(nPointsBeforeEvent) ' before event and ' num2str(takeNPointsAfterEvent) ' from each unit, trial by trial']);
% take only trials where event occurred
mergedData=mergedData(didEventOccur==1,:);
% fill in nans with zeros, because no spiking
% mergedData(isnan(mergedData))=0;

mergedAllData=[mergedAllData; mergedData];
labels=[labels; ones(size(mergedData,1),1)];

% drop columns that are all nan
% mergedAllData=mergedAllData(:,~all(isnan(mergedAllData),1));
% % drop rows that are all nan
% allnanrows=all(isnan(mergedAllData),2);
% mergedAllData=mergedAllData(~allnanrows,:);
% labels=labels(~allnanrows);

% close all;

figure();
imagesc([mergedAllData repmat(labels,1,10).*max(mergedAllData,[],'all','omitnan')]);
xlim([0 size(mergedAllData,2)+10]);

% figure(); imagesc(isnan(mergedAllData));

[tens,labels]=reshapeIntoTensor(mergedAllData,labels,nPointsBeforeEvent,takeNPointsAfterEvent);
save([saveDir '_tensor.mat'],'tens');
save([saveDir '_labels.mat'],'labels');

end

function Response=makeUnitsUnique(Response)

uSess=unique(Response.fromWhichSess_forTrials);
sessUOffset=0;
triOffset=0;
backupSess=Response.fromWhichSess_forTrials;
backupUs=Response.fromWhichUnit;
backupTrials=Response.fromWhichTrial;
for i=1:length(uSess)
    currSess=uSess(i);
    whichUs=unique(Response.fromWhichUnit(backupSess==currSess));
    for j=1:length(whichUs)
        currU=whichUs(j);
        Response.fromWhichUnit(backupSess==currSess & backupUs==currU)=currU+sessUOffset;
    end
    currTris=unique(backupTrials(backupSess==currSess));
    for j=1:length(currTris)
        currT=currTris(j);
        Response.fromWhichTrial(backupSess==currSess & backupTrials==currT)=currT+triOffset;
    end
    sessUOffset=sessUOffset+max(whichUs,[],'all','omitnan')+1;
    triOffset=triOffset+max(currTris,[],'all','omitnan')+1;
end

end

function LDA_analysis_sub(mergedAllData,mergedDataLabels,typMod)

% each row of mergedAllData is a different observation (e.g., trial)
% and each column is a predictor variable, e.g., firing rate at a different
% time point aligned to some event
% labels are the group labels for each row

testFraction=0.3;
trainFraction=1-testFraction;
nModels=100;

% takes, for example, physiology data points from post-outcome period
% and tries to classify according to labels

% set up the tables
statsTable=table;
statsTable.mouseID=strings(0,1);
statsTable.condition=zeros(0,1);
statsTable.channel=zeros(0,1);
statsTable.channelSwapped=false(0,1);
statsTable.analysisMode=strings(0,1);
statsTable.modelType=strings(0,1);
statsTable.randomShuffle=false(0,1);
statsTable.testFraction=zeros(0,1);
statsTable.trainFraction=zeros(0,1);
statsTable.LDA_nModels=zeros(0,1);
statsTable.LDA_trainPred_avg=zeros(0,1);
statsTable.LDA_trainPred_sd=zeros(0,1);
statsTable.LDA_testPred_avg=zeros(0,1);
statsTable.LDA_testPred_sd=zeros(0,1);
statsTable.LDA_hyp_trainPred=zeros(0,1);
statsTable.LDA_hyp_testPred=zeros(0,1);

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
nConditions=length(unique(mergedDataLabels));
for conditionCounter=1:length(conditionSets)

    if numel(mergedAllData)>0
        % LDA analysis
        trainPred=zeros(nModels, 1);
        testPred=zeros(nModels, 1);
        testPredHyp=zeros(nModels, 1);

        testMSE_Train=zeros(1, nModels);
        testMSE_Test=zeros(1, nModels);

        trainMoments=zeros(1, 4, nModels);

        % linear and z-score
        dataSize=size(mergedAllData);
        mergedAllData=reshape(mergedAllData, 1, numel(mergedAllData)); % linearize
        mergedAllData=normalize(mergedAllData); % z-score
        mergedAllData=reshape(mergedAllData, dataSize);

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
                mergedTestData=[];
                mergedTrainData=[];
                mergedTestLabels=[];
                mergedTrainLabels=[];

                allMinSizes_train=nan(1,nConditions);
                allMinSizes_test=nan(1,nConditions);
                for sCounter=1:nConditions
                    howManyTrials=sum(mergedDataLabels==sCounter,'all','omitnan');
                    allMinSizes_train(sCounter)=min(floor(trainFraction * howManyTrials));
                    allMinSizes_test(sCounter)=min(floor(testFraction * howManyTrials));
                end
                minTrainSize=min(allMinSizes_train,[],'all');
                minTestSize=min(allMinSizes_test,[],'all');
                if randCounter==1 && mCounter==1
                    disp(['using ' num2str(minTrainSize) ' trials for train']);
                    disp(['using ' num2str(minTestSize) ' trials for test']);
                end

                for sCounter=1:nConditions
%                     howManyTrials=sum(mergedDataLabels==sCounter,'all','omitnan');
%                     minTrainSize=min(floor(trainFraction * howManyTrials));
%                     minTestSize=min(floor(testFraction * howManyTrials));

                    sIndices=find(mergedDataLabels==sCounter);

                    nIndices=length(sIndices);
                    testII=randperm(nIndices, minTestSize);
                    testI=sort(sIndices(testII));
                    
                    mergedTestData=[mergedTestData; mergedAllData(testI,:)];
                    mergedTestLabels=[mergedTestLabels; mergedDataLabels(testI)];

                    nonTestIndices=setdiff(sIndices, testI);
                    trainII=randperm(length(nonTestIndices), minTrainSize);
                    trainSets=sort(nonTestIndices(trainII));

                    mergedTrainData=[mergedTrainData; mergedAllData(trainSets,:)];
                    mergedTrainLabels=[mergedTrainLabels; mergedDataLabels(trainSets)];
                end

                % randomize labels and data on this balanced set
                originalTestLabels=mergedTestLabels;
                if randomShuffle % shuffle the labels maintaining balance
                    mergedTestLabels=mergedTestLabels(randperm(length(mergedTestLabels), length(mergedTestLabels)));
                    mergedTrainLabels=mergedTrainLabels(randperm(length(mergedTrainLabels), length(mergedTrainLabels)));
                end

                % if all nan col, then need to drop from both training and
                % test set
                allNanCol=all(isnan(mergedTrainData(mergedTrainLabels==1,:)),1) | all(isnan(mergedTrainData(mergedTrainLabels==2,:)),1);
                mergedTrainData=mergedTrainData(:,~allNanCol);
                mergedTestData=mergedTestData(:,~allNanCol);
                mergedTrainData(isnan(mergedTrainData))=0; 
                mergedTestData(isnan(mergedTestData))=0; 

                switch typMod
                    case 'LDA'
                        % run LDA model
                        % with default - may overfit and reduce fit quality on test data
                        %                 try
                        %                     modelType='quadratic'; % linear may be better for the means (vs. all)
                        %                     Mdl = fitcdiscr(mergedTrainData, mergedTrainLabels, 'DiscrimType', modelType);
                        %                 catch
                        try
                            modelType='pseudoQuadratic'; % linear may be better for the means (vs. all)
                            stilltakenonan=~any(isnan(mergedTrainData),2);
                            mergedTrainData=mergedTrainData(stilltakenonan,:);
                            mergedTrainLabels=mergedTrainLabels(stilltakenonan);
                            if size(mergedTrainData,1)<2
                                continue
                            else
                                Mdl = fitcdiscr(mergedTrainData, mergedTrainLabels, 'DiscrimType', modelType);
                            end
                        catch
                            continue
                        end

                        yy=predict(Mdl, mergedTrainData);
                        trainPred(mCounter)=mean((mergedTrainLabels==yy));
                        yy=predict(Mdl, mergedTestData);
                        testPred(mCounter)=mean((mergedTestLabels==yy));
                    case 'linearRegression'
                        % do linear regression instead
                        Mdl=fitlm(mergedTrainData,mergedTrainLabels,'constant');
                    case 'SVM'
                        % SVM
                        X=mergedTrainData;
                        Y=mergedTrainLabels;
                        Mdl=fitcsvm(X,Y,'KernelScale','auto','Standardize',true,'OutlierFraction',0.05);
                        sv=Mdl.SupportVectors;
%                         figure;
%                         gscatter(X(:,1),X(:,3),Y,'br','xo');
                        predict_trainData=predict(Mdl,X);
%                         figure;
%                         gscatter(X(:,1),X(:,3),predictions,'rb','ox');
%                         title('Prediction of Training Data View 1');
                        predictions=predict(Mdl,mergedTestData);
%                         figure;
%                         gscatter(X(:,1),X(:,3),predictions,'rb','ox');
%                         title('Prediction of Test Data View 1');
                        modelType='SVM';
                end

                % comparisons to means and calc moments
                for sCounter=1:nConditions
                    cTrainData=mergedTrainData(mergedTrainLabels==sCounter, :);
                    cTestData=mergedTestData(originalTestLabels==sCounter, :); %mergedTestLabels==sCounter,:);

                    switch typMod
                        case 'LDA'
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
                        case 'SVM'
                            testMSE_Train(sCounter, mCounter)=...
                                sum(predict_trainData==mergedTrainLabels)/length(mergedTrainLabels); % fraction correct

                            testMSE_Test(sCounter, mCounter)=...
                                sum(predictions==mergedTestLabels)/length(mergedTestLabels); % fraction correct
                            % moments
                            cTrainData=reshape(cTrainData, 1, numel(cTrainData));
                            trainMoments(sCounter, 1, mCounter)=mean(cTrainData);
                            trainMoments(sCounter, 2, mCounter)=var(cTrainData);
                            trainMoments(sCounter, 3, mCounter)=skewness(cTrainData);
                            trainMoments(sCounter, 4, mCounter)=kurtosis(cTrainData);
                    end
                end
            end

            switch typMod
                    case 'LDA'
                        %   with hyper parameter sweep and cross-validation - better fits
                        Mdl = fitcdiscr(mergedTrainData,mergedTrainLabels, ...
                            'OptimizeHyperparameters', 'auto',...
                            'discrimType','pseudoLinear',...
                            'HyperparameterOptimizationOptions', ...
                            struct('ShowPlots', false, 'Verbose', 0, ...
                            'Repartition', true, ...
                            'AcquisitionFunctionName','expected-improvement-plus'));
                        %    close all
                        yy=predict(Mdl, mergedTestData);
                        testPredHyp=mean((mergedTestLabels==yy));
                        yy=predict(Mdl, mergedTrainData);
                        trainPredHyp=mean((mergedTrainLabels==yy));

%                         disp(['OVERALL hyp test set accuracy: '  num2str(testPredHyp) ' train ' num2str(trainPredHyp)])
%                         statsTable.LDA_hyp_trainPred(tableCounter)=trainPredHyp;
%                         statsTable.LDA_hyp_testPred(tableCounter)=testPredHyp;
                case 'SVM'
                    disp('for SVM, look at, for example, MSE_train_cond1_avg or MSE_test_cond1_avg for average accuray');
            end


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

            
            for sCounter=1:nConditions
%                 statsTable.(['cond' num2str(sCounter)])(tableCounter)=conditionStrings{sCounter};

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

disp(statsTable);

% calcalate dPrimes
dPrimeTable=table;
% dPrimeTable.mouseID=statsTable.mouseID(1:2:end);
% dPrimeTable.condition=statsTable.condition(1:2:end);
% dPrimeTable.channel=statsTable.channel(1:2:end);

colNames=statsTable.Properties.VariableNames;
avgColumns=find(contains(colNames, '_avg'));

% for cCounter=1:nConditions
%     dPrimeTable.(['cond' num2str(cCounter)])=statsTable.(['cond' num2str(cCounter)])(1:2:end);
% end

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
%     dPrimeNames{end+1}=removeDash(cNameStub);
end

disp(dPrimeTable);

end