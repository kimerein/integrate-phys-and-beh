function bootstrapFitDataToReachingRTModel(data,dataForModel,fits,behaviorEvent,ITI,saveDir)

nRuns=50;
as=nan(1,nRuns);
bs=nan(1,nRuns);
lses=nan(1,nRuns);
errs=nan(1,nRuns);
for i=1:nRuns
    disp(i);
    [~,a,b,lse,err,n,best_predict,differ]=compareReachingRTModelAndData(data,dataForModel,fits,behaviorEvent,ITI,1);
    if i==1
        sum_predictions=zeros(size(best_predict));
        sum_differs=zeros(size(differ));
    end
    sum_predictions=sum_predictions+best_predict;
    sum_differs=sum_differs+differ;
    as(i)=a;
    bs(i)=b;
    lses(i)=lse;
    errs(i)=err;
end
sum_predictions=sum_predictions./nRuns;
sum_differs=sum_differs./nRuns;

% saveDir=[saveDir '\fitsToData\' data.event_name];
saveDir=[saveDir '\fitsToData_excludeFastRTSlowing\' data.event_name];
mkdir(saveDir);
RPEAndRateFits.realData=n;
RPEAndRateFits.meanPredictions=sum_predictions;
RPEAndRateFits.meanDiffers=sum_differs;
RPEAndRateFits.rateCoefficient=as;
RPEAndRateFits.rpeCoefficient=bs;
RPEAndRateFits.lse=lses;
RPEAndRateFits.error=errs;
save([saveDir '\RPEAndRateFits.mat'],'RPEAndRateFits');