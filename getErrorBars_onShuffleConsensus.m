function all_accs=getErrorBars_onShuffleConsensus(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,nRuns)

countunits=10:10:200;
all_accs=nan(length(countunits),nRuns);

for unitsCounter=1:length(countunits)
for countRuns=1:nRuns

    shuffle_consensus=cued_success_Response.consensus_idx;
    f=find(~isnan(shuffle_consensus));
    r=randperm(length(f));
    shuffle_consensus(f)=shuffle_consensus(f(r));

    binsForTuning{1}=[-10 -9 10]; binsForTuning{2}=[-10 -9 10];
    tuningOutput=plotUnitSummariesAfterTCAlabels(shuffle_consensus,cued_success_Response.cXfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
    [allgp1_cuedfailFR,allgp1_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[1 5]);
    [allgp2_cuedfailFR,allgp2_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',2,[1 5]);
    plotTuningOutputScatter(tuningOutput,'grp1_succ_uncue','grp1_fail_uncue',2,[1 5]);
    [allgp1_cuedsuccFR,allgp1_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',2,[1 5]);
    [allgp2_cuedsuccFR,allgp2_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',2,[1 5]);
    plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_fail',2,[1 5]);
    close all;

    figure(); nBoot=100; nUnits=countunits(unitsCounter);
    takeThese_gp1=nan(nBoot,nUnits); takeThese_gp2=nan(nBoot,nUnits);
    for i=1:nBoot
        takeThese_gp1(i,:)=randsample(length(allgp1_cuedsuccFR),nUnits); takeThese_gp2(i,:)=randsample(length(allgp2_cuedsuccFR),nUnits);
    end
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_cuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedsuccFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[temp2-temp1 temp2+temp1]; ylabels=[ones(size(temp1,1),1)];
    scatter(temp2-temp1,temp2+temp1,[],'g'); hold on; scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'g','filled'); cuedsuccmeanx=nanmean(temp2-temp1); cuedsuccmeany=nanmean(temp2+temp1);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_cuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedfailFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 2*ones(size(temp1,1),1)];
    scatter(temp2-temp1,temp2+temp1,[],'r'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'r','filled'); cuedfailmeanx=nanmean(temp2-temp1); cuedfailmeany=nanmean(temp2+temp1);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_uncuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedsuccFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 3*ones(size(temp1,1),1)];
    scatter(temp2-temp1,temp2+temp1,[],'b'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'b','filled'); uncuedsuccmeanx=nanmean(temp2-temp1); uncuedsuccmeany=nanmean(temp2+temp1);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_uncuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedfailFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 4*ones(size(temp1,1),1)];
    scatter(temp2-temp1,temp2+temp1,[],'y'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'y','filled'); uncuedfailmeanx=nanmean(temp2-temp1); uncuedfailmeany=nanmean(temp2+temp1);
    scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
    xlabel('Gp 2 minus gp 1 average unit firing rate'); ylabel('Gp 2 plus gp 1 average unit firing rate');
    % LDA
    ldaModel=fitcdiscr(Xmatrix,ylabels);
    predictedY=predict(ldaModel,Xmatrix);
    accuracy=sum(predictedY==ylabels)/length(ylabels);
    disp(['Accuracy of LDA on training set: ', num2str(accuracy * 100), '%']);
    all_accs(unitsCounter,countRuns)=accuracy * 100;
    close all;

end
end

end