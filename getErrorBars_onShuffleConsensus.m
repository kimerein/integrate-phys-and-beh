function [all_accs,all_accs_3way]=getErrorBars_onShuffleConsensus(cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,nRuns)

countunits=10:10:200;
all_accs=nan(length(countunits),nRuns);
all_accs_3way=nan(length(countunits),nRuns);

for unitsCounter=1:length(countunits)
for countRuns=1:nRuns

    % just shuffle gp1 vs gp2
    shuffle_consensus=cued_success_Response.consensus_idx;
    f=find(~isnan(shuffle_consensus));
    r=randperm(length(f));
    shuffle_consensus(f)=shuffle_consensus(f(r));

%     shuffle_consensus=cued_success_Response.consensus_idx;

    binsForTuning{1}=[-10 -9 10]; binsForTuning{2}=[-10 -9 10];
    tuningOutput=plotUnitSummariesAfterTCAlabels(shuffle_consensus,cued_success_Response.cXfail_sus1to5sec,cued_success_Response,cued_failure_Response,uncued_success_Response,uncued_failure_Response,[],'uncuedOverCued','tuning',binsForTuning);
    [allgp1_cuedfailFR,allgp1_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp1_fail','grp1_fail_uncue',2,[1 5]);
    [allgp2_cuedfailFR,allgp2_uncuedfailFR]=plotTuningOutputScatter(tuningOutput,'grp2_fail','grp2_fail_uncue',2,[1 5]);
    plotTuningOutputScatter(tuningOutput,'grp1_succ_uncue','grp1_fail_uncue',2,[1 5]);
    [allgp1_cuedsuccFR,allgp1_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp1_succ','grp1_succ_uncue',2,[1 5]);
    [allgp2_cuedsuccFR,allgp2_uncuedsuccFR]=plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_succ_uncue',2,[1 5]);
    plotTuningOutputScatter(tuningOutput,'grp2_succ','grp2_fail',2,[1 5]);
    close all;

%     % shuffle trial labels
%     all1=[allgp1_cuedfailFR; allgp1_uncuedfailFR; allgp1_cuedsuccFR; allgp1_uncuedsuccFR];
%     all2=[allgp2_cuedfailFR; allgp2_uncuedfailFR; allgp2_cuedsuccFR; allgp2_uncuedsuccFR];
%     rall1=randperm(length(all1)); rall2=randperm(length(all2));
%     all1=all1(rall1); all2=all2(rall2);
%     inds1=[ones(size(allgp1_cuedfailFR)); 2*ones(size(allgp1_uncuedfailFR)); 3*ones(size(allgp1_cuedsuccFR)); 4*ones(size(allgp1_uncuedsuccFR))];
%     inds2=[ones(size(allgp2_cuedfailFR)); 2*ones(size(allgp2_uncuedfailFR)); 3*ones(size(allgp2_cuedsuccFR)); 4*ones(size(allgp2_uncuedsuccFR))];
%     allgp1_cuedfailFR=all1(inds1==1); allgp1_uncuedfailFR=all1(inds1==2); allgp1_cuedsuccFR=all1(inds1==3); allgp1_uncuedsuccFR=all1(inds1==4);
%     allgp2_cuedfailFR=all2(inds2==1); allgp2_uncuedfailFR=all2(inds2==2); allgp2_cuedsuccFR=all2(inds2==3); allgp2_uncuedsuccFR=all2(inds2==4);

    figure(); nBoot=100; nUnits=countunits(unitsCounter);
    takeThese_gp1=nan(nBoot,nUnits); takeThese_gp2=nan(nBoot,nUnits);
    for i=1:nBoot
        takeThese_gp1(i,:)=randsample(length(allgp1_cuedsuccFR),nUnits,true); takeThese_gp2(i,:)=randsample(length(allgp2_cuedsuccFR),nUnits,true);
    end
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_cuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedsuccFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[temp1 temp2]; ylabels=[ones(size(temp1,1),1)];
    scatter(temp1,temp2,[],'g'); hold on; scatter(nanmean(temp1),nanmean(temp2),[],'g','filled'); cuedsuccmeanx=nanmean(temp1); cuedsuccmeany=nanmean(temp2);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_cuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedfailFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp1 temp2]]; ylabels=[ylabels; 2*ones(size(temp1,1),1)];
    scatter(temp1,temp2,[],'r'); scatter(nanmean(temp1),nanmean(temp2),[],'r','filled'); cuedfailmeanx=nanmean(temp1); cuedfailmeany=nanmean(temp2);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_uncuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedsuccFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp1 temp2]]; ylabels=[ylabels; 3*ones(size(temp1,1),1)];
    scatter(temp1,temp2,[],'b'); scatter(nanmean(temp1),nanmean(temp2),[],'b','filled'); uncuedsuccmeanx=nanmean(temp1); uncuedsuccmeany=nanmean(temp2);
    temp1=nan(nBoot,1); temp2=nan(nBoot,1);
    for i=1:nBoot
        temp1(i)=nanmean(allgp1_uncuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedfailFR(takeThese_gp2(i,:)));
    end
    Xmatrix=[Xmatrix; [temp1 temp2]]; ylabels=[ylabels; 4*ones(size(temp1,1),1)];
    scatter(temp1,temp2,[],'y'); scatter(nanmean(temp1),nanmean(temp2),[],'y','filled'); uncuedfailmeanx=nanmean(temp1); uncuedfailmeany=nanmean(temp2);
    scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
    xlabel('Gp 1 average unit firing rate'); ylabel('Gp 2 average unit firing rate');

%     takeThese_gp1=nan(nBoot,nUnits); takeThese_gp2=nan(nBoot,nUnits);
%     for i=1:nBoot
%         takeThese_gp1(i,:)=randsample(length(allgp1_cuedsuccFR),nUnits,true); takeThese_gp2(i,:)=randsample(length(allgp2_cuedsuccFR),nUnits,true);
%     end
%     temp1=nan(nBoot,1); temp2=nan(nBoot,1);
%     for i=1:nBoot
%         temp1(i)=nanmean(allgp1_cuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedsuccFR(takeThese_gp2(i,:)));
%     end
%     Xmatrix=[temp2-temp1 temp2+temp1]; ylabels=[ones(size(temp1,1),1)];
%     scatter(temp2-temp1,temp2+temp1,[],'g'); hold on; scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'g','filled'); cuedsuccmeanx=nanmean(temp2-temp1); cuedsuccmeany=nanmean(temp2+temp1);
%     temp1=nan(nBoot,1); temp2=nan(nBoot,1);
%     for i=1:nBoot
%         temp1(i)=nanmean(allgp1_cuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_cuedfailFR(takeThese_gp2(i,:)));
%     end
%     Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 2*ones(size(temp1,1),1)];
%     scatter(temp2-temp1,temp2+temp1,[],'r'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'r','filled'); cuedfailmeanx=nanmean(temp2-temp1); cuedfailmeany=nanmean(temp2+temp1);
%     temp1=nan(nBoot,1); temp2=nan(nBoot,1);
%     for i=1:nBoot
%         temp1(i)=nanmean(allgp1_uncuedsuccFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedsuccFR(takeThese_gp2(i,:)));
%     end
%     Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 3*ones(size(temp1,1),1)];
%     scatter(temp2-temp1,temp2+temp1,[],'b'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'b','filled'); uncuedsuccmeanx=nanmean(temp2-temp1); uncuedsuccmeany=nanmean(temp2+temp1);
%     temp1=nan(nBoot,1); temp2=nan(nBoot,1);
%     for i=1:nBoot
%         temp1(i)=nanmean(allgp1_uncuedfailFR(takeThese_gp1(i,:))); temp2(i)=nanmean(allgp2_uncuedfailFR(takeThese_gp2(i,:)));
%     end
%     Xmatrix=[Xmatrix; [temp2-temp1 temp2+temp1]]; ylabels=[ylabels; 4*ones(size(temp1,1),1)];
%     scatter(temp2-temp1,temp2+temp1,[],'y'); scatter(nanmean(temp2-temp1),nanmean(temp2+temp1),[],'y','filled'); uncuedfailmeanx=nanmean(temp2-temp1); uncuedfailmeany=nanmean(temp2+temp1);
%     scatter((cuedsuccmeanx+cuedfailmeanx+uncuedsuccmeanx+uncuedfailmeanx)/4,(cuedsuccmeany+cuedfailmeany+uncuedsuccmeany+uncuedfailmeany)/4,[],'k','filled');
%     xlabel('Gp 2 minus gp 1 average unit firing rate'); ylabel('Gp 2 plus gp 1 average unit firing rate');
    
    % LDA
    ldaModel=fitcdiscr(Xmatrix,ylabels);
    predictedY=predict(ldaModel,Xmatrix);
    accuracy=sum(predictedY==ylabels)/length(ylabels);
    disp(['Accuracy of LDA on training set: ', num2str(accuracy * 100), '%']);
    all_accs(unitsCounter,countRuns)=accuracy * 100;
    close all;

    % 3-way classification
    ylabels(ylabels==4)=2;
    ldaModel=fitcdiscr(Xmatrix,ylabels);
    predictedY=predict(ldaModel,Xmatrix);
    accuracy=sum(predictedY==ylabels)/length(ylabels);
    disp(['Accuracy of LDA on training set: ', num2str(accuracy * 100), '%']);
    all_accs_3way(unitsCounter,countRuns)=accuracy * 100;
    close all;

end
end

end