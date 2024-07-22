function [all_accs,all_accs_3way]=getErrorBars_onShuffleConsensus_trialbytrial(a,b,cued_success_Response,nRuns)

countunits=10:10:200; %10:10:200;
all_accs=nan(length(countunits),nRuns);
all_accs_3way=nan(length(countunits),nRuns);

temper_backup=[a.fortbytclass.unitfr_success; a.fortbytclass.unitfr_failure; b.fortbytclass.unitfr_success; b.fortbytclass.unitfr_failure];

for unitsCounter=1:length(countunits)
for countRuns=1:nRuns

    shuffle_consensus=cued_success_Response.consensus_idx;
    f=find(~isnan(shuffle_consensus));
    r=randperm(length(f));
    shuffle_consensus(f)=shuffle_consensus(f(r));

%     shuffle_consensus=cued_success_Response.consensus_idx;
    
    if countunits(unitsCounter)<75
        nUnitsNow=100;
    elseif countunits(unitsCounter)>75 && countunits(unitsCounter)<125
        nUnitsNow=150;
    elseif countunits(unitsCounter)>125 && countunits(unitsCounter)<175
        nUnitsNow=200;
    else
        nUnitsNow=250;
    end

    % shuffle before
%     temper=[temper_backup];
%     temper_ids=[ones(size(a.fortbytclass.unitfr_success)); 2*ones(size(a.fortbytclass.unitfr_failure)); 3*ones(size(b.fortbytclass.unitfr_success)); 4*ones(size(b.fortbytclass.unitfr_failure))];
%     temper=temper(randperm(length(temper)));
%     a.fortbytclass.unitfr_success=temper(temper_ids==1);
%     a.fortbytclass.unitfr_failure=temper(temper_ids==2);
%     b.fortbytclass.unitfr_success=temper(temper_ids==3);
%     b.fortbytclass.unitfr_failure=temper(temper_ids==4);

    out_decode=decodeTrialByTrialType(a.fortbytclass,b.fortbytclass,shuffle_consensus,100,nUnitsNow,countunits(unitsCounter),true,false,false,false,false); % nBoots,nUnits,nTrials,withReplacement,addThirdAxis,nanAllZeros,justBoostrapTrials,collapseWithinUnit
    
    % shuffle trial labels
%     all1=[out_decode.cuedsucc_temp1; out_decode.cuedfail_temp1; out_decode.uncuedsucc_temp1; out_decode.uncuedfail_temp1];
%     all2=[out_decode.cuedsucc_temp2; out_decode.cuedfail_temp2; out_decode.uncuedsucc_temp2; out_decode.uncuedfail_temp2];
%     rall1=randperm(length(all1)); rall2=randperm(length(all2));
%     all1=all1(rall1); all2=all2(rall2);
%     inds1=[ones(size(out_decode.cuedsucc_temp1)); 2*ones(size(out_decode.cuedfail_temp1)); 3*ones(size(out_decode.uncuedsucc_temp1)); 4*ones(size(out_decode.uncuedfail_temp1))];
%     inds2=[ones(size(out_decode.cuedsucc_temp2)); 2*ones(size(out_decode.cuedfail_temp2)); 3*ones(size(out_decode.uncuedsucc_temp2)); 4*ones(size(out_decode.uncuedfail_temp2))];
%     out_decode.cuedsucc_temp1=all1(inds1==1); out_decode.cuedfail_temp1=all1(inds1==2); out_decode.uncuedsucc_temp1=all1(inds1==3); out_decode.uncuedfail_temp1=all1(inds1==4);
%     out_decode.cuedsucc_temp2=all2(inds2==1); out_decode.cuedfail_temp2=all2(inds2==2); out_decode.uncuedsucc_temp2=all2(inds2==3); out_decode.uncuedfail_temp2=all2(inds2==4);
    
%     % one mapping
%     figure();
%     scatter(out_decode.cuedsucc_temp2-out_decode.cuedsucc_temp1,out_decode.cuedsucc_temp2+out_decode.cuedsucc_temp1,[],'g'); hold on;
%     Xmatrix=[out_decode.cuedsucc_temp2-out_decode.cuedsucc_temp1 out_decode.cuedsucc_temp2+out_decode.cuedsucc_temp1]; ylabels=[ones(size(out_decode.cuedsucc_temp1,1),1)];
%     scatter(out_decode.cuedfail_temp2-out_decode.cuedfail_temp1,out_decode.cuedfail_temp2+out_decode.cuedfail_temp1,[],'r');
%     Xmatrix=[Xmatrix; [out_decode.cuedfail_temp2-out_decode.cuedfail_temp1 out_decode.cuedfail_temp2+out_decode.cuedfail_temp1]]; ylabels=[ylabels; 2*ones(size(out_decode.cuedsucc_temp1,1),1)];
%     scatter(out_decode.uncuedsucc_temp2-out_decode.uncuedsucc_temp1,out_decode.uncuedsucc_temp2+out_decode.uncuedsucc_temp1,[],'b');
%     Xmatrix=[Xmatrix; [out_decode.uncuedsucc_temp2-out_decode.uncuedsucc_temp1 out_decode.uncuedsucc_temp2+out_decode.uncuedsucc_temp1]]; ylabels=[ylabels; 3*ones(size(out_decode.cuedsucc_temp1,1),1)];
%     scatter(out_decode.uncuedfail_temp2-out_decode.uncuedfail_temp1,out_decode.uncuedfail_temp2+out_decode.uncuedfail_temp1,[],'y');
%     Xmatrix=[Xmatrix; [out_decode.uncuedfail_temp2-out_decode.uncuedfail_temp1 out_decode.uncuedfail_temp2+out_decode.uncuedfail_temp1]]; ylabels=[ylabels; 4*ones(size(out_decode.cuedsucc_temp1,1),1)];
%     xlabel('gp2 minus gp1'); ylabel('gp2 plus gp1');

    % gp1 v gp2 mapping
    figure();
    scatter(out_decode.cuedsucc_temp1,out_decode.cuedsucc_temp2,[],'g'); hold on;
    Xmatrix=[out_decode.cuedsucc_temp1 out_decode.cuedsucc_temp2]; ylabels=[ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.cuedfail_temp1,out_decode.cuedfail_temp2,[],'r');
    Xmatrix=[Xmatrix; [out_decode.cuedfail_temp1 out_decode.cuedfail_temp2]]; ylabels=[ylabels; 2*ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.uncuedsucc_temp1,out_decode.uncuedsucc_temp2,[],'b');
    Xmatrix=[Xmatrix; [out_decode.uncuedsucc_temp1 out_decode.uncuedsucc_temp2]]; ylabels=[ylabels; 3*ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.uncuedfail_temp1,out_decode.uncuedfail_temp2,[],'y');
    Xmatrix=[Xmatrix; [out_decode.uncuedfail_temp1 out_decode.uncuedfail_temp2]]; ylabels=[ylabels; 4*ones(size(out_decode.cuedsucc_temp1,1),1)];
    xlabel('gp1'); ylabel('gp2');


%     ldaModel=fitcdiscr(Xmatrix,ylabels);
%     predictedY=predict(ldaModel,Xmatrix);
%     accuracy=sum(predictedY==ylabels)/length(ylabels);
%     disp(['Accuracy of LDA on training set: ', num2str(accuracy * 100), '%']);
%     all_accs(unitsCounter,countRuns)=accuracy * 100;
%     close all;

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