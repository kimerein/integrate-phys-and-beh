function all_accs=getErrorBars_onShuffleConsensus_trialbytrial(a,b,cued_success_Response,nRuns)

countunits=10:10:170;
all_accs=nan(length(countunits),nRuns);

for unitsCounter=1:length(countunits)
for countRuns=1:nRuns

    shuffle_consensus=cued_success_Response.consensus_idx;
    f=find(~isnan(shuffle_consensus));
    r=randperm(length(f));
    shuffle_consensus(f)=shuffle_consensus(f(r));

    out_decode=decodeTrialByTrialType(a.fortbytclass,b.fortbytclass,shuffle_consensus,100,150,countunits(unitsCounter),false,false,false,false,false); % nBoots,nUnits,nTrials,withReplacement,addThirdAxis,nanAllZeros,justBoostrapTrials,collapseWithinUnit
    % one mapping
    figure();
    scatter(out_decode.cuedsucc_temp2-out_decode.cuedsucc_temp1,out_decode.cuedsucc_temp2+out_decode.cuedsucc_temp1,[],'g'); hold on;
    Xmatrix=[out_decode.cuedsucc_temp2-out_decode.cuedsucc_temp1 out_decode.cuedsucc_temp2+out_decode.cuedsucc_temp1]; ylabels=[ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.cuedfail_temp2-out_decode.cuedfail_temp1,out_decode.cuedfail_temp2+out_decode.cuedfail_temp1,[],'r');
    Xmatrix=[Xmatrix; [out_decode.cuedfail_temp2-out_decode.cuedfail_temp1 out_decode.cuedfail_temp2+out_decode.cuedfail_temp1]]; ylabels=[ylabels; 2*ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.uncuedsucc_temp2-out_decode.uncuedsucc_temp1,out_decode.uncuedsucc_temp2+out_decode.uncuedsucc_temp1,[],'b');
    Xmatrix=[Xmatrix; [out_decode.uncuedsucc_temp2-out_decode.uncuedsucc_temp1 out_decode.uncuedsucc_temp2+out_decode.uncuedsucc_temp1]]; ylabels=[ylabels; 3*ones(size(out_decode.cuedsucc_temp1,1),1)];
    scatter(out_decode.uncuedfail_temp2-out_decode.uncuedfail_temp1,out_decode.uncuedfail_temp2+out_decode.uncuedfail_temp1,[],'y');
    Xmatrix=[Xmatrix; [out_decode.uncuedfail_temp2-out_decode.uncuedfail_temp1 out_decode.uncuedfail_temp2+out_decode.uncuedfail_temp1]]; ylabels=[ylabels; 4*ones(size(out_decode.cuedsucc_temp1,1),1)];
    xlabel('gp2 minus gp1'); ylabel('gp2 plus gp1');
    ldaModel=fitcdiscr(Xmatrix,ylabels);
    predictedY=predict(ldaModel,Xmatrix);
    accuracy=sum(predictedY==ylabels)/length(ylabels);
    disp(['Accuracy of LDA on training set: ', num2str(accuracy * 100), '%']);
    all_accs(unitsCounter,countRuns)=accuracy * 100;
    close all;

end
end

end