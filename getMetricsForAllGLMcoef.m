function metrics=getMetricsForAllGLMcoef(glm_coef,glm_intercept,feature_names,timestep,nShiftsBefore)

metrics.preCueAmp=nan(size(glm_coef,1),1);
metrics.postCueAmp_over1sec=nan(size(glm_coef,1),1);
metrics.postCueAmp_at1sec=nan(size(glm_coef,1),1);
metrics.allDrop_sustained=nan(size(glm_coef,1),1);
metrics.cXdrop_sustained=nan(size(glm_coef,1),1);
metrics.allSucc_sustained=nan(size(glm_coef,1),1);
metrics.cXsucc_sustained=nan(size(glm_coef,1),1);
metrics.allMiss_sustained=nan(size(glm_coef,1),1);
metrics.cXmiss_sustained=nan(size(glm_coef,1),1);
metrics.allFail_sustained=nan(size(glm_coef,1),1);
metrics.cXfail_sustained=nan(size(glm_coef,1),1);
for i=1:size(glm_coef,1)
    [~,~,met]=plotGLMcoef(glm_coef(i,:),glm_intercept,feature_names,timestep,nShiftsBefore,'mean',true,[]);
    metrics.preCueAmp(i)=met.preCueAmp;
    metrics.postCueAmp_over1sec(i)=met.postCueAmp_over1sec;
    metrics.postCueAmp_at1sec(i)=met.postCueAmp_at1sec;
    metrics.allDrop_sustained(i)=met.allDrop_sustained;
    metrics.cXdrop_sustained(i)=met.cXdrop_sustained;
    metrics.allSucc_sustained(i)=met.allSucc_sustained;
    metrics.cXsucc_sustained(i)=met.cXsucc_sustained;
    metrics.allMiss_sustained(i)=met.allMiss_sustained;
    metrics.cXmiss_sustained(i)=met.cXmiss_sustained;
    metrics.allFail_sustained(i)=(met.allDrop_sustained+met.allMiss_sustained)/2;
    metrics.cXfail_sustained(i)=(met.cXsucc_sustained+met.cXmiss_sustained)/2;
end

end