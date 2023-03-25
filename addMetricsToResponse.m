function cued_success_Response=addMetricsToResponse(cued_success_Response,metrics,all_glm_coef,indexGLMcellsIntoUnitNames,whichGLMinds)

f=fieldnames(metrics);
for i=1:length(f)
    cued_success_Response=addToR(cued_success_Response,metrics.(f{i}),f{i},indexGLMcellsIntoUnitNames);
end

for i=1:length(whichGLMinds)
    cued_success_Response=addToR(cued_success_Response,all_glm_coef(:,whichGLMinds(i)),['glmcoef_index' num2str(whichGLMinds(i))],indexGLMcellsIntoUnitNames);
end

end

function cued_success_Response=addToR(cued_success_Response,thingtoadd,fieldcalled,indexGLMcellsIntoUnitNames)

preCueAmp=nan(size(cued_success_Response.unitbyunit_x,1),1);
preCueAmp(indexGLMcellsIntoUnitNames(~isnan(indexGLMcellsIntoUnitNames)))=thingtoadd(~isnan(indexGLMcellsIntoUnitNames));
cued_success_Response.(fieldcalled)=preCueAmp;

end