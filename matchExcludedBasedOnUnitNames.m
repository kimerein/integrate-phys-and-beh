function cuedReach_Response=matchExcludedBasedOnUnitNames(unitbyunitnames_cuedreach,unitbyunit_names,cuedReach_Response)

% index unitbyunitnames_cuedreach into unitbyunit_names

indexCuedreachcellsIntoUnitNames=getNamesIndexIntoNamesList(unitbyunitnames_cuedreach.names,unitbyunit_names);
f=find(unitbyunit_names.excluded==0); indsIntoExcluded=f(indexCuedreachcellsIntoUnitNames); 
cuedReach_Response.excluded=zeros(size(unitbyunit_names.excluded)); cuedReach_Response.excluded(~ismember(1:length(cuedReach_Response.excluded),indsIntoExcluded))=1;

end