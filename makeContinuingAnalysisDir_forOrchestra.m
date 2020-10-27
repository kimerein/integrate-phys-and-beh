function outDir=makeContinuingAnalysisDir_forOrchestra(addAtFront,continuingAnalysisDir,cutAtWord)

for i=1:length(continuingAnalysisDir)
    temp=continuingAnalysisDir{i};
    f=regexp(temp,cutAtWord,'ONCE');
    outDir{i}=[addAtFront temp(f+length(cutAtWord)+1:end)];
end

