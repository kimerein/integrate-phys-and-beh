function outDir=makeContinuingAnalysisDir_forOrchestra(addAtFront,continuingAnalysisDir,cutAtWord,goingToO2)

for i=1:length(continuingAnalysisDir)
    temp=continuingAnalysisDir{i};
    f=regexp(temp,cutAtWord,'ONCE');
    outDir{i}=[addAtFront temp(f+length(cutAtWord)+1:end)];
end

% If going to O2, replace all w /
if goingToO2==true
    for i=1:length(outDir)
        temp=outDir{i};
        f=regexp(temp,'\');
        temp(f)='/';
        outDir{i}=temp;
    end
end

