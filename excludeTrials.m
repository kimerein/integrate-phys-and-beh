function [alltbt,metadata,out]=excludeTrials(alltbt,metadata,out,metadata_nanField)

excludeTheseTrials=isnan(metadata.(metadata_nanField));
disp('excluding sess ids');
disp(unique(metadata.sessid(excludeTheseTrials)));

if any(excludeTheseTrials==1)==0
    return
end

f=fieldnames(alltbt);
for i=1:length(f)
    temp=alltbt.(f{i});
    alltbt.(f{i})=temp(~excludeTheseTrials,:);   
end
f=fieldnames(out);
for i=1:length(f)
    temp=out.(f{i});
    out.(f{i})=temp(~excludeTheseTrials);
end
f=fieldnames(metadata);
for i=1:length(f)
    temp=metadata.(f{i});
    metadata.(f{i})=temp(~excludeTheseTrials);
end

% fix sessid
sessid_un=unique(metadata.sessid);
for i=1:length(sessid_un)
    metadata.sessid(metadata.sessid==sessid_un(i))=i;
end

end