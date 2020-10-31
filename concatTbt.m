function [alltbt,out,metadata]=concatTbt(tbt1,tbt2,out1,out2,metadata1,metadata2)

f=fieldnames(tbt1);
for i=1:length(f)
    alltbt.(f{i})=[tbt1.(f{i}); tbt2.(f{i})];
end

f=fieldnames(out1);
for i=1:length(f)
    out.(f{i})=[out1.(f{i}); out2.(f{i})];
end

f=fieldnames(metadata1);
for i=1:length(f)
    metadata.(f{i})=[metadata1.(f{i}); metadata2.(f{i})];
end