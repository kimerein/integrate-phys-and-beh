function [newtbt,newout,newmetadata]=filtTbt(alltbt,out,sortField,range_values,metadata)

fieldForSort=out.(sortField);

f=fieldnames(alltbt);
for i=1:length(f)
    temp=alltbt.(f{i});
    if size(temp,1)~=length(fieldForSort)
        continue
    end
    if islogical(range_values)
        newtbt.(f{i})=temp(fieldForSort==true,:);
    elseif length(range_values)>1
        newtbt.(f{i})=temp(fieldForSort>=range_values(1) & fieldForSort<=range_values(2),:);
    end
end

f=fieldnames(out);
for i=1:length(f)
    temp=out.(f{i});
    if islogical(range_values)
        newout.(f{i})=temp(fieldForSort==true);
    elseif length(range_values)>1
        newout.(f{i})=temp(fieldForSort>=range_values(1) & fieldForSort<=range_values(2));
    end
end

f=fieldnames(metadata);
for i=1:length(f)
    temp=metadata.(f{i});
    if islogical(range_values)
        newmetadata.(f{i})=temp(fieldForSort==true);
    elseif length(range_values)>1
        newmetadata.(f{i})=temp(fieldForSort>=range_values(1) & fieldForSort<=range_values(2));
    end
end
