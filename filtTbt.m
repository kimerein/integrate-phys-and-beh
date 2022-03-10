function [newtbt,newout,newmetadata,takeThese]=filtTbt(alltbt,out,sortField,range_values,metadata,displayProgress)

fieldForSort=out.(sortField);

f=fieldnames(alltbt);
if length(range_values)==1
    error('Need at least two range values for filtTbt');
elseif length(range_values)==2
    disp(['Taking fieldforsort BETWEEN two values']);
    takeThese=fieldForSort>=range_values(1) & fieldForSort<=range_values(2);
elseif length(range_values)>2
    takeThese=ismember(fieldForSort,range_values);
end
for i=1:length(f)
    if displayProgress==true && i~=length(f)
        disp([num2str((i/length(f))*100,2) '% of the way through filtering']);
    end
    temp=alltbt.(f{i});
    if size(temp,1)~=length(fieldForSort)
        continue
    end
    if islogical(range_values)
        newtbt.(f{i})=temp(fieldForSort==true,:);
    elseif length(range_values)>1
        newtbt.(f{i})=temp(takeThese,:);
    end
end

f=fieldnames(out);
for i=1:length(f)
    temp=out.(f{i});
    if islogical(range_values)
        newout.(f{i})=temp(fieldForSort==true);
    elseif length(range_values)>1
        newout.(f{i})=temp(takeThese);
    end
end

if displayProgress==true
    disp('99% of the way through filtering');
end

f=fieldnames(metadata);
for i=1:length(f)
    temp=metadata.(f{i});
    if islogical(range_values)
        newmetadata.(f{i})=temp(fieldForSort==true);
    elseif length(range_values)>1
        newmetadata.(f{i})=temp(takeThese);
    end
end
