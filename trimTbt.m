function tbt=trimTbt(tbt,fieldsOfSize1,toSize1,fieldsOfSize2,toSize2)

f=fieldnames(tbt);
for i=1:length(f)
    temp=tbt.(f{i});
    if size(temp,2)==fieldsOfSize1
        temp=temp(:,1:toSize1);
    elseif size(temp,2)==fieldsOfSize2
        temp=temp(:,1:toSize2);
    end
    tbt.(f{i})=temp;
end

end