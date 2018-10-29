function data=padAtEnd(data)

for i=1:size(data,1)
    f=find(isnan(data(i,:)),1,'first');
    if f~=1
        data(i,f:end)=data(i,f-1);
    end
end
