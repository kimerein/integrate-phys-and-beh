function data=smoothMatrix(data,smooby)

for i=1:size(data,1)
    data(i,:)=smooth(data(i,:),smooby);
end

end