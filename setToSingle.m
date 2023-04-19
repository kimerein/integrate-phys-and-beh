function alltbt=setToSingle(alltbt,fie,thresh)

for i=1:length(fie)
    temp=alltbt.(fie{i});
    temp(temp>=thresh)=1;
    temp(temp<thresh)=0;
    alltbt.(fie{i})=single(temp);
end

end