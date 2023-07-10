function newvec=downSampSum(vec,n)

if size(vec,1)>1
    newvec=zeros(length(1:n:length(vec)),1);
else 
    newvec=zeros(1,length(1:n:length(vec)));
end
stepInds=1:n:length(vec);
for i=1:length(stepInds)
    if i==length(stepInds)
        newvec(i)=sum(vec(stepInds(i):length(vec)),'all','omitnan');
    else
        newvec(i)=sum(vec(stepInds(i):stepInds(i+1)),'all','omitnan');
    end
end