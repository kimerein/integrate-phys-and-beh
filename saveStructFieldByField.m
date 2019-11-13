function saveStructFieldByField(s,saveDir)

f=fieldnames(s); 
for i=1:length(f)
    disp(i);
    temp=s.(f{i});
    save([saveDir '\' f{i} '.mat'],'temp');
end