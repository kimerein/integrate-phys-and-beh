function saveStructFieldByField(s,saveDir)

if isstruct(s) && length(s)==1
    f=fieldnames(s);
    for i=1:length(f)
        disp(i);
        temp=s.(f{i});
        save([saveDir '\' f{i} '.mat'],'temp');
    end
elseif isstruct(s) && length(s)>1
    for i=1:length(s)
        temp=s(i);
        save([saveDir '\scell' num2str(i) '.mat'],'temp');
    end
elseif iscell(s)
    for i=1:length(s)
        temp=s{i};
        save([saveDir '\cell' num2str(i) '.mat'],'temp');
    end
end