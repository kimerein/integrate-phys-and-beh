function correctedDistributions=load_corrected_distributions(datadir)

ls=dir(datadir);
for i=3:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if thisisdir==1
        continue
    end
    temp=regexp(thisname,'.mat','once');
    if ~isempty(temp)
        a=load([datadir '\' thisname]);
        correctedDistributions.(thisname(1:temp-1))=a.temp;
    end
end