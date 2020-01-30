function out=loadStructFieldByField(datadir)

doingcellarray=false;

ls=dir(datadir);
for i=3:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if thisisdir==1
        continue
    end
    temp=regexp(thisname,'.mat','once');
    if ~isempty(regexp(thisname,'scell','once'))
        doingcellarray=true;
        if ~isempty(temp)
            a=load([datadir '\' thisname]);
            indtostr=regexp(thisname,'scell','once');
            cellind=str2num(thisname(indtostr+5:temp-1));
            outcell(cellind)=a.temp;
        end
    elseif ~isempty(regexp(thisname,'cell','once'))
        doingcellarray=true;
        if ~isempty(temp)
            a=load([datadir '\' thisname]);
            indtostr=regexp(thisname,'cell','once');
            cellind=str2num(thisname(indtostr+4:temp-1));
            outcell{cellind}=a.temp;
        end        
    else
        if ~isempty(temp)
            a=load([datadir '\' thisname]);
            if length(fieldnames(a))==0
            else
                out.(thisname(1:temp-1))=a.temp;
            end
        end
    end
end

if doingcellarray==true
    out=outcell;
end