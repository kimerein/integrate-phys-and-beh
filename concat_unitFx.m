function unitFx=concat_unitFx(expt_dir)

ls=dir(expt_dir);
unitFx=[];
for i=1:length(ls)
    thisname=ls(i).name;
    if ~isempty(strfind(thisname,'unitFx'))
        a=load([expt_dir '\' thisname]);
        currUnitFx=a.unitFx;
        if isempty(unitFx)
            unitFx=currUnitFx;
        else
            unitFx(length(unitFx)+1:length(unitFx)+1+length(currUnitFx)-1)=currUnitFx;
        end
    end
end