function concat_unitFx(expt_dir)

ls=dir(expt_dir);
j=1;
unitFx=[];
for i=1:length(ls)
    thisname=ls(i).name;
    if ~isempty(regexp(thisname,'unitFx'))
        a=load([expt_dir '\' thisname]);
        currUnitFx=a.unitFx;
        j=j+1;
    end
end

if isempty(tbt)
    disp('No tbt data saved in this directory');
    return
end