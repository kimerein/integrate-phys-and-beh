function [all_glm_coef,unitnames,fromWhichSess,all_glm_pva]=readInGLMCoefs(datadir,settings)

maxUnitsPerSess=settings.maxUnitsPerSess;
if iscell(datadir)
    dd=datadir;
    s=dd{1};
else
    dd=1;
    s=datadir;
end

% Is this matlab glm or python glm?
r=regexp(s,'matglm', 'once');
if ~isempty(r)
    isMatGLM=true;
else
    isMatGLM=false;
    r=regexp(s,'forglm', 'once');
    if isempty(r)
        error('directories need to point to matglm or forglm folders in readInGLMCoefs.m');
    end
end

units_count=1;
unitnames=cell(maxUnitsPerSess*length(dd),1);
fromWhichSess=nan(maxUnitsPerSess*length(dd),1);
unitscount_coef=1;
all_glm_coef=[];
all_glm_pva=[];
for j=1:length(dd)
    if iscell(dd)
        datadir=dd{j};
    end
    if isMatGLM==false
        if exist([datadir sep 'unitnames.mat'],'file')
            a=load([datadir sep 'unitnames.mat']);
            una=a.unitnames;
            unitnames(units_count:units_count+length(una)-1)=una;
            fromWhichSess(units_count:units_count+length(una)-1)=j;
            units_count=units_count+length(una);
            datadir=[datadir sep 'output'];
        else
            continue
        end
    end
    disp(['reading in from ' datadir]);
    ls=dir(datadir);
    for i=3:length(ls)
        if strcmp(ls(i).name,'r2scores')
            continue
        end
        switch isMatGLM
            case true
                if ~contains(ls(i).name,'coef')
                    if contains(ls(i).name,'unitnames')
                        a=load([ls(i).folder sep ls(i).name]);
                        una=a.unitnames;
                        unitnames(units_count:units_count+length(una)-1)=una;
                        fromWhichSess(units_count:units_count+length(una)-1)=j;
                        units_count=units_count+length(una);
                    end
                    continue
                end 
                a=load([ls(i).folder sep ls(i).name]);
                glm_coef=a.coef;
                rp=regexp(ls(i).name,'coef.mat');
                a=load([ls(i).folder sep ls(i).name(1:rp-1) 'p.mat']);
                glm_pva=a.pva;
            case false
                if contains(ls(i).name,'metadata')
                    continue
                end
                a=load([ls(i).folder sep ls(i).name]);
                glm_coef=a.glm_coef;
                glm_pva=nan(size(glm_coef));
        end
        if isempty(all_glm_coef)
            all_glm_coef=nan(maxUnitsPerSess*length(dd),length(glm_coef));
        end
        all_glm_coef(unitscount_coef,:)=glm_coef;
        all_glm_pva(unitscount_coef,:)=glm_pva;
        unitscount_coef=unitscount_coef+1;
    end
end
unitscount_coef=unitscount_coef-1;
units_count=units_count-1;
all_glm_coef=all_glm_coef(1:unitscount_coef,:);
unitnames=unitnames(1:unitscount_coef);
fromWhichSess=fromWhichSess(1:unitscount_coef);
