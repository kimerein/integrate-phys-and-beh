function [ex,fname]=SU_QC_file_exists(filedir, unit, onCh)

dd=dir(filedir);
ex=false;
fname=[];
for i=1:length(dd)
    na=dd(i).name;
    if length(onCh)==1
        if ~isempty(regexp(na,['unit' num2str(unit) 'onCh' num2str(onCh)]))
            ex=true;
            fname=na;
            break
        end
    else
        % look for name on any of these channels
        for j=1:length(onCh)
            if ~isempty(regexp(na,['unit' num2str(unit) 'onCh' num2str(onCh(j))]))
                ex=true;
                fname=na;
                break
            end
        end
        if ex==true
            break
        end
    end
end

end