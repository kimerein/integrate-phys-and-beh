function [ex,fname]=cue_alignment_file_exists(filedir, unit, onCh)

dd=dir(filedir);
ex=false;
fname=[];
for i=1:length(dd)
    na=dd(i).name;
    if ~isempty(regexp(na,['unit' num2str(unit) 'onCh' num2str(onCh)]))
        ex=true;
        fname=na;
        break
    end
end

end