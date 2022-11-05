function out=getCriteriaForUnitsToPlot(varargin)

persistent crit 

if isempty(varargin)
    if isempty(crit)
        error('failed to initialize persistent variable crit in getCriteriaForUnitsToPlot.m');
    end
    out=crit;
elseif length(varargin)==1
    disp('setting value of persistent variable crit in getCriteriaForUnitsToPlot.m');
    crit=varargin{1};
    out=crit;
end

end