function [parmhat]=fitAndPlot(varargin)

if length(varargin)<5
    disp('Need more inputs to fitAndPlot');
elseif length(varargin)==5
    hax=varargin{1};
    data=varargin{2};
    fitt=varargin{3};
    binsForHist=varargin{4};
    c=varargin{5};
    parmhat_input=[];
elseif length(varargin)>5
    hax=varargin{1};
    data=varargin{2};
    fitt=varargin{3};
    binsForHist=varargin{4};
    c=varargin{5};
    parmhat_input=varargin{6};
end

axes(hax);

% plot histogram of data
[n,x]=histcounts(data,binsForHist);
[n,x]=cityscape_hist(n,x);
plot(x,n./nansum(n),'Color',c);
hold on;

switch fitt
    case 'lognfit'
        if isempty(parmhat_input)
            [parmhat,parmci]=lognfit(data(~isnan(data)));
            %parmhat=mle(data(~isnan(data)),'distribution','logn');
            % plot fit
            Y=lognpdf(nanmin(binsForHist):0.01:nanmax(binsForHist),parmhat(1),parmhat(2));
            plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y*(binsForHist(2)-binsForHist(1))*0.55,'Color',c);
        else
            if isempty(parmhat_input(1)) || isempty(parmhat_input(2))
                error('pass in all parameters for fit, or pass in no parameters for fit');
            end
            Y=lognpdf(nanmin(binsForHist):0.01:nanmax(binsForHist),parmhat_input(1),parmhat_input(2));
            plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y*(binsForHist(2)-binsForHist(1))*0.6,'Color',c);
        end
    case 'sum_2_lognfit'
        if isempty(parmhat_input)
            error('sum_2_lognfit requires user to pass in parameters for one of the 2 logn terms');
        end
        % fit function 
        % a*logn(step size, sqrt(2)*sqrt(parmhat(2)) + b*logn(0,sqrt(2)*sqrt(parmhat(2))
        % to data
        enoughPointsToFit=1;
        if enoughPointsToFit==1
            x=nanmean([binsForHist(1:end-1); binsForHist(2:end)],1);
            y=histcounts(data(~isnan(data)),binsForHist);
            thisFunc=@(a,b,fitmu,sigmaSq,noChangeMu,x) a*lognpdf(x,fitmu,sigmaSq) + b*lognpdf(x,noChangeMu,sigmaSq);
            ft=fittype(thisFunc,'problem',{'sigmaSq','noChangeMu'},'independent','x','dependent','y','coefficients',{'a','b','fitmu'});
            myfit=fit(x',y',ft,'problem',{parmhat_input(2),parmhat_input(1)},'StartPoint',[0.8,1,-0.5],'Lower',[0,0,-1.5],'Upper',[10,10,1.5]);
            disp(myfit);
            xpoints=nanmin(binsForHist):0.01:nanmax(binsForHist);
            fittedy=thisFunc(myfit.a,myfit.b,myfit.fitmu,parmhat_input(2),parmhat_input(1),xpoints);
            scaleFac=(nanmax(y./nansum(y))/2)./nanmax(fittedy);
            scaleFac=scaleFac*0.85;
            plot(xpoints,fittedy*scaleFac,'Color',c);
            plot(xpoints,myfit.a*lognpdf(xpoints,myfit.fitmu,parmhat_input(2)).*scaleFac,'Color',[0.5 0.5 0.5]);
            plot(xpoints,myfit.b*lognpdf(xpoints,parmhat_input(1),parmhat_input(2)).*scaleFac,'Color',[0.5 0.5 0.5]);
            disp(['fit mean is ' num2str(myfit.fitmu)]);
        else
            Y=lognpdf(nanmin(binsForHist):0.01:nanmax(binsForHist),parmhat_input(1),parmhat_input(2));
            plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y*(binsForHist(2)-binsForHist(1))*0.6,'Color',c);
            Y2=lognpdf(nanmin(binsForHist):0.01:nanmax(binsForHist),-0.5,parmhat_input(2));
            plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y2*(binsForHist(2)-binsForHist(1))*0.3,'Color',c);
            plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y*(binsForHist(2)-binsForHist(1))*0.6+Y2*(binsForHist(2)-binsForHist(1))*0.3,'Color','b');
        end
    otherwise
        disp('do not recognize this fittype');
        return
end


