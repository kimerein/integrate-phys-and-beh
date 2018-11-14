function [parmhat]=fitAndPlot(hax,data,fittype,binsForHist,c)

axes(hax);

% plot histogram of data
[n,x]=histcounts(data,binsForHist);
[n,x]=cityscape_hist(n,x);
plot(x,n./nansum(n),'Color',c);
hold on;

switch fittype
    case 'lognfit'
        [parmhat,parmci]=lognfit(data(~isnan(data)));
        %parmhat=mle(data(~isnan(data)),'distribution','logn');
        % plot fit
        Y=lognpdf(nanmin(binsForHist):0.01:nanmax(binsForHist),parmhat(1),parmhat(2));
        plot(nanmin(binsForHist):0.01:nanmax(binsForHist),Y*(binsForHist(2)-binsForHist(1))*0.6,'Color',c);
        pause
    otherwise
        disp('do not recognize this fittype');
        return
end


