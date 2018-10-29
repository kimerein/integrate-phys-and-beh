function makeCityscapeAvs(xdata1,ydata1,xdata2,ydata2)

[y1,x1]=cityscape_hist(nanmean(ydata1,1),nanmean(xdata1,1));
[y1_plus_se,~]=cityscape_hist(nanmean(ydata1,1)+nanstd(ydata1,[],1)./sqrt(size(ydata1,1)),nanmean(xdata1,1));
[y1_minus_se,~]=cityscape_hist(nanmean(ydata1,1)-nanstd(ydata1,[],1)./sqrt(size(ydata1,1)),nanmean(xdata1,1));

if ~isempty(xdata2)
    [y2,x2]=cityscape_hist(nanmean(ydata2,1),nanmean(xdata2,1));
    [y2_plus_se,~]=cityscape_hist(nanmean(ydata2,1)+nanstd(ydata2,[],1)./sqrt(size(ydata2,1)),nanmean(xdata2,1));
    [y2_minus_se,~]=cityscape_hist(nanmean(ydata2,1)-nanstd(ydata2,[],1)./sqrt(size(ydata2,1)),nanmean(xdata2,1));
end

figure();
for i=1:length(x1)-1
    h=fill([x1(i) x1(i+1) x1(i+1) x1(i)],[y1_plus_se(i) y1_plus_se(i+1) y1_minus_se(i+1) y1_minus_se(i)],[0.5 0.5 0.5]);
    hold on;
    set(h,'EdgeColor','none');
end
hold on;
plot(x1,y1,'Color','k');

if ~isempty(xdata2)
    for i=1:length(x2)-1
        h=fill([x2(i) x2(i+1) x2(i+1) x2(i)],[y2_plus_se(i) y2_plus_se(i+1) y2_minus_se(i+1) y2_minus_se(i)],[0.9 0.6 0.6]);
        hold on;
        set(h,'EdgeColor','none');
    end
    hold on;
    plot(x2,y2,'Color','r');
end