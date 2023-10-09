function testWhereCyanPurpleDifference(x,y,label,cmap,ratiocutoffs)

xrange(1)=nanmin(x);
xrange(2)=nanmax(x);
yrange(1)=nanmin(y);
yrange(2)=nanmax(y);
% xbins=xrange(1):0.015:xrange(2)+0.0001;
% ybins=-1.05:0.15:1.05; %yrange(1):0.15:yrange(2);
xbins=-0.066:0.09:xrange(2); %xrange(1):0.015:xrange(2)+0.0001;
xbins(end)=xbins(end)+0.0001;
ybins=-1:1:1; %yrange(1):0.15:yrange(2);
ybins(end)=ybins(end)+0.0001;

label1_for_bin=nan(length(xbins)-1,length(ybins)-1);
label2_for_bin=nan(length(xbins)-1,length(ybins)-1);
for i=1:length(xbins)-1
    for j=1:length(ybins)-1
        label1count=nansum(x(label==1)>=xbins(i) & x(label==1)<xbins(i+1) & y(label==1)>=ybins(j) & y(label==1)<ybins(j+1));
        label2count=nansum(x(label==2)>=xbins(i) & x(label==2)<xbins(i+1) & y(label==2)>=ybins(j) & y(label==2)<ybins(j+1));
        label1_for_bin(i,j)=label1count;
        label2_for_bin(i,j)=label2count;
    end
end

label1_for_bin(label1_for_bin==0)=0.01;
label2_for_bin(label2_for_bin==0)=0.01;
ratio_1to2=label1_for_bin./label2_for_bin;
ratio_1to2(label1_for_bin+label2_for_bin<2.03)=nan; 

figure(); 
s=scatter(x(label==1),y(label==1),100,cmap(1,:),'filled'); 
s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; hold on;
s=scatter(x(label==2),y(label==2),100,cmap(2,:),'filled'); 
s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;

figure(); 
for i=1:length(xbins)-1
    for j=1:length(ybins)-1
        if ratio_1to2(i,j)<ratiocutoffs(1) % 2 preferred
            patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(2,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on;
        elseif ratio_1to2(i,j)>ratiocutoffs(2) % 1 preferred
            patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(1,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on;
        end
    end
end
s=scatter(x(label==2),y(label==2),100,cmap(2,:),'filled'); 
s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; hold on;
s=scatter(x(label==1),y(label==1),100,cmap(1,:),'filled'); 
s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; 

end