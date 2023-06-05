function plotReachTrajectoriesOnTop(noLEDx,noLEDy,noLEDz,LEDx,LEDy,LEDz,reachTrajTimes)

plotSE=true;
ds=10;

figure(); 
colorsUpTo=size(LEDx,2);
cmap=colormap(cool(colorsUpTo));
scatter3(nanmean(LEDx,1),nanmean(LEDy,1),nanmean(LEDz,1),30,cmap,'LineWidth',1);
hold on;
if plotSE==true
    ds_LEDx=downSampMatrix(LEDx,ds);
    ds_LEDy=downSampMatrix(LEDy,ds);
    ds_LEDz=downSampMatrix(LEDz,ds);
    for i=1:size(ds_LEDx,2)
        plot3([nanmean(ds_LEDx(:,i),1)-nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1)) nanmean(ds_LEDx(:,i),1)+nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1))],...
            [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
            [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
        plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
            [nanmean(ds_LEDy(:,i),1)-nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1)) nanmean(ds_LEDy(:,i),1)+nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1))],...
            [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
        plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
            [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
            [nanmean(ds_LEDz(:,i),1)-nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1)) nanmean(ds_LEDz(:,i),1)+nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1))],'Color','k');
    end
end


scatter3(nanmean(LEDx,1),nanmean(LEDy,1),nanmean(LEDz,1),10,'filled','MarkerFaceColor','r');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Average');
 
colorsUpTo=size(noLEDx,2);
cmap=colormap(cool(colorsUpTo));
scatter3(nanmean(noLEDx,1),nanmean(noLEDy,1),nanmean(noLEDz,1),30,cmap);
if plotSE==true

end


end