function plotReachTrajectoriesOnTop(noLEDx,noLEDy,noLEDz,LEDx,LEDy,LEDz,reachTrajTimes)

plotSE=false; %true;
doOtherPoints=true;
flipZ=false; %true;
ds=1; %5;
startAtSamePoint=false; %true;

if startAtSamePoint==true
    noLEDx=noLEDx-repmat(nanmean(nanmean(noLEDx(:,1:500),2),1),size(noLEDx,1),size(noLEDx,2));
    noLEDy=noLEDy-repmat(nanmean(nanmean(noLEDy(:,1:500),2),1),size(noLEDy,1),size(noLEDy,2));
    noLEDz=noLEDz-repmat(nanmean(nanmean(noLEDz(:,1:500),2),1),size(noLEDz,1),size(noLEDz,2));

    LEDx=LEDx-repmat(nanmean(nanmean(LEDx(:,1:500),2),1),size(LEDx,1),size(LEDx,2));
    LEDy=LEDy-repmat(nanmean(nanmean(LEDy(:,1:500),2),1),size(LEDy,1),size(LEDy,2));
    LEDz=LEDz-repmat(nanmean(nanmean(LEDz(:,1:500),2),1),size(LEDz,1),size(LEDz,2));
end

if flipZ==true
    noLEDz=-noLEDz;
    LEDz=-LEDz;
end

% first show all
figure();
togX=[noLEDx; LEDx];
togY=[noLEDy; LEDy];
togZ=[noLEDz; LEDz];
togLabels=[zeros(size(noLEDx,1),1); ones(size(LEDx,1),1)];
rearr=randperm(size(togX,1));
togX=togX(rearr,:);
togY=togY(rearr,:);
togZ=togZ(rearr,:);
togLabels=togLabels(rearr);
for i=1:size(togLabels,1)
    if togLabels(i)==0
        scatter3(togX(i,:),togY(i,:),togZ(i,:),10,'k');
    else
        scatter3(togX(i,:),togY(i,:),togZ(i,:),10,'r');
    end
    hold on; 
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('All');
view(33.6808,-22.8782);

% then show summary
figure(); 
colorsUpTo=size(noLEDx,2);
cmap=colormap(cool(colorsUpTo));
% scatter3(nanmean(noLEDx,1),nanmean(noLEDy,1),nanmean(noLEDz,1),30,cmap);
if plotSE==true
    ds_LEDx=downSampMatrix(noLEDx,ds);
    ds_LEDy=downSampMatrix(noLEDy,ds);
    ds_LEDz=downSampMatrix(noLEDz,ds);
    for i=10:size(ds_LEDx,2)-80
        if i==10
            lastLoc=[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)];
        else
            if euclDistance(lastLoc,[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)])<1
                lastLoc=[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)];
                continue
            end
        end
        filledOval(nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1)),nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1)),-nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1)),...
            0.1,nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1),'k'); hold on;
%         plot3([nanmean(ds_LEDx(:,i),1)-nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1)) nanmean(ds_LEDx(:,i),1)+nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1))],...
%             [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
%             [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
%         plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
%             [nanmean(ds_LEDy(:,i),1)-nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1)) nanmean(ds_LEDy(:,i),1)+nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1))],...
%             [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
%         plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
%             [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
%             [nanmean(ds_LEDz(:,i),1)-nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1)) nanmean(ds_LEDz(:,i),1)+nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1))],'Color','k');
    end
end
if doOtherPoints==true
    scatter3(nanmean(noLEDx,1),nanmean(noLEDy,1),nanmean(noLEDz,1),60,cmap); 
    hold on;
    scatter3(nanmean(noLEDx,1),nanmean(noLEDy,1),nanmean(noLEDz,1),30,'filled','MarkerFaceColor','k');
end

% Then red
colorsUpTo=size(LEDx,2);
cmap=colormap(cool(colorsUpTo));
if plotSE==true
    ds_LEDx=downSampMatrix(LEDx,ds);
    ds_LEDy=downSampMatrix(LEDy,ds);
    ds_LEDz=downSampMatrix(LEDz,ds);
    for i=10:size(ds_LEDx,2)-80
        if i==10
            lastLoc=[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)];
        else
            if euclDistance(lastLoc,[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)])<1
                lastLoc=[nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1)];
                continue
            end
        end
        filledOval(nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1)),nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1)),-nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1)),...
            0.1,nanmean(ds_LEDx(:,i),1),nanmean(ds_LEDy(:,i),1),nanmean(ds_LEDz(:,i),1),'r'); hold on;
%         plot3([nanmean(ds_LEDx(:,i),1)-nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1)) nanmean(ds_LEDx(:,i),1)+nanstd(ds_LEDx(:,i),[],1)./sqrt(size(ds_LEDx(:,i),1))],...
%             [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
%             [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
%         plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
%             [nanmean(ds_LEDy(:,i),1)-nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1)) nanmean(ds_LEDy(:,i),1)+nanstd(ds_LEDy(:,i),[],1)./sqrt(size(ds_LEDy(:,i),1))],...
%             [nanmean(ds_LEDz(:,i),1) nanmean(ds_LEDz(:,i),1)],'Color','k');
%         plot3([nanmean(ds_LEDx(:,i),1) nanmean(ds_LEDx(:,i),1)],...
%             [nanmean(ds_LEDy(:,i),1) nanmean(ds_LEDy(:,i),1)],...
%             [nanmean(ds_LEDz(:,i),1)-nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1)) nanmean(ds_LEDz(:,i),1)+nanstd(ds_LEDz(:,i),[],1)./sqrt(size(ds_LEDz(:,i),1))],'Color','k');
    end
end
if doOtherPoints==true
    scatter3(nanmean(LEDx,1),nanmean(LEDy,1),nanmean(LEDz,1),60,cmap,'LineWidth',1);
    hold on;
    scatter3(nanmean(LEDx,1),nanmean(LEDy,1),nanmean(LEDz,1),30,'filled','MarkerFaceColor','r');
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Average');

view(33.6808,-22.8782);
[caz,cel]=view();
disp('Current azimuth');
disp(caz);
disp('Current elevation');
disp(cel);
% view(-116.4279,56.0009);
% view(-162.2125,9.8783);
% view(-330.8644,17.7600);

end

function d=euclDistance(firstpoint,secondpoint)

d=sqrt((firstpoint(1)-secondpoint(1)).^2 + (firstpoint(2)-secondpoint(2)).^2 + (firstpoint(3)-secondpoint(3)).^2);

end

function filledOval(radiusX,radiusY,radiusZ,alph,Xoffset,Yoffset,Zoffset,c)

% Specify the dimensions of the spherical shape
% radius = 5; % Radius of the sphere

% Create a grid of points for the sphere
theta = linspace(0, 2*pi, 100);
phi = linspace(0, pi, 50);
[theta, phi] = meshgrid(theta, phi);
x = radiusX * sin(phi) .* cos(theta);
y = radiusY * sin(phi) .* sin(theta);
z = radiusZ * cos(phi);
x = x+Xoffset;
y = y+Yoffset;
z = z+Zoffset;

% Create a figure and plot the spherical shape in 3D
% figure;
surf(x, y, z, 'FaceColor', c, 'EdgeColor', 'none');
alpha(alph);

% % Set the axis labels
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% 
% % Set the aspect ratio of the plot
% daspect([1, 1, 1]);
% 
% % Add a title to the plot
% title('Filled Spherical Shape in 3D');
% 
% % Set the view angle
% view(3); % Adjust the view angle as desired
% 
% % Adjust the plot limits based on the dimensions of the sphere
% xlim([-radius, radius]);
% ylim([-radius, radius]);
% zlim([-radius, radius]);

end