function [is1preferred,is2preferred,xbins,ybins,ratio_1to2,zbins]=testWhereCyanPurpleDifference(varargin)

zbins=[];
if length(varargin)==5
    % 2D boxes
    x=varargin{1};
    y=varargin{2};
    label=varargin{3};
    cmap=varargin{4};
    ratiocutoffs=varargin{5};

    xrange(1)=nanmin(x);
    xrange(2)=nanmax(x);
    yrange(1)=nanmin(y);
    yrange(2)=nanmax(y);
    % xbins=xrange(1):0.015:xrange(2)+0.0001;
    % ybins=-1.05:0.15:1.05; %yrange(1):0.15:yrange(2);
    % xbins=-0.066:0.09:xrange(2); %xrange(1):0.015:xrange(2)+0.0001;
    xbins=-0.08:0.09:xrange(2); %xrange(1):0.015:xrange(2)+0.0001;
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

    is2preferred=zeros(size(ratio_1to2));
    is1preferred=zeros(size(ratio_1to2));
    figure();
    for i=1:length(xbins)-1
        for j=1:length(ybins)-1
            if ratio_1to2(i,j)<ratiocutoffs(1) % 2 preferred
                patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(2,:),'FaceAlpha',0.4); hold on;
                is2preferred(i,j)=1;
            elseif ratio_1to2(i,j)>ratiocutoffs(2) % 1 preferred
                patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(1,:),'FaceAlpha',0.4); hold on;
                is1preferred(i,j)=1;
            end
        end
    end
    s=scatter(x(label==2),y(label==2),100,cmap(2,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; hold on;
    s=scatter(x(label==1),y(label==1),100,cmap(1,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
elseif length(varargin)==6
    % 3D boxes

    x=varargin{1};
    y=varargin{2};
    label=varargin{3};
    cmap=varargin{4};
    ratiocutoffs=varargin{5};
    z=varargin{6};

    xrange(1)=nanmin(x);
    xrange(2)=nanmax(x);
    yrange(1)=nanmin(y);
    yrange(2)=nanmax(y);
    zrange(1)=nanmin(z);
    zrange(2)=nanmax(z);
    % xbins=xrange(1):0.015:xrange(2)+0.0001;
    % ybins=-1.05:0.15:1.05; %yrange(1):0.15:yrange(2);
    % xbins=-0.066:0.09:xrange(2); %xrange(1):0.015:xrange(2)+0.0001;
    
    xbins=xrange(1):0.01:xrange(2);
    xbins(end)=xbins(end)+0.0001;
    ybins=-1:0.2:1; %yrange(1):0.15:yrange(2);
    ybins(end)=ybins(end)+0.0001;
    zbins=-1:0.2:1;
    zbins(end)=zbins(end)+0.0001;

    label1_for_bin=nan(length(xbins)-1,length(ybins)-1,length(zbins)-1);
    label2_for_bin=nan(length(xbins)-1,length(ybins)-1,length(zbins)-1);
    for i=1:length(xbins)-1
        for j=1:length(ybins)-1
            for k=1:length(zbins)-1
                label1count=nansum(x(label==1)>=xbins(i) & x(label==1)<xbins(i+1) & y(label==1)>=ybins(j) & y(label==1)<ybins(j+1) & z(label==1)>=zbins(k) & z(label==1)<zbins(k+1));
                label2count=nansum(x(label==2)>=xbins(i) & x(label==2)<xbins(i+1) & y(label==2)>=ybins(j) & y(label==2)<ybins(j+1) & z(label==2)>=zbins(k) & z(label==2)<zbins(k+1));
                label1_for_bin(i,j,k)=label1count;
                label2_for_bin(i,j,k)=label2count;
            end
        end
    end

    label1_for_bin(label1_for_bin==0)=0.01;
    label2_for_bin(label2_for_bin==0)=0.01;
    ratio_1to2=label1_for_bin./label2_for_bin;
%     ratio_1to2(label1_for_bin+label2_for_bin<2.03)=nan;

    figure();
    s=scatter3(x(label==1),y(label==1),z(label==1),100,cmap(1,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; hold on;
    s=scatter3(x(label==2),y(label==2),z(label==2),100,cmap(2,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;

    is2preferred=zeros(size(ratio_1to2));
    is1preferred=zeros(size(ratio_1to2));
    figure();
    for i=1:length(xbins)-1
        for j=1:length(ybins)-1
            for k=1:length(zbins)-1
                if ratio_1to2(i,j,k)<ratiocutoffs(1) % 2 preferred
%                     patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(2,:),'FaceAlpha',0.4); hold on;
                    plot_cube(xbins(i), ybins(j), zbins(k), xbins(i+1)-xbins(i),ybins(j+1)-ybins(j),zbins(k+1)-zbins(k),cmap(2,:)); hold on;
                    is2preferred(i,j,k)=1;
                elseif ratio_1to2(i,j,k)>ratiocutoffs(2) % 1 preferred
%                     patch([xbins(i) xbins(i) xbins(i+1) xbins(i+1)],[ybins(j) ybins(j+1) ybins(j+1) ybins(j)],cmap(1,:),'FaceAlpha',0.4); hold on;
                    plot_cube(xbins(i), ybins(j), zbins(k), xbins(i+1)-xbins(i),ybins(j+1)-ybins(j),zbins(k+1)-zbins(k),cmap(1,:)); hold on;
                    is1preferred(i,j,k)=1;
                end
            end
        end
    end
    xlabel('X'); ylabel('Y'); zlabel('Z');
%     figure();
    s=scatter3(x(label==2),y(label==2),z(label==2),100,cmap(2,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4; hold on;
    s=scatter3(x(label==1),y(label==1),z(label==1),100,cmap(1,:),'filled');
    s.MarkerFaceAlpha=0.4; s.MarkerEdgeAlpha=0.4;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

end

function plot_cube(x0, y0, z0, ax,ay,az,c)
    % Define the cube vertices
    vertices = [x0, y0, z0;
                x0+ax, y0, z0;
                x0+ax, y0+ay, z0;
                x0, y0+ay, z0;
                x0, y0, z0+az;
                x0+ax, y0, z0+az;
                x0+ax, y0+ay, z0+az;
                x0, y0+ay, z0+az];

    % Define the faces of the cube
    faces = [1 2 3 4;  2 6 7 3; 3 7 8 4; 1 5 8 4; 1 2 6 5; 5 6 7 8];

    % Plot the cube using patch
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', c, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % The 'FaceAlpha' property makes the cube transparent

%     hold on;
% 
%     % Add a 3D scatter point
%     scatter3(0.5*(2*x0+a), 0.5*(2*y0+a), 0.5*(2*z0+a), 'bo'); % this plots a point at the center of the cube for demonstration
    
    % Set plot limits
%     xlim([x0-0.5*a, x0+1.5*a]);
%     ylim([y0-0.5*a, y0+1.5*a]);
%     zlim([z0-0.5*a, z0+1.5*a]);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    view(3); % sets the default 3D view
    hold off;
end
