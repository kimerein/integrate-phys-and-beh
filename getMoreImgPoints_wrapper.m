function [side_image_points,under_image_points]=getMoreImgPoints_wrapper(varargin)

framesPerChunk=1500;
ds=10;
maxPoints=20;

% Get file name of video with reaches
filename=varargin{1};
discardFirstNFrames=varargin{2};

% Read beginning of movie and discard unused frames
videoFReader = vision.VideoFileReader(filename,'PlayCount',1);

n=discardFirstNFrames; % How many frames to read initially
if discardFirstNFrames>0
    for i=1:n
        [~,EOF]=step(videoFReader);
        if EOF==true
            break
        end
    end
end

for i=1:framesPerChunk % How many frames to read initially
    [frame,EOF]=step(videoFReader);
    if EOF==true
        allframes=allframes(:,:,:,1:i-1);
        break
    end
    if i==1
        allframes=nan([size(frame,1) size(frame,2) size(frame,3) framesPerChunk]);
    end
    allframes(:,:,:,i)=frame;
end

% X is across the frame horizontally
% Y is down the frame vertically
% X and Y zeros are in upper left corner

% Have user find vertices for every ds_th frame
% More image points
disp('Click side-view then under-view. These points all need to match!');
all_vs=nan(length(1:ds:framesPerChunk),6);
j=1;
for i=1:size(allframes,4)
    if mod(i,ds)==0
        all_vs(j,:)=getVertices(allframes(:,:,:,i),'paw positions','paw pos');
        j=j+1;
        if j>maxPoints
            all_vs=all_vs(1:maxPoints,:);
            break
        end
    end
end

% Convert object points to real points in 3D space in units of millimeters
% Keep track of relative image points for side view and under view
% X is paw forward-back from perch to pellet (+ is closer to perch and away
% from pellet)
% Y is side-to-side paw movements (+ is from reaching paw, i.e., right paw, resting position
% inward toward nose)
% Z is up-and-down paw movements (+ is downward further and further below
% the presenter wheel)
side_image_points=[];
under_image_points=[];
% Other points
for i=1:size(all_vs,1)
    side_image_points=[side_image_points; [all_vs(i,1), all_vs(i,2)]];
    under_image_points=[under_image_points; [all_vs(i,3), all_vs(i,4)]];
end

% Change mirrors to "virtual cameras"
% Need to flip both mirror points horizontally
side_image_points=flipImagePoints(side_image_points,size(frame,2),size(frame,1),'horizontal');
under_image_points=flipImagePoints(under_image_points,size(frame,2),size(frame,1),'horizontal');

end

function objPoint=convertImagePointToObjectPoint(side_x,side_y,under_x,under_y,reference)

% get position in X-Y plane from bottom view
% project vector from origin to [under_x,under_y] onto X axis, then Y axis,
% get magnitudes
x_pos=reference.under_scaleTomm*normKeepSign(vecProject([under_x-reference.under_origin_x,under_y-reference.under_origin_y],reference.under_X_vec,true),1);
y_pos=reference.under_scaleTomm*normKeepSign(vecProject([under_x-reference.under_origin_x,under_y-reference.under_origin_y],reference.under_Y_vec,true),2);

% get position in X-Z plane from side view
% project vector from origin to [side_x,side_y] onto X axis, then Z axis,
% get magnitudes
x_pos2=reference.side_scaleTomm*normKeepSign(vecProject([side_x-reference.side_origin_x,side_y-reference.side_origin_y],reference.side_X_vec,true),1);
z_pos=reference.side_scaleTomm*normKeepSign(vecProject([side_x-reference.side_origin_x,side_y-reference.side_origin_y],reference.side_Z_vec,true),2);
% after conversion into real space, x_pos and x_pos2 should be approximately the same
disp([x_pos x_pos2]);

objPoint=[x_pos,y_pos,z_pos];

end

function vecnorm=normKeepSign(vec,signFromInd)

if signFromInd==1
    if vec(1)<0
        vecnorm=-norm(vec);
    else
        vecnorm=norm(vec);
    end
elseif signFromInd==2
    if vec(2)<0
        vecnorm=-norm(vec);
    else
        vecnorm=norm(vec);
    end
else
    error('normKeepSign assumes a 2D vector');
end
    
end

function C=vecProject(A,B,makeB_unitVector)

% A,B are row vectors
if makeB_unitVector==true
    B=B./norm(B);
end
C=(sum(A.*B)/norm(B)^2)*B;

end

function showPointsOnImage(frame,side_image_points,under_image_points,object_points)

if ~isempty(frame)
    figure();
    imagesc(frame);
    colormap gray
end
hold on;
cmap=colormap('cool');
stepcmap=floor(size(cmap,1)./size(side_image_points,1));
indsIntoCmap=1:stepcmap:size(cmap,1);
indsIntoCmap=indsIntoCmap(1:size(side_image_points,1));
scatter(side_image_points(:,1),side_image_points(:,2),[],cmap(indsIntoCmap,:));
scatter(under_image_points(:,1),under_image_points(:,2),[],cmap(indsIntoCmap,:));

figure();
% need to flip Z axis for this display
scatter3(object_points(:,1),-object_points(:,2),-object_points(:,3),[],cmap(indsIntoCmap,:),'filled');
daspect([1 1 1]);

end

function showReferencesOnImage(frame,origin_x,origin_y,slope_Xaxis,slope_Yaxis,b_Xaxis,b_Yaxis)

xChange=100;
if ~isempty(frame)
    figure();
    imagesc(frame);
    colormap gray
end
hold on;
line([origin_x origin_x+xChange],[slope_Xaxis*origin_x+b_Xaxis slope_Xaxis*(origin_x+xChange)+b_Xaxis],'Color','g');
line([origin_x origin_x+xChange],[slope_Yaxis*origin_x+b_Yaxis slope_Yaxis*(origin_x+xChange)+b_Yaxis],'Color','y');

end

function [x0,y0] = intersectPoints(m1,m2,b1,b2)
% Insersection point of two lines with known slope and constant
% parameters.
% [x0 y0] = intersectPoints(m1,m2,b1,b1)
% where m's are slope, and b's are constants.
x0 = (b2-b1)/(m1-m2); % find the x point
y0 = m1*x0+b1;
end

function perpangle=getPerpAngle(corner_points) 

perpangle=getVectorAngle(corner_points([3 4 1 2]),corner_points([3 4 5 6]));

end

function ang=getVectorAngle(vec_points1,vec_points2)

% vec_points should be x_1, y_1, x_2, y_2
u=[vec_points1(3)-vec_points1(1) vec_points1(4)-vec_points1(2)];
v=[vec_points2(3)-vec_points2(1) vec_points2(4)-vec_points2(2)];
CosTheta=max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
ang=real(acosd(CosTheta));

end

function vs=getVertices(frame,nameOfFig,whichInput)

global vertexViews
global continueAnalysis

makeCheckerboardOnVid(frame,nameOfFig,whichInput);
pause;
temp=vertexViews;
temp=temp';
vs=temp(1:end); % x_v1, y_v1, x_v2, y_v2, x_v3, y_v3

if continueAnalysis==1
    disp('Succesfully defined');
else
    disp('Failed to define');
end

end