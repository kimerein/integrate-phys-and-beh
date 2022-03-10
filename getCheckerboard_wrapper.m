function [object_points,side_image_points,under_image_points,mirror_vs]=getCheckerboard_wrapper(varargin)

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

% Play movie until at good frame for defining axes
fig=implay(allframes,30);
fig.Parent.Position=[100 100 800 800];
pause;
if ~isempty(regexp(version,'2017b','once')) || ~isempty(regexp(version,'2018b','once')) ||  ~isempty(regexp(version,'2020b','once')) || ~isempty(regexp(version,'2021a','once'))
    currentFrameNumber=fig.DataSource.Controls.CurrentFrame;
else
    currentFrameNumber=fig.data.Controls.CurrentFrame;
end

% Instructions to user
continuebutton=questdlg('Choose frame where pellet stopped in front of mouse.','Instructions','Yes','Cancel','Cancel');
switch continuebutton
    case 'Yes'
    case 'Cancel'
        return
end

frame=allframes(:,:,:,currentFrameNumber);

% X is across the frame horizontally
% Y is down the frame vertically
% X and Y zeros are in upper left corner

% Get "camera" (actually mirror) centers
disp('Click side-view mirror center then under-view mirror center.');
mirror_vs=getVertices(frame,'mirror centers','mirror centers');

% Draw perch line onto both views
disp('Click side-view then under-view perch line.');
perch_vs=getVertices(frame,'perch line','perch line');
% Defines Y axis and Y origin
% Get image angle for perch lines
sideview_perchangle=getVectorAngle([0 0 1 0],perch_vs(1:4));
underview_perchangle=getVectorAngle([0 0 1 0],perch_vs(5:8));

if sign((perch_vs(4)-perch_vs(2))/(perch_vs(3)-perch_vs(1)))==-1
    % slope is positive, for viewer of image
    sideview_perchangle=abs(sideview_perchangle);
    side_Y_m=abs((perch_vs(4)-perch_vs(2))/(perch_vs(3)-perch_vs(1)));
else
    % slope is negative
    sideview_perchangle=-abs(sideview_perchangle);
    side_Y_m=-abs((perch_vs(4)-perch_vs(2))/(perch_vs(3)-perch_vs(1)));
end
if sign((perch_vs(8)-perch_vs(6))/(perch_vs(7)-perch_vs(5)))==-1
    % slope is positive
    underview_perchangle=abs(underview_perchangle);
    under_Y_m=abs((perch_vs(8)-perch_vs(6))/(perch_vs(7)-perch_vs(5)));
else
    % slope is negative
    underview_perchangle=-abs(underview_perchangle);
    under_Y_m=-abs((perch_vs(8)-perch_vs(6))/(perch_vs(7)-perch_vs(5)));
end

% Draw wheel cutout onto both views
disp('Click side-view then under-view wheel cutout.');
wheel_vs=getVertices(frame,'wheel cutout','wheel cutout');
% Defines perpendicular and X axis
% Get image angle for lines that are perpendicular in real space
underview_perpangle=getPerpAngle(wheel_vs(7:12));
under_Xangle=underview_perchangle+underview_perpangle;
under_X_m=tand(under_Xangle);
% flip slope signs, because Y axis is inverted by imagesc
under_X_m=-under_X_m;
under_Y_m=-under_Y_m;
side_Y_m=-side_Y_m;

% Get position of stopped pellet
disp('Click side-view then under-view stopped pellet.');
wheel_vs=getVertices(frame,'stopped pellet','stopped pellet');
% Defines X,Y origin
under_X_b=wheel_vs(4)-under_X_m*wheel_vs(3);
under_Y_b=wheel_vs(4)-under_Y_m*wheel_vs(3);
%[under_perchorigin_x,under_perchorigin_y]=intersectPoints(under_Y_m,under_X_m,under_Y_b,under_X_b);
under_origin_x=wheel_vs(3);
under_origin_y=wheel_vs(4);
side_origin_x=wheel_vs(1);
side_origin_y=wheel_vs(2);
showReferencesOnImage(frame,under_origin_x,under_origin_y,under_X_m,under_Y_m,under_X_b,under_Y_b);
reference.under_origin_x=under_origin_x; % in under image, x position of origin
reference.under_origin_y=under_origin_y; % in under image, y position of origin
reference.under_X_m=under_X_m; % in under image, slope of X axis
reference.under_Y_m=under_Y_m; % in under image, slope of Y axis
reference.under_X_b=under_X_b; % in under image, y-intercept/offset of X axis
reference.under_Y_b=under_Y_b; % in under image, y-intercept/offset of Y axis
reference.under_X_vec=[1 reference.under_X_m*1]; % direction of vector from origin
reference.under_Y_vec=[1 reference.under_Y_m*1]; % direction of vector from origin

% Get line in front of wheel cutout
disp('Click side-view then under-view cutout front edge for Z.');
cutoutZ_vs=getVertices(frame,'cutout front edge','cutout front edge');
% Defines Z axis and origin
side_Z_m=(cutoutZ_vs(4)-cutoutZ_vs(2))/(cutoutZ_vs(3)-cutoutZ_vs(1));
side_Z_b=side_origin_y-side_Z_m*side_origin_x;
side_Y_b=side_origin_y-side_Y_m*side_origin_x;
showReferencesOnImage([],side_origin_x,side_origin_y,side_Y_m,side_Z_m,side_Y_b,side_Z_b);
reference.side_origin_x=side_origin_x; % in side image, x position of origin
reference.side_origin_y=side_origin_y; % in side image, y position of origin
reference.side_Z_m=side_Z_m; % in side image, slope of Z axis
reference.side_Y_m=side_Y_m; % in side image, slope of Y axis
reference.side_Z_b=side_Z_b; % in side image, y-intercept/offset of Z axis
reference.side_Y_b=side_Y_b; % in side image, y-intercept/offset of Y axis
reference.side_X_vec=[1 reference.side_Y_m*1]; % direction of vector from origin, note that am using perch angle (actually Y axis) as approximation of X axis, because view is side-on
reference.side_Z_vec=[1 reference.side_Z_m*1]; % direction of vector from origin

% Get paw length when max outstretched
fig=implay(allframes,30);
fig.Parent.Position=[100 100 800 800];
pause;
if ~isempty(regexp(version,'2017b','once')) || ~isempty(regexp(version,'2018b','once')) ||  ~isempty(regexp(version,'2020b','once')) || ~isempty(regexp(version,'2021a','once'))
    currentFrameNumber=fig.DataSource.Controls.CurrentFrame;
else
    currentFrameNumber=fig.data.Controls.CurrentFrame;
end
% Instructions to user
continuebutton=questdlg('Choose frame where paw fully outstretched in front of mouse.','Instructions','Yes','Cancel','Cancel');
switch continuebutton
    case 'Yes'
    case 'Cancel'
        return
end
disp('Click side-view then under-view outstretched paw length.');
frame=allframes(:,:,:,currentFrameNumber);
pawlength_vs=getVertices(frame,'paw length','paw length');
% Unit conversion to real space
reference.side_pawlength=sqrt((pawlength_vs(3)-pawlength_vs(1))^2+(pawlength_vs(4)-pawlength_vs(2))^2); % in pixels
reference.under_pawlength=sqrt((pawlength_vs(7)-pawlength_vs(5))^2+(pawlength_vs(8)-pawlength_vs(6))^2); % in pixels
% From
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226774,
% fig. 1A
% Average length of mouse paw outstretched is 8.5 mm
reference.side_scaleTomm=8.5/reference.side_pawlength;
reference.under_scaleTomm=8.5/reference.under_pawlength;

% Have user find vertices for every ds_th frame
% More object and image points
disp('Click side-view then under-view then top-view. These points all need to match!');
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
% Some simplifying approximations:
% 1. Bottom view gives movement of paw in X-Y plane
% 2. Side view gives movement of paw in X-Z plane
% Here are the points:
% A. Pellet stopped position (0,0,0) -- origin
% B. Nose
% C. Paw positions during reach
% D. Top corners of pellet presenter wheel
% Origin
object_points=[0,0,0];
side_image_points=[reference.side_origin_x,reference.side_origin_y];
under_image_points=[reference.under_origin_x,reference.under_origin_y];
% Other points
for i=1:size(all_vs,1)
    temp=convertImagePointToObjectPoint(all_vs(i,1),all_vs(i,2),all_vs(i,3),all_vs(i,4),reference);
    object_points=[object_points; temp];
    side_image_points=[side_image_points; [all_vs(i,1), all_vs(i,2)]];
    under_image_points=[under_image_points; [all_vs(i,3), all_vs(i,4)]];
end
showPointsOnImage([],side_image_points,under_image_points,object_points);

% Change mirrors to "virtual cameras"
% Need to flip both mirror points horizontally
side_image_points=flipImagePoints(side_image_points,size(frame,1),size(frame,2),'horizontal');
under_image_points=flipImagePoints(under_image_points,size(frame,1),size(frame,2),'horizontal');
mirror_vs(1:2)=flipImagePoints(mirror_vs(1:2),size(frame,1),size(frame,2),'horizontal');
mirror_vs(3:4)=flipImagePoints(mirror_vs(3:4),size(frame,1),size(frame,2),'horizontal');

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