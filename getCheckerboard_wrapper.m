function all_vs=getCheckerboard_wrapper(varargin)

framesPerChunk=1500;
ds=10;

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
sideview_perpangle=getPerpAngle(wheel_vs(1:6));
underview_perpangle=getPerpAngle(wheel_vs(7:12));
side_Xangle=sideview_perchangle-sideview_perpangle;
under_Xangle=underview_perchangle+underview_perpangle;
side_X_m=tand(side_Xangle);
under_X_m=tand(under_Xangle);
% flip slope signs, because Y axis is inverted by imagesc
side_X_m=-side_X_m;
under_X_m=-under_X_m;
side_Y_m=-side_Y_m;
under_Y_m=-under_Y_m;
side_Y_b=perch_vs(2)-side_Y_m*perch_vs(1);
under_Y_b=perch_vs(6)-under_Y_m*perch_vs(5);

% Get position of stopped pellet
disp('Click side-view then under-view stopped pellet.');
wheel_vs=getVertices(frame,'stopped pellet','stopped pellet');
% Defines X,Y origin
side_X_b=wheel_vs(2)-side_X_m*wheel_vs(1);
under_X_b=wheel_vs(4)-under_X_m*wheel_vs(3);
[side_origin_x,side_origin_y]=intersectPoints(side_Y_m,side_X_m,side_Y_b,side_X_b);
[under_origin_x,under_origin_y]=intersectPoints(under_Y_m,under_X_m,under_Y_b,under_X_b);

showReferencesOnImage(frame,side_origin_x,side_origin_y,side_X_m,side_Y_m,side_X_b,side_Y_b);
showReferencesOnImage([],under_origin_x,under_origin_y,under_X_m,under_Y_m,under_X_b,under_Y_b);

% Get line in front of wheel cutout
disp('Click side-view then under-view cutout front edge for Z.');
cutoutZ_vs=getVertices(frame,'cutout front edge','cutout front edge');
% Defines Z axis and origin

% Get paw length when max outstretched
disp('Click side-view then under-view outstretched paw length.');
pawlength_vs=getVertices(frame,'paw length','paw length');
% Unit conversion to real space

% Have user find vertices for every ds_th frame
% More object and image points
disp('Click side-view then under-view then top-view.');
all_vs=nan(length(1:ds:framesPerChunk),6);
j=1;
for i=1:size(allframes,4)
    if mod(i,ds)==0
        all_vs(j,:)=getVertices(allframes(:,:,:,i),'paw positions','paw pos');
        j=j+1;
    end
end

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