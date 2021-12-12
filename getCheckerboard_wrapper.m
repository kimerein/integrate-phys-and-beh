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

% Draw perch line onto both views
perch_vs=getPerchLine(frame);

return

% Have user find vertices for every ds_th frame
all_vs=nan(length(1:ds:framesPerChunk),6);
j=1;
for i=1:size(allframes,4)
    if mod(i,ds)==0
        all_vs(j,:)=getVertices(allframes(:,:,:,i));
        j=j+1;
    end
end

end

function vs=getPerchLine(frame)

global vertexViews
global continueAnalysis

makeCheckerboardOnVid(frame,'Perch line','perch line');
disp('Click side-view perch line then under-view perch line.');
pause;
temp=vertexViews;
temp=temp';
vs=temp(1:end); % x_v1, y_v1, x_v2, y_v2, x_v3, y_v3, x_v4, y_v4

if continueAnalysis==1
    disp('Succesfully defined');
else
    disp('Failed to define');
end

end

function vs=getVertices(frame)

global vertexViews
global continueAnalysis

makeCheckerboardOnVid(frame,'Paw positions','paw pos');
disp('Click side-view then under-view then top-view.');
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