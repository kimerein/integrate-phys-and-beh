function plotEpilinesOnImage(varargin)

framesPerChunk=1500;
ds=10;
maxPoints=20;

% Get file name of video with reaches
filename=varargin{1};
discardFirstNFrames=varargin{2};
points=varargin{3};
epilines=varargin{4};

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

frame=allframes(:,:,:,currentFrameNumber);

% flip image
frame=flip(frame,2);

figure();
imagesc(frame);
colormap gray
hold on;
cmap=colormap('jet');
stepcmap=floor(size(cmap,1)./size(points,1));
indsIntoCmap=1:stepcmap:size(cmap,1);
indsIntoCmap=indsIntoCmap(1:size(points,1));

xrange=80:0.01:170;
for i=1:size(points,1)
    scatter(points(i,1),points(i,2),[],cmap(indsIntoCmap(i),:));
    a=epilines(i,1);
    b=epilines(i,2);
    c=epilines(i,3);
    y=(-a/b)*xrange-(c/b);
    plot(xrange,y,'Color',cmap(indsIntoCmap(i),:));
end


