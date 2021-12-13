function varargout = makeCheckerboardOnVid(varargin)
% makeCheckerboardOnVid MATLAB code for makeCheckerboardOnVid.fig
%      makeCheckerboardOnVid, by itself, creates a new makeCheckerboardOnVid or raises the existing
%      singleton*.
%
%      H = makeCheckerboardOnVid returns the handle to a new makeCheckerboardOnVid or the handle to
%      the existing singleton*.
%
%      makeCheckerboardOnVid('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PERCHZONEGUI.M with the given input arguments.
%
%      makeCheckerboardOnVid('Property','Value',...) creates a new makeCheckerboardOnVid or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before perchZoneGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to perchZoneGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help perchZoneGUI

% Last Modified by GUIDE v2.5 09-Dec-2021 17:27:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @makeCheckerboardOnVid_OpeningFcn, ...
                   'gui_OutputFcn',  @makeCheckerboardOnVid_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before makeCheckerboardOnVid is made visible.
function makeCheckerboardOnVid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to makeCheckerboardOnVid (see VARARGIN)

% Choose default command line output for makeCheckerboardOnVid
handles.output = hObject;
handles.firstVertex=[];
handles.secondVertex=[];
handles.thirdVertex=[];
handles.fourthVertex=[];
handles.fifthVertex=[];
handles.sixthVertex=[];

slice=varargin{1};
guititle=varargin{2};
whichInput=varargin{3}; % 'perch line','wheel cutout','stopped pellet','cutout front edge','paw pos','paw length'
h=imagesc(slice);
colormap gray
set(h.Parent.Parent,'Position',[0 0 560 350]);
set(h.Parent.Parent,'Name',guititle);
handles.h=h;
handles.slice=slice;
handles.whichInput=whichInput;
handles.vertices=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes makeCheckerboardOnVid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = makeCheckerboardOnVid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.whichInput
    case 'paw length'
        targetSize=5;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
            line([handles.firstVertex(1) handles.secondVertex(1)],[handles.firstVertex(2) handles.secondVertex(2)],'Color','y');
        elseif ~isfield(handles,'isThirdPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.thirdVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','g');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','g');
        elseif ~isfield(handles,'isFourthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.fourthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','r');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','r');
            line([handles.thirdVertex(1) handles.fourthVertex(1)],[handles.thirdVertex(2) handles.fourthVertex(2)],'Color','y');
        end
    case 'paw pos'
        targetSize=20;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
        elseif ~isfield(handles,'isThirdPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.thirdVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','y');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','y');
        end
    case 'perch line'
        targetSize=10;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
            line([handles.firstVertex(1) handles.secondVertex(1)],[handles.firstVertex(2) handles.secondVertex(2)],'Color','y');
        elseif ~isfield(handles,'isThirdPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.thirdVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','g');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','g');
        elseif ~isfield(handles,'isFourthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.fourthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','r');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','r');
            line([handles.thirdVertex(1) handles.fourthVertex(1)],[handles.thirdVertex(2) handles.fourthVertex(2)],'Color','y');
        end
    case 'cutout front edge'
        targetSize=5;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
            line([handles.firstVertex(1) handles.secondVertex(1)],[handles.firstVertex(2) handles.secondVertex(2)],'Color','y');
        elseif ~isfield(handles,'isThirdPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.thirdVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','g');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','g');
        elseif ~isfield(handles,'isFourthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.fourthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','r');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','r');
            line([handles.thirdVertex(1) handles.fourthVertex(1)],[handles.thirdVertex(2) handles.fourthVertex(2)],'Color','y');
        end
    case 'wheel cutout'
        targetSize=5;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
            line([handles.firstVertex(1) handles.secondVertex(1)],[handles.firstVertex(2) handles.secondVertex(2)],'Color','y');
        elseif ~isfield(handles,'isThirdPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.thirdVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','g');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','g');
            line([handles.secondVertex(1) handles.thirdVertex(1)],[handles.secondVertex(2) handles.thirdVertex(2)],'Color','y');
        elseif ~isfield(handles,'isFourthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.fourthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isFifthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=false;
            handles.isFifthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.fifthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
            line([handles.fourthVertex(1) handles.fifthVertex(1)],[handles.fourthVertex(2) handles.fifthVertex(2)],'Color','y');
        elseif ~isfield(handles,'isSixthPress')
            handles.isFirstPress=false;
            handles.isSecondPress=false;
            handles.isThirdPress=false;
            handles.isFourthPress=false;
            handles.isFifthPress=false;
            handles.isSixthPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.sixthVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','g');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','g');
            line([handles.fifthVertex(1) handles.sixthVertex(1)],[handles.fifthVertex(2) handles.sixthVertex(2)],'Color','y');
        end
    case 'stopped pellet'
        targetSize=20;
        if ~isfield(handles,'isFirstPress')
            handles.isFirstPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.firstVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','m');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','m');
        elseif ~isfield(handles,'isSecondPress')
            handles.isFirstPress=false;
            handles.isSecondPress=true;
            [currVertex_x,currVertex_y]=ginput(1);
            handles.secondVertex=[currVertex_x,currVertex_y];
            line([nanmax([currVertex_x-targetSize 0]) currVertex_x+targetSize],[currVertex_y currVertex_y],'Color','c');
            line([currVertex_x currVertex_x],[nanmax([currVertex_y-targetSize 0]) currVertex_y+targetSize],'Color','c');
        end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global continueAnalysis vertexViews

% Hint: delete(hObject) closes the figure
% disp(handles.firstVertex);
% disp(handles.secondVertex);
% disp(handles.thirdVertex);
continueAnalysis=1;
vertexViews=[];
switch handles.whichInput
    case 'paw length'
        vertexViews=[handles.firstVertex; handles.secondVertex; handles.thirdVertex; handles.fourthVertex];
    case 'paw pos'
        vertexViews=[handles.firstVertex; handles.secondVertex; handles.thirdVertex];
    case 'perch line'
        vertexViews=[handles.firstVertex; handles.secondVertex; handles.thirdVertex; handles.fourthVertex];
    case 'cutout front edge'
        vertexViews=[handles.firstVertex; handles.secondVertex; handles.thirdVertex; handles.fourthVertex];
    case 'wheel cutout'
        vertexViews=[handles.firstVertex; handles.secondVertex; handles.thirdVertex; handles.fourthVertex; handles.fifthVertex; handles.sixthVertex];
    case 'stopped pellet'
        vertexViews=[handles.firstVertex; handles.secondVertex];
end
delete(hObject);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
