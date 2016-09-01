function varargout = scintillate(varargin)
% SCINTILLATE MATLAB code for scintillate.fig
%      SCINTILLATE, by itself, creates a new SCINTILLATE or raises the existing
%      singleton*.
%
%      H = SCINTILLATE returns the handle to a new SCINTILLATE or the handle to
%      the existing singleton*.
%
%      SCINTILLATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCINTILLATE.M with the given input arguments.
%
%      SCINTILLATE('Property','Value',...) creates a new SCINTILLATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scintillate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scintillate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% See readme.MD for dependency information.

% Edit the above text to modify the response to help scintillate

% Last Modified by GUIDE v2.5 07-Jun-2016 12:03:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @scintillate_OpeningFcn, ...
    'gui_OutputFcn',  @scintillate_OutputFcn, ...
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

% ------------------------------------------------------------------------
% --- Executes just before scintillate is made visible.
function scintillate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scintillate (see VARARGIN)

% Choose default command line output for scintillate

handles.output = hObject;

% Add scintillate working directory and subdirectories to path
p = genpath(pwd);
addpath (p);
% Update handles structure
guidata(hObject, handles);

%scintillateData.filename = [];
%guidata(hObject, scintillateData);
% UIWAIT makes scintillate wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% ------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = scintillate_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------------------------------------------------------------------
% --- Executes on button press in openbutton.
function [Filename]=openbutton_Callback(hObject, eventdata, handles)
tsStack = [];
loadedstring = char('Scintillate; Dublon et al. 2016 Kemisk Ekologi, Sveriges lantbruksuniversitet, Alnarp.');
set(handles.statustext,'String', loadedstring);

% Read in tiff data
[Filename,PathName] = uigetfile( ...
    {'*.tif;*.tiff;', 'Tiff File (*.tif,*.tiff)';
    ...
    '*.*', 'All Files (*.*)'}, ...
    'Open tiff file');

if isequal(Filename,0)
    helpdlg('Load cancelled','Load');
    loadedstring = char('Load cancelled, previous data if present remains resident in memory.');
    set(handles.statustext,'String', loadedstring);
    return
else
    loadedstring = char(sprintf('Loaded: %s' , Filename));
    short_loadedstring = char(sprintf('%s' , Filename));
    set(handles.statustext,'String', loadedstring);
    set(handles.filename_string,'String',short_loadedstring);
    [pathstr, name, ext] = fileparts(Filename);
end
%tic
tsStack = TIFFStack([PathName,Filename]);
no_of_frames = size(tsStack,3);
%set handles info to enable accessibility by other functions
handles.dat = tsStack;
handles.filename = Filename;
handles.pathstr = PathName;
handles.name = name;
handles.ext = ext;
handles.no_of_frames = no_of_frames;
loadedstring = char('Plotting all frames...');
set(handles.statustext,'String', loadedstring);
%use guidata to enable accessibility by other functions
guidata(hObject,handles);

% Enable buttons now we have data
set(handles.view_movie,'Enable','on');
set(handles.delta_up,'Enable','on');
set(handles.background_subtraction,'Enable','on');
set(handles.threshold_button,'Enable','on');
set(handles.pre_stim_mean_button,'Enable','on');
set(handles.ica_button,'Enable','on');
set(handles.impixelregion,'Enable','on');
set(handles.col_filter_button,'Enable','on');
set(handles.crop_button,'Enable','on');

% Calibrate the slider making the maximum number of frames equal to
% no_of_frames and not the default 250 (changeable in GUIDE) then enable the slider.
set(handles.slider_frame,'Max',no_of_frames);
set(handles.slider_frame,'Enable','on');

% Define stimulus delivery frame
dlg_title = 'Stimulus';
num_lines = 1;
prompt{1} = 'Stimulus delivered at frame: ';
stim_defaults = {'11'};  % This default can be varied here
options.Resize='off';
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,stim_defaults);

stim = str2double(answer{1});

% Produce reference value for basal fluorescence: sumImage
limit =stim -1;

for k = 1 : limit
    thisImage = tsStack(:,:,k);
    if k == 1
        sumImage = thisImage;
    else
        sumImage = sumImage + thisImage;
    end
end
sumImage = sumImage / limit;

% Thresholding options
pixel_filter = str2double(get(handles.edit_text_pixel_filter,'String'));
threshold = str2double(get(handles.edit_text_threshold,'String'));
if threshold == 0
    level = graythresh(sumImage); % global image threshold using Otsu's method
    bw = im2bw(sumImage,level);
    threshold_msg = sprintf('Using global image threshold using Otsu method. Returned as %d',level);
    replacement_threshold_string = char(sprintf('%s', level));
    set(handles.edit_text_threshold,'String', replacement_threshold_string);
    set(handles.statustext,'String', threshold_msg);
else
    bw = im2bw(sumImage,threshold);
end
bw_cut = bwareaopen(bw, pixel_filter); %remove small objects fewer than pixel_filter pixels from image to cut down computational time
colorbar;
colval = get(handles.col_listbox,'Value');
switch colval;
    case 1
        colormap bone
    case 2
        colormap jet
    case 3
        colormap cool
    case 4
        colormap hsv
    case 5
        colormap parula
end

% Plot all the frames and arrange according to size of stack. Colour
% according to colourmap

for sp=1:no_of_frames;
    val=ceil(sqrt(no_of_frames));
    ssp=subplot(val,val,sp);
    imagesc(tsStack(:, :, sp));
    grid('off');
    axis('off');
    colval = get(handles.col_listbox,'Value');
    switch colval;
        case 1
            colormap bone
        case 2
            colormap jet
        case 3
            colormap cool
        case 4
            colormap hsv
        case 5
            colormap parula
    end
    msg = sprintf('Frames up to %d',no_of_frames);
    set(handles.statustext,'String', msg)
    
    % define figure handle so that the figure can later be destroyed
    figHandles = findall(0,'Type','axes');
    handles.subplotaxes=figHandles;
    guidata(gcbo, handles);
end
% Plot the pre-stimulus reference frame (sumImage) in its own window
figure('Name', 'Pre-stimulus reference average','NumberTitle','off');
imagesc(sumImage);
switch colval;
    case 1
        colormap bone
    case 2
        colormap jet
    case 3
        colormap cool
    case 4
        colormap hsv
    case 5
        colormap parula
end
axis('off');
colorbar;
%toc

% ------------------------------------------------------------------------
% --- Executes on button press in view_movie.
function view_movie_Callback(~, ~, handles)
tsStack = handles.dat;
no_of_frames = handles.no_of_frames;

% Delete existing figures in the main window but retain UI elements
subplots = handles.subplotaxes;
try
    delete(subplots);
end
viewer =0;

% Go through the stack displaying each frame and pausing by the user defined
% pause time (ptime) variable
for viewer = 1:no_of_frames;
    movfigure = imagesc(tsStack(:, :, viewer));
    colorbar;
    colval = get(handles.col_listbox,'Value');
    switch colval;
        case 1
            colormap bone
        case 2
            colormap jet
        case 3
            colormap cool
        case 4
            colormap hsv
        case 5
            colormap parula
            
    end
    
    drawnow;
    m = [num2str(viewer)];
    set(handles.statustext,'String', sprintf('Frame %s' ,m));
    set(handles.slider_frame, 'value',viewer);
    ptime = str2double(get(handles.pausetime,'String'));
    pause(ptime);
end


% ------------------------------------------------------------------------
% --- Executes on button press in background_subtraction.
function background_subtraction_Callback(~, ~, handles)

tsStack = handles.dat;
% Delete existing figures in the main window but retain UI elements
subplots=handles.subplotaxes;
no_of_frames = handles.no_of_frames;

tsStack = handles.dat;
subplots=handles.subplotaxes;
deltaF = [char(916) 'F '];
try
    delete(subplots);
end
dlg_title = 'Average Pre Stimulus background';
num_lines = 1;
prompt{1} = 'Stimulus delivered at frame: ';
stim_defaults = {'11'};
options.Resize='off';
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,stim_defaults);
stim = str2double(answer{1});
limit =stim -1;
stim_end = size(tsStack,3);

for k = 1 : limit
    thisImage = tsStack(:,:,k);
    if k == 1
        sumImage = thisImage;
    else
        sumImage = sumImage + thisImage;
    end
end
sumImage = sumImage / limit;
post_frame =[];
tic
for post_frame = stim:stim_end
    sub_image(:,:,post_frame) = tsStack(:,:,post_frame) - sumImage;
    sub_image2(:,:,post_frame) = sub_image(:,:,post_frame)./sumImage; %this allows for division by the pre stimulus average
end
toc

for viewer = stim:stim_end
    imagesc(sub_image(:,:,viewer)) % Use sub_image2 instead of sub_image for division by the pre stimulus average
    colorbar;
    colval = get(handles.col_listbox,'Value');
    m = [deltaF num2str(viewer)];
    set(handles.statustext,'String', m);
    ptime = str2double(get(handles.pausetime,'String'));
    set(handles.slider_frame, 'value',viewer);
    pause(ptime);
end

% ------------------------------------------------------------------------
% --- Executes on button press in delta_up.
function delta_up_Callback(hObject, ~, handles)

tsStack = handles.dat;
% Delete existing figures in the main window but retain UI elements
subplots=handles.subplotaxes;
delta = char(916);
no_of_frames = handles.no_of_frames;
try
    delete(subplots);
end

button_click = get(hObject,'Value');
if button_click ==1;
    set(hObject, 'String' , 'frame delta up');
    
    % Use imabsdiff to see absolute matrix by matrix differences
    
    for absdiff = 1:no_of_frames-1
        ni = imabsdiff(tsStack(:,:,absdiff),tsStack(:,:,absdiff+1));
        imagesc(ni);
        drawnow;
        set(handles.statustext,'String', sprintf('Frame %s up %s' ,delta, num2str(absdiff)));
        set(handles.slider_frame, 'Value',(absdiff));
        ptime = str2double(get(handles.pausetime,'String'));
        pause(ptime);
    end
    
else
    
    set(hObject, 'String' , 'frame delta down')
    for k = no_of_frames-1:-1:1
        d(:, :, k) = imabsdiff(tsStack(:, :, k), tsStack(:, :, k+1));
        imagesc(d(:, :, k));
        drawnow
        m = [num2str(k)];
        set(handles.statustext,'String', sprintf('Frame %s down %s' ,delta, num2str(m)));
        set(handles.slider_frame, 'value',k);
        ptime = str2double(get(handles.pausetime,'String'));
        pause(ptime);
    end
end
% ------------------------------------------------------------------------
function pausetime_Callback(hObject, ~, handles)
% Hints: get(hObject,'String') returns contents of pausetime as text
% UI element for pause time
ptime = str2double(get(hObject,'String'));
set(handles.statustext,'String', num2str(ptime));
% ------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function pausetime_CreateFcn(hObject, ~, ~)
% hObject    handle to pausetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% ------------------------------------------------------------------------
% --- Executes on button press in col_filter_button.
function col_filter_button_Callback(~, ~, handles)
tsStack = handles.dat;
% Delete existing figures in the main window but retain UI elements
subplots=handles.subplotaxes;
no_of_frames = handles.no_of_frames;
try
    delete(subplots);
end
set(handles.statustext,'String', 'thinking...')
dlg_title = 'Column filter dimensions';
num_lines = 1;
prompt{1} = 'Filter size: ';
smooth_defaults = {'3'};
options.Resize='off';
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,smooth_defaults);
smooth = str2double(answer{1});
% column filter
for colfilter = 1:no_of_frames
    filtered_tsStack = uint16(colfilt(tsStack(:,:,colfilter),[smooth,smooth],'sliding',@mean));
    imagesc(filtered_tsStack)
    drawnow;
    set(handles.statustext,'String', sprintf('Frame %s' ,num2str(colfilter)));
    set(handles.slider_frame, 'value', colfilter);
end
% ------------------------------------------------------------------------
% --- Executes on button press in crop_button.
function crop_button_Callback(hObject, ~, handles)
% hObject    handle to crop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tsStack = handles.dat;
no_of_frames = handles.no_of_frames;
% Delete existing figures in the main window but retain UI elements
subplots = handles.subplotaxes;
try
    delete(subplots);
end
viewer =0;

tss= imagesc(tsStack(:, :, 1));
drawnow expose;
set(handles.statustext,'String', 'Select area, secondary click and choose Crop Image');
[sel rect] = imcrop(tss);
[file,path] = uiputfile('*.tif','Save cropped output',sprintf('cropped_%s', handles.name));
fp = [path file];
% Cancel gracefully if the user cancels
if isequal(file,0)
    helpdlg('Save as .tif cancelled','Save')
end

divider = 10^(floor(log10(no_of_frames))-1);
tic
% Crop using save as TIFF and update user as to progress
for i=1:no_of_frames
    cropped_tsStack(:, :, i)= imcrop(tsStack(:, :, i),[rect]);
    saveastiff (uint16(cropped_tsStack), fp)
    
    if (round(i/divider)==i/divider)
        progress = fprintf('Frame %d written in %.0f seconds, %2d percent complete, time left=%.0f seconds \n', ...
            i, toc, i/no_of_frames*100, (no_of_frames - i)/(i/toc));
        
    end
end
% Graceful exit if required
if isequal(file,0);
    set(handles.statustext,'String', 'Crop save is cancelled');
else
    set(handles.statustext,'String', sprintf('Saved as %s ' ,[path file]));
end
% ------------------------------------------------------------------------
% --- Executes on button press in ica_button.
function ica_button_Callback(hObject, eventdata, handles)
% hObject    handle to ica_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global tsStack

tsStack = handles.dat;
% Generate prompt for ICA variables
prompt = {'subtract reference (1 is yes) : ','subtract global reference (1 is yes): ', ...
    'use one colour only (1 is yes): ', 'emphasize colours (1 is yes): ','ind reference: ', ...
    'number of ICA components : ' ,'background weighting : ' , ...
    'font size for graphs : ', 'downsample to resolution (px) : ', ...
    'false colour coding scheme : ','frames to skip at the start : ','frames to skip at the end : '};
dlg_title = 'ICA Parameters';
options.Resize='off';
num_lines = 1;
% Here default values may be replaced.
def = {'1','1','0','0','1:11','3','0.1','12','150','@hsv','1','1'};
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,def);
% Graceful exit if user cancels.
if isempty(answer) ;
    return;
else
    conf.do_subtract_ref = str2num(answer{1});
    conf.do_subtract_global = str2num(answer{2});
    conf.do_color_one_only = str2num(answer{3});
    conf.do_emphasize_col = str2num(answer{4});
    conf.ind_ref = str2num(answer{5});
    conf.n_ica_comp = str2num(answer{6});
    conf.w_bg = str2num(answer{7});
    conf.my_fontsize = str2num(answer{8});
    conf.lowres = str2num(answer{9});
    conf.f_col = str2func(answer{10});
    skip_start = str2num(answer{11});
    conf.ind_skip_start = [1:skip_start]; %We always skip a minimum of 1 frame at the start
    skip = str2num(answer{12})-1; % We always skip a minimum of 1 frame at the end
    skip_frame = size(tsStack,3)-skip;
    conf.ind_skip_end = [skip_frame:size(tsStack,3)];
    % Option to crop here
    i_range = conf.ind_skip_start:size(tsStack,1); % cropping in height...
    j_range = conf.ind_skip_end:size(tsStack,2); % and here cropping in length
end
plot_ica_analysis(tsStack, conf, i_range, j_range);
% ------------------------------------------------------------------------
% --- Executes on selection change in col_listbox.
function col_listbox_Callback(hObject, ~, ~)
% hObject    handle to col_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns col_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from col_listbox
contents = cellstr(get(hObject,'String'));
colbar_type = contents{get(hObject,'Value')};
% ------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function col_listbox_CreateFcn(hObject, ~, ~)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% ------------------------------------------------------------------------
% --- Executes on button press in pre_stim_mean_button.
function pre_stim_mean_button_Callback(~, ~, handles)
% hObject    handle to pre_stim_mean_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tsStack = handles.dat;
subplots = handles.subplotaxes;
no_of_frames = handles.no_of_frames;
try
    delete(subplots);
end
dlg_title = 'Average Pre Stimulus background';
num_lines = 1;
prompt{1} = 'Stimulus delivered at frame: ';
stim_defaults = {'11'};
options.Resize='off';
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,stim_defaults);
stim = str2double(answer{1});
limit =stim -1;

for k = 1 : limit
    thisImage = tsStack(:,:,k);
    if k == 1
        sumImage = thisImage;
    else
        sumImage = sumImage + thisImage;
    end
end
sumImage = sumImage / limit;
% Here pixel filter value is used
pixel_filter = str2double(get(handles.edit_text_pixel_filter,'String'));
threshold = str2double(get(handles.edit_text_threshold,'String'));
if threshold == 0
    level = graythresh(sumImage); % global image threshold using Otsu's method
    bw = im2bw(sumImage,level);
    threshold_msg = sprintf('Using global image threshold using Otsu method. Returned as %d',level);
    replacement_threshold_string = char(sprintf('%s', level));
    set(handles.edit_text_threshold,'String', replacement_threshold_string);
    set(handles.statustext,'String', threshold_msg);
else
    bw = im2bw(sumImage,threshold);
end
bw_cut = bwareaopen(bw, pixel_filter); %remove small objects fewer than pixel_filter pixels from image to cut down computational time
cc = bwconncomp(bw_cut, 6); %create connected regions;

numPixels = cellfun(@numel,cc.PixelIdxList);
%[biggest,idx] = max(numPixels);
%bw(cc.PixelIdxList{idx}) = 0;
imagesc(sumImage);
colorbar;
colval = get(handles.col_listbox,'Value');
msg = sprintf('Displaying average background summation for frames 1 to %d',limit);
set(handles.statustext,'String', msg)

titlestring=['Frames 1 - ' num2str(no_of_frames)];

newFrameByFramefig=figure('Name', titlestring,'NumberTitle','off');
for sp=1:no_of_frames;
    val=ceil(sqrt(no_of_frames));
    
    %ssp=subplot(val,val,sp);
    subplot('Position',[(mod(sp-1,8))/8 1-(ceil(sp/8))/8 1/8 1/8])
    imagesc(tsStack(:, :, sp));
    grid('off');
    axis('off');
    
    msg = sprintf('Frames up to %d',no_of_frames);
    %set(handles.statustext,'String', msg)
end

% ------------------------------------------------------------------------
function edit_text_threshold_Callback(~, ~, ~)
% hObject    handle to edit_text_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_text_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_text_threshold as a double

% ------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function edit_text_threshold_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_text_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
filter_threshold = str2double(get(hObject,'String'));
update = num2str(filter_threshold);
try
    set(handles.statustext,'String', update);
end

% ------------------------------------------------------------------------
function edit_text_pixel_filter_Callback(~, ~, ~)
% hObject    handle to edit_text_pixel_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_text_pixel_filter as text
%        str2double(get(hObject,'String')) returns contents of edit_text_pixel_filter as a double

% ------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function edit_text_pixel_filter_CreateFcn(hObject, ~, handles)
% hObject    handle to edit_text_pixel_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%global tsStack
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

pixel_filter = str2double(get(hObject,'String'));
update = num2str(pixel_filter);
try
    set(handles.statustext,'String', update);
end
% ------------------------------------------------------------------------
% --- Executes when figure1 is resized.
function figure1_ResizeFcn(~, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ------------------------------------------------------------------------
% --- Executes on button press in threshold_button.
function threshold_button_Callback(~, ~, handles)
% hObject    handle to threshold_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global tsStack
tsStack = handles.dat;
subplots=handles.subplotaxes;
no_of_frames = handles.no_of_frames;

try
    delete(subplots);
end

dlg_title = 'Average Pre Stimulus background';
num_lines = 1;
prompt{1} = 'Stimulus delivered at frame: ';
stim_defaults = {'11'};
options.Resize='off';
answer = inputdlg_no_cancel(prompt,dlg_title,num_lines,stim_defaults);
stim = str2double(answer{1});
limit =stim -1;
pixel_filter = str2double(get(handles.edit_text_pixel_filter,'String'));
threshold = str2double(get(handles.edit_text_threshold,'String'));
if threshold == 0
    level = graythresh(tsStack(:,:,1)); % global image threshold using Otsu's method
    threshold_msg = sprintf('Using global image threshold using Otsu method. Returned as %d',level);
    replacement_threshold_string = char(sprintf('%s', level));
    set(handles.edit_text_threshold,'String', replacement_threshold_string);
    set(handles.statustext,'String', threshold_msg);
    no_of_frames = handles.no_of_frames;
else
    level = threshold;
end
bw = [];
bw_cut = [];
for viewer = 1:no_of_frames;
    bw(:, :, viewer) = im2bw(tsStack(:,:,viewer),level);
    bw_cut(:, :, viewer)  = bwareaopen(bw(:, :, viewer), pixel_filter);
    cc(viewer) = bwconncomp(bw_cut(:, :, viewer), 6);
    [numPixels] = cellfun(@numel,cc(viewer).PixelIdxList);
    drawme = imagesc(bw_cut(:, :, viewer));
    drawnow;
    m = [num2str(viewer)];
    set(handles.statustext,'String', sprintf('Frame %s' ,m));
    ptime = str2double(get(handles.pausetime,'String'));
    set(handles.slider_frame, 'value',viewer);
    pause(ptime);
end

% ------------------------------------------------------------------------
% --- Executes on button press in impixelregion.
function impixelregion_Callback(hObject,~,handles)
% hObject    handle to impixelregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tsStack = handles.dat;
impixelregion;

% ------------------------------------------------------------------------
% --- Executes on slider movement.
function slider_frame_Callback(hObject, ~, handles)
% hObject    handle to slider_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
tsStack = handles.dat;
no_of_frames = handles.no_of_frames;
min = 1;
set(hObject,'SliderStep',[1/no_of_frames, 10/no_of_frames],'Max',no_of_frames);
val = round(get(hObject,'Value'));

if val == 0;
    val =1;
end

colval = get(handles.col_listbox,'Value');
switch colval;
    case 1
        colormap bone
    case 2
        colormap jet
    case 3
        colormap cool
    case 4
        colormap hsv
    case 5
        colormap parula
end
imagesc(tsStack(:,:,val));
m = [num2str(val)];
set(handles.statustext,'String', m);
drawnow;
% ------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function slider_frame_CreateFcn(hObject, ~, handles)
% hObject    handle to slider_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% ------------------------------------------------------------------------
function main_menu_Callback(~, ~, ~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
function about_menuitem_Callback(~, ~, ~)
msgbox({'Scintillate','v.1.0.0 (b20160607)',...
    'A MATLAB based GUI graphic visualiser for multi-stacked TIFF files.','',...
    'Ian Dublon, Markus Nilsson, Anna Balkenius, Peter Anderson and Mattias Larsson','2016',...,'',...,'',...
    'This program uses some 3rd party freeware components; please see License menu item for more information.','',...
    'ian.dublon@slu.se | ian.dublon@gmail.com'},'About');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
function quit_menuitem_Callback(~, ~, ~)
close all
clc
quit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
function license_menuitem_Callback(~, ~, ~)
ai = char(228);
FastICAauthors = ['Hugo G' ai 'vert, Jarmo Hurri, Jaako S' ai 'rel' ai ', and Aapo Hyv' ai 'rinen'];
msgbox({'This product is provided "as is", with absolutely no warranty.','',...
    'Please cite this software if you use it.',...
    'Licensed under the terms of the Creative Commons Attribution 3.0 license.',...
    'http://creativecommons.org/licenses/by/3.0/','',...
    'We gratefully acknowledge the following authors for their freeware functions, used herein:','',...
    'TIFFStack.m (C) Dylan Muir (http://www.mathworks.com/matlabcentral/fileexchange/32025-tiffstack/content/@TIFFStack/TIFFStack.m)',''...
    'FastICA package for MATLAB (C) 1996-2005',...
    FastICAauthors 'Laboratory of Information and Computer Science, Helsinki University of Technology, Finland (http://research.ics.aalto.fi/ica/software.shtml)',''...
    'saveastiff.m (C) YoonOh Tak (http://www.mathworks.com/matlabcentral/fileexchange/35684-save-and-load-data-as-multi-frame-tiff-format/content//saveastiff.m)',''...
    'fastsmooth.m (C) Tom O''Haver (http://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function)'},'Licensing');
% --------------------------------------------------------------------
function getting_started_menuitem_Callback(~, ~, ~)
% hObject    handle to getting_started_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whatplatform = ispc;
if whatplatform>0
    winopen('getting_started.pdf');
else
    open('getting_started.pdf');
end
