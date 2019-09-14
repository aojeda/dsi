function SourceConnectivity
fig = findall(0,'Name','Source Connectivity');
if ~isempty(fig)
    figure(fig);
    return;
end
fig = figure('Name','Source Connectivity','NumberTitle',false,'Menu','none','Toolbar','none','Units','normalized');
fig.Position(3:4) = [0.3766    0.5157];
fig.UserData = struct('M',[],'Y',[],'C',[],'L',[]);

p = fileparts(which('headModel.m'));
p = fullfile(p,'+vis','icons');
playIcon = imread([p filesep 'Gnome-media-playback-start.svg.png']);
helpIcon = imread([p filesep 'Gnome-help-browser.svg.png']);
plotIcon = imread([p filesep 'Gnome-emblem-photos.svg.png']);
                              
pnt = uipanel(fig,'Units','normalized','Position',[0 0.9015 1 0.1000]);
uicontrol(pnt,'style','text','String','Data directory','Units','Normalized','Position',[0.0057    0.5    0.2075    0.375]);
uicontrol(pnt,'style','edit','Units','Normalized',                          'Position',[0.2095    0.5    0.6188    0.375],'tag','Data directory');
uicontrol(pnt,'style','pushbutton','String','Select','Units','Normalized',  'Position',[0.8302    0.5    0.1472    0.375],'callback',@callback_select1);

uicontrol(pnt,'style','text','String','Results directory','Units','Normalized','Position',[0.0057    0.0542    0.2075    0.375]);
uicontrol(pnt,'style','edit','Units','Normalized',                          'Position',[0.2095    0.1015    0.6188    0.375],'tag','Results directory');
uicontrol(pnt,'style','pushbutton','String','Select','Units','Normalized',  'Position',[0.8302    0.0859    0.1472    0.375],'callback',@callback_select2);

ax = axes(fig, 'Units','normalized');
ax.Position = [0.1300    0.2    0.75 0.65];
pn = uipanel(fig,'Units','normalized');
pn.Position = [0      0    1 0.0529];
lb = uicontrol(pn,'style','text','string','Metric','Units','normalized');
lb.Position = [0.0059         0    0.06    0.8000];
tx = uicontrol(pn,'style','popupmenu','string',{'Mutual Information','Transfer Entropy','Correlation'},'Units','normalized');
tx.Position = [0.07    0.0400    0.2    0.8800];

lb = uicontrol(pn,'style','text','string','WinSize (sec)','Units','normalized');
lb.Position = [0.2860         0    0.1322    0.8000];
tx = uicontrol(pn,'style','edit','string','0.5','Units','normalized','tag','winSize');
tx.Position = [0.4178    0.0800    0.0789    0.8000];

lb = uicontrol(pn,'style','text','string','Overlap (%)','Units','normalized');
lb.Position = [0.5078         0    0.1322    0.8000];
tx = uicontrol(pn,'style','edit','string','50','Units','normalized','tag','overlap');
tx.Position = [0.6335    0.0800    0.0626    0.8000];

tx = uicontrol(pn,'style','radiobutton','string','Normalize','Units','normalized','tag','normalize');
tx.Position = [0.7192    0.0800    0.1215    0.8000];

bt1 = uicontrol(pn,'style','pushbutton','CData',playIcon,'TooltipString','Run','Units','normalized','callback',@callback_run, 'Units','normalized');
bt1.Position = [0.8624         0    0.0436    1.0000];
bt2 = uicontrol(pn,'style','pushbutton','CData',plotIcon,'TooltipString','Vualize','callback',@callback_plot,'Units','normalized');
bt2.Position = [0.9084         0    0.0436    1.0000];
bt3 = uicontrol(pn,'style','pushbutton','CData',helpIcon,'TooltipString','Help','callback',@callback_help,'Units','normalized');
bt3.Position = [0.9530         0    0.0436    1.0000];
end

%%
function callback_select1(src,evnt)
fig = src.Parent;
selpath = uigetdir(pwd, 'Select data folder');
if selpath
    set(findall(fig,'tag','Data directory'),'String',selpath);
end
end

function callback_select2(src,evnt)
fig = src.Parent;
selpath = uigetdir(pwd, 'Select data folder');
if selpath
    set(findall(fig,'tag','Results directory'),'String',selpath);
end
end

function callback_run(src,evnt)
fig = src.Parent.Parent;
src = get(findall(fig,'tag','Data directory'), 'String');
if ~exist(src,'dir')
    return;
end
dest = get(findall(fig,'tag','Results directory'), 'String');
h = findall(fig,'style','popupmenu');
metric = lower(h.String{h.Value});
winSize = str2double(get(findall(fig,'tag','winSize'), 'String'));
overlap = str2double(get(findall(fig,'tag','overlap'), 'String'))/100;
normalize = get(findall(fig,'tag','normalize'), 'Value');
dest = fullfile(dest,metric);
mkdir(dest)

if isempty(winSize), error('Invalid value of parameter WinSize');end
if isempty(overlap), error('Invalid value of parameter Overlap');end

files = pickfiles(src, '.set');
n = size(files,1);
jFrame = get(handle(fig),'JavaFrame');
try
    jRootPane = jFrame.fHG1Client.getWindow;
catch
    jRootPane = jFrame.fHG2Client.getWindow;
end     
jStatusBar = com.mathworks.mwswing.MJStatusBar;
javaObjectEDT(jStatusBar);
jProgressBar = javax.swing.JProgressBar;
javaObjectEDT(jProgressBar);
set(jProgressBar, 'Minimum',0, 'Maximum',100, 'Value',0);
jStatusBar.add(jProgressBar,'West');
jRootPane.setStatusBar(jStatusBar);
jStatusBar.setVisible(1);
jStatusBar.setText('Running...');
drawnow
for subject=1:n
    file = deblank(files(subject,:));
    [~,fileName] = fileparts(file);
    EEG = pop_loadset(file);
    EEG = srcConnectivityCont(EEG, metric, winSize, overlap, normalize);
    pop_saveset(EEG,'filepath',dest,'filename',['conn_' fileName],'savemode','onefile');
    jProgressBar.setValue(100*subject/n);
end
jStatusBar.setText('');
jStatusBar.setVisible(0);
jProgressBar.setValue(0);
end

function callback_plot(src,evnt)
fig = src.Parent.Parent;
dest = get(findall(fig,'tag','Results directory'), 'String');
files = pickfiles(dest, '.set');
if ~exist(dest,'dir')
    return;
end
n = size(files,1);
C = [];
for subject=1:n
    file = deblank(files(subject,:));
    EEG = pop_loadset(file);
    C = cat(4,C,EEG.etc.conn.C);
end
h = findall(fig,'style','popupmenu');
metric = lower(h.String{h.Value});

ax = findall(fig,'type','axes');
axes(ax)

T = mean(C,4);
T = bsxfun(@minus,T, mean(T,3));
mx = prctile(abs(T(:)),95);
cnn = Connectivity(T(:,:,1), EEG.etc.src.roi, metric);
cnn.image(ax)

for k=2:length(EEG.etc.conn.times)
    cnn.C = T(:,:,k);
    cnn.image(ax);
    xlabel([num2str(EEG.etc.conn.times(k)) ' ms'])
    caxis ([-mx mx]);
    pause(0.25)
end
end

function callback_help(src,evnt)
web('https://github.com/aojeda/dsi/wiki/How-to-add-new-artifact-scalp-projections-to-the-existing-artifact-dictionary','-browser');
end