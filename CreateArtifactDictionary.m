function fig = CreateArtifactDictionary
fig = findall(0,'Name','Create Artifact Dictionary');
if ~isempty(fig)
    figure(fig);
    return;
end
fig = figure('Name','Create Artifact Dictionary','NumberTitle',false,'Menu','none','Toolbar','none','Units','normalized');
fig.Position(3:4) = [0.3536    0.5157];
t = uitoolbar(fig);
fig.UserData = struct('M',[],'Y',[],'C',[],'L',[]);

p = fileparts(which('headModel.m'));
p = fullfile(p,'+vis','icons');
playIcon = imread([p filesep 'Gnome-media-playback-start.svg.png']);
saveIcon = imread([p filesep 'Dialog-apply.svg.png']);
helpIcon = imread([p filesep 'Gnome-help-browser.svg.png']);
plotIcon = imread([p filesep 'Gnome-emblem-photos.svg.png']);
                     
uitoggletool(t,'CData',playIcon,'TooltipString','Run','OnCallback',@callback_run,'OffCallback',@callback_run);
uitoggletool(t,'CData',plotIcon,'TooltipString','Vualize','OnCallback',@callback_plot,'OffCallback',@callback_plot);
uitoggletool(t,'CData',saveIcon,'TooltipString','Commit','OnCallback',@callback_save,'OffCallback',@callback_save);         
uitoggletool(t,'CData',helpIcon,'TooltipString','Help','OnCallback',@callback_help,'OffCallback',@callback_help);
            
pnt = uipanel(fig,'Units','normalized','Position',[0.0000    0.9428    1.0000    0.0570]);
uicontrol(pnt,'style','text','String','Data directory','Units','Normalized','Position',[0.0057    0.0542    0.2075    0.7001]);
uicontrol(pnt,'style','edit','Units','Normalized',                          'Position',[0.2095    0.1015    0.6188    0.8235]);
uicontrol(pnt,'style','pushbutton','String','Select','Units','Normalized',  'Position',[0.8302    0.0710    0.1472    0.8235],'callback',@callback_select);

ax = axes(fig, 'Units','normalized');
ax.Position = [0.1300    0.2    0.75 0.65];
pn = uipanel(fig,'Units','normalized');
pn.Position = [0      0    1 0.0529];
lb = uicontrol(pn,'style','text','string','Number of clusters: ','Units','normalized');
lb.Position = [0.0236    0.0000    0.1985    0.8000];
tx = uicontrol(pn,'style','edit','string','8','Units','normalized');
tx.Position = [0.2102    0.0400    0.0844    0.8800];
bt = uicontrol(pn,'style','pushbutton','string','K-MEANS','Units','normalized','callback',@callback_clustering);
bt.Position = [0.3252         0    0.1388    1.0000];
sel = uicontrol(pn,'style','text','string','Selection: ','Units','normalized');
sel.Position = [0.5946    0.0    0.092    0.8000];
sel_tx = uicontrol(pn,'style','edit','string','','Units','normalized');
sel_tx.Position = [0.6865    0.0400    0.1973    0.8800];
end

function callback_select(src,evnt)
fig = src.Parent;
selpath = uigetdir(pwd, 'Select data folder');
if selpath
    set(findall(fig,'style','edit'),'String',selpath);
end
end

function callback_run(src,evnt)
fig = src.Parent.Parent;
h = findall(fig,'style','edit');
dataFolder = get(h(end),'String');
if ~exist(dataFolder,'dir')
    return;
end
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
[M, Y, C, L] = buildArtifactDictionary({dataFolder}, jProgressBar);
jStatusBar.setText('');
jStatusBar.setVisible(0);
jProgressBar.setValue(0);
fig.UserData.M = M;
fig.UserData.C = C;
fig.UserData.L = L;
fig.UserData.Y = Y;
end

function callback_plot(src,evnt)
fig = src.Parent.Parent;
if isempty(fig.UserData.M)
    return
end
bt = findall(fig,'string','K-MEANS');
callback_clustering(bt);
end

function callback_clustering(src,evnt)
fig = src.Parent.Parent;
tx = findall(findall(fig,'type','uipanel'),'style','edit');
nc = str2double(get(tx(end-1),'String'));
[L, C] = kmeans(fig.UserData.M',nc,'distance','sqeuclidean', 'Replicates',9);
C = C';
fig.UserData.C = C;
fig.UserData.L = L;
Y = tsne(fig.UserData.M','Distance','correlation','NumDimensions',2,'Verbose',1,'Options',struct('MaxIter',500,'OutputFcn',[],'TolFun',1e-10),'Perplexity',50);
fig.UserData.Y = Y;
M = fig.UserData.M;
color = jet(nc);
ax = findall(fig,'type','axes');
if length(ax)>1
    delete(ax(1:end-1));
    ax = ax(end);
end
cla(ax)
scatter(ax,fig.UserData.Y(:,1),fig.UserData.Y(:,2),15,color(fig.UserData.L,:),'filled','Marker','o','MarkerFaceAlpha',0.5);
colormap(color)
cb = colorbar(ax);
cb.Ticks = linspace(0,1,nc);cb.TickLabels = num2cell(1:nc);
ang = linspace(0,360*(nc-1)/nc,nc)*pi/180;%+5*pi/6;
template = headModel.loadDefault;
chanlocs = template.makeChanlocs;
AXi = [];
for i=1:nc
    % [ax_x,ax_y] = pol2cart(ang(i),r(i));
    [ax_x,ax_y] = pol2cart(ang(i),0.85);
    axi = axes(fig,'Position',[ax_x/2+0.5-0.05,ax_y/2+0.5-0.05,0.1,0.1]);
    axes(axi); %#ok
    topoplot(mean(M(:,L==i),2),chanlocs,'electrodes','off');
    title(num2str(i))
    AXi = [AXi axi];
    axis(axi,'on')
    set(axi,'XTickLabel',[],'YTickLabel',[],'Box','on');
    axi.YAxis.Color = color(i,:);
    axi.XAxis.Color = color(i,:);
end
for i=1:nc
    colormap(AXi(i),bipolar(256,0.8));
end
end

function callback_save(src,evnt)
fig = src.Parent.Parent;
C = fig.UserData.C;
tx = findall(fig,'style','edit');
ind = str2num(tx(1).String); %#ok
Ai = C(:,ind);
load('Artifact_dictionary.mat'); %#ok
save(which('Artifact_dictionary.mat'),'A','Ai','templateName'); %#ok
end

function callback_help(src,evnt)
web('https://github.com/aojeda/dsi/wiki/How-to-add-new-artifact-scalp-projections-to-the-existing-artifact-dictionary','-browser');
end