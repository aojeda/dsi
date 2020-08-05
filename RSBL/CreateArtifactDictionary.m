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
sel.Position = [0.4865         0    0.0920    0.8000];
sel_tx = uicontrol(pn,'style','edit','string','','Units','normalized');
sel_tx.Position = [0.5872    0.0400    0.1973    0.8800];
chk = uicontrol(pn,'style','checkbox','Units','normalized', 'string','Overwrite');
chk.Position = [0.8124    0.04    0.1973    0.8800];
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
[M, Y, C, L, ICAact, ICApsd] = buildArtifactDictionary({dataFolder}, jProgressBar);
jStatusBar.setText('');
jStatusBar.setVisible(0);
jProgressBar.setValue(0);
fig.UserData.M = M;
fig.UserData.C = C;
fig.UserData.L = L;
fig.UserData.Y = Y;
fig.UserData.ICAact = ICAact;
fig.UserData.ICApsd = ICApsd;
fig.UserData.Cact = [];
fig.UserData.Cpsd = [];
end

function callback_plot(src,evnt)
fig = src.Parent.Parent;
if isempty(fig.UserData.M)
    return
end
bt = findall(fig,'string','K-MEANS');
callback_clustering(bt);
end

function [cact, cpsd] = select_centroid_act(M, L, C, ICAact, ICApsd)
cact = cell(size(C,2),1);
cpsd = cell(size(C,2),1);
for k=1:size(C,2)
    sel = find(L==k);
    [~, loc_mn] = min(sum((M(:, sel) - C(:,k)).^2, 1));
    %cact{k} = ICAact{sel(loc_mn)};
    %cpsd{k} = ICApsd{sel(loc_mn)};
    cact{k} = ICAact(sel);
    cpsd{k} = ICApsd(sel);
end
end

function plot_ic(src, ~)
figure; 
subplot(211);
for i=1:length(src.UserData)
    plot(src.UserData(i).ic.time, src.UserData(i).ic.data);
    if i==1
        hold on;
    end
end
xlabel('Time (sec)')
ylabel('IC')
grid()
subplot(212);
for i=1:length(src.UserData)
    loglog(src.UserData(i).psd.freq, src.UserData(i).psd.data);
    if i==1
        hold on;
    end
end
xlabel('Frequency (Hz)')
ylabel('dB/Hz')
grid()
end

function callback_clustering(src,evnt)
fig = src.Parent.Parent;
fig.Pointer='watch';
tx = findall(findall(fig,'type','uipanel'),'style','edit');
nc = str2double(get(tx(end-1),'String'));
[L, C, e] = kmeans(fig.UserData.M',nc,'distance','sqeuclidean', 'Replicates',9);
C = C';
[fig.UserData.cact, fig.UserData.cpsd] = select_centroid_act(fig.UserData.M,L,C,fig.UserData.ICAact, fig.UserData.ICApsd);
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
    h = topoplot(mean(M(:,L==i),2),chanlocs,'electrodes','off');
    set(h, 'ButtonDownFcn', @plot_ic, 'UserData', struct('ic', fig.UserData.cact{i}, 'psd', fig.UserData.cpsd{i}))
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
fig.Pointer='arrow';
end

function callback_save(src,evnt)
fig = src.Parent.Parent;
C = fig.UserData.C;
tx = findall(fig,'style','edit');
overwrite = findall(fig,'style','checkbox');
ind = str2num(tx(1).String); %#ok
load('Artifact_dictionary.mat'); %#ok
if overwrite.Value
    [p,n] = fileparts(which('Artifact_dictionary.mat'));
    if ~exist(fullfile(p,[n '_orig.mat']), 'file')
        save(fullfile(p,[n '_orig.mat']),'A','templateName'); %#ok
    end
    A = C(:,ind);
    save(which('Artifact_dictionary.mat'),'A','templateName'); 
else
    Ai = C(:,ind);
    save(which('Artifact_dictionary.mat'),'A','Ai','templateName');
end 
end

function callback_help(src,evnt)
web('https://github.com/aojeda/dsi/wiki/How-to-add-new-artifact-scalp-projections-to-the-existing-artifact-dictionary','-browser');
end