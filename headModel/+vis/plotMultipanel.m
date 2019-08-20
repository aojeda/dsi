function fig = plotMultipanel(hm, X)

[hm.fvLeft,hm.fvRight, hm.leftH, hm.rightH] = geometricTools.splitBrainHemispheres(hm.cortex);

fig = figure('Color',[1 1 1]);
fig.Position(3:4) = [642   755];
mx = prctile(abs(X(:)), 95);
%%

ax = subplot(421);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[-mx mx]);
view([-90 90])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(ax,bipolar(512, 0.99))
%xlabel('Superior','FontWeight','bold')
ax.Position([1 3 4]) = [0.1 0.3991 0.2037];


ax = subplot(422);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[-mx mx]);
view([90 -90])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(ax,bipolar(512, 0.99))
% xlabel('Inferior','FontWeight','bold')
ax.Position([1 3 4]) = [0.52 0.3991 0.2037];


ax = subplot(423);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[-mx mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
%title('Left')
%zlabel('Lateral','FontWeight','bold')

ax = subplot(424);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

mx = prctile(abs(X(:)), 95);
set(ax,'Clim',[-mx mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(ax,bipolar(512, 0.99))
%title('Right')


ax = subplot(425);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

set(ax,'Clim',[-mx mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
% zlabel('Medial','FontWeight','bold')

ax = subplot(426);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,1),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[-mx mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))

ax = subplot(427);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X,...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

set(ax,'Clim',[-mx mx]);
view([90 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
ax.Position([1 3 4]) = [0.18 0.2350    0.1024];

ax = subplot(428);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X,...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[-mx mx]);
view([-90 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
ax.Position([1 3 4]) = [0.63 0.2350    0.1024];
set(findall(fig,'type','axes'),'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','off', 'visible','off')
fig.Position(3:4) = [405   755];