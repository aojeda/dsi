classdef Connectivity < handle
    properties
        C
        roiName
        measureName
    end
    properties(Hidden)
        DT
    end
    methods
        function obj = Connectivity(C, roiName, measureName)
            if nargin < 3
                measureName = 'Connectivity';
            end
            obj.C = C;
            obj.roiName = roiName;
            obj.measureName = measureName;
        end
        function t = rankConnections(obj, prcAlpha)
            if nargin < 2
                indnz = find(obj.C);
            elseif isempty(prcAlpha)
                indnz = find(obj.C);
            else
                th = prctile(nonzeros(obj.C),prcAlpha);
                indnz = find(obj.C<th(1) | obj.C>th(2));
            end
            [i,j] = ind2sub(size(obj.C),indnz);
            [~,sorting] = sort(obj.C(indnz),'descend');
            Conn = obj.C(indnz(sorting));
            I = obj.roiName(i(sorting));
            J = obj.roiName(j(sorting));
            t = table(Conn, I, J);
        end
        function ax = image(obj, ax, scaleType)
            if nargin < 2
                hFigure = figure;
                ax = axes;
            else
                hFigure = get(ax,'Parent');
                cla(ax);
            end
            if nargin < 3, scaleType = 'bipolar';end
            Cs = obj.C;
            mx = prctile(abs(nonzeros(Cs)),95);
            mx(mx==0) = max(abs(Cs(:)));
            if strcmp(scaleType,'bipolar')
                scale = [-mx mx];
                color = bipolar(256,0.99);
            else
                scale = [0 mx];
                color = bipolar(256,0.99);
                color = color(129:end,:);
            end
            imagesc(Cs,scale);
            colormap(color);
            colorbar;
            xlabel('Sources');
            ylabel('Sources');
            title(obj.measureName);
            dcmHandle = datacursormode(hFigure);
            dcmHandle.SnapToDataVertex = 'off';
            set(dcmHandle,'UpdateFcn',@(src,event)showConnection(obj,event));
            dcmHandle.Enable = 'off';
            set(hFigure,'UserData', dcmHandle);
        end
        function fig = plotGraph(obj, hm, nodePower, prcAlpha, ax)
            if nargin < 3
                nodePower = [];
            end
            if nargin < 4
                prcAlpha = [];
            end
            if nargin < 5
                fig = figure;
                ax = axes(fig);
            else
                fig = ax.Parent;
            end
            nodePower = nodePower(:);
            hcolor = vis.human(3);
            h = patch(hm.cortex,'FaceColor',hcolor(2,:),'FaceLighting','phong','LineStyle','none','FaceAlpha',0.1);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            material(h,'dull')
            camlight
            axis vis3d equal tight
            view([-90 90]);
            dcmHandle = datacursormode(fig);
            dcmHandle.SnapToDataVertex = 'off';
            set(dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            dcmHandle.Enable = 'off';
            set(fig,'UserData', dcmHandle);
            fig.UserData = hm;
            
            hold(ax,'on');
            [indi,indj] = ind2sub(size(obj.C),find(obj.C~=0));
            src  = hm.getCentroidROI(hm.atlas.label(indi));
            dest = hm.getCentroidROI(hm.atlas.label(indj));
            n = length(indi);
            c = zeros(n,1);
            for k=1:n
                c(k) = obj.C(indi(k),indj(k));
            end
            if ~isempty(prcAlpha)
                th = prctile(c,prcAlpha);
                if length(th) == 2
                    rmind = c>th(1) & c<th(2);
                else
                    rmind = c<th;
                end
                c(rmind) = [];
                src(rmind,:) = [];
                dest(rmind,:) = [];
                indi(rmind) = [];
                indj(rmind) = [];
            end
            rmind = contains(hm.atlas.label(indi),'fusiform') | contains(hm.atlas.label(indj),'fusiform')...
                | contains(hm.atlas.label(indi),'temporalpole') | contains(hm.atlas.label(indj),'temporalpole')...
                | contains(hm.atlas.label(indi),'entorhinal') & contains(hm.atlas.label(indj),'entorhinal');
            
            c(rmind) = [];
            src(rmind,:) = [];
            dest(rmind,:) = [];
            indi(rmind) = [];
            indj(rmind) = [];
            n = length(indi);
            
            %c = sign(c);
            mx = prctile(abs(c),95);
            clim = [-mx mx];
            
            lw = griddedInterpolant(linspace(0,max(clim),256)',logspace(-0.1,0.2,256)','linear','nearest');
            cm = griddedInterpolant(linspace(clim(1),clim(2),256)',1:256,'nearest','nearest');
            cmap = bipolar(256,0.8);
            lineColor = cmap(cm(c),:);
            lineWidth = lw(abs(c));
            for k=1:n
                hl = line([src(k,1) dest(k,1)],[src(k,2) dest(k,2)],[src(k,3) dest(k,3)],'linewidth',lineWidth(k),'color',lineColor(k,:));
                hl.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
            if isempty(nodePower) || sum(abs(nodePower))<eps
                croi = flipud(jet(length(hm.atlas.label)));
                center = hm.getCentroidROI(hm.atlas.label);
                [~,sorting] = sort(center(:,1),'descend');
                center = center(sorting,:);
                atlasLabel = hm.atlas.label(sorting);
                roi = unique(hm.atlas.label([indi;indj]),'stable');
                loc = find(ismember(atlasLabel,roi));
                roi = atlasLabel(loc);
                center = center(loc,:);
                croi = croi(loc,:);
                r = 0.005*ones(length(roi),1);
            else
                roi = hm.atlas.label;
                center = hm.getCentroidROI(hm.atlas.label);
                mx = prctile(abs(nonzeros(nodePower)),95);
                mx(mx==0) = max([abs(nodePower);1]);
                cm = griddedInterpolant(linspace(-mx,mx,256)',1:256,'nearest','nearest');
                croi = cmap(cm(nodePower),:);
                rm = griddedInterpolant(linspace(0,mx,256)',linspace(2.5e-3,7.5e-3,256),'linear','nearest');
                r = rm((nodePower));
%                 roi = unique(hm.atlas.label([indi;indj]),'stable');
%                 loc = ismember(hm.atlas.label,roi);
%                 center = hm.getCentroidROI(hm.atlas.label(loc));
%                 nodePower = nodePower(loc);
%                 mx = prctile(abs(nodePower),95);
%                 mx(mx==0) = max([abs(nodePower);1]);
%                 cm = griddedInterpolant(linspace(-mx,mx,256)',1:256,'nearest','nearest');
%                 croi = cmap(cm(nodePower),:);
%                 rm = griddedInterpolant(linspace(0,mx,256)',linspace(2.5e-3,7.5e-3,256),'linear','nearest');
%                 r = rm(abs(nodePower));
            end
            hs = [];
            for k=1:length(roi)
                [sx,sy,sz] = ellipsoid(center(k,1),center(k,2),center(k,3),r(k),r(k),r(k));
                hs(end+1) = surf(sx,sy,sz,'LineStyle','None','FaceColor',croi(k,:));
            end
            hold(ax,'off')
            axis off;
        end
        function fig = plotDirectedGraph(obj, hm, prcAlpha, ax, invertColorScale)
            if nargin < 3
                prcAlpha = [];
            end
            if nargin < 4
                fig = figure('Name',obj.measureName);
                ax = axes(fig);
            elseif isempty(ax)
                fig = figure('Name',obj.measureName);
                ax = axes(fig);
            else
                fig = ax.Parent;
            end
            if nargin < 5
                invertColorScale = false;
            end
            hcolor = vis.human(3);
            h = patch(hm.cortex,'FaceColor',hcolor(2,:),'FaceLighting','phong','LineStyle','none','FaceAlpha',0.1);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            material(h,'dull')
            camlight
            axis vis3d equal tight
            view([-90 90]);
            dcmHandle = datacursormode(fig);
            dcmHandle.SnapToDataVertex = 'off';
            set(dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            dcmHandle.Enable = 'off';
            set(fig,'UserData', dcmHandle);
            fig.UserData = hm;
            
            hold(ax,'on');
            Ca = triu(obj.C-obj.C');
            S = sign(Ca);
            Ca = abs(Ca);
            [indi,indj] = ind2sub(size(obj.C),find(Ca~=0));
            dest  = hm.getCentroidROI(hm.atlas.label(indi));
            src = hm.getCentroidROI(hm.atlas.label(indj));
            n = length(indi);
            c = zeros(n,1);
            for k=1:n
                c(k) = Ca(indi(k),indj(k));
            end
            if ~isempty(prcAlpha)
                th = prctile(c,prcAlpha);
                if length(th) == 2
                    rmind = c>th(1) & c<th(2);
                else
                    rmind = c<th;
                end
                c(rmind) = [];
                src(rmind,:) = [];
                dest(rmind,:) = [];
                indi(rmind) = [];
                indj(rmind) = [];
            end
            rmind = contains(hm.atlas.label(indi),'fusiform') | contains(hm.atlas.label(indj),'fusiform')...
                | contains(hm.atlas.label(indi),'temporalpole') | contains(hm.atlas.label(indj),'temporalpole')...
                | contains(hm.atlas.label(indi),'frontalpole') | contains(hm.atlas.label(indj),'frontalpole')...
                | contains(hm.atlas.label(indi),'entorhinal') & contains(hm.atlas.label(indj),'entorhinal');
            
            c(rmind) = [];
            src(rmind,:) = [];
            dest(rmind,:) = [];
            indi(rmind) = [];
            indj(rmind) = [];
            n = length(c);
            
            c = sign(c);
            mx = prctile(abs(c),95);
            clim = [0 mx];
            
            lw = griddedInterpolant(linspace(0,max(clim),256)',logspace(-0.1,0.2,256)','linear','nearest');
            cmap = bipolar(256,0.8);
            if invertColorScale
                cmap = flipud(cmap(1:128,:));
            else
                cmap = cmap(129:end,:);
            end
            cm = griddedInterpolant(linspace(clim(1),clim(2),128)',1:128,'nearest','nearest');
            
            lineWidth = lw(abs(c));
            colormap(ax,cmap);
            cn = 50;
            hl = [];
            for k=1:n
                cij = [c(k),0];
                if S(indi(k),indj(k))>=0    
                    hl(k) = clinep(linspace(src(k,1),dest(k,1),cn)',linspace(src(k,2),dest(k,2),cn)',linspace(src(k,3),dest(k,3),cn)',cm(linspace(cij(1),cij(2),cn)),lineWidth(k));
                else
                    hl(k) = clinep(linspace(dest(k,1),src(k,1),cn)',linspace(dest(k,2),src(k,2),cn)',linspace(dest(k,3),src(k,3),cn)',cm(linspace(cij(1),cij(2),cn)),lineWidth(k));
                end
            end
            %nodePower = abs(sum(Ca,1)');
            nodePower = abs(sum(sign(obj.C),1)');
            %nodePower = sum(obj.C,2);
            roi = hm.atlas.label;
            center = hm.getCentroidROI(hm.atlas.label);
            mx = prctile(abs(nonzeros(nodePower)),95);
            mx(mx==0) = max([abs(nodePower);1]);
            cm = griddedInterpolant(linspace(0,mx,128)',1:128,'nearest','nearest');
            croi = cmap(cm(nodePower),:);
            rm = griddedInterpolant(linspace(0,mx,256)',linspace(1e-3,5e-3,256),'linear','nearest');
            r = rm((nodePower));
            hs = [];
            for k=1:length(roi)
                [sx,sy,sz] = ellipsoid(center(k,1),center(k,2),center(k,3),r(k),r(k),r(k));
                hs(end+1) = surf(sx,sy,sz,'LineStyle','None','FaceColor',croi(k,:));
            end
            hold(ax,'off')
            axis off;
        end
        
    end
    methods(Hidden)
        function drawConnection(obj, src, dest, cij)
            d = sqrt(sum((src-dest).^2));
            r = linspace(cij(1),cij(2));
            th = linspace(0,2*pi);
            [R,T] = meshgrid(r,th);
            R = flipud(R);
            X = R.*cos(T) ;
            Y = R.*sin(T) ;
            Z = d*fliplr(R);
            v = dest-src;
            v = v/(norm(v)+eps);
            z = [0, 0, 1];
            u = cross(z,v);
            u = u/(norm(u)+eps);
            ang = acos(dot(z, v));
            A = zeros(3);
            A(1, 2) = -u(3);
            A(1, 3) = u(2);
            A(2, 1) = u(3);
            A(2, 3) = -u(1);
            A(3, 1) = -u(2);
            A(3, 2) = u(1);
            R = eye(3)+sin(ang)*A+(1-cos(ang))*A^2;
            XYZ = R*[X(:)';Y(:)';Z(:)'];
            X = reshape(XYZ(1,:),size(X)) + src(1);
            Y = reshape(XYZ(2,:),size(Y)) + src(2);
            Z = reshape(XYZ(3,:),size(Z)) + src(3);
            surf(X,Y,Z);xlabel('X');ylabel('Y');zlabel('Z')
        end
        function output_txt = showConnection(obj,eventObj)
            pos = get(eventObj,'Position');
            try
                roi_i = obj.roiName{pos(2)};
                roi_j = obj.roiName{pos(1)};
                output_txt = ['(' roi_i ', ' roi_j ', ' num2str(eventObj.Target.CData(pos(2),pos(1))) ')'];
            catch
                output_txt = 'No labeled';
            end
        end
        function output_txt = showLabel(obj,event_obj)
            try
                hm = get(gcf,'UserData');
            catch
                return
            end
            if isempty(obj.DT)
                vertices = hm.cortex.vertices;
                obj.DT = delaunayTriangulation(vertices(:,1),vertices(:,2),vertices(:,3));
            end
            pos = get(event_obj,'Position');
            loc = nearestNeighbor(obj.DT, pos);
            try
                output_txt = hm.atlas.label{hm.atlas.colorTable(loc)};
            catch
                output_txt = 'No labeled';
            end
        end
    end
end