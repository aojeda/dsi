classdef Connectivity < handle
    properties
        C
        roiName
        measureName
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
            mx = prctile(nonzeros(Cs),95);
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
            set(dcmHandle,'UpdateFcn',@(src,event)showLabel(obj,event));
            dcmHandle.Enable = 'off';
            set(hFigure,'UserData', dcmHandle);
        end
    end
    methods(Hidden)
        function output_txt = showLabel(obj,eventObj)
            pos = get(eventObj,'Position');
            try
                roi_i = obj.roiName{pos(2)};
                roi_j = obj.roiName{pos(1)};
                output_txt = ['(' roi_i ', ' roi_j ', ' num2str(eventObj.Target.CData(pos(2),pos(1))) ')'];
            catch
                output_txt = 'No labeled';
            end
        end
    end
end