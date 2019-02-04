classdef EEGSource < handle
    properties
        x
    end
    properties(GetAccess=private)
        H
        indG
        indV
        gamma
        indGamma
        roi
        P
    end
    properties(Dependent)
        y
        g
        g_roi
        v
    end
    methods
        function y = get.y(obj)
            y = obj.H(:,obj.indG)*obj.x(obj.indG,:,:);
        end
        function g = get.g(obj)
            g = obj.x(obj.indG,:,:);
        end
        function g = get_source_trial(obj,trial)
            g = obj.x(obj.indG,:,trial);
        end
        function g_roi = get.g_roi(obj)
            dim = size(obj.x);
            if length(dim) < 3
                g_roi = obj.P*obj.x(obj.indG,:);
            else
                g_roi = reshape(obj.P*reshape(obj.x(obj.indG,:,:),length(obj.indG),[]),[length(obj.roi), dim(2:3)]);
            end
        end
        function v = get.v(obj)
            v = obj.x(obj.indV,:,:);
        end
        function obj = EEGSource(EEG)
            obj.H = EEG.etc.src.H;
            obj.x = EEG.etc.src.actFull;
            obj.indG = EEG.etc.src.indG;
            obj.indV = EEG.etc.src.indV;
            obj.gamma = EEG.etc.src.gamma;
            obj.indGamma = EEG.etc.src.indGamma;
            hm = headModel.loadFromFile(EEG.etc.src.hmfile);
            obj.P = double(hm.indices4Structure(hm.atlas.label));
            obj.P = sparse(bsxfun(@rdivide,obj.P, sum(obj.P)))';
            if length(obj.indG) == size(hm.cortex.vertices,1)*3
                obj.P = kron(obj.P,[1 1 1]);
            end
        end
    end
end