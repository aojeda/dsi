classdef LargeTensor < handle
    properties(SetAccess=protected)
        mmf
    end
    properties(Hidden,SetAccess=protected)
        file
    end
    methods
        function obj = LargeTensor(dims, filename)
            if nargin < 1, error('Tensor size is the first argument! Example:  a = LargeTensor([2,3,4]);');end
            if nargin < 2
                filename = tempname;
            end
            [fpath, fname] = fileparts(filename);
            free = java.io.File(fpath).getFreeSpace();
            req = prod(dims)*8;
            if free > req
                obj.file = filename;
            elseif java.io.File(getHomeDir).getFreeSpace() > req
                obj.file = fullfile(getHomeDir, fname);
            else
                error('No space left on device');
            end
            if ~exist(obj.file,'file')
                try
                    rt = system(['fallocate -l ' num2str(prod(dims)*8) ' ' obj.file]);
                    if rt ~=0, error('fallocate is not installed.');end
                catch
                    disp('Creating mmf file...')
                    fid = fopen(obj.file,'w');
                    z = zeros(prod(dims(setdiff(1:length(dims),2))),1);
                    for k=1:dims(2)
                        fwrite(fid, z, 'double');
                    end
                    fclose(fid);
                    disp('done.')
                end
            end
            obj.mmf = memmapfile(obj.file,'Format',{'double' dims 'x'},'Writable',true);
        end
        function delete(obj)
            if exist(obj.file,'file')
                delete(obj.file);
                [p,n] = fileparts(obj.file);
                hdr = fullfile(p,[n '.hdr']);
                bin = fullfile(p,[n '.bin']);
                if exist(hdr,'file'), delete(hdr);end
                if exist(bin,'file'), delete(bin);end
            end
        end
        function slice = subsref(obj,s) %#ok
            ind = '';
            n = length(s.subs);
            for k=1:n
                if ischar(s.subs{k})
                    ind = [ind s.subs{k}];
                else
                    ind = cat(2,ind,['[' num2str(s.subs{k}(:)') ']']);
                end
                if k<n
                    ind(end+1) = ',';
                end
            end
            cmd = ['slice=obj.mmf.Data.x(' ind ');'];
            eval(cmd);
        end
        function slice = subsasgn(obj,s,value) %#ok
            ind = '';
            for k=1:length(s.subs)
                if ischar(s.subs{k})
                    ind = [ind ',' s.subs{k}];
                else
                    ind = [ind ',[' num2str(s.subs{k}) ']'];
                end
            end
            ind(ind(1) == ',') = []; 
            cmd = ['obj.mmf.Data.x(' ind ') = value;'];
            eval(cmd);
            slice = obj;
            % slice = subsasgn(obj.mmf.Data.x,s,value);
        end
        function saveToFile(obj,filename)
            [p,n] = fileparts(filename);
            dims = size(obj);
            hdr = fullfile(p,[n '.hdr']);
            fid = fopen(hdr,'w');
            fprintf(fid,'[%s]',num2str(dims));
            fclose(fid);
            bin = fullfile(p,[n '.bin']);
            copyfile(obj.mmf.Filename,bin);
        end
        %%
        function dims = size(obj,d)
            if nargin <2, d = [];end
            dims = obj.mmf.Format{2};
            if ~isempty(d), dims = dims(d);end
        end
        function obj = reshape(obj,dims)
            obj.mmf.Format{2} = dims;
        end
        function obj = minus(obj,value)
            obj.mmf.Data.x = obj.mmf.Data.x - value.mmf.Data.x;
        end
        function obj = plus(obj,value)
            obj.mmf.Data.x = obj.mmf.Data.x + value.mmf.Data.x;
        end
        function obj = times(obj,value)
            obj.mmf.Data.x = obj.mmf.Data.x.*value.mmf.Data.x;
        end
        function obj = power(obj,value)
            obj.mmf.Data.x = obj.mmf.Data.x.^value;
        end
        function obj = sqrt(obj)
            obj.mmf.Data.x = sqrt(obj.mmf.Data.x);
        end
        function obj = mtimes(obj,value)
            % Handle the multiplication by a matrix on the left
            if isa(value,'LargeTensor') && ~isa(obj,'LargeTensor')
                tmp = obj;
                obj = value;
                value = tmp;
                clear tmp;
                obj = value*obj.mmf.Data.x;
                return
            end
            if ~isa(value,class(obj))
                obj.mmf.Data.x = obj.mmf.Data.x*value;
            else
                obj.mmf.Data.x = obj.mmf.Data.x.*value.mmf.Data.x;
            end
        end
        function obj = mrdivide(obj,value)
            if ~isa(value,class(obj))
                obj.mmf.Data.x = obj.mmf.Data.x/value;
            else
                obj.mmf.Data.x = obj.mmf.Data.x/value.mmf.Data.x;
            end
        end
    end
    methods(Static)
        function obj = loadFromFile(filename)
            [p,n] = fileparts(filename);
            hdr = fullfile(p,[n '.hdr']);
            bin = fullfile(p,[n '.bin']);
            fid = fopen(hdr,'r');
            dims = eval(fgets(fid));
            obj = LargeTensor(dims,bin);
        end
    end
end

function homeDir = getHomeDir
if ispc
    homeDir= getenv('USERPROFILE');
else
    homeDir = getenv('HOME');
end
end