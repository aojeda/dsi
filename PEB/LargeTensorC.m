classdef LargeTensorC
    properties(Hidden,SetAccess=protected)
        Real
        Imag
    end
    methods 
        function obj = LargeTensorC(dims)
            obj.Real = LargeTensor(dims);
            obj.Imag = LargeTensor(dims);
        end
        function delete(obj)
            delete(obj.Real);
            delete(obj.Imag);
        end
        function slice = subsref(obj,s)
            slice = subsref(obj.Real,s) + 1i*subsref(obj.Imag,s);
        end
        function slice = subsasgn(obj,s,value)
            obj.Real = subsasgn(obj.Real,s,real(value));
            obj.Imag = subsasgn(obj.Imag,s,imag(value));
            slice = obj;
        end
        function dims = size(obj,d)
            dims = size(obj.Real,d);
        end
        function obj = reshape(obj,dims)
            obj.Real = reshape(obj.Real,dims);
            obj.Imag = reshape(obj.Imag,dims);
        end
        function obj = minus(obj,value)
            obj.Real = minus(obj.Real,real(value));
            obj.Imag = minus(obj.Imag,imag(value));
        end
        function obj = plus(obj,value)
            obj.Real = plus(obj.Real,real(value));
            obj.Imag = plus(obj.Imag,imag(value));
        end
    end
end