classdef Scanner < handle
    % Scanner
    % Wrapper around Pulseq system limits struct returned by mr.opts().
    %
    % Use:
    %   scanner = Scanner.prodModel("Philips.Ingenia1.5T.Ambition");

    properties (SetAccess = private)
        prodId (1,1) string = ""
        sys struct = struct()     % pulseq sys struct from mr.opts
        tol (1,1) double = 1e-12
    end

    methods
        function obj = Scanner(prodId)
            % Constructor is intended to be called by factory functions.
            if nargin >= 1
                obj.prodId = string(prodId); 
                prodMap = Scanner.ProdModelMap.getInstance();
                obj.sys = prodMap.get(prodId);
            end
        end

        function p = getMaxGrad(obj)
            p = obj.sys.maxGrad;
        end

        function p = getMaxSlew(obj)
            p = obj.sys.maxSlew;
        end

        function p = getRfRingdownTime(obj)
            p = obj.sys.rfRingdownTime;
        end

        function p = getRfDeadTime(obj)
            p = obj.sys.rfDeadTime;
        end

        function p = getAdcDeadTime(obj)
            p = obj.sys.adcDeadTime;
        end

        function p = getRfRasterTime(obj)
            p = obj.sys.rfRasterTime;
        end

        function p = getGradRasterTime(obj)
            p = obj.sys.gradRasterTime;
        end

        function p = getAdcRasterTime(obj)
            p = obj.sys.adcRasterTime;
        end

        function p = getBaseRasterTime(obj)
            p = obj.sys.baseRasterTime;
        end

        function p = getMaxB1(obj)
            p = obj.sys.maxB1;
        end

        function txt = describe(obj)
            txt = sprintf("%s | maxGrad=%s Hz/m, maxSlew=%s Hz/m/s, rfRaster=%s s, gradRaster=%s s, adcRaster=%s s, baseRaster=%s s, maxB1=%s", ...
                obj.prodId, ...
                obj.val2str(obj.getMaxGrad()), ...
                obj.val2str(obj.getMaxSlew()), ...
                obj.val2str(obj.getRfRasterTime()), ...
                obj.val2str(obj.getGradRasterTime()), ...
                obj.val2str(obj.getAdcRasterTime()), ...
                obj.val2str(obj.getBaseRasterTime()), ...
                obj.val2str(obj.getMaxB1()));
        end
    end
    methods (Access = private)
        function s = val2str(~, v)
            if isempty(v)
                s = "[]";
            elseif isnumeric(v) && isscalar(v)
                s = char(string(num2str(double(v), "%.6g")));
            else
                s = "<non-scalar>";
            end
        end
    end
end

