classdef Scanner < handle
    % Scanner
    % Wrapper around Pulseq system limits struct returned by mr.opts().
    %
    % Use:
    %   scanner = Scanner.prodModel("Philips.Ingenia1.5T.Ambition");

    properties (SetAccess = private)
        prodId (1,1) string = ""
        vendor (1,1) string = ""
        series (1,1) string = ""
        fieldStrengthT (1,1) double = NaN
        model (1,1) string = ""

        sys struct = struct()     % pulseq sys struct from mr.opts
        tol (1,1) double = 1e-12
    end

    methods
        function obj = Scanner(prodId, sysStruct, parts)
            % Constructor is intended to be called by factory functions.
            if nargin >= 1, obj.prodId = string(prodId); end
            if nargin >= 2, obj.sys = sysStruct; end
            if nargin >= 3 && isstruct(parts)
                obj.vendor = string(parts.vendor);
                obj.series = string(parts.series);
                obj.fieldStrengthT = double(parts.fieldStrengthT);
                obj.model  = string(parts.model);
            end
        end

        function s = getSys(obj)
            s = obj.sys;
        end

        function spec = getSpec(obj)
            % Return normalized specs (unit assumptions depend on how you created sys)
            spec = struct();
            spec.prodId = obj.prodId;
            spec.vendor = obj.vendor;
            spec.series = obj.series;
            spec.fieldStrengthT = obj.fieldStrengthT;
            spec.model = obj.model;

            % Common Pulseq sys fields (may or may not exist depending on mr.opts version)
            spec.maxGrad = getFieldOr(obj.sys, "maxGrad", []);
            spec.maxSlew = getFieldOr(obj.sys, "maxSlew", []);
            spec.rfRasterTime   = getFieldOr(obj.sys, "rfRasterTime", []);
            spec.gradRasterTime = getFieldOr(obj.sys, "gradRasterTime", []);
            spec.adcRasterTime  = getFieldOr(obj.sys, "adcRasterTime", []);
            spec.rfDeadTime     = getFieldOr(obj.sys, "rfDeadTime", []);
            spec.rfRingdownTime = getFieldOr(obj.sys, "rfRingdownTime", []);
            spec.adcDeadTime    = getFieldOr(obj.sys, "adcDeadTime", []);

            % Optional: if you define these in your model db
            spec.maxB1 = getFieldOr(obj.sys, "maxB1", []);
        end

        function txt = describe(obj)
            spec = obj.getSpec();
            txt = sprintf("%s | maxGrad=%s, maxSlew=%s, rfRaster=%s, adcRaster=%s", ...
                obj.prodId, val2str(spec.maxGrad), val2str(spec.maxSlew), ...
                val2str(spec.rfRasterTime), val2str(spec.adcRasterTime));
        end
    end

    methods (Static)
        function obj = prodModel(prodId)
            % Convenience static wrapper: Scanner.prodModel(...)
            obj = Scanner.prodModel(prodId); %#ok<MCNPN>
        end
    end
end

% ---- local helpers ----
function v = getFieldOr(s, name, defaultV)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = defaultV;
end
end

function s = val2str(v)
if isempty(v)
    s = "[]";
elseif isnumeric(v) && isscalar(v)
    s = string(num2str(v, "%.4g"));
else
    s = "<non-scalar>";
end
end