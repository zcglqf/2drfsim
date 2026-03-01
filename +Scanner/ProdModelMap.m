classdef ProdModelMap < handle
    % ProdModelMap
    % Case-insensitive map from product model string -> Pulseq sys struct

    properties (Access = private)
        map containers.Map
        tol (1,1) double = 1e-12
    end

    methods (Access = private)
        function obj = ProdModelMap()
            obj.map = containers.Map('KeyType','char','ValueType','any');
            obj.initialize();
        end
    end

    methods (Access = private)
        function initialize(obj)
            % Add product entries here
            % "any" baseline model
            sys = mr.opts( ...
                'MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 1e-5, ...
                'adcRasterTime', 1e-7);
            sys.baseRasterTime = 1e-7;   % make sure all raster times are integer multiple of baseRaster.
            obj.add("any", sys);

            % Philips Ingenia 1.5T Ambition, reuse baseline for now
            sys = mr.opts( ...
                'MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 1e-5, ...
                'adcRasterTime', 1e-7);
            sys.baseRasterTime = 1e-7;   % make sure all raster times are integer multiple of baseRaster.
            obj.add("Philips.Ingenia1.5T.Ambition", sys);

            % Siemens Prisma 3.0T
            sys = mr.opts( ...
                'MaxGrad', 80, 'GradUnit', 'mT/m', ...
                'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 30e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 10e-6, ...
                'adcRasterTime', 0.1e-6);
            sys.baseRasterTime = 1e-7;
            obj.add("Siemens.Prisma3.0T", sys);

            % GE Signa 3.0T Hero
            sys = mr.opts( ...
                'MaxGrad', 50, 'GradUnit', 'mT/m', ...
                'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 10e-6, ...
                'adcRasterTime', 0.1e-6);
            sys.baseRasterTime = 1e-7;
            obj.add("GE.Signa3.0T.Hero", sys);
        end
    end

    methods
        function add(obj, prodId, sys)
            % Validate sys.baseRasterTime and raster multiples
            obj.validateSysRaster(sys, prodId);

            key = obj.normalizeKey(prodId);
            obj.map(key) = sys;
        end

        function sys = get(obj, prodId)
            key = obj.normalizeKey(prodId);

            if ~isKey(obj.map, key)
                error("Scanner:UnknownProdModel", ...
                    "Unknown product model: %s", prodId);
            end

            sys = obj.map(key);
        end
    end

    methods (Static)
        function instance = getInstance()
            persistent singleton
            if isempty(singleton)
                singleton = Scanner.ProdModelMap();
            end
            instance = singleton;
        end
    end

    methods (Access = private)
        function key = normalizeKey(~, prodId)
            key = char(lower(strtrim(string(prodId))));
        end

        function validateSysRaster(obj, sys, prodId)
            if ~isstruct(sys)
                error("Scanner:BadSys", "sys must be a struct for %s.", string(prodId));
            end

            if ~isfield(sys, "baseRasterTime") || isempty(sys.baseRasterTime)
                error("Scanner:MissingBaseRasterTime", ...
                    "sys.baseRasterTime missing for %s.", string(prodId));
            end

            base = double(sys.baseRasterTime);
            if ~(isfinite(base) && base > 0)
                error("Scanner:BadBaseRasterTime", ...
                    "sys.baseRasterTime must be finite and >0 for %s.", string(prodId));
            end

            % Required raster fields
            required = ["rfRasterTime","gradRasterTime","adcRasterTime"];
            for i = 1:numel(required)
                f = required(i);
                if ~isfield(sys, f) || isempty(sys.(f))
                    error("Scanner:MissingRasterField", ...
                        "sys.%s missing for %s.", f, string(prodId));
                end
            end

            % Check integer multiple within tolerance
            obj.assertMultiple(sys.rfRasterTime, base, "rfRasterTime", prodId);
            obj.assertMultiple(sys.gradRasterTime, base, "gradRasterTime", prodId);
            obj.assertMultiple(sys.adcRasterTime, base, "adcRasterTime", prodId);
        end

        function assertMultiple(obj, value, base, fieldName, prodId)
            ratio = double(value) / double(base);
            ratioRounded = round(ratio);
            if abs(ratio - ratioRounded) > obj.tol
                error("Scanner:RasterNotMultiple", ...
                    "sys.%s (%.15g) is not an integer multiple of baseRasterTime (%.15g) for %s.", ...
                    fieldName, double(value), double(base), string(prodId));
            end
        end
    end
end
