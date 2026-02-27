% File: +ObjStretcher/Seg.m
%
% A lightweight value-class representing one constant "meta data" segment
% for Bloch simulation.
%
% Fields: rf, gx, gy, gz are assumed constant over [tStart, tStart+dur).
%
% ADC handling:
%   adc is a uint8 bit-field.
%     bit0 (1): ADC active (inside dwell integration window)
%     bit1 (2): Dwell start (first adcRasterTime inside a dwell)
%     bit2 (4): Dwell end   (last  adcRasterTime inside a dwell)
%
% Notes:
% - Keep this object CPU-side. For GPU simulation, compile numeric arrays
%   (structure-of-arrays) and move those arrays to gpuArray.
% - Uses single precision by default for GPU-friendliness.
%
classdef Seg
    properties
        % Timing
        tStart (1,1) double = 0          % segment start time (s)
        dur    (1,1) double = 0          % duration (s)

        % Constant controls over this segment
        rf (1,1) single = complex(single(0), single(0)) % B1 (in Hz)
        gx (1,1) single = single(0)      % gradient x (your chosen units)
        gy (1,1) single = single(0)      % gradient y
        gz (1,1) single = single(0)      % gradient z

        % bit0 (1): ADC active
        % bit1 (2): dwell start
        % bit2 (4): dwell end
        adc (1,1) uint8 = uint8(0)
    end

    properties (Dependent)
        tEnd
        tCenter
    end

    methods
        function obj = Seg(varargin)
            % Constructor options:
            %   Seg()                                        -> empty
            %   Seg('dur',..., 'rf',..., ...)                -> name/value
            %   Seg(tStart, dur, rf, gx, gy, gz, adc)  -> positional
            %
            if nargin == 0
                return;
            end

            if nargin >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                % Positional signature
                obj.tStart = double(varargin{1});
                obj.dur    = double(varargin{2});

                if nargin >= 3 && ~isempty(varargin{3}), obj.rf  = castComplexSingle(varargin{3}); end
                if nargin >= 4 && ~isempty(varargin{4}), obj.gx  = single(varargin{4}); end
                if nargin >= 5 && ~isempty(varargin{5}), obj.gy  = single(varargin{5}); end
                if nargin >= 6 && ~isempty(varargin{6}), obj.gz  = single(varargin{6}); end
                if nargin >= 7 && ~isempty(varargin{7}), obj.adc = varargin{7}; end

                obj.validateBasic();
                return;
            end

            % Name/value
            if mod(nargin,2) ~= 0
                error('ObjStretcher:Seg:badArgs', 'Name/value inputs must come in pairs.');
            end

            for i = 1:2:nargin
                name = varargin{i};
                val  = varargin{i+1};
                switch lower(string(name))
                    case "tstart"
                        obj.tStart = double(val);
                    case "dur"
                        obj.dur = double(val);
                    case "rf"
                        obj.rf = castComplexSingle(val);
                    case "gx"
                        obj.gx = single(val);
                    case "gy"
                        obj.gy = single(val);
                    case "gz"
                        obj.gz = single(val);
                    case "adc"
                        obj.adc = uint8(val);
                    otherwise
                        error('ObjStretcher:Seg:unknownField', 'Unknown field: %s', name);
                end
            end

            obj.validateBasic();
        end

        function v = get.tEnd(obj)
            v = obj.tStart + obj.dur;
        end

        function v = get.tCenter(obj)
            v = obj.tStart + 0.5 * obj.dur;
        end

        function obj = set.adc(obj, v)

            v = uint8(v);

            % Allowed values:
            % 0  = no ADC
            % 1  = ADC active (mid)
            % 3  = ADC active + start
            % 5  = ADC active + end
            allowed = uint8([0 1 3 5]);

            if ~ismember(v, allowed)
                error('ObjStretcher:Seg:invalidAdc', ...
                    'adc must be one of {0, 1, 3, 5}.');
            end

            obj.adc = v;
        end

        function tf = isAdcActive(obj)
            tf = bitand(obj.adc, uint8(1)) ~= 0;
        end

        function tf = isAdcStart(obj)
            tf = bitand(obj.adc, uint8(2)) ~= 0;
        end

        function tf = isAdcEnd(obj)
            tf = bitand(obj.adc, uint8(4)) ~= 0;
        end

        function validateBasic(obj)
            if ~(isfinite(obj.tStart) && isfinite(obj.dur))
                error('ObjStretcher:Seg:badTime', ...
                    'tStart and dur must be finite.');
            end
            if obj.dur <= 0
                error('ObjStretcher:Seg:badDur', ...
                    'dur must be positive.');
            end

            hasRF  = abs(obj.rf) > 0;
            hasADC = obj.isAdcActive();

            if hasRF && hasADC
                error('ObjStretcher:Seg:rfAdcOverlap', ...
                    'RF and ADC cannot coexist in same segment.');
            end

            % If ADC is not active, start/end bits must not be set.
            if ~hasADC && obj.adc ~= 0
                error('ObjStretcher:Seg:adcInvalid', ...
                    'adc has bits set while ADC is inactive.');
            end
        end

        function validateAgainstSystem(obj, sys, adcDwell, tol)

            if nargin < 4 || isempty(tol)
                tol = 1e-12;
            end

            obj.validateBasic();

            hasRF  = abs(obj.rf) > 0;
            hasADC = obj.isAdcActive();

            % RF alignment
            if hasRF
                if ~isfield(sys, 'rfRasterTime')
                    error('ObjStretcher:Seg:missingRfRaster', ...
                        'sys.rfRasterTime missing.');
                end

                rfRaster = double(sys.rfRasterTime);

                if abs(obj.dur - rfRaster) > tol
                    error('ObjStretcher:Seg:rfRasterMisalign', ...
                        'RF segment duration not aligned to rfRasterTime.');
                end

                % check RF amplitude
                if isfield(sys, 'maxB1') && abs(obj.rf) > sys.maxB1
                    error('ObjStretcher:Seg:rfAmplitudeOverFlow', ...
                        'RF amplitude exceed system limit.');
                end
            end

            % ADC alignment
            if hasADC
                if isempty(adcDwell)
                    error('ObjStretcher:Seg:missingAdcDwell', ...
                        'ADC dwell required for validation.');
                else
                    if ~isfield(sys, 'adcRasterTime')
                        error('ObjStretcher:Seg:missingAdcRaster', ...
                            'sys.adcRasterTime missing.');
                    end

                    % dwell time must be multiple of adcRasterTime
                    adcRaster = double(sys.adcRasterTime);

                    ratio = adcDwell / adcRaster;
                    ratioRounded = round(ratio);

                    if abs(ratio - ratioRounded) > tol
                        error('ObjStretcher:Seg:adcDwellRasterMisalign', ...
                            'ADC dwell must be an integer multiple of adcRasterTime.');
                    end
                end

                % due must be the same as adcRasterTime
                if abs(obj.dur - adcRaster) > tol
                    error('ObjStretcher:Seg:adcDwellMisalign', ...
                        'ADC segment duration not aligned to dwell.');
                end
            end

            % Check gradient value
            if ~isfield(sys, 'maxGrad')
                error('ObjStretcher:Seg:missingMaxGrad', ...
                    'sys.maxGrad missing.');
            end

            if abs(obj.gx) > sys.maxGrad || abs(obj.gy) > sys.maxGrad || abs(obj.gz) > sys.maxGrad
                error('ObjStretcher:Seg:gradOverFlow', ...
                    'Gradient amplitude above system limit.');
            end

        end

        function s = toStruct(obj)
            % Convert to plain struct (useful before compiling to buffers)
            s = struct( ...
                'tStart',  obj.tStart, ...
                'dur',     obj.dur, ...
                'rf',      obj.rf, ...
                'gx',      obj.gx, ...
                'gy',      obj.gy, ...
                'gz',      obj.gz, ...
                'adc',     obj.adc);
        end
    end
end

% ---- local helper ----
function z = castComplexSingle(v)
% Accept real/complex, cast to complex single explicitly.
v = single(v);
if ~isreal(v)
    z = complex(real(v), imag(v));
else
    z = complex(v, single(0));
end
end