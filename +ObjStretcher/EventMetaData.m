classdef EventMetaData
    % Represents a small segment of an MR event
    % (RF, Gradient, ADC, etc.)
    %
    % Units:
    %   Time         : seconds (double)
    %   Value        : single complex (Hz or Hz/m depending on type)
    %
    % EventType:
    %   ObjStretcher.EventType enum
    properties
        type (1,1) ObjStretcher.EventType = ObjStretcher.EventType.RF

        tStart    (1,1) double = 0
        vStart   (1,1) single = complex(single(0))

        tEnd      (1,1) double = 0
        vEnd     (1,1) single = complex(single(0))

        vAve (1,1) single = complex(single(0))

    end

    properties (Dependent)
        dur
    end

    methods
        function obj = EventMetaData(varargin)
            % Name/value constructor

            if nargin == 0
                return;
            end

            if mod(nargin,2) ~= 0
                error('ObjStretcher:EventMetaData:badArgs', ...
                      'Name/value inputs must come in pairs.');
            end

            for i = 1:2:nargin
                name = lower(string(varargin{i}));
                val  = varargin{i+1};

                switch name
                    case "type"
                        obj.type = val;

                    case "tstart"
                        obj.tStart = double(val);

                    case "vstart"
                        obj.vStart = castComplexSingle(val);

                    case "tend"
                        obj.tEnd = double(val);

                    case "vend"
                        obj.vEnd = castComplexSingle(val);

                    case "vave"
                        obj.vAve = castComplexSingle(val);

                    case "dur"
                        obj.dur = double(val);

                    otherwise
                        error('ObjStretcher:EventMetaData:unknownField', ...
                              'Unknown field: %s', name);
                end
            end
            validate(obj);
        end

        function v = get.dur(obj)
            v = obj.tEnd - obj.tStart;
        end

        function validate(obj)
            if obj.tEnd < obj.tStart
                error('ObjStretcher:EventMetaData:badTimeOrder', ...
                      'tEnd must be >= tStart.');
            end
        end
    end
end


% ---- local helper ----
function z = castComplexSingle(v)
v = single(v);
if ~isreal(v)
    z = complex(real(v), imag(v));
else
    z = complex(v, single(0));
end
end