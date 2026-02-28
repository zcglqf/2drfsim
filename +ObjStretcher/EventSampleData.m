classdef EventSampleData
    % Represents a time-point sample extracted from EventMetaData.
    %
    % Fields:
    %   t        : time point in seconds (double)
    %   startOrEnd  : logical (true = start, false = end)
    %   type        : ObjStretcher.EventType
    %   v       : waveform value at that time (complex single)
    %
    % Units:
    %   t  -> seconds
    %   v -> Hz or Hz/m depending on Type

    properties
        t       (1,1) double = 0
        startOrEnd (1,1) logical = false
        type       (1,1) ObjStretcher.EventType = ObjStretcher.EventType.RF
        v      (1,1) single = complex(single(0))
    end

    methods
        function obj = EventSampleData(varargin)
            % Name/value constructor

            if nargin == 0
                return;
            end

            if mod(nargin,2) ~= 0
                error('ObjStretcher:EventSampleData:badArgs', ...
                      'Name/value inputs must come in pairs.');
            end

            for i = 1:2:nargin
                name = lower(string(varargin{i}));
                val  = varargin{i+1};

                switch name
                    case "t"
                        obj.t = double(val);

                    case "startorend"
                        obj.startOrEnd = logical(val);

                    case "type"
                        obj.type = val;

                    case "value"
                        obj.v = castComplexSingle(val);

                    otherwise
                        error('ObjStretcher:EventSampleData:unknownField', ...
                              'Unknown field: %s', name);
                end
            end
            validate(obj);
        end

        function validate(obj)
            if ~isfinite(obj.t)
                error('ObjStretcher:EventSampleData:badTime', ...
                      'Time must be finite.');
            end
        end

        function tf = isStart(obj)
            tf = obj.startOrEnd;
        end

        function tf = isEnd(obj)
            tf = ~obj.startOrEnd;
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