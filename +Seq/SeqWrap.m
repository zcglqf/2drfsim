classdef SeqWrap < handle
    % SeqWrap
    % Wrapper around Pulseq Sequence (*.seq)
    %
    % Responsibilities:
    % - load Pulseq sequence
    % - extract waveform/time data
    % - compute dwell time
    % - expose raster times
    % Usage:
    % - create an empty object
    %      wrap = Seq.SeqWrap(); 
    % - construct from a sequence
    %      wrap = Seg.SeqWrap(seq); % seq is a pulseq object created using
    %      seq=mr.Sequence();
    % - construct from load a '*.seq' file
    %      wrap = Seq.SeqWrap();
    %      wrap.loadFromFile(filePath)

    properties
        seq            % Pulseq sequence object

        wfRf           % cached waveforms (values OR [t,value] depending on Pulseq)
        wfGx
        wfGy
        wfGz
        wfAdc          % ADC time points (recommended) or ADC struct returned by Pulseq

        tol (1,1) double = 1e-12
    end

    methods (Access = private)
        function extWaveform(obj)
            if isempty(obj.seq)
                error("SeqWrap:NoSeq", "No sequence object loaded.");
            end

            [wfg, ~, ~, obj.wfAdc] = obj.seq.waveforms_and_times(true);

            % wfg is typically a cell: {gx, gy, gz, rf}
            obj.wfGx = wfg{1};
            obj.wfGy = wfg{2};
            obj.wfGz = wfg{3};
            obj.wfRf = wfg{4};
        end
    end

    methods
        function obj = SeqWrap(seqObj)
            if nargin >= 1
                obj.seq = seqObj;
                obj.extWaveform();
            else
                obj.seq = [];
            end            
        end

        function obj = loadFromFile(obj, fileName)
            s = mr.Sequence();
            s.read(fileName);
            obj.seq = s;

            obj.extWaveform();
        end

        function waveform = getWaveform(obj, name, refresh)
            % getWaveform(obj, name, refresh)
            % name: 'rf','gx','gy','gz','adc'
            if nargin < 3 || isempty(refresh)
                refresh = false;
            end

            if isempty(obj.seq)
                error("SeqWrap:NoSeq", "No sequence loaded.");
            end

            if isempty(obj.wfRf) || refresh
                obj.extWaveform();
            end

            name = lower(string(name));

            switch name
                case "rf"
                    waveform = obj.wfRf;
                case "gx"
                    waveform = obj.wfGx;
                case "gy"
                    waveform = obj.wfGy;
                case "gz"
                    waveform = obj.wfGz;
                case "adc"
                    waveform = obj.wfAdc;
                otherwise
                    error("SeqWrap:UnknownChannel", "Unknown channel name: %s", name);
            end
        end

        function dwell = getAdcDwell(obj)
            % getAdcDwell
            if isempty(obj.wfAdc)
                dwell = [];
                return;
            end

            dt = diff(double(obj.wfAdc(:)));
            dt = dt(dt > 0);

            if isempty(dt)
                dwell = [];
                return;
            end

            % Check if all dwell intervals are equal within tolerance
            allEqual = max(abs(dt(1:10) - dt(1))) < obj.tol;
            if allEqual
                dwell = dt(1);
            else
                % If there are gaps between readouts, dt may not be constant.
                % A robust alternative is median(dt(dt < some_threshold)).
                dwell = [];
            end
        end
    end
end