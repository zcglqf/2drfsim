function segs = samplesToSegs(samplesIn, sys, adcDwell, minDuration, tol)
%ObjStretcher.samplesToSegs Convert EventSampleData[] to Seg[].
%
% This follows the Python logic:
%   - sort by Time
%   - maintain Depth(type) and Value(type) stacks
%   - for each consecutive time interval:
%       if duration < minDuration: just update record and continue
%       else: emit one Seg with current accumulated RF/Grad values and duration
%
% Inputs:
%   samplesIn   : ObjStretcher.EventSampleData[]
%   sys         : pulseq mr.opts struct (used only if you choose to validate)
%   adcDwell    : dwell time in seconds (optional; passed to validateAgainstSystem if used)
%   minDuration : default 1e-7
%   tol         : default 1e-12
%
% Output:
%   segs        : ObjStretcher.Seg[]
if nargin < 4 || isempty(minDuration)
    minDuration = 1e-7;
end
if nargin < 5 || isempty(tol)
    tol = 1e-12;
end

if isempty(samplesIn)
    segs = ObjStretcher.Seg.empty(0,1);
    return;
end

% ---- Sort by Time, tie-breaker: process starts before ends at same time ----
t = arrayfun(@(s) s.t, samplesIn);
k = arrayfun(@(s) s.startOrEnd, samplesIn);   % true=start, false=end
% sortrows: primary Time asc, secondary StartOrEnd desc (start first)
[~, order] = sortrows([t(:), -double(k(:))], [1 2]);
samples = samplesIn(order);

% ---- Setup Depth/Value stacks over EventType enum values ----
allTypes = enumeration('ObjStretcher.EventType');
nTypes = numel(allTypes);

Depth = zeros(nTypes, 1, 'int16');
Value = complex(zeros(nTypes, 1, 'single'));

% Map enum -> 1-based index (assumes enum underlying values are 0..)
typeToIdx = @(tp) double(uint8(tp)) + 1;

% Helper to apply one sample record
    function updateRecord(sp)
        idx = typeToIdx(sp.Type);
        v = single(sp.Value);

        if sp.StartOrEnd
            Depth(idx) = Depth(idx) + 1;
            Value(idx) = Value(idx) + v;
        else
            Depth(idx) = Depth(idx) - 1;
            Value(idx) = Value(idx) - v;
        end
    end

% Indices for commonly used types
iRF  = typeToIdx(ObjStretcher.EventType.RFAMPLITUDE);
iGM  = typeToIdx(ObjStretcher.EventType.GRADIENTM);
iGP  = typeToIdx(ObjStretcher.EventType.GRADIENTP);
iGS  = typeToIdx(ObjStretcher.EventType.GRADIENTS);
iACQ = typeToIdx(ObjStretcher.EventType.ACQ);

% Initialize
updateRecord(samples(1));

% Preallocate worst-case number of segments
segsTmp(numel(samples)-1,1) = ObjStretcher.Seg();
segCount = 0;

for i = 2:numel(samples)
    dt = samples(i).Time - samples(i-1).Time;

    if dt < minDuration
        updateRecord(samples(i));
        continue;
    end

    rfComplex = Value(iRF);
    gx = real(Value(iGM));
    gy = real(Value(iGP));
    gz = real(Value(iGS));

    % Acquisition active if ACQ depth > 0
    doAcq = Depth(iACQ) > 0;

    % adc flag: 0 or 1 only here (start/end flags handled in later expansion step)
    adcFlag = uint8(0);
    if doAcq
        adcFlag = uint8(1);
    end

    segCount = segCount + 1;

    segsTmp(segCount) = ObjStretcher.Seg( ...
        'tStart', samples(i-1).Time, ...
        'dur', dt, ...
        'rf', rfComplex, ...
        'gx', gx, ...
        'gy', gy, ...
        'gz', gz, ...
        'adc', adcFlag);

    % Optional validation (enable if desired)
    % segsTmp(segCount).validateAgainstSystem(sys, adcDwell, tol);

    updateRecord(samples(i));
end

segs = segsTmp(1:segCount);

end