function grads = CreateRephasingGradientWaveformByArea(mpsarea, sys, slope)
% CalcRephasingGradientWaveformByArea
% Calculate rephasing trapezoid gradients given per-axis areas.
%
% Inputs
%   mpsarea  : vector length 1..3, area in 1/m
%   sys      : (optional) Pulseq system struct (mr.opts).
%   slope  : (optional) ramp time in s. If omitted, choose minimal ramp based on maxSlew/maxGrad.
%
% Output
%   grads : 1xN cell array of Pulseq gradient objects (N = numel(mpsarea))

    if nargin < 2 || isempty(sys)
        sys = mr.opts();
    end
    if nargin < 3 || isempty(slope)
        slope = 0;
    end

    mpsarea = double(mpsarea(:)).'; % 1xN
    n = numel(mpsarea);
    if n < 1 || n > 3
        error('mpsarea must have 1, 2, or 3 elements.');
    end

    RephAreaMax = max(abs(mpsarea));

    % Determine SlopDur and FlatDur (seconds)
    if slope == 0
        % Triangle or trapezoid with slew constraint
        if RephAreaMax * sys.maxSlew < sys.maxGrad * sys.maxGrad
            SlopDur = sqrt(RephAreaMax / sys.maxSlew);
            FlatDur = 0.0;
        else
            SlopDur = sys.maxGrad / sys.maxSlew;
            FlatDur = RephAreaMax / sys.maxGrad - sys.maxGrad / sys.maxSlew;
        end
    else
        % Enforce amplitude limit given fixed ramp
        if (RephAreaMax / slope) < sys.maxGrad
            SlopDur = slope;
            FlatDur = 0.0;
        else
            SlopDur = slope;
            FlatDur = RephAreaMax / sys.maxGrad - slope;
        end
    end

    % Rasterize times to gradRasterTime if available
    SlopDur = ceil(SlopDur / sys.gradRasterTime) * sys.gradRasterTime;
    FlatDur = ceil(FlatDur / sys.gradRasterTime) * sys.gradRasterTime;

    % Gradient strengths in SI: G = Area / (SlopDur + FlatDur)
    denom = SlopDur + FlatDur;
    if denom <= 0
        error('Invalid timing (SlopDur+FlatDur <= 0). Check inputs/limits.');
    end
    gradStr = mpsarea / denom; % T/m

    chans = {'x','y','z'};
    grads = cell(1,n);
    for i = 1:n
        grads{i} = mr.makeTrapezoid(chans{i}, sys, ...
            'amplitude', gradStr(i), ...
            'riseTime',  SlopDur, ...
            'flatTime',  FlatDur, ...
            'fallTime',  SlopDur);
    end
end

