function grads = CalcRephasingGradientWaveformByArea(mpsarea, sys, slope)
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
        slope = sys.maxGrad / sys.maxSlew;
    end
    slope = ceil(slope / sys.gradRasterTime) * sys.gradRasterTime;

    mpsarea = double(mpsarea(:)).'; % 1xN
    n = numel(mpsarea);
    if n < 1 || n > 3
        error('mpsarea must have 1, 2, or 3 elements.');
    end

    RephAreaMax = max(abs(mpsarea));

    % Determine SlopDur and FlatDur (seconds)
    if slope == 0
        % Triangle or trapezoid with slew constraint
        if RephAreaMax * maxslew_Tms < maxgrad_Tm * maxgrad_Tm
            SlopDur = sqrt(RephAreaMax / maxslew_Tms);
            FlatDur = 0.0;
        else
            SlopDur = maxgrad_Tm / maxslew_Tms;
            FlatDur = RephAreaMax / maxgrad_Tm - maxgrad_Tm / maxslew_Tms;
        end
    else
        ramp_s = double(slope) * 1e-3;
        % Enforce amplitude limit given fixed ramp
        if (RephAreaMax / ramp_s) < maxgrad_Tm
            SlopDur = ramp_s;
            FlatDur = 0.0;
        else
            SlopDur = ramp_s;
            FlatDur = RephAreaMax / maxgrad_Tm - ramp_s;
        end
    end

    % Rasterize times to gradRasterTime if available
    if isfield(sys,'gradRasterTime') && ~isempty(sys.gradRasterTime) && sys.gradRasterTime > 0
        SlopDur = ceil(SlopDur / sys.gradRasterTime) * sys.gradRasterTime;
        FlatDur = ceil(FlatDur / sys.gradRasterTime) * sys.gradRasterTime;
    end

    % Gradient strengths in SI: G = Area / (SlopDur + FlatDur)
    denom = SlopDur + FlatDur;
    if denom <= 0
        error('Invalid timing (SlopDur+FlatDur <= 0). Check inputs/limits.');
    end
    gradStr_Tm = areas_Tm_s / denom; % T/m

    % Convert amplitudes into sys.gradUnit for mr.makeTrapezoid
    amps = physicalToSysAmplitude(gradStr_Tm, sys);

    chans = {'x','y','z'};
    grads = cell(1,n);
    for i = 1:n
        grads{i} = mr.makeTrapezoid(chans{i}, sys, ...
            'amplitude', amps(i), ...
            'riseTime',  SlopDur, ...
            'flatTime',  FlatDur, ...
            'fallTime',  SlopDur);
    end
end

