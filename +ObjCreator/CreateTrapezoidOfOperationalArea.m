function grad = CreateTrapezoidOfOperationalArea(signedOperationalArea, sys, slope, blankSlopeRatio, cha)
% CreateTrapezoidOfOperationalArea
%
% Inputs
%   signedOperationalArea : scalar numeric, opeartional area under the
%   trapezoid gradient
%   sys : (optional) Pulseq system struct from mr.opts(...)
%   slope : (optional) user-requested ramp duration (s)
%   blankSlopRatio : (optional) double, ratio of ramp time that is not being used for
%   encoding.
%   cha : (optional) 'x'/'y'/'z'
%
% Outputs
%   grad : Pulseq gradient object

    % ---- Defaults for optional inputs ----
    % sys
    if nargin < 2 || isempty(sys)
        sys = mr.opts();   % Pulseq default system limits
    end
    
    % slope 
    % If not provided, derive something reasonable from sys.
    % NOTE: sys.maxGrad is typically in Hz/m and sys.maxSlew in Hz/m/s in MATLAB Pulseq,
    % so slope ~= maxGrad/maxSlew gives seconds.
    if nargin < 3 || isempty(slope)
        if isfield(sys,'maxGrad') && isfield(sys,'maxSlew') && sys.maxSlew > 0
            slope = sys.maxGrad / sys.maxSlew;   % seconds
        else
            error('Cannot infer slope because sys.maxGrad/maxSlew is missing or invalid.');
        end
    end
    slope = ceil(slope / sys.gradRasterTime) * sys.gradRasterTime;
    
    % blankSlopeRatio, default is 1.0
    if nargin < 4 || isempty(blankSlopeRatio)
        blankSlopeRatio = 1.0;
    end

    % channel, default is 'x'
    if nargin < 5 || isempty(cha)
        cha = 'x';
    end

    % slopFactor = 1 - blankSlopRatio^2
    slopFactor = 1.0 - blankSlopeRatio * blankSlopeRatio;

    % return a dummy gradient if requested area is zero. but note pulseq
    % will round the rise/fall time to gradient raster time.
    if abs(signedOperationalArea) < eps
        grad = mr.makeTrapezoid(cha, sys, 'amplitude', 0, 'riseTime', 0, 'flatTime', 0, 'fallTime', 0);
        return;
    end

    % minimumRampDurationMaxAmplitude = (maxStrength/(maxSlewRate*grDwellTime))*grDwellTime;
    % simplifies to maxStrength/maxSlewRate, but keep the dwellTime form for clarity
    minimumRampDurationMaxAmplitude = sys.maxGrad / sys.maxSlew;
    minimumRampDurationMaxAmplitude = ceil(minimumRampDurationMaxAmplitude/ sys.gradRasterTime) * sys.gradRasterTime;

    % if user input ramp is reasonably large enough, then use it.
    if minimumRampDurationMaxAmplitude > slope
        slope = minimumRampDurationMaxAmplitude;
    end

    % Determine plateauDuration based on triangle vs trapezoid condition
    if abs(signedOperationalArea) < slopFactor * sys.maxGrad * slope
        plateauDuration = 0.0; % triangle
    else
        plateauDuration = abs(signedOperationalArea) / sys.maxGrad - slopFactor * slope;
        plateauDuration = ceil(plateauDuration/sys.gradRasterTime) * sys.gradRasterTime;
        if plateauDuration < 0
            plateauDuration = 0.0; % numeric guard
        end
    end

    % effectiveDuration = plateauDuration + slopFactor * minimumRampDurationMaxAmplitude
    effectiveDuration = plateauDuration + slopFactor * slope;

    amplitude = signedOperationalArea / effectiveDuration;

    % Return trapezoid
    grad = mr.makeTrapezoid(cha, sys, 'amplitude', amplitude, 'riseTime', slope, 'flatTime', plateauDuration, 'fallTime', slope);
end
