function [gradsMPS, aux] = CreateMPSTrapezoidOfOperationalArea(signedOperationalAreas, sys, slope, blankSlopeRatio)
% CreateMPSTrapezoidOfOperationalArea
%
% Inputs
%   signedOperationalAreas : 1x3 numeric, operational areas for M/P/S
%   sys                   : (optional) Pulseq sys from mr.opts()
%   blankSlopeRatio        : (optional) slope blank ratio (0..1)
%
% Outputs
%   gradsMPS : 1x3 cell array {gM, gP, gS} (Pulseq gradient events)
%   aux      : struct of intermediate values (for debugging/tests)

    if nargin < 2 || isempty(sys)
        sys = mr.opts();
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

    if nargin < 4 || isempty(blankSlopeRatio)
        blankSlopeRatio = 1.0; % typical default; set to 1.0 if you want "flat-only"
    end
    blankSlopeRatio = min(max(blankSlopeRatio, 0), 1);

    signedOperationalAreas = double(signedOperationalAreas(:)).'; % force 1x3
    n = numel(signedOperationalAreas);
    if n < 1 || n > 3
        error('CreateMPSTrapezoidOfOperationalArea:invalidLength', 'signedOperationalAreas must have 1, 2, or 3 elements.');
    end
    chans = {'x','y','z'};
    gradsMPS = cell(1, n);

    % maxAbsArea = max over axes of abs(area)
    maxAbsArea = max(abs(signedOperationalAreas));

    % slopFactor = 1 - blankSlopeRatio^2
    slopFactor = 1.0 - blankSlopeRatio * blankSlopeRatio;

    % CloseEnough(maxAbsArea, 0, minPossibleAbsoluteArea)
    % Pulseq doesn't expose minPossibleAbsoluteArea directly; use eps-based threshold.
    if maxAbsArea < eps
        for i = 1:n
            gradsMPS{i} = mr.makeTrapezoid(chans{i}, sys, 'amplitude', 0, 'riseTime', 0, 'flatTime', 0, 'fallTime', 0);
        end
        aux = packAux(maxAbsArea, slopFactor, 0, 0, 0, signedOperationalAreas);
        return;
    end

    % minimumRampDurationMaxAmplitude = maxStrength/maxSlewRate (rasterized)
    % (the C++ expression with grDwellTime simplifies to this)
    minimumRampDurationMaxAmplitude = sys.maxGrad / sys.maxSlew; % seconds
    minimumRampDurationMaxAmplitude = ceil(minimumRampDurationMaxAmplitude / sys.gradRasterTime) * sys.gradRasterTime;

    % if user input ramp is reasonably large enough, then use it.
    if minimumRampDurationMaxAmplitude > slope
        slope = minimumRampDurationMaxAmplitude;
    end
    
    % Determine plateauDuration based on triangle vs trapezoid condition
    if maxAbsArea < slopFactor * sys.maxGrad * slope
        plateauDuration = 0.0; % triangle
    else
        plateauDuration = maxAbsArea / sys.maxGrad - slopFactor * slope;
        plateauDuration = ceil(plateauDuration / sys.gradRasterTime) * sys.gradRasterTime;
        plateauDuration = max(plateauDuration, 0);
    end

    effectiveDuration = plateauDuration + slopFactor * slope;
    if effectiveDuration <= 0
        error('effectiveDuration <= 0. Check blankSlopeRatio and system limits.');
    end

    amps = signedOperationalAreas(:).' / effectiveDuration;  % 1xN

    for i = 1:n
        gradsMPS{i} = mr.makeTrapezoid(chans{i}, sys, ...
            'amplitude', amps(i), ...
            'riseTime',  slope, ...
            'flatTime',  plateauDuration, ...
            'fallTime',  slope);
    end

    aux = packAux(maxAbsArea, slopFactor, plateauDuration, minimumRampDurationMaxAmplitude, effectiveDuration, signedOperationalAreas);
end

function aux = packAux(maxAbsArea, slopFactor, plateauDuration, rampTime, effectiveDuration, signedOperationalAreas)
    aux = struct();
    aux.maxAbsArea = maxAbsArea;
    aux.slopFactor = slopFactor;
    aux.plateauDuration = plateauDuration;
    aux.rampTime = rampTime;
    aux.effectiveDuration = effectiveDuration;
    aux.signedOperationalAreas = signedOperationalAreas;
end
