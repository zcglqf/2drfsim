function [LeadRephGradsMPS, ExcitationGradsMPS, TailRephGradsMPS] = CalcExcitationGradientWaveforms( ...
    kstart, kend, sys)
% CalcExcitationGradientWaveforms
% Returns Pulseq gradient objects (mr.makeTrapezoid) for:
%   - Excitation selection gradient
%   - Tail rephase gradient
%   - Lead rephase gradient
%
% Inputs:
%   kstart, kend : 1x3 vectors in m^-1 (M/P/S)
%   sys          : Pulseq system struct (optional). If empty, will create.
%
% Outputs:
%   LeadRephGradsMPS, ExcitationGradsMPS, TailRephGradsMPS : 1x3 cell arrays of grad objects

    if nargin < 3 || isempty(sys)
        sys = mr.opts();
        sys.slopeRatio = 0;
    end

    kstart = double(kstart(:)).';   % 1x3
    kend   = double(kend(:)).';     % 1x3

    if numel(kstart) ~= numel(kend)
        error('kstart and kend must have the same length.');
    end

    % kstride in m^-1
    kstride = kend - kstart;

    % --- Excitation gradient durations/areas ---
    % Python:
    % ExciFlatDur = 1e3 * max(|kstride|) / (maxgrad * gamma)
    % Interpretation: convert mm^-1 -> m^-1 by 1e3. Dur in seconds.
    ExciFlatDur = max(abs(kstride)) / (sys.maxGrad * gamma);  % s

    % ExciFlatArea = 1e3 * kstride / gamma  (units: T/m * s)
    ExciFlatArea = (1e3 * kstride) / gamma;

    % ExciGradStr = ExciFlatArea / ExciFlatDur  (T/m)
    ExciGradStr = ExciFlatArea / ExciFlatDur;

    % ExciSlopDur = max(|ExciGradStr| / maxslew)  (s)
    ExciSlopDur = max(abs(ExciGradStr) ./ maxslew);

    % ExciSlopArea = 0.5 * ExciGradStr * ExciSlopDur  (T/m*s)
    ExciSlopArea = 0.5 * ExciGradStr * ExciSlopDur;

    % Reference time: your code uses axis index 1 (Python) which is MATLAB 2
    % ref = (0 - kstart[1])/(kend[1]-kstart[1]) * ExciFlatDur + ExciSlopDur
    if abs(kend(2) - kstart(2)) < eps
        warning('kend(2) == kstart(2); reference time set to ExciSlopDur.');
        ref = ExciSlopDur;
    else
        ref = (0.0 - kstart(2)) / (kend(2) - kstart(2)) * ExciFlatDur + ExciSlopDur;
    end

    % --- Tail rephase gradient (after excitation): cancel kend and excitation ramps ---
    TailAreaMax = max(abs(-(1e3 * kend) / gamma - ExciSlopArea));

    if TailAreaMax * maxslew < maxgrad * maxgrad
        TailSlopDur = sqrt(TailAreaMax / maxslew);
        TailFlatDur = 0.0;
    else
        TailSlopDur = maxgrad / maxslew;
        TailFlatDur = TailAreaMax / maxgrad - maxgrad / maxslew;
    end

    TailArea    = -(1e3 * kend) / gamma - ExciSlopArea;
    TailGradStr = TailArea / (TailSlopDur + TailFlatDur);

    % --- Lead rephase gradient (before excitation): prewind to start kstart etc. ---
    LeadAreaMax = max(abs((1e3 * kend) / gamma - ExciSlopArea - ExciFlatArea));

    if LeadAreaMax * maxslew < maxgrad * maxgrad
        LeadSlopDur = sqrt(LeadAreaMax / maxslew);
        LeadFlatDur = 0.0;
    else
        LeadSlopDur = maxgrad / maxslew;
        LeadFlatDur = LeadAreaMax / maxgrad - maxgrad / maxslew;
    end

    LeadArea    = (1e3 * kend) / gamma - ExciSlopArea - ExciFlatArea;
    LeadGradStr = LeadArea / (LeadSlopDur + LeadFlatDur);

    % --- Build Pulseq gradient objects for x/y/z (map M/P/S to x/y/z) ---
    chans = {'x','y','z'};  % choose your mapping M->x, P->y, S->z

    LeadRephGradsMPS      = cell(1,3);
    ExcitationGradsMPS    = cell(1,3);
    TailRephGradsMPS      = cell(1,3);

    for ori = 1:3
        ch = chans{ori};

        % Excitation trapezoid: rise=ExciSlopDur, flat=ExciFlatDur, fall=ExciSlopDur
        ExcitationGradsMPS{ori} = mr.makeTrapezoid(ch, sys, ...
            'amplitude', ExciGradStr(ori), ...
            'riseTime',  ExciSlopDur, ...
            'flatTime',  ExciFlatDur, ...
            'fallTime',  ExciSlopDur);

        % Tail rephase trapezoid
        TailRephGradsMPS{ori} = mr.makeTrapezoid(ch, sys, ...
            'amplitude', TailGradStr(ori), ...
            'riseTime',  TailSlopDur, ...
            'flatTime',  TailFlatDur, ...
            'fallTime',  TailSlopDur);

        % Lead rephase trapezoid
        LeadRephGradsMPS{ori} = mr.makeTrapezoid(ch, sys, ...
            'amplitude', LeadGradStr(ori), ...
            'riseTime',  LeadSlopDur, ...
            'flatTime',  LeadFlatDur, ...
            'fallTime',  LeadSlopDur);
    end

    % NOTE:
    % Pulseq grad objects do not carry an internal "reference time" like your
    % TrapezoidStaticAmplitude class does. To emulate your ref/ref1/ref2 timing,
    % you typically place these grads into different blocks (mr.makeBlockPulse/seq.addBlock)
    % with appropriate delays, or you pad with mr.makeDelay().
end
