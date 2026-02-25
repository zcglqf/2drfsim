function [LeadRephGradsMPS, ExcitationGradsMPS, TailRephGradsMPS] = CreateExcitationGradientWaveforms( ...
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
    end

    kstart = double(kstart(:)).';   % 1x3
    kend   = double(kend(:)).';     % 1x3

    if numel(kstart) ~= numel(kend)
        error('kstart and kend must have the same length.');
    end
    nel = numel(kstart);

    % kstride in m^-1
    kstride = kend - kstart;

    % --- Excitation gradient durations/areas ---
    % Python:
    % ExciFlatDur = 1e3 * max(|kstride|) / (sys.maxGrad * gamma)
    % Interpretation: convert mm^-1 -> m^-1 by 1e3. Dur in seconds.
    ExciFlatDur = max(abs(kstride)) / sys.maxGrad;  % s
    ExciFlatDur = ceil(ExciFlatDur/sys.gradRasterTime) * sys.gradRasterTime;

    % ExciFlatArea = 1e3 * kstride / gamma  (units: T/m * s)
    ExciFlatArea = kstride;

    % ExciGradStr = ExciFlatArea / ExciFlatDur  (T/m)
    ExciGradStr = ExciFlatArea / ExciFlatDur;

    % ExciSlopDur = max(|ExciGradStr| / sys.maxSlew)  (s)
    ExciSlopDur = max(abs(ExciGradStr) / sys.maxSlew);
    ExciSlopDur = ceil(ExciSlopDur / sys.gradRasterTime) * sys.gradRasterTime;

    % ExciSlopArea = 0.5 * ExciGradStr * ExciSlopDur  (T/m*s)
    ExciSlopArea = 0.5 * ExciGradStr * ExciSlopDur;

    % --- Tail rephase gradient (after excitation): cancel kend and excitation ramps ---
    TailAreaMax = max(abs(-kend - ExciSlopArea));
    TailArea    = -kend - ExciSlopArea;
    [TailSlopDur, TailFlatDur] = trapTimingFromArea(TailAreaMax, sys);
    TailGradStr = TailArea / (TailSlopDur + TailFlatDur);

    % --- Lead rephase gradient (before excitation): prewind to start kstart etc. ---
    LeadAreaMax = max(abs(kend - ExciSlopArea - ExciFlatArea));
    LeadArea    = kend - ExciSlopArea - ExciFlatArea;
    [LeadSlopDur, LeadFlatDur] = trapTimingFromArea(LeadAreaMax, sys);
    LeadGradStr = LeadArea / (LeadSlopDur + LeadFlatDur);

    % --- Build Pulseq gradient objects for x/y/z (map M/P/S to x/y/z) ---
    chans = {'x','y','z'};  % choose your mapping M->x, P->y, S->z

    LeadRephGradsMPS      = cell(1,nel);
    ExcitationGradsMPS    = cell(1,nel);
    TailRephGradsMPS      = cell(1,nel);

    for ori = 1:numel(kstart)
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

end


function [slopDur, flatDur] = trapTimingFromArea(areaMax, sys)
    % areaMax in 1/m, sys.maxGrad in Hz/m, sys.maxSlew in Hz/m/s
    if areaMax < eps
        slopDur = 0; flatDur = 0; return;
    end

    if areaMax * sys.maxSlew < sys.maxGrad^2
        slopDur = sqrt(areaMax / sys.maxSlew);
        flatDur = 0;
    else
        slopDur = sys.maxGrad / sys.maxSlew;
        flatDur = areaMax / sys.maxGrad - sys.maxGrad / sys.maxSlew;
    end

    slopDur = ceil(slopDur / sys.gradRasterTime) * sys.gradRasterTime;
    flatDur = ceil(flatDur / sys.gradRasterTime) * sys.gradRasterTime;
end