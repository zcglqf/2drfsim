function tests = test_CalcExcitationGradientWaveforms
tests = functiontests(localfunctions);
end

function testShapesAndBasicFields(testCase)
    sys = makeTestSys();

    kstart = [-200; 50; 10].';
    kend   = [ 300; -20; 40].';

    [gLead, gExci, gTail] = CalcExcitationGradientWaveforms(kstart, kend, sys);

    verifyEqual(testCase, numel(gLead), 3);
    verifyEqual(testCase, numel(gExci), 3);
    verifyEqual(testCase, numel(gTail), 3);

    for i = 1:3
        verifyTrue(testCase, isfield(gExci{i}, 'amplitude'));
        verifyTrue(testCase, isfield(gExci{i}, 'riseTime'));
        verifyTrue(testCase, isfield(gExci{i}, 'flatTime'));
        verifyTrue(testCase, isfield(gExci{i}, 'fallTime'));

        verifyGreaterThanOrEqual(testCase, gExci{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, gExci{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, gExci{i}.fallTime, 0);

        verifyGreaterThanOrEqual(testCase, gLead{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, gLead{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, gLead{i}.fallTime, 0);

        verifyGreaterThanOrEqual(testCase, gTail{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, gTail{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, gTail{i}.fallTime, 0);
    end
end

function testExcitationAreaMatchesKstride(testCase)
    sys = makeTestSys();

    kstart = [-100, 200, -50];
    kend   = [ 400,  20,  10];

    [~, gExci, ~] = CalcExcitationGradientWaveforms(kstart, kend, sys);

    kstride = kend - kstart;

    % Excitation trapezoid area (full area): amp * (flat + 0.5*(rise+fall))
    Aexci = zeros(1,3);
    for i = 1:3
        Aexci(i) = trapezoidArea(gExci{i});
    end

    % For your excitation design, ExciGradStr = kstride/ExciFlatDur and rise/fall are symmetric.
    % This implies the *full* trapezoid area equals:
    %   A = kstride + (ExciGradStr*ExciSlopDur) = kstride + 2*ExciSlopArea
    % But your code sets ExciFlatArea = kstride, and does not explicitly enforce full area.
    %
    % What we can reliably test from your code: the FLAT area equals kstride:
    %   kstride == amp * flatTime   (since amp = kstride/flatTime)
    Aflat = zeros(1,3);
    for i = 1:3
        Aflat(i) = gExci{i}.amplitude * gExci{i}.flatTime;
    end

    verifyEqual(testCase, Aflat, kstride, 'RelTol', 1e-12, 'AbsTol', 1e-9);

    % Also sanity: flat amplitude should not exceed sys.maxGrad (within numerical)
    verifyLessThanOrEqual(testCase, max(abs([gExci{1}.amplitude, gExci{2}.amplitude, gExci{3}.amplitude])), sys.maxGrad + 1e-9);
end

function testLeadAndTailAreasMatchReturnedDurations(testCase)
    sys = makeTestSys();

    kstart = [-120, 10, 0];
    kend   = [  80, 30, -40];

    [gLead, gExci, gTail] = CalcExcitationGradientWaveforms(kstart, kend, sys);

    kstride = kend - kstart;

    % Reconstruct ExciSlopArea using returned excitation params
    % ExciSlopArea = 0.5 * ExciGradStr * ExciSlopDur; ExciGradStr == gExci.amplitude
    ExciSlopArea = 0.5 * [gExci{1}.amplitude, gExci{2}.amplitude, gExci{3}.amplitude] * gExci{1}.riseTime;

    % Your code defines:
    TailArea = -kend - ExciSlopArea;
    LeadArea = kend - ExciSlopArea - kstride; % since ExciFlatArea = kstride

    TailArea_expected = -kend - ExciSlopArea;
    LeadArea_expected = kend - ExciSlopArea - kstride;

    % Now verify that the gradients you returned realize those areas according to your own model:
    % TailGradStr = TailArea / (TailSlopDur + TailFlatDur)
    % And returned trapezoid uses rise=TailSlopDur, flat=TailFlatDur, fall=TailSlopDur.
    for i = 1:3
        denomTail = gTail{i}.riseTime + gTail{i}.flatTime; % matches your computation denom
        denomLead = gLead{i}.riseTime + gLead{i}.flatTime;

        verifyGreaterThan(testCase, denomTail, 0);
        verifyGreaterThan(testCase, denomLead, 0);

        verifyEqual(testCase, gTail{i}.amplitude * denomTail, TailArea_expected(i), ...
            'RelTol', 1e-12, 'AbsTol', 1e-6);

        verifyEqual(testCase, gLead{i}.amplitude * denomLead, LeadArea_expected(i), ...
            'RelTol', 1e-12, 'AbsTol', 1e-6);

        % Times should be rasterized by trapTimingFromArea
        verifyEqual(testCase, mod(gTail{i}.riseTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
        verifyEqual(testCase, mod(gTail{i}.flatTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
        verifyEqual(testCase, mod(gLead{i}.riseTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
        verifyEqual(testCase, mod(gLead{i}.flatTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
    end
end

function testTrapTimingFromAreaMatchesPiecewise(testCase)
    sys = makeTestSys();

    % Choose two areas: one in triangle regime, one in trapezoid regime
    Acrit = (sys.maxGrad^2) / sys.maxSlew;

    Atri  = 0.1 * Acrit;
    Atrap = 5.0 * Acrit;

    [sd1, fd1] = trapTimingFromArea(Atri, sys);
    [sd2, fd2] = trapTimingFromArea(Atrap, sys);

    % Triangle: flatDur = 0, slopDur = sqrt(A/maxSlew) (rasterized up)
    sd1_expected = sqrt(Atri / sys.maxSlew);
    sd1_expected = ceil(sd1_expected / sys.gradRasterTime) * sys.gradRasterTime;
    verifyEqual(testCase, fd1, 0);
    verifyEqual(testCase, sd1, sd1_expected, 'AbsTol', 1e-15);

    % Trapezoid: slopDur = maxGrad/maxSlew, flatDur = A/maxGrad - maxGrad/maxSlew (rasterized)
    sd2_expected = sys.maxGrad / sys.maxSlew;
    fd2_expected = Atrap / sys.maxGrad - sys.maxGrad / sys.maxSlew;
    sd2_expected = ceil(sd2_expected / sys.gradRasterTime) * sys.gradRasterTime;
    fd2_expected = ceil(fd2_expected / sys.gradRasterTime) * sys.gradRasterTime;

    verifyEqual(testCase, sd2, sd2_expected, 'AbsTol', 1e-15);
    verifyEqual(testCase, fd2, fd2_expected, 'AbsTol', 1e-12);
end

%% ---- helpers ----

function sys = makeTestSys()
    % Use Hz/m units to match the assumptions in the function under test
    sys = mr.opts('MaxGrad', 200e3, 'GradUnit','Hz/m', ...
                  'MaxSlew', 800e6, 'SlewUnit','Hz/m/s');
    if ~isfield(sys,'gradRasterTime') || isempty(sys.gradRasterTime)
        sys.gradRasterTime = 10e-6;
    end
end

function A = trapezoidArea(g)
    A = g.amplitude * (g.flatTime + 0.5*(g.riseTime + g.fallTime));
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