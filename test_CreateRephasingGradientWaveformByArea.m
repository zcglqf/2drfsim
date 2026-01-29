function tests = test_CreateRephasingGradientWaveformByArea
% Unit tests for CreateRephasingGradientWaveformByArea.m
tests = functiontests(localfunctions);
end

function testReturnsCorrectCount(testCase)
    sys = makeTestSys();

    g1 = CreateRephasingGradientWaveformByArea(10, sys, 0);
    g2 = CreateRephasingGradientWaveformByArea([10 20], sys, 0);
    g3 = CreateRephasingGradientWaveformByArea([10 20 30], sys, 0);

    verifyEqual(testCase, numel(g1), 1);
    verifyEqual(testCase, numel(g2), 2);
    verifyEqual(testCase, numel(g3), 3);
end

function testRejectsInvalidLength(testCase)
    sys = makeTestSys();

    didThrow = false;
    try
        CreateRephasingGradientWaveformByArea([], sys, 0);
    catch
        didThrow = true;
    end
    verifyTrue(testCase, didThrow);
end

function testRasterization(testCase)
    sys = makeTestSys();

    areas = [100 50 25];
    grads = CreateRephasingGradientWaveformByArea(areas, sys, 0);

    % Check rise/flat are multiples of gradRasterTime
    for i = 1:numel(grads)
        verifyEqual(testCase, mod(grads{i}.riseTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
        verifyEqual(testCase, mod(grads{i}.flatTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
        verifyEqual(testCase, grads{i}.fallTime, grads{i}.riseTime);
    end
end

function testTriangleBranchSlope0(testCase)
    sys = makeTestSys();

    % Choose a small area so RephAreaMax * maxSlew < maxGrad^2
    % i.e., RephAreaMax < maxGrad^2 / maxSlew
    Acrit = (sys.maxGrad^2) / sys.maxSlew;
    areas = 0.1 * Acrit; % definitely triangle

    grads = CreateRephasingGradientWaveformByArea(areas, sys, 0);

    % In triangle branch, FlatDur should be 0 (after rasterization it should stay 0)
    verifyEqual(testCase, grads{1}.flatTime, 0);

    % Area consistency with your model: denom = SlopDur + FlatDur = SlopDur
    % amplitude = area/denom => amplitude*denom == area
    denom = grads{1}.riseTime + grads{1}.flatTime;
    verifyGreaterThan(testCase, denom, 0);
    verifyEqual(testCase, grads{1}.amplitude * denom, areas, 'RelTol', 1e-12, 'AbsTol', 1e-9);
end

function testFixedSlopeNoFlat(testCase)
    sys = makeTestSys();

    slope = 200e-6; % 0.2 ms
    % Pick area such that area/slope < maxGrad so FlatDur=0
    areas = 0.5 * sys.maxGrad * slope;

    grads = CreateRephasingGradientWaveformByArea(areas, sys, slope);

    verifyEqual(testCase, grads{1}.riseTime, ceil(slope/sys.gradRasterTime)*sys.gradRasterTime);
    verifyEqual(testCase, grads{1}.flatTime, 0);

    denom = grads{1}.riseTime + grads{1}.flatTime;
    verifyEqual(testCase, grads{1}.amplitude * denom, areas, 'RelTol', 1e-12, 'AbsTol', 1e-9);
end

function testFixedSlopeWithFlat(testCase)
    sys = makeTestSys();

    slope = 200e-6; % 0.2 ms
    % Pick area such that area/slope >= maxGrad so FlatDur>0
    areas = 2.0 * sys.maxGrad * slope; % forces FlatDur positive

    grads = CreateRephasingGradientWaveformByArea(areas, sys, slope);

    verifyGreaterThan(testCase, grads{1}.flatTime, 0);

    denom = grads{1}.riseTime + grads{1}.flatTime;
    verifyEqual(testCase, grads{1}.amplitude * denom, areas, 'RelTol', 1e-12, 'AbsTol', 1e-6);
end

function testLargeAreaSlope0HitsTrapezoidBranch(testCase)
    % This will FAIL until you fix the bug sys.sys.maxGrad -> sys.maxGrad
    sys = makeTestSys();

    % Large area so RephAreaMax*maxSlew >= maxGrad^2
    Acrit = (sys.maxGrad^2) / sys.maxSlew;
    areas = 5.0 * Acrit;

    didThrow = false;
    try
        grads = CreateRephasingGradientWaveformByArea(areas, sys, 0);
    catch ME
        didThrow = true;
        % Helpful message if your current typo triggers:
        disp(ME.message);
    end

    verifyFalse(testCase, didThrow, ...
        'Function threw an error in slope==0 trapezoid branch. Fix sys.sys.maxGrad typo.');

    % After you fix the typo, you can also assert flatTime > 0
    grads = CreateRephasingGradientWaveformByArea(areas, sys, 0);
    verifyGreaterThan(testCase, grads{1}.flatTime, 0);
end

%% ---- Helper ----
function sys = makeTestSys()
    % Build sys in units accepted by your Pulseq install.
    sys = mr.opts('MaxGrad', 200, 'GradUnit','Hz/m', ...
                  'MaxSlew', 800, 'SlewUnit','Hz/m/s');
    % Ensure gradRasterTime exists; mr.opts() usually provides it.
    if ~isfield(sys,'gradRasterTime') || isempty(sys.gradRasterTime)
        sys.gradRasterTime = 10e-6;
    end
end
