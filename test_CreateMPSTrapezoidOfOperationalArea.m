function tests = test_CreateMPSTrapezoidOfOperationalArea
% Unit tests for CreateMPSTrapezoidOfOperationalArea.m
tests = functiontests(localfunctions);
end

%% ---------- Tests ----------

function testReturnsCorrectCount(testCase)
    sys = makeTestSys();
    [g1, ~] = CreateMPSTrapezoidOfOperationalArea(10, sys, [], 0.5);
    [g2, ~] = CreateMPSTrapezoidOfOperationalArea([10 20], sys, [], 0.5);
    [g3, ~] = CreateMPSTrapezoidOfOperationalArea([10 20 30], sys, [], 0.5);

    verifyEqual(testCase, numel(g1), 1);
    verifyEqual(testCase, numel(g2), 2);
    verifyEqual(testCase, numel(g3), 3);
end

function testGradientAmplitudeRight(testCase)
    sys = mr.opts('MaxGrad', 33, 'GradUnit','mT/m', ...
               'MaxSlew', 100, 'SlewUnit','T/m/s');

    g = CreateMPSTrapezoidOfOperationalArea(1000, sys);
    g = g{1};
    verifyEqual(testCase, g.flatArea, 1000, 'AbsTol', 1e-15);
    verifyEqual(testCase, g.amplitude * g.flatTime, 1000, 'AbsTol', 1e-15);

    % After times are rounded to gradient raster time, the gradient
    % amplitude may deviate from the system max.
    verifyEqual(testCase, g.amplitude / sys.gamma * 1e3, 33, 'RelTol', 2e-2);    
    
end

function testRejectsInvalidLength(testCase)
    sys = makeTestSys();
    verifyError(testCase, @() CreateMPSTrapezoidOfOperationalArea([], sys, [], 0.5), ...
        'CreateMPSTrapezoidOfOperationalArea:invalidLength');
end

function testZeroAreaBranch(testCase)
    sys = makeTestSys();
    [grads, aux] = CreateMPSTrapezoidOfOperationalArea([0 0], sys, [], 0.5);

    verifyEqual(testCase, numel(grads), 2);
    for i = 1:2
        verifyEqual(testCase, grads{i}.amplitude, 0);
        % Note: Pulseq may rasterize 0 times to one tick; we don't require strict 0.
        verifyGreaterThanOrEqual(testCase, grads{i}.riseTime, sys.gradRasterTime);
        verifyGreaterThanOrEqual(testCase, grads{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, grads{i}.fallTime, sys.gradRasterTime);
    end

    verifyEqual(testCase, aux.maxAbsArea, 0);
    verifyEqual(testCase, aux.plateauDuration, 0);
    verifyEqual(testCase, aux.effectiveDuration, 0);
end

function testRampTimeRasterization(testCase)
    sys = makeTestSys();
    A = [100 50 0];
    [~, aux] = CreateMPSTrapezoidOfOperationalArea(A, sys, [], 0.3);

    % rampTime must be a multiple of gradRasterTime
    verifyEqual(testCase, mod(aux.rampTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);

    % rampTime is ceil((maxGrad/maxSlew)/raster)*raster
    rampExpected = sys.maxGrad / sys.maxSlew;
    rampExpected = ceil(rampExpected / sys.gradRasterTime) * sys.gradRasterTime;

    verifyEqual(testCase, aux.rampTime, rampExpected, 'AbsTol', 1e-15);
end

function testTriangleVsTrapezoidBranches(testCase)
    sys = makeTestSys();

    % Choose blankSlopeRatio so slopFactor > 0
    blankSlopeRatio = 0.0; % slopFactor = 1
    slopFactor = 1 - blankSlopeRatio^2;

    % Compute rampTime per code
    rampTime = sys.maxGrad/sys.maxSlew;
    rampTime = ceil(rampTime/sys.gradRasterTime)*sys.gradRasterTime;

    % --- Triangle case: maxAbsArea < slopFactor*maxGrad*rampTime ---
    Atri = 0.5 * (slopFactor * sys.maxGrad * rampTime);
    [~, auxTri] = CreateMPSTrapezoidOfOperationalArea([Atri], sys, [], blankSlopeRatio);
    verifyEqual(testCase, auxTri.plateauDuration, 0);

    % --- Trapezoid case: maxAbsArea large ---
    Atrap = 5.0 * (slopFactor * sys.maxGrad * rampTime);
    [~, auxTrap] = CreateMPSTrapezoidOfOperationalArea([Atrap], sys, [], blankSlopeRatio);
    verifyGreaterThanOrEqual(testCase, auxTrap.plateauDuration, 0);
end

function testAmplitudeMatchesDefinition(testCase)
    sys = makeTestSys();
    blankSlopeRatio = 0.2;
    slopFactor = 1 - blankSlopeRatio^2;

    signedOperationalAreas = [120, -60, 30];
    [grads, aux] = CreateMPSTrapezoidOfOperationalArea(signedOperationalAreas, sys, [], blankSlopeRatio);

    % The code sets: amps = areas / effectiveDuration
    for i = 1:numel(signedOperationalAreas)
        verifyEqual(testCase, grads{i}.amplitude, signedOperationalAreas(i)/aux.effectiveDuration, ...
            'RelTol', 1e-12, 'AbsTol', 1e-9);
    end

    % Also verify effectiveDuration formula in aux
    verifyEqual(testCase, aux.effectiveDuration, aux.plateauDuration + slopFactor*aux.rampTime, ...
        'AbsTol', 1e-15);
end

function testPlateauDurationFormulaMatchesCode(testCase)
    % This test recomputes plateauDuration using the SAME formula as your code,
    % including sys.gamma, and checks it matches aux.plateauDuration.
    sys = makeTestSys();

    blankSlopeRatio = 0.3;
    slopFactor = 1 - blankSlopeRatio^2;

    A = [200, 0]; % ensure trapezoid branch likely (depends on sys)
    maxAbsArea = max(abs(A));

    % rampTime per code
    rampTime = sys.maxGrad/sys.maxSlew;
    rampTime = ceil(rampTime/sys.gradRasterTime)*sys.gradRasterTime;

    % Run function
    [~, aux] = CreateMPSTrapezoidOfOperationalArea(A, sys, [], blankSlopeRatio);

    if maxAbsArea < slopFactor*sys.maxGrad*rampTime
        plateauExpected = 0;
    else
        plateauExpected = maxAbsArea / (sys.maxGrad/sys.gamma) - slopFactor*rampTime;
        plateauExpected = ceil(plateauExpected/sys.gradRasterTime)*sys.gradRasterTime;
        plateauExpected = max(plateauExpected, 0);
    end

    verifyEqual(testCase, aux.plateauDuration, plateauExpected, 'AbsTol', 1e-12);
end

%% ---------- Helpers ----------

function sys = makeTestSys()
    % Create a sys struct that works with your function, and includes sys.gamma.
    % Use units that mr.opts accepts in your install.
    sys = mr.opts('MaxGrad', 32, 'GradUnit','mT/m', ...
                  'MaxSlew', 130, 'SlewUnit','T/m/s', ...
                  'rfDeadTime', 100e-6, ...
                  'adcDeadTime', 10e-6);

    % Inject gamma field because your function uses sys.gamma.
    % Use Hz/T for proton.
    sys.gamma = 42.576e6;
end
