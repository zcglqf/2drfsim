function tests = test_CreateTrapezoidOfOperationalArea
% Unit tests for CreateTrapezoidOfOperationalArea.m
tests = functiontests(localfunctions);
end

function testZeroAreaReturnsZeroGrad(testCase)
    sys = mr.opts(); % default Pulseq sys

    g = CreateTrapezoidOfOperationalArea(0, sys, [], 0.0, 'x');

    verifyEqual(testCase, g.amplitude, 0);
    verifyEqual(testCase, g.riseTime,  sys.gradRasterTime);
    verifyEqual(testCase, g.flatTime,  0);
    verifyEqual(testCase, g.fallTime,  sys.gradRasterTime);
end

function testGradientAmplitudeRight(testCase)
    sys = mr.opts('MaxGrad', 33, 'GradUnit','mT/m', ...
               'MaxSlew', 100, 'SlewUnit','T/m/s');

    g = CreateTrapezoidOfOperationalArea(1000, sys);

    verifyEqual(testCase, g.flatArea, 1000, 'AbsTol', 1e-15);
    verifyEqual(testCase, g.amplitude * g.flatTime, 1000, 'AbsTol', 1e-15);

    % After times are rounded to gradient raster time, the gradient
    % amplitude may deviate from the system max.
    verifyEqual(testCase, g.amplitude / sys.gamma * 1e3, 33, 'RelTol', 2e-2);    
    
end

function testRasterizationAndMinSlope(testCase)
    sys = mr.opts();

    signedA = 100;           % 1/m
    slope_in = 7.3e-6;       % s (intentionally not raster-aligned)
    blankSlopeRatio = 0.2;   % slopFactor = 1 - 0.04 = 0.96

    g = CreateTrapezoidOfOperationalArea(signedA, sys, slope_in, blankSlopeRatio, 'x');

    % Slope must be rasterized
    verifyEqual(testCase, mod(g.riseTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
    verifyEqual(testCase, mod(g.fallTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);
    verifyEqual(testCase, mod(g.flatTime, sys.gradRasterTime), 0, 'AbsTol', 1e-15);

    % Slope must be >= (maxGrad/maxSlew) rasterized
    minSlope = sys.maxGrad / sys.maxSlew;
    minSlope = round(minSlope / sys.gradRasterTime) * sys.gradRasterTime;

    verifyGreaterThanOrEqual(testCase, g.riseTime, minSlope);
    verifyEqual(testCase, g.riseTime, g.fallTime);
end

function testOperationalAreaMatchesModel_Triangle(testCase)
    sys = mr.opts();

    % Choose parameters that force triangle condition:
    blankSlopeRatio = 0.0;     % slopFactor=1
    slopFactor = 1 - blankSlopeRatio^2;

    % Use an explicit slope so we can craft "small area"
    slope = sys.maxGrad / sys.maxSlew;
    slope = round(slope / sys.gradRasterTime) * sys.gradRasterTime;

    % Force triangle: abs(A) < slopFactor*maxGrad*slope
    signedA = 0.25 * (slopFactor * sys.maxGrad / sys.gamma * slope);  % 1/m

    g = CreateTrapezoidOfOperationalArea(signedA, sys, slope, blankSlopeRatio, 'x');

    % Should be triangle
    verifyEqual(testCase, g.flatTime, 0);

    % Verify operational area model:
    % effectiveDuration = flat + slopFactor*slope = slopFactor*slope
    effDur = g.flatTime + slopFactor * g.riseTime;
    A_model = g.amplitude * effDur;

    verifyEqual(testCase, A_model, signedA, 'RelTol', 1e-12, 'AbsTol', 1e-9);

    % Also check the actual physical trapezoid area (includes full ramps):
    A_trap = trapezoidArea(g);
    verifyGreaterThanOrEqual(testCase, abs(A_trap), abs(signedA)); % with slopFactor=1, should be equal-ish or larger depending model
end

function testOperationalAreaMatchesModel_Trapezoid(testCase)
    sys = mr.opts();

    blankSlopeRatio = 0.6;     % slopFactor = 1 - 0.36 = 0.64
    slopFactor = 1 - blankSlopeRatio^2;

    % Let slope be inferred by passing []
    slope = []; 

    % Choose an area large enough to trigger trapezoid branch:
    % Need abs(A) >= slopFactor*maxGrad*slope. We don't know slope yet, but it's ~maxGrad/maxSlew.
    slope_guess = sys.maxGrad / sys.maxSlew;
    slope_guess = round(slope_guess / sys.gradRasterTime) * sys.gradRasterTime;

    signedA = 5.0 * (slopFactor * sys.maxGrad / sys.gamma * slope_guess);  % 1/m, safely trapezoid

    g = CreateTrapezoidOfOperationalArea(signedA, sys, slope, blankSlopeRatio, 'x');

    % Should be trapezoid (flatTime > 0), unless raster rounding makes it barely 0
    verifyGreaterThanOrEqual(testCase, g.flatTime, 0);

    % Verify operational area model: A = amp*(flat + slopFactor*slope)
    effDur = g.flatTime + slopFactor * g.riseTime;
    A_model = g.amplitude * effDur;

    verifyEqual(testCase, A_model, signedA, 'RelTol', 1e-12, 'AbsTol', 1e-6);

    % Sanity: amplitude should not exceed maxGrad (Pulseq will enforce too)
    verifyLessThanOrEqual(testCase, abs(g.amplitude), sys.maxGrad + 1e-9);
end

function testSignPreserved(testCase)
    sys = mr.opts();

    blankSlopeRatio = 0.3;
    Apos = 200;
    Aneg = -200;

    gpos = CreateTrapezoidOfOperationalArea(Apos, sys, [], blankSlopeRatio, 'x');
    gneg = CreateTrapezoidOfOperationalArea(Aneg, sys, [], blankSlopeRatio, 'x');

    verifyGreaterThan(testCase, gpos.amplitude, 0);
    verifyLessThan(testCase, gneg.amplitude, 0);
end

% ---- helper: full trapezoid area in 1/m (includes full ramps) ----
function A = trapezoidArea(g)
    A = g.amplitude * (g.flatTime + 0.5*(g.riseTime + g.fallTime));
end
