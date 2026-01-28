function tests = test_CalcExcitationGradientWaveforms
% Unit tests for CalcExcitationGradientWaveforms.m (Pulseq-native units)
% Expects Pulseq gradients with amplitude in Hz/m.

tests = functiontests(localfunctions);
end

function testAreasMatchKStride(testCase)
    gamma = 42.576e6; % Hz/T (proton), used only for converting physical limits

    % ----- Define physical limits (for readability) -----
    maxgrad_Tm  = 40e-3;   % T/m
    maxslew_Tms = 150;     % T/m/s

    % ----- Convert to Pulseq MATLAB units -----
    maxgrad_Hzm  = gamma * maxgrad_Tm;    % Hz/m
    maxslew_Hzms = gamma * maxslew_Tms;   % Hz/m/s

    sys = mr.opts('MaxGrad', maxgrad_Hzm,  'GradUnit','T/m', ...
                  'MaxSlew', maxslew_Hzms, 'SlewUnit','T/m/s');

    % Non-degenerate case (avoid kend(2)==kstart(2))
    kstart = [-2.0, -3.0,  0.5];   % m^-1
    kend   = [ 4.0,  1.0, -1.0];   % m^-1

    % IMPORTANT:
    % This test assumes CalcExcitationGradientWaveforms expects maxgrad/maxslew in T/m and T/m/s
    % (physical units) AND internally converts to Hz/m for Pulseq gradients.
    [glead, gexci, gtail, ref] = CalcExcitationGradientWaveforms( ...
        kstart, kend, maxgrad_Tm, maxslew_Tms, gamma, sys);

    % Basic checks
    verifyEqual(testCase, numel(gexci), 3);
    verifyEqual(testCase, numel(gtail), 3);
    verifyEqual(testCase, numel(glead), 3);
    verifyGreaterThanOrEqual(testCase, ref, 0);

    % Trapezoid area in Pulseq units:
    % area = amplitude(Hz/m) * [flat + 0.5*(rise+fall)](s)  => 1/m
    trapArea = @(g) g.amplitude * (g.flatTime + 0.5*(g.riseTime + g.fallTime));

    Aexci = zeros(1,3);
    Atail = zeros(1,3);
    Alead = zeros(1,3);

    for i = 1:3
        % Sanity: durations are nonnegative
        verifyGreaterThanOrEqual(testCase, gexci{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, gexci{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, gexci{i}.fallTime, 0);

        verifyGreaterThanOrEqual(testCase, gtail{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, gtail{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, gtail{i}.fallTime, 0);

        verifyGreaterThanOrEqual(testCase, glead{i}.riseTime, 0);
        verifyGreaterThanOrEqual(testCase, glead{i}.flatTime, 0);
        verifyGreaterThanOrEqual(testCase, glead{i}.fallTime, 0);

        Aexci(i) = trapArea(gexci{i}); % 1/m
        Atail(i) = trapArea(gtail{i}); % 1/m
        Alead(i) = trapArea(glead{i}); % 1/m
    end

    % Convert Pulseq area (1/m) to k (mm^-1)
    k_from_area_mmInv = @(A) 1e-3 * A;

    % === 1) Excitation achieves kstride ===
    kstride_expected = kend - kstart;             % mm^-1
    kstride_actual   = k_from_area_mmInv(Aexci);  % mm^-1

    tol_k = 1e-8 + 1e-6 * max(1, max(abs(kstride_expected)));
    verifyLessThanOrEqual(testCase, max(abs(kstride_actual - kstride_expected)), tol_k);

    % === 2) Tail design invariant (Pulseq-native form) ===
    % Original python (physical units):
    %   TailArea_T = -(1e3*kend)/gamma - ExciSlopArea_single_T
    %
    % Convert both sides to Hz/m*s (i.e., 1/m) by multiplying by gamma:
    %   TailArea_Hz = -(1e3*kend) - ExciSlopArea_single_Hz
    %
    % We can compute ExciSlopArea_single_Hz from returned excitation trapezoids:
    ExciAmp_Hzpm = cellfun(@(g) g.amplitude, gexci); % Hz/m
    ExciSlopDur  = gexci{1}.riseTime;                % s (common)

    ExciSlopArea_single_Hz = 0.5 * ExciAmp_Hzpm * ExciSlopDur; % (Hz/m)*s = 1/m

    Atail_expected = -(1e3 * kend) - ExciSlopArea_single_Hz;   % 1/m

    tol_area = 1e-4 + 1e-6 * max(1, max(abs(Atail_expected)));
    verifyLessThanOrEqual(testCase, max(abs(Atail - Atail_expected)), tol_area);

    % === 3) Lead design invariant (Pulseq-native form) ===
    % Original python:
    %   LeadArea_T = +(1e3*kend)/gamma - ExciSlopArea_single_T - ExciFlatArea_T
    % Multiply by gamma:
    %   LeadArea_Hz = +(1e3*kend) - ExciSlopArea_single_Hz - ExciFlatArea_Hz
    %
    % Compute ExciFlatArea_Hz from Aexci = ExciFlatArea + (two ramps):
    ExciSlopArea_total_Hz = ExciAmp_Hzpm * ExciSlopDur; % two ramps together
    ExciFlatArea_Hz       = Aexci - ExciSlopArea_total_Hz;

    Alead_expected = (1e3 * kend) - ExciSlopArea_single_Hz - ExciFlatArea_Hz;

    tol_area2 = 1e-4 + 1e-6 * max(1, max(abs(Alead_expected)));
    verifyLessThanOrEqual(testCase, max(abs(Alead - Alead_expected)), tol_area2);
end

function testHandlesZeroStrideAxis(testCase)
    gamma = 42.576e6;

    maxgrad_Tm  = 40e-3;
    maxslew_Tms = 150;

    maxgrad_Hzm  = gamma * maxgrad_Tm;
    maxslew_Hzms = gamma * maxslew_Tms;

    sys = mr.opts('MaxGrad', maxgrad_Hzm,  'GradUnit','Hz/m', ...
                  'MaxSlew', maxslew_Hzms, 'SlewUnit','Hz/m/s');

    kstart = [0, -1, 0];
    kend   = [0,  2, 0];

    [~, gexci, ~, ref] = CalcExcitationGradientWaveforms( ...
        kstart, kend, maxgrad_Tm, maxslew_Tms, gamma, sys);

    verifyEqual(testCase, numel(gexci), 3);
    verifyGreaterThanOrEqual(testCase, ref, 0);
end
