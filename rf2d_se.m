%% demo_excitation_epi_like.m
% Requires:
%   - Pulseq MATLAB toolbox (mr.*)
%   - Your functions on path:
%       CalcExcitationGradientWaveforms.m
%       CalcRephasingGradientWaveformByArea.m  (optional, if you use inter-line rephasing)

clear; 
close all;
debugShowPlot = true;

%% ---------------------------
% 0) System / constants
% ----------------------------
% Use Hz/m units to match your current gradient math in 1/m.
sys = mr.opts('MaxGrad', 33, 'GradUnit','mT/m', ...
              'MaxSlew', 130, 'SlewUnit','T/m/s');
seq = mr.Sequence(sys);

% Indices: P=1, S=2
PORI = 1; SORI = 2;

% "2D excitation" example: M is "unlimited" so we only use P and S
FovM = 50e-3;
FovP = 50e-3;   % meters (40 mm)
FovS = 10e-3;    % meters (5 mm)
FOV  = [FovP, FovS];
Offc = [0.02, 0.02]'; % meters
flipAngle = pi/2;   % 90 deg excitation
tiltAngle = 0;
% norm vector of the k-space lines
klineNormVec = [0, 1]';
% unit vector along k-space lines
klineProcVec = [1, 0]';
TiltRot = [cos(tiltAngle) -sin(tiltAngle); sin(tiltAngle)  cos(tiltAngle)];
klineNormVec = TiltRot * klineNormVec;
klineProcVec = TiltRot * klineProcVec;
fprintf('K-Space line norm vector:\n');
fprintf('%g ', klineNormVec);   % transpose so fprintf prints row-wise
fprintf('\n');

% Off-center (meters) – used to add phase ramp across k
Offc = [20e-3, 20e-3]';   % (M,P,S) in meters

% Choose EPI style (between lines)
EPIType = "BiPolar";  % "FlyBack" or "BiPolar"

%% ---------------------------
% 1) RF 1D shapes (example)
% ----------------------------
% Normalize (optional but typical)
thisFile = mfilename('fullpath');          % full path to this .m file (without .m)
thisDir  = fileparts(thisFile);            % folder containing this .m file
matPath  = fullfile(thisDir, 'rfshapes', 'sincgauss21.mat');  % subfolder "data"
rf_shape = load(matPath);
amShape = rf_shape.shape.am;
amShape = amShape / max(amShape);
nrfs = numel(amShape);

% TBW + ref ratios (from your python object)
tbwP = rf_shape.shape.tbp;
tbwS = rf_shape.shape.tbp;
refP = rf_shape.shape.ref;
refS = rf_shape.shape.ref;

%% ---------------------------
% 2) Compute excitation k-bounds (1/m)
% ----------------------------
% k-size per axis: ksz = TBW / FOV
% (FOV in meters -> k in 1/m)
ksz = [tbwP/FovP, tbwS/FovS];
kmin = -ksz .* [refP, refS];
kmax = kmin + ksz;

fprintf('kmin(P,S) = [%g, %g] 1/m\n', kmin(PORI), kmin(SORI));
fprintf('kmax(P,S) = [%g, %g] 1/m\n', kmax(PORI), kmax(SORI));

% now we have rf waveform and k-space supports in P and S
kcoorP = linspace(kmin(PORI), kmax(PORI), nrfs);
kcoorS = linspace(kmin(SORI), kmax(SORI), nrfs);
rfwFuncP = griddedInterpolant(kcoorP, amShape, 'linear', 'nearest');
rfwFuncS = griddedInterpolant(kcoorS, amShape, 'linear', 'nearest');

%% ---------------------------
% 3) Build RF in excitation k-space: RfVal(P,S) = rfP(kP)*rfS(kS)
% ----------------------------
% Create 1D k grids for P and S matching the RF sample count
assert(numel(amShape) == nrfs);
t01 = linspace(0,1,nrfs).'; % normalized time index
kP = kmin(PORI) + (kmax(PORI)-kmin(PORI))*t01;
kS = kmin(SORI) + (kmax(SORI)-kmin(SORI))*t01;

% Pad endpoints like Python (optional)
kPpad = [2*kP(1)-kP(2); kP; 2*kP(end)-kP(end-1)];
kSpad = [2*kS(1)-kS(2); kS; 2*kS(end)-kS(end-1)];
rfPpad = [amShape(1); amShape; amShape(end)];
rfSpad = [amShape(1); amShape; amShape(end)];

% Build 2D grid values
[KP, KS] = ndgrid(kPpad, kSpad);
[RP, RS] = ndgrid(rfPpad, rfSpad);
RfVal2D = RP .* RS;

if debugShowPlot
    figure;
    surf(KP, KS, RfVal2D, 'EdgeColor','none');
    xlabel('kP (1/m)'); ylabel('kS (1/m)'); zlabel('RF weight');
    title('2D RF in excitation k-space = rfP(kP)*rfS(kS)');
    view(45,30); grid off;
end

%% ---------------------------
% 4) Generate 2D k-lines (P direction, stepped in S)
% ----------------------------
% Choose FOE step in meters for S (just like your python uses FOE and 1/FOE)
FOE_S = 200e-3;              % meters
dK_S  = 1/FOE_S;            % 1/m
% ------------------------------------


    % ---- corners of rectangle in (P,S) ----
    coordTR = [kmax(PORI), kmax(SORI)];
    coordBR = [kmax(PORI), kmin(SORI)];
    coordBL = [kmin(PORI), kmin(SORI)];
    coordTL = [kmin(PORI), kmax(SORI)];

    vBL2TR = coordTR - coordBL;
    vBR2TL = coordTL - coordBR;

    % ---- dk and safe margin ----
    % dk = 1.0 / fieldOfExcitation;

    kspaceAspectRatio = (kmax(SORI) - kmin(SORI)) / ...
                        (kmax(PORI) - kmin(PORI));

    % C++: safeMargin = sin(2*atan(tan(tilt)*aspectRatio)) * dk;
    safeMargin = abs(sin(2.0 * atan(tan(tiltAngle) * kspaceAspectRatio)) * abs(dK_S));
    % safeMargin = 0.2 * dK_S;

    % ---- pick diagonal that spans more in n-direction ----
    proj1 = dot(klineNormVec, vBL2TR);
    proj2 = dot(klineNormVec, vBR2TL);

    if abs(proj1) > abs(proj2)
        kProjLowEnd = dot(klineNormVec, coordBL) + sign(proj1) * safeMargin;   % projection of lower end corner to the normal vector
        KProjHighEnd = dot(klineNormVec, coordTR) - sign(proj1) * safeMargin;  % projection of higher end corner to the normal vector
        % dK_S = sign(proj1) * (dK_S);
    else
        kProjLowEnd = dot(klineNormVec, coordBR) + sign(proj2) * safeMargin;
        KProjHighEnd = dot(klineNormVec, coordTL) - sign(proj2) * safeMargin;
        % dK_S = sign(proj2) * (dK_S);
    end

    % ---- convert kmin/kmax to line indices ----
    numStart = int32(sign(kProjLowEnd) * floor(abs(kProjLowEnd) / abs(dK_S)));
    numStop  = int32(sign(KProjHighEnd) * floor(abs(KProjHighEnd) / abs(dK_S)));
    % numStart = round(kProjLowEnd / abs(dK_S));
    % numStop  = round(KProjHighEnd / abs(dK_S));

    if numStart > numStop 
        [numStart, numStop] = deal(numStop, numStart);
    end

    fprintf('K-Space acquisition starts from line no %d to line no %d\n', numStart, numStop);
    fprintf('Delta k = %f\n', dK_S);



% line equations of the 4 edges, represented by norm vectors and
% interceptions.
edgeNormVecs = [1, 0; 0, 1];
edgeNormVecs = kron(edgeNormVecs, [1, 1]');
edgeDists = [kmin(1), kmax(1), kmin(2), kmax(2)]';

% For each k-space line, its norm is n = (P, S) = (a11, a12), which is (0, 1) 
% without tilt. Interception (distance from origin to the line) is d1 = n * dK_S, 
% n belong to [numStart, numStop]; with rotation, n -> R * n;
% the line can be represented as 
%       (a11, a12) x (x, y)^T = d1
% On the other hand, the norm of an edge is, e.g, n = (1, 0), interception is e.g. kmin
% , the edge line can be represented as
%       (a21, a22) x (x, y)^T = d2
% the intersection between the the k-line and edge can be solved using
%       |a11, a12| |x|     d1
%       |a21, a22| |y|  =  d2
% mix all line and all 4 edges, can create a larger matrix Ax = b, and
% solve for the intersection points. suppose there are N lines, they have
% same norm vectors but different interceptions - 1 norm vector, N
% interceptions - leading to 1 row and N colums for the linear equation. 
% For the 4 edges, it will be 2 norm vectors and 4 interceptions - leading to 4 rows for the linear equaitons.

% Each k-line equation, combined with each of the 4 k-space edges, forms a linear system A*x = b.
% Together, A will be a 8 x 8 matrix, unknow marix X will be a 8 x 4N block matrix. and b is 8 x 4N too. 
% Solving these systems gives the intersection points between each k-line and the k-space edges.
% From the two valid intersection points (endpoints) of each k-line, we can compute the gradient
% and RF waveforms
%
% To do this efficiently for many lines, we build a large block-diagonal (partitioned) matrix A
% and a corresponding stacked vector/matrix b, then solve for all intersection points at once.
% The solution matrix X is organized in the same partitioned way. In the sketch below, there are
% 5 k-lines:
%
%        _         A        _   _                 X                   _     _                 B                  _
%       |  a a <--edge norm  | |  x x x x x                              |   |  b b b b b  <-- edge interc             |
%       |  a a <--line norm  | |  x x x x x                              |   |  b b b b b  <-- line interc             |
%       |      a a           | |            x x x x x                    |   |            b b b b b                    |
%       |      a a           | |            x x x x x                    | = |            b b b b b                    |
%       |          a a       | |                      x x x x x          |   |                      b b b b b          |
%       |          a a       | |                      x x x x x          |   |                      b b b b b          |
%       |              a a   | |                                x x x x x|   |                                b b b b b|
%       |_             a a  _| |_                               x x x x x_|   |_                               b b b b b_|

A = zeros(8, 2);
A(1:2:end) = edgeNormVecs;
A(2:2:end) = repmat(klineNormVec', [4, 1]);
A = blkdiag(A(1:2,:), A(3:4,:), A(5:6,:), A(7:8,:));

lineDists = (numStart : sign(numStop - numStart) : numStop).' * dK_S; % S-planes
nLines = numel(lineDists);

B = zeros(8, nLines);
B(1:2:end) = repmat(edgeDists, 1, nLines);
B(2:2:end) = repmat(lineDists', 4, 1);
B = blkdiag(B(1:2,:), B(3:4,:), B(5:6,:), B(7:8,:));

if rcond(A) < 1e-12
    X = pinv(A) * B;          % SVD-based pseudoinverse (robust, slower)
else
    X = A\B;
end

% Organize 8x4N block matrix into 8xN matrix 
endPoints = [X(1:2, 1:nLines); ...
    X(3:4,nLines+1:2*nLines); ...
    X(5:6,2*nLines+1:3*nLines); ...
    X(7:8,3*nLines+1:4*nLines)];
endPoints = reshape(endPoints, [2, 4, nLines]); % coordinates, edge, line

% plot before pruning
figure; hold on; axis equal; grid on;
C = turbo(nLines);
bx = [kmin(PORI) kmax(PORI) kmax(PORI) kmin(PORI) kmin(PORI)];
by = [kmin(SORI) kmin(SORI) kmax(SORI) kmax(SORI) kmin(SORI)];
hFirst = []; hLast = [];
for i = 1:nLines
    h = plot(squeeze(endPoints(1,:,i)), squeeze(endPoints(2,:,i)), '-o', 'Color', C(i,:), 'LineWidth', 1.5);
    if i == 1, hFirst = h; end
    if i == nLines, hLast  = h; end
end
legend([hFirst hLast], {'First line (cool)', 'Last line (warm)'}, 'Location','best');
plot(bx, by, 'k-', 'LineWidth', 1.5); hold off; title('k-space lines before prunning');

pntOnEdge = isOnRectEdge(squeeze(endPoints(1,:,:)), squeeze(endPoints(2,:,:)), kmin(PORI), kmax(PORI), kmin(SORI), kmax(SORI), 1e-3*abs(dK_S));
assert(all(sum(pntOnEdge) == 2), 'Each k-space line has two intersection points with edges of k-space.');
endPoints = endPoints(:, pntOnEdge);
endPoints = reshape(endPoints, [2, 2, nLines]); % coordinates, 2 edges expected, line

% sort the endPoints order, so that they align with line proceed direction
% endPoints: 2x2xN
%   endPoints(:,1,i) = point1 [x;y]
%   endPoints(:,2,i) = point2 [x;y]
dpnt = squeeze(endPoints(:,2,:) - endPoints(:,1,:));   % 2xN, coordiantes, line
dotv = klineProcVec.' * dpnt;                          % 1xN

swapMask = dotv < 0.0;

% swap point1 and point2 for the lines that are opposite to u
endPoints(:, [1 2], swapMask) = endPoints(:, [2 1], swapMask);



x = squeeze(endPoints(1,:,:));   % 2x30
y = squeeze(endPoints(2,:,:));   % 2x30
figure; hold on; axis equal; grid on;
hFirst = []; hLast = [];
for i = 1:nLines
    % h = plot(x(:,i), y(:,i), '-o', 'Color', C(i,:), 'LineWidth', 1.5);
    h = quiver(endPoints(1,1,i), endPoints(2,1,i), ...
        endPoints(1,2,i) - endPoints(1,1,i), endPoints(2,2,i) - endPoints(2,1,i), ...
        0, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'Color', C(i,:));
    h.AutoScale = 'off';
    if i == 1, hFirst = h; end
    if i == nLines, hLast  = h; end
end

plot(bx, by, 'k-', 'LineWidth', 1.5);

legend([hFirst hLast], {'First line (cool)', 'Last line (warm)'}, 'Location','best');

hold off;


% %% ---------------------------
% % 5) Build sequence: add 3 blocks per line (lead -> exci(+RF) -> tail)
% % ----------------------------
prevTail = [];  
prevExci = [];
prevLead = [];
Nsamp = 1e4; 

% estimate total k-space line length and the gradient object that can
% achieve it.
diffLines = squeeze(endPoints(:,2,:) - endPoints(:,1,:));  % 2xN
sumDiff = sum(diffLines, 2);  
[~,dummyExci,~] = CalcExcitationGradientWaveforms([0,0]', sumDiff, sys);
dummytotflatdur = dummyExci{1}.flatTime;

% estimated total rf dur. The length of each line is different. We scale
% the no of samples of each rf pulse against the center line, with the RF 
% pulse of the center line having 256 samples.
ilineCen = (0 - numStart + 1);
klineLength = squeeze(vecnorm(endPoints(:,2,:) - endPoints(:,1,:), 2, 1)); % 2-norm, 1-first column
klineRatio = klineLength / klineLength(ilineCen);
RfShapes  = cell(1, nLines);
for li = 1:nLines
    kstart = squeeze(endPoints(:,1,li));
    kend   = squeeze(endPoints(:,2,li));
    nsamples = round(klineRatio(li) * Nsamp);
    % --- Build RF waveform along this k-line by sampling RfVal2D ---
    t = linspace(0, 1, nsamples);               % 1xN
    coord = kstart + t.*(kend - kstart);        % 2xN, N points from A to B
    am = rfwFuncP(coord(PORI,:)) .* rfwFuncS(coord(SORI,:));
    RfShapes{li} = am;
end
% compare the Rf total dur and Gradient flattop total dur
dummyrfshape = [RfShapes{:}];
dummyrf = mr.makeArbitraryRf(dummyrfshape, flipAngle, 'system', sys);
dummyrfdur = dummyrf.shape_dur;
sys.maxB1 = 20e-6 * sys.gamma; % do this for pulseq v1.4.2
% adjust number of samples, therefore dur, to fully utilize max B1 
% dummyrf.signal has absolute rf signal amplitudes in Hz
rfScaleFactor = max(abs(dummyrf.signal)) / (sys.maxB1 * 0.98);
dummyrfdurMaxB1 = dummyrfdur * rfScaleFactor;
% final dur for both rf and gradients flat top.
dummyFinalDur = max(dummyrfdurMaxB1, dummytotflatdur);
rfScaleFactor = dummyFinalDur / dummyrfdur;
grScaleFactor = dummyFinalDur / dummytotflatdur;
% adjust samples for rf pulses, this is the only way to change their
% duratioon
Nsamp = Nsamp * rfScaleFactor;
% generate the RF again with the right no of samples.
for li = 1:nLines
    kstart = squeeze(endPoints(:,1,li));
    kend   = squeeze(endPoints(:,2,li));
    nsamples = round(klineRatio(li) * Nsamp);
    % --- Build RF waveform along this k-line by sampling RfVal2D ---
    t = linspace(0, 1, nsamples);               % 1xN
    coord = kstart + t.*(kend - kstart);        % 2xN, N points from A to B
    am = rfwFuncP(coord(PORI,:)) .* rfwFuncS(coord(SORI,:));
    RfShapes{li} = am;
end
dummyrfshape = [RfShapes{:}];
dummyrf = mr.makeArbitraryRf(dummyrfshape, flipAngle, 'system', sys);
dummyrfdur = dummyrf.shape_dur;
% stores the signal back into the cell array
startIdx = 1;
endIdx = length(RfShapes{1});
for li = 1:nLines
    RfShapes{li} = dummyrf.signal(startIdx:endIdx)';
    startIdx = endIdx+1;
    if li < nLines
        endIdx = endIdx + length(RfShapes{li+1});
    end
end

% generate a new system object with new gradient constraint.
sysg = sys;
sysg.maxGrad = sys.maxGrad / grScaleFactor;

% now we can calcualte the gradients and rf pulses.
% LeadGrads = cell(1, nLines);
% ExciGrads = cell(1, nLines);
% TailGrads = cell(1, nLines);
% RfShapes  = zeros(Nsamp, nLines);
anchorTime = 0;
kspCenTime = 0; % k-space center time
for li = 1:nLines
     bFlip = (mod(li,2)==0) && (EPIType=="BiPolar");
     if bFlip
         kstart = squeeze(endPoints(:,2,li));
         kend   = squeeze(endPoints(:,1,li));
     else
         kstart = squeeze(endPoints(:,1,li));
         kend   = squeeze(endPoints(:,2,li));
     end
     [gLead, gExci, gTail] = CalcExcitationGradientWaveforms(kstart, kend, sysg);

    % --- Build RF waveform along this k-line by sampling RfVal2D ---
    nsamples = length(RfShapes{li});
    t = linspace(0, 1, nsamples);               % 1xN
    coord = kstart + t.*(kend - kstart);     % 2xN, N points from A to B

    ph = 2*pi * Offc.' * coord;
    b1_signal = RfShapes{li} .* exp(1j*ph);
    subFA_rad = 2*pi * sum(RfShapes{li}) * sys.rfRasterTime;
    % --- Optional: merge prev tail + current lead into a rephasing block (like Python) ---
    if li > 1
        % Compute areas in 1/m from trapezoids (full area)
        totalRephArea = cellfun(@(s) s.area, prevTail) + cellfun(@(s) s.area, gLead);
        % If you have CalcRephasingGradientWaveformByArea (expects 1/m), use it:
        gReph = CreateRephasingGradientWaveformByArea(totalRephArea, sys, 0); % slope=0 -> auto        
        gReph{1}.channel = 'x';
        gReph{2}.channel = 'y';
        % Add rephasing as its own block
        seq.addBlock(gReph{1}, gReph{2});
        anchorTime = anchorTime + mr.calcDuration(gReph{1});
    end
%     % Add off-center phase term: phase = -2*pi * k·Offc (k in 1/m, Offc in m)
%     kLineSamples = [zeros(Nsamp,1), kPline, kSline];
%     ph = -2*pi * (kLineSamples * Offc(:));   % radians
%     b1 = rfw .* exp(1j*ph);
% 
%     % Resample b1 to RF raster length matching excitation plateau
%     Tflat = gExci{1}.flatTime;
%     Nr = max(1, round(Tflat / sys.rfRasterTime));
%     b1r = interp1(linspace(0,1,Nsamp), b1, linspace(0,1,Nr), 'linear', 0).';
% 
    rf = mr.makeArbitraryRf(RfShapes{li}, subFA_rad, 'use','excitation', 'system', sys);
    % rf.signal = abs(b1_signal);
    % rf.t = (0:nsamples-1).' * sys.rfRasterTime; 
    rf.delay = gExci{1}.riseTime;   % start RF on plateau
% 
    % Add excitation gradients + RF in the same block
    gExci{1}.channel = 'x';
    gExci{2}.channel = 'y';
    seq.addBlock(gExci{1}, gExci{2}, rf);
    
    if li == ilineCen
        t = linspace(0, 1, 4096);
        coord = kstart + t.*(kend - kstart);
        [~, idxc] =  min(vecnorm(coord, 2, 1));
        kspCenTime = anchorTime + gExci{1}.riseTime + idxc / 4096 * gExci{1}.flatTime;
        fprintf('k-space center time is %f\n', kspCenTime);
    end

    anchorTime = anchorTime + mr.calcDuration(gExci{1});
% 
%     % Tail as its own block
%     addGradBlock(seq, gTail);
% 
    % store for next iteration
    prevTail = gTail;
    prevExci = gExci;
    prevLead = gLead;
end
gTail{1}.channel = 'x';
gTail{2}.channel = 'y';
seq.addBlock(gTail{1}, gTail{2});

exciFinTime = seq.duration();

%% Define 180 refocusing part and EPI acquisition.
fov = 256e-3;
Nx = 32; Ny = 32;
% FOVM = 200e-3;
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 3.2e-4;
gx = mr.makeTrapezoid('x',sys,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,sys,'Duration',gx.flatTime,'Delay',gx.riseTime);

% Pre-phasing gradients
preTime=8e-4;
%gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2-deltak/2,'Duration',preTime);
% gzReph = mr.makeTrapezoid('z',sys,'Area',-gz.area/2,'Duration',preTime);
%gyPre = mr.makeTrapezoid('y',sys,'Area',-Ny/2*deltak,'Duration',preTime);
% we need no minus for in-plane prephasers because of the spin-echo (position reflection in k-space)
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2 + deltak,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',sys,'Area',-Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/sys.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',sys,'Area',deltak,'Duration',dur);

% Refocusing pulse with spoiling gradients
% rf180 = mr.makeBlockPulse(pi,sys,'Duration',500e-6,'use','refocusing');
[rf180, gz180] = mr.makeSincPulse(pi,sys,'Duration', 4000e-6,...
    'SliceThickness',FovS*1.3,'apodization',0.5,'timeBwProduct',4,...
    'use','refocusing');
gz180.channel = 'y';
gzSpoil = mr.makeTrapezoid('y',sys,'Area', 2*6/FovS,'Duration',3*preTime);

% minTE1, minTE2
durationToCenter = (Ny/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy); % start of EPI train to echo center
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
exciFinTime = round(exciFinTime / sys.gradRasterTime) * sys.gradRasterTime;
kspCenTime = round(kspCenTime / sys.gradRasterTime) * sys.gradRasterTime;
minTE1 = (exciFinTime - kspCenTime)  + mr.calcDuration(gzSpoil) + (rf180centerInclDelay);
minTE2 = (rf180centerInclDelay) +  mr.calcDuration(gzSpoil) + durationToCenter + mr.calcDuration(gxPre);
minTEHalf = max(minTE1, minTE2);
minTEHalf = round(minTEHalf / sys.gradRasterTime) * sys.gradRasterTime;
% Calculate delay time %% MZ: I thisk this is very wrong!
% minTEhalf = (exciFinTime - kspCenTime) + mr.calcDuration(gzSpoil) + (rf180centerInclDelay);
TE = minTEHalf * 2; % 120e-3; %minTEhalf * 2;


% rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);

% delayTE1=TE/2 - (mr.calcDuration(gz) - rfCenterInclDelay) - preTime - mr.calcDuration(gzSpoil) - rf180centerInclDelay;

delayTE1=TE/2 - (exciFinTime - kspCenTime) - preTime - mr.calcDuration(gzSpoil) - rf180centerInclDelay;
delayTE2=TE/2 - (mr.calcDuration(rf180) - rf180centerInclDelay) - mr.calcDuration(gzSpoil) - durationToCenter;
delayTE1 = round(delayTE1 / sys.gradRasterTime) * sys.gradRasterTime;
delayTE2 = round(delayTE2 / sys.gradRasterTime) * sys.gradRasterTime;

% Define sequence blocks
% seq.addBlock(rf,gz);
if delayTE1 > sys.gradRasterTime
    seq.addBlock(mr.makeDelay(delayTE1));
end
seq.addBlock(gzSpoil);
seq.addBlock(rf180, gz180);
seq.addBlock(gzSpoil);
if delayTE2 > sys.gradRasterTime
    seq.addBlock(mr.makeDelay(delayTE2));
end
seq.addBlock(gxPre,gyPre); 
for i=1:Ny
    seq.addBlock(gx,adc);           % Read one line of k-space
    seq.addBlock(gy);               % Phase blip
    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
end
seq.addBlock(mr.makeDelay(1e-4));

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% export and visualization
seq.setDefinition('FOV', [fov fov FovS]);
seq.setDefinition('Name', 'epise');
seq.write('epi_se.seq');   % Output sequence for scanner
seq.plot();             % Plot sequence waveforms

%% calculate trajectory 
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

%% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display

%% sanity checks
TE_check=(t_refocusing(1)-t_excitation(1))*2;
fprintf('intended TE=%.03f ms, actual spin echo TE=%.03fms\n', TE*1e3, TE_check*1e3); 
% % 
% %% ---------------------------
% % 6) Plot the sequence
% % ----------------------------
seq.plot('timeDisp','ms');







function onEdge = isOnRectEdge(x, y, xmin, xmax, ymin, ymax, tol)
    if nargin < 7, tol = 1e-12; end

    % inside = (x >= xmin-tol) & (x <= xmax+tol) & (y >= ymin-tol) & (y <= ymax+tol);
    % onBoundaryValue = (abs(x-xmin) <= tol) | (abs(x-xmax) <= tol) | ...
    %                   (abs(y-ymin) <= tol) | (abs(y-ymax) <= tol);

    onBoundaryValue = ((abs(x-xmin) <= tol) & ((y >= ymin-tol) & (y <= ymax+tol))) | ...   % left edge
                  ((abs(x-xmax) <= tol) & ((y >= ymin-tol) & (y <= ymax+tol))) | ...   % right edge
                  ((abs(y-ymin) <= tol) & ((x >= xmin-tol) & (x <= xmax+tol))) | ...   % bottom edge
                  ((abs(y-ymax) <= tol) & ((x >= xmin-tol) & (x <= xmax+tol)));        % top edge
    onEdge = onBoundaryValue;
end



