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
              'MaxSlew', 100, 'SlewUnit','T/m/s');
seq = mr.Sequence(sys);

% Indices: P=1, S=2
PORI = 1; SORI = 2;

% "2D excitation" example: M is "unlimited" so we only use P and S
FovM = inf;
FovP = 40e-3;   % meters (40 mm)
FovS = 8e-3;    % meters (5 mm)
FOV  = [FovP, FovS];
Offc = [0.02, 0.02]'; % meters
flipAngle = pi/6;   % 90 deg excitation
tiltAngle = pi/3;
% norm vector of the k-space lines
klineNormVec = [0, 1]';
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
EPIType = "FlyBack";  % "FlyBack" or "BiPolar"

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
N = numel(amShape);
t01 = linspace(0,1,N).'; % normalized time index
kP = kmin(1) + (kmax(1)-kmin(1))*t01;
kS = kmin(2) + (kmax(2)-kmin(2))*t01;

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
FOE_S = 50e-3;              % meters
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



% length of traspass across the k-space region by norm vector
edgeNormVecs = [1, 0; 0, 1];
edgeNormVecs = kron(edgeNormVecs, [1, 1]');
edgeDists = [kmin(1), kmax(1), kmin(2), kmax(2)]';

% For each k-space line, its norm is n = (P, S) = (a11, a12) = (0, 1), distance to origin is
% d1 = n * dK_S, n belong to [numStart, numStop]; with rotation, n -> R * n;
% the line can be represented as 
%       (a11, a12) x (x, y)^T = d1
% On the other hand, the norm of an edge is n = (1, 0), distance is (e.g.) kmin
% , the edge line can be represented as
%       (a21, a22) x (x, y)^T = d2
% the intersection between the the k-line and edge can be solved using
%       |a11, a12| |x|     d1
%       |a21, a22| |y|  =  d2
% mix all line and all 4 edges, can create a larger matrix Ax = b, and
% solve for the intersection points.

% Each k-line equation, combined with each of the 4 k-space edges, forms a linear system A*x = b.
% Solving these systems gives the intersection points between each k-line and the k-space edges.
% From the two valid intersection points (endpoints) of each k-line, we can compute the gradient
% and RF waveforms.
%
% To do this efficiently for many lines, we build a large block-diagonal (partitioned) matrix A
% and a corresponding stacked vector/matrix b, then solve for all intersection points at once.
% The solution matrix X is organized in the same partitioned way. In the sketch below, there are
% 5 k-lines:
%
%        _         A        _   _                 X                   _     _                 B                  _
%       |  a a               | |  x x x x x                              |   |  b b b b b                              |
%       |  a a               | |  x x x x x                              |   |  b b b b b                              |
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

endPoints = [X(1:2, 1:nLines); ...
    X(3:4,nLines+1:2*nLines); ...
    X(5:6,2*nLines+1:3*nLines); ...
    X(7:8,3*nLines+1:4*nLines)];
endPoints = reshape(endPoints, [2, 4, nLines]); % coordinates, line, edge

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
plot(bx, by, 'k-', 'LineWidth', 1.5); hold off;

pntOnEdge = isOnRectEdge(squeeze(endPoints(1,:,:)), squeeze(endPoints(2,:,:)), kmin(1), kmax(1), kmin(2), kmax(2), 1e-3*abs(dK_S));
assert(all(sum(pntOnEdge) == 2), 'Each k-space line has two intersection points with edges of k-space.');
endPoints = endPoints(:, pntOnEdge);
endPoints = reshape(endPoints, [2, 2, nLines]);

% sort that endPoints order, so that they align with line proceed direction
% endPoints: 2x2xN
%   endPoints(:,1,i) = point1 [x;y]
%   endPoints(:,2,i) = point2 [x;y]
dpnt = squeeze(endPoints(:,2,:) - endPoints(:,1,:));   % 2xN, direction point1->point2
dotv = klineProcVec.' * dpnt;                          % 1xN

swapMask = dotv < 0.0;

% swap point1 and point2 for the lines that are opposite to u
endPoints(:, [1 2], swapMask) = endPoints(:, [2 1], swapMask);


% plot lines and box
% c0 = [0 1 0];   % green
% c1 = [0 0 1];   % blue
% t  = linspace(0,1,nLines).';
% C  = (1-t).*c0 + t.*c1;   % Nx3 colors
% C = turbo(nLines);

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
Nsamp = 256; 

% estimate total k-space line length
diffLines = squeeze(endPoints(:,2,:) - endPoints(:,1,:));  % 2x18
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
    kstart = squeeze(endPoints(:,2,li));
    kend   = squeeze(endPoints(:,1,li));
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
% adjust number of samples, therefore dur to fully utilize max B1 
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
    kstart = squeeze(endPoints(:,2,li));
    kend   = squeeze(endPoints(:,1,li));
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
     % LeadGrads{li} = gLead;
     % ExciGrads{li} = gExci;
     % TailGrads{li} = gTail;
     % gLead{1}.delay = anchorTime;
     % gLead{2}.delay = anchorTime;
     % gExci{1}.delay = anchorTime + mr.calcDuration(gLead{1});
     % gExci{2}.delay = anchorTime + mr.calcDuration(gLead{2});
     

    % --- Build RF waveform along this k-line by sampling RfVal2D ---
    nsamples = length(RfShapes{li});
    t = linspace(0, 1, nsamples);               % 1xN
    coord = kstart + t.*(kend - kstart);     % 2xN, N points from A to B
    % am = rfwFuncP(coord(PORI,:)) .* rfwFuncS(coord(SORI,:));
    % RfShapes(:, li) = am;
    ph = 2*pi * Offc.' * coord;
    b1_signal = RfShapes{li} .* exp(1j*ph);
    % --- Optional: merge prev tail + current lead into a rephasing block (like Python) ---
    if li > 1
        % Compute areas in 1/m from trapezoids (full area)
        totalRephArea = cellfun(@(s) s.area, prevTail) + cellfun(@(s) s.area, gLead);
        % If you have CalcRephasingGradientWaveformByArea (expects 1/m), use it:
        gReph = CreateRephasingGradientWaveformByArea(totalRephArea, sys, 0); % slope=0 -> auto        
        gReph{1}.channel = 'y';
        gReph{2}.channel = 'z';
        % Add rephasing as its own block
        seq.addBlock(gReph{1}, gReph{2});
    else
        % First line: add lead as its own block
        % seq.addBlock(seq, gLead{1}, gLead{2});
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
    rf = mr.makeArbitraryRf(ones(1, nsamples), pi/100, 'system', sys);
    rf.signal = b1_signal;
    rf.t = (0:nsamples-1).' * sys.rfRasterTime; 
    rf.delay = gExci{1}.riseTime;   % start RF on plateau
% 
    % Add excitation gradients + RF in the same block
    gExci{1}.channel = 'y';
    gExci{2}.channel = 'z';
    seq.addBlock(gExci{1}, gExci{2}, rf);
% 
%     % Tail as its own block
%     addGradBlock(seq, gTail);
% 
    % store for next iteration
    prevTail = gTail;
    prevExci = gExci;
    prevLead = gLead;
end
gTail{1}.channel = 'y';
gTail{2}.channel = 'z';
seq.addBlock(gTail{1}, gTail{2});

% % 
% %% ---------------------------
% % 6) Plot the sequence
% % ----------------------------
seq.plot('timeDisp','ms');



% %% ---------------------------
% % Helper functions
% % ----------------------------
% function addGradBlock(seq, gcell)
%     % Add up to 3 gradients as a single block
%     if isempty(gcell), return; end
%     if numel(gcell)==3
%         seq.addBlock(gcell{1}, gcell{2}, gcell{3});
%     elseif numel(gcell)==2
%         seq.addBlock(gcell{1}, gcell{2});
%     else
%         seq.addBlock(gcell{1});
%     end
% end

% function A = trapArea(g)
%     % Full trapezoid area in Pulseq gradient units (Hz/m*s = 1/m)
%     A = g.amplitude * (g.flatTime + 0.5*(g.riseTime + g.fallTime));
% end


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



