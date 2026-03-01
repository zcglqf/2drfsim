function scanner = prodModel(prodId)
%Scanner.prodModel Create a Scanner object from a product model id string.
%
% Example:
%   scanner = Scanner.prodModel("Philips.Ingenia1.5T.Ambition");

prodId = string(prodId);

parts = parseProdId(prodId);

db = modelDb();
key = parts.vendor + "." + parts.series + parts.strengthToken + "." + parts.model;

if ~isfield(db, matlab.lang.makeValidName(key))
    error("Scanner:UnknownModel", "Unknown product model: %s", prodId);
end

entry = db.(matlab.lang.makeValidName(key));

% Build mr.opts sys from entry
sys = mr.opts( ...
    'MaxGrad', entry.maxGrad, 'GradUnit', entry.gradUnit, ...
    'MaxSlew', entry.maxSlew, 'SlewUnit', entry.slewUnit, ...
    'rfRingdownTime', entry.rfRingdownTime, ...
    'rfDeadTime', entry.rfDeadTime, ...
    'adcDeadTime', entry.adcDeadTime, ...
    'rfRasterTime', entry.rfRasterTime, ...
    'gradRasterTime', entry.gradRasterTime, ...
    'adcRasterTime', entry.adcRasterTime);

% Optional fields you may want
if isfield(entry, "maxB1")
    sys.maxB1 = entry.maxB1;
end

scanner = Scanner.Scanner(prodId, sys, parts);

end

% ---- helpers ----

function parts = parseProdId(prodId)
% Expect: Vendor.SeriesStrength.Model
% Example: Philips.Ingenia1.5T.Ambition

tokens = split(prodId, ".");
if numel(tokens) ~= 3
    error("Scanner:BadProdId", ...
        "prodId must have format Vendor.SeriesStrength.Model, got: %s", prodId);
end

vendor = tokens(1);
seriesStrength = tokens(2);
model = tokens(3);

% Parse seriesStrength like "Ingenia1.5T"
m = regexp(seriesStrength, "^(.*?)(\d+(\.\d+)?)T$", "tokens", "once");
if isempty(m)
    error("Scanner:BadProdId", ...
        "Cannot parse seriesStrength (expected like Ingenia1.5T): %s", seriesStrength);
end

series = string(m{1});
fieldStrengthT = str2double(m{2});
strengthToken = string(m{2}) + "T";

parts = struct( ...
    'vendor', vendor, ...
    'series', series, ...
    'fieldStrengthT', fieldStrengthT, ...
    'strengthToken', strengthToken, ...
    'model', model);
end

function db = modelDb()
% Minimal built-in database.
% You should expand this as you add models.

db = struct();

key = matlab.lang.makeValidName("Philips.Ingenia1.5T.Ambition");
db.(key) = struct( ...
    'maxGrad', 40, 'gradUnit', "mT/m", ...
    'maxSlew', 120, 'slewUnit', "T/m/s", ...
    'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6, ...
    'rfRasterTime', 1e-6, ...
    'gradRasterTime', 10e-6, ...
    'adcRasterTime', 0.1e-6, ...
    'maxB1', []); % optional, set if you want

end