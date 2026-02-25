function results = run_all_tests()
%RUN_ALL_TESTS Run all unit tests in this repository.

repoRoot = fileparts(mfilename('fullpath'));

% Make sure your package code is on the path (adjust if needed)
addpath(genpath(repoRoot));

% Prefer tests folder
testsFolder = fullfile(repoRoot, "tests");

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TestReportPlugin
import matlab.unittest.plugins.CodeCoveragePlugin

% Build suite
if isfolder(testsFolder)
    suite = TestSuite.fromFolder(testsFolder, "IncludingSubfolders", true);
else
    % Fallback: run any tests discoverable on path
    suite = TestSuite.fromFolder(repoRoot, "IncludingSubfolders", true);
end

% Runner + report
runner = TestRunner.withTextOutput("Verbosity", 2);

% Optional: HTML report
reportsDir = fullfile(repoRoot, "test_reports");
if ~isfolder(reportsDir), mkdir(reportsDir); end
runner.addPlugin(TestReportPlugin.producingHTML(reportsDir));

% Optional: coverage for your package folder(s)
pkgFolder = fullfile(repoRoot, "+mypkg");  % change to your package folder
if isfolder(pkgFolder)
    runner.addPlugin(CodeCoveragePlugin.forFolder(pkgFolder));
end

% Run
results = runner.run(suite);

% Make CI-friendly: fail the script if any test fails
assert(all([results.Passed]), "Some unit tests failed.");
end