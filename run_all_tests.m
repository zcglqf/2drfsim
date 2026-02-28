function results = run_all_tests(testFile)
%RUN_ALL_TESTS Run all unit tests, or run a single specified test file.
%
% Usage:
%   run_all_tests()                                   % run all tests under ./tests
%   run_all_tests("tests/ObjStretcher/test_Seg.m")    % run one test file
%   run_all_tests("test_EventMetaData.m")             % run one test file (searched under ./tests)
%
% Returns:
%   results : matlab.unittest.TestResult array

repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(repoRoot));

testsFolder = fullfile(repoRoot, "tests");

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TestReportPlugin

% Build suite
if nargin < 1 || isempty(testFile)
    if isfolder(testsFolder)
        suite = TestSuite.fromFolder(testsFolder, "IncludingSubfolders", true);
    else
        suite = TestSuite.fromFolder(repoRoot, "IncludingSubfolders", true);
    end
else
    tf = string(testFile);

    % If they pass only a filename, search under tests/
    if ~isfile(tf)
        candidates = dir(fullfile(testsFolder, "**", tf));
        if isempty(candidates)
            error("run_all_tests:TestFileNotFound", ...
                "Test file not found: %s", tf);
        end
        tf = fullfile(candidates(1).folder, candidates(1).name);
    end

    suite = TestSuite.fromFile(char(tf));
end

% Runner + HTML report
runner = TestRunner.withTextOutput("Verbosity", 2);

reportsDir = fullfile(repoRoot, "test_reports");
if ~isfolder(reportsDir), mkdir(reportsDir); end
runner.addPlugin(TestReportPlugin.producingHTML(reportsDir));

% Run
results = runner.run(suite);

% Fail if any test failed
assert(all([results.Passed]), "Some unit tests failed.");

end
