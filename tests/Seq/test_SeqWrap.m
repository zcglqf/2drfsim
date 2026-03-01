classdef test_SeqWrap < matlab.unittest.TestCase
    % Unit tests for Seq.SeqWrap using tests/ObjStretcher/res/epi.seq

    properties
        epiFile (1,1) string
        wrap
    end

    methods (TestMethodSetup)
        function setupEach(testCase)
            thisFile = string(mfilename("fullpath"));
            thisDir = fileparts(thisFile);

            testCase.epiFile = fullfile(thisDir, "res", "epi.seq");
            testCase.wrap = Seq.SeqWrap(); % no seq loaded initially
        end
    end

    methods (Test)
        function testConstructorWithoutSeq(testCase)

            testCase.verifyEmpty(testCase.wrap.seq);
            testCase.verifyEmpty(testCase.wrap.wfRf);
            testCase.verifyEmpty(testCase.wrap.wfGx);
            testCase.verifyEmpty(testCase.wrap.wfGy);
            testCase.verifyEmpty(testCase.wrap.wfGz);
            testCase.verifyEmpty(testCase.wrap.wfAdc);
        end

        function testConstructorWithSeqLoadsWaveforms(testCase)

            testCase.assumeTrue(isfile(testCase.epiFile), ...
                "epi.seq not found: " + testCase.epiFile);

            % create pulseq object
            s = mr.Sequence();
            s.read(testCase.epiFile);

            testCase.wrap = Seq.SeqWrap(s);

            testCase.verifyNotEmpty(testCase.wrap.seq);

            % wavefroms should have been filled automatically
            testCase.verifyEqual(length(testCase.wrap.wfRf), 3002);
            testCase.verifyEqual(length(testCase.wrap.wfGx), 259);
            testCase.verifyEqual(length(testCase.wrap.wfGy), 196);
            testCase.verifyEqual(length(testCase.wrap.wfGz), 8);
            testCase.verifyEqual(length(testCase.wrap.wfAdc), 4096);
            testCase.verifyEqual(testCase.wrap.getAdcDwell(), 4e-6, 'AbsTol', 1e-12);
        end

        function testLoadFromFileSetsSeqAndClearsCache(testCase)
            testCase.assumeTrue(isfile(testCase.epiFile), ...
                "epi.seq not found at: " + testCase.epiFile);

            % Pre-fill cache fields to make sure load clears them
            testCase.wrap.wfRf  = 1;
            testCase.wrap.wfGx  = 1;
            testCase.wrap.wfGy  = 1;
            testCase.wrap.wfGz  = 1;
            testCase.wrap.wfAdc = 1;

            testCase.verifyWarningFree(@() testCase.wrap.loadFromFile(testCase.epiFile));

            testCase.verifyFalse(isempty(testCase.wrap.seq));
            testCase.verifyEqual(length(testCase.wrap.wfRf), 3002);
            testCase.verifyEqual(length(testCase.wrap.wfGx), 259);
            testCase.verifyEqual(length(testCase.wrap.wfGy), 196);
            testCase.verifyEqual(length(testCase.wrap.wfGz), 8);
            testCase.verifyEqual(length(testCase.wrap.wfAdc), 4096);
            testCase.verifyEqual(testCase.wrap.getAdcDwell(), 4e-6, 'AbsTol', 1e-12);
        end

        function testGetWaveformUnknownChannelThrows(testCase)
            testCase.assumeTrue(isfile(testCase.epiFile), ...
                "epi.seq not found at: " + testCase.epiFile);

            testCase.wrap.loadFromFile(testCase.epiFile);

            testCase.verifyError(@() testCase.wrap.getWaveform("banana", false), ...
                "SeqWrap:UnknownChannel");
        end

    end
end