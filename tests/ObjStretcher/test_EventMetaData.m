classdef test_EventMetaData < matlab.unittest.TestCase
    % Unit tests for ObjStretcher.EventMetaData

    methods (Test)

        function testDefaultCtor(testCase)
            m = ObjStretcher.EventMetaData();

            testCase.verifyEqual(m.type, ObjStretcher.EventType.RF);
            testCase.verifyEqual(m.tStart, 0);
            testCase.verifyEqual(m.tEnd, 0);
            testCase.verifyEqual(m.dur, 0);
            testCase.verifyEqual(m.vStart, complex(single(0)));
            testCase.verifyEqual(m.vEnd, complex(single(0)));
            testCase.verifyEqual(m.vAve, complex(single(0)));
        end

        function testNameValueCtor(testCase)
            m = ObjStretcher.EventMetaData( ...
                'type', ObjStretcher.EventType.GX, ...
                'tStart', 1e-3, ...
                'tEnd', 1.01e-3, ...
                'vStart', 2, ...
                'vEnd', 2, ...
                'vAve', 2);

            testCase.verifyEqual(m.type, ObjStretcher.EventType.GX);
            testCase.verifyEqual(m.tStart, 1e-3);
            testCase.verifyEqual(m.tEnd, 1.01e-3);
            testCase.verifyEqual(m.dur, 1e-5, 'AbsTol', 1e-12);

            testCase.verifyClass(m.vStart, 'single');
            testCase.verifyEqual(real(m.vStart), single(2));
            testCase.verifyEqual(imag(m.vStart), single(0));
        end

        function testComplexSingleCasting(testCase)
            m = ObjStretcher.EventMetaData( ...
                'tStart', 0, ...
                'tEnd', 1e-6, ...
                'vStart', complex(1,2), ...
                'vEnd', complex(1,2), ...
                'vAve', complex(1,2));

            testCase.verifyClass(m.vStart, 'single');
            testCase.verifyEqual(real(m.vStart), single(1));
            testCase.verifyEqual(imag(m.vStart), single(2));
        end

        function testEndBeforeStartThrows(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventMetaData( ...
                    'tStart', 2e-6, ...
                    'tEnd', 1e-6), ...
                'ObjStretcher:EventMetaData:badTimeOrder');
        end

        function testNameValueOddCountThrows(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventMetaData('tStart', 0, 'tEnd'), ...
                'ObjStretcher:EventMetaData:badArgs');
        end

        function testUnknownFieldThrows(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventMetaData('banana', 1), ...
                'ObjStretcher:EventMetaData:unknownField');
        end

    end
end