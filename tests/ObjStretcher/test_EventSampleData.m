classdef test_EventSampleData < matlab.unittest.TestCase
    % Unit tests for ObjStretcher.EventSampleData

    methods (Test)

        function testDefaultCtor(testCase)
            s = ObjStretcher.EventSampleData();

            testCase.verifyEqual(s.t, 0);
            testCase.verifyEqual(s.startOrEnd, false);
            testCase.verifyEqual(s.type, ObjStretcher.EventType.RF);
            testCase.verifyEqual(s.v, complex(single(0)));
        end

        function testNameValueCtorAndHelpers(testCase)
            s = ObjStretcher.EventSampleData( ...
                't', 1e-3, ...
                'startOrEnd', true, ...
                'type', ObjStretcher.EventType.GX, ...
                'value', complex(1, 2));

            testCase.verifyEqual(s.t, 1e-3, "AbsTol", 0);
            testCase.verifyEqual(s.startOrEnd, true);
            testCase.verifyEqual(s.type, ObjStretcher.EventType.GX);

            % value cast to complex single
            testCase.verifyClass(s.v, 'single');
            testCase.verifyTrue(~isreal(s.v));
            testCase.verifyEqual(real(s.v), single(1), "AbsTol", single(0));
            testCase.verifyEqual(imag(s.v), single(2), "AbsTol", single(0));

            % helpers
            testCase.verifyTrue(s.isStart());
            testCase.verifyFalse(s.isEnd());
        end

        function testIsEndWhenStartOrEndFalse(testCase)
            s = ObjStretcher.EventSampleData( ...
                't', 2e-3, ...
                'startOrEnd', false, ...
                'type', ObjStretcher.EventType.RF, ...
                'value', 3);

            testCase.verifyFalse(s.isStart());
            testCase.verifyTrue(s.isEnd());
        end

        function testBadArgsOddNameValueCountThrows(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventSampleData('t', 0, 'startOrEnd'), ...
                'ObjStretcher:EventSampleData:badArgs');
        end

        function testUnknownFieldThrows(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventSampleData('banana', 1), ...
                'ObjStretcher:EventSampleData:unknownField');
        end

        function testValidateRejectsNaNTime(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventSampleData('t', NaN), ...
                'ObjStretcher:EventSampleData:badTime');
        end

        function testValidateRejectsInfTime(testCase)
            testCase.verifyError(@() ...
                ObjStretcher.EventSampleData('t', Inf), ...
                'ObjStretcher:EventSampleData:badTime');
        end

    end
end