classdef test_ProdModelMap < matlab.unittest.TestCase
    % Unit tests for Scanner.ProdModelMap

    properties
        map
    end

    methods (TestMethodSetup)
        function setupEach(testCase)
            testCase.map = Scanner.ProdModelMap.getInstance();
        end
    end

    methods (Test)

        function testSingletonReturnsSameInstance(testCase)
            m1 = Scanner.ProdModelMap.getInstance();
            m2 = Scanner.ProdModelMap.getInstance();
            testCase.verifyTrue(m1 == m2);
        end

        function testCaseInsensitiveLookup(testCase)
            sys1 = testCase.map.get("Philips.Ingenia1.5T.Ambition");
            sys2 = testCase.map.get("philips.ingenia1.5t.ambition");
            sys3 = testCase.map.get("  PHILIPS.INGENIA1.5T.AMBITION  ");

            testCase.verifyTrue(isstruct(sys1));
            testCase.verifyTrue(isstruct(sys2));
            testCase.verifyTrue(isstruct(sys3));

            testCase.verifyTrue(isequal(sys1, sys2));
        end

        function testUnknownModelThrows(testCase)
            testCase.verifyError(@() testCase.map.get("Unknown.Vendor.Model"), ...
                "Scanner:UnknownProdModel");
        end

        function testAnyModelExists(testCase)
            sys = testCase.map.get("any");

            testCase.verifyTrue(isstruct(sys));
            testCase.verifyTrue(isfield(sys, "baseRasterTime"));
            testCase.verifyGreaterThan(double(sys.baseRasterTime), 0);
        end

        function testBaseRasterMultipleForAllKnownEntries(testCase)
            prodIds = [
                "any"
                "Philips.Ingenia1.5T.Ambition"
                "Siemens.Prisma3.0T"
                "GE.Signa3.0T.Hero"
            ];

            tol = 1e-12;

            for i = 1:numel(prodIds)
                sys = testCase.map.get(prodIds(i));

                base = double(sys.baseRasterTime);
                testCase.verifyGreaterThan(base, 0);

                testCase.verifyTrue(isMultiple(sys.rfRasterTime, base, tol));
                testCase.verifyTrue(isMultiple(sys.gradRasterTime, base, tol));
                testCase.verifyTrue(isMultiple(sys.adcRasterTime, base, tol));
            end
        end

        function testAddRejectsMissingBaseRasterTime(testCase)
            sys = mr.opts( ...
                'MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 1e-5, ...
                'adcRasterTime', 1e-7);

            % Do NOT set sys.baseRasterTime
            testCase.verifyError(@() testCase.map.add("Test.NoBase", sys), ...
                "Scanner:MissingBaseRasterTime");
        end

        function testAddRejectsRasterNotMultiple(testCase)
            sys = mr.opts( ...
                'MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 1e-5, ...
                'adcRasterTime', 1e-7);

            sys.baseRasterTime = 3e-7; % makes 1e-6 not an integer multiple

            testCase.verifyError(@() testCase.map.add("Test.BadMultiple", sys), ...
                "Scanner:RasterNotMultiple");
        end

        function testAddAcceptsValidRasterMultiple(testCase)
            sys = mr.opts( ...
                'MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, ...
                'gradRasterTime', 1e-5, ...
                'adcRasterTime', 1e-7);

            sys.baseRasterTime = 1e-7;

            testCase.verifyWarningFree(@() testCase.map.add("Test.Valid", sys));

            sys2 = testCase.map.get("test.valid"); % case-insensitive key
            testCase.verifyEqual(sys2.baseRasterTime, sys.baseRasterTime);
            testCase.verifyEqual(sys2.rfRasterTime, sys.rfRasterTime);
        end

    end
end

function tf = isMultiple(value, base, tol)
ratio = double(value) / double(base);
tf = abs(ratio - round(ratio)) <= tol;
end