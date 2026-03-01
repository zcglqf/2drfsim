classdef test_Scanner < matlab.unittest.TestCase
    % Unit tests for Scanner.Scanner

    methods (Test)

        function testConstructorLoadsSys(testCase)
            sc = Scanner.Scanner("Philips.Ingenia1.5T.Ambition");

            testCase.verifyEqual(sc.prodId, "Philips.Ingenia1.5T.Ambition");
            testCase.verifyTrue(isstruct(sc.sys));
            testCase.verifyTrue(isfield(sc.sys, "maxGrad"));
            testCase.verifyTrue(isfield(sc.sys, "rfRasterTime"));
        end

        function testConstructorIsCaseInsensitive(testCase)
            sc1 = Scanner.Scanner("Philips.Ingenia1.5T.Ambition");
            sc2 = Scanner.Scanner("philips.ingenia1.5t.ambition");

            % prodId preserves input string, so compare sys contents instead
            testCase.verifyTrue(isequaln(sc1.sys, sc2.sys));
        end

        function testGettersMatchSysFields(testCase)
            sc = Scanner.Scanner("Siemens.Prisma3.0T");

            testCase.verifyEqual(sc.getMaxGrad(), sc.sys.maxGrad);
            testCase.verifyEqual(sc.getMaxSlew(), sc.sys.maxSlew);
            testCase.verifyEqual(sc.getRfRingdownTime(), sc.sys.rfRingdownTime);
            testCase.verifyEqual(sc.getRfDeadTime(), sc.sys.rfDeadTime);
            testCase.verifyEqual(sc.getAdcDeadTime(), sc.sys.adcDeadTime);
            testCase.verifyEqual(sc.getRfRasterTime(), sc.sys.rfRasterTime);
            testCase.verifyEqual(sc.getGradRasterTime(), sc.sys.gradRasterTime);
            testCase.verifyEqual(sc.getAdcRasterTime(), sc.sys.adcRasterTime);
            testCase.verifyEqual(sc.getBaseRasterTime(), sc.sys.baseRasterTime);
            testCase.verifyEqual(sc.getMaxB1(), sc.sys.maxB1);
        end

        function testDescribeReturnsString(testCase)
            sc = Scanner.Scanner("GE.Signa3.0T.Hero");

            txt = sc.describe();

            testCase.verifyTrue(ischar(txt) || isstring(txt));
            testCase.verifyNotEmpty(txt);

            % should contain product id and some keywords
            testCase.verifyTrue(contains(string(txt), "GE.Signa3.0T.Hero"));
            testCase.verifyTrue(contains(string(txt), "maxGrad="));
            testCase.verifyTrue(contains(string(txt), "rfRaster="));
        end

        function testUnknownModelThrows(testCase)
            testCase.verifyError(@() Scanner.Scanner("Unknown.Model"), ...
                "Scanner:UnknownProdModel");
        end

    end
end