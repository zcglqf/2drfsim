classdef test_Seg < matlab.unittest.TestCase
    % Unit tests for ObjStretcher.Seg
    methods (Test)

        function testDefaultCtor(testCase)
            s = ObjStretcher.Seg();
            testCase.verifyEqual(s.tStart, 0);
            testCase.verifyEqual(s.dur, 0);
            testCase.verifyEqual(s.gx, single(0));
            testCase.verifyEqual(s.gy, single(0));
            testCase.verifyEqual(s.gz, single(0));
            testCase.verifyEqual(s.adc, uint8(0));
        end

        function testCtorPositionalAndDerivedTimes(testCase)
            t0 = 1.25;
            dur = 2e-6;
            rf  = complex(single(1.5), single(-0.25));
            gx  = 3e-3;
            gy  = 0;
            gz  = -4e-3;
            adc = 0;

            s = ObjStretcher.Seg(t0, dur, rf, gx, gy, gz, adc);

            testCase.verifyEqual(s.tStart, t0);
            testCase.verifyEqual(s.dur, dur);

            testCase.verifyEqual(s.tEnd, t0 + dur, "AbsTol", 0);
            testCase.verifyEqual(s.tCenter, t0 + 0.5*dur, "AbsTol", 0);

            testCase.verifyEqual(s.gx, single(gx));
            testCase.verifyEqual(s.gy, single(gy));
            testCase.verifyEqual(s.gz, single(gz));
            testCase.verifyEqual(s.adc, uint8(adc));

            % rf stored as complex single
            testCase.verifyClass(s.rf, 'single');
            testCase.verifyTrue(~isreal(s.rf));
            testCase.verifyEqual(real(s.rf), real(rf), "AbsTol", single(0));
            testCase.verifyEqual(imag(s.rf), imag(rf), "AbsTol", single(0));
        end

        function testCtorNameValue(testCase)
            s = ObjStretcher.Seg( ...
                'tStart', 0.1, ...
                'dur', 1e-6, ...
                'rf', complex(2, 3), ...
                'gx', 1, ...
                'gy', 2, ...
                'gz', 3, ...
                'adc', 0);

            testCase.verifyEqual(s.tStart, 0.1);
            testCase.verifyEqual(s.dur, 1e-6);
            testCase.verifyEqual(s.gx, single(1));
            testCase.verifyEqual(s.gy, single(2));
            testCase.verifyEqual(s.gz, single(3));
            testCase.verifyEqual(s.adc, uint8(0));

            testCase.verifyClass(s.rf, 'single');
            testCase.verifyTrue(~isreal(s.rf));
            testCase.verifyEqual(real(s.rf), single(2), "AbsTol", single(0));
            testCase.verifyEqual(imag(s.rf), single(3), "AbsTol", single(0));
        end

        function testRfSetterForcesComplexSingleEvenIfReal(testCase)
            s = ObjStretcher.Seg('tStart', 0, 'dur', 1e-6, 'rf', 5); % real input
            testCase.verifyClass(s.rf, 'single');
            testCase.verifyTrue(~isreal(s.rf)); % should become complex with imag=0
            testCase.verifyEqual(real(s.rf), single(5), "AbsTol", single(0));
            testCase.verifyEqual(imag(s.rf), single(0), "AbsTol", single(0));
        end

        function testToStruct(testCase)
            s = ObjStretcher.Seg(0, 1e-6, complex(1,2), 3,4,5,0);
            st = s.toStruct();

            testCase.verifyTrue(isstruct(st));
            testCase.verifyEqual(st.tStart, s.tStart);
            testCase.verifyEqual(st.dur, s.dur);
            testCase.verifyEqual(st.rf, s.rf);
            testCase.verifyEqual(st.gx, s.gx);
            testCase.verifyEqual(st.gy, s.gy);
            testCase.verifyEqual(st.gz, s.gz);
            testCase.verifyEqual(st.adc, s.adc);
        end

        function testValidateRejectsNonPositiveDuration(testCase)
            testCase.verifyError(@() ObjStretcher.Seg(0, 0), 'ObjStretcher:Seg:badDur');
            testCase.verifyError(@() ObjStretcher.Seg(0, -1e-6), 'ObjStretcher:Seg:badDur');
        end

        function testValidateRejectsNonFinite(testCase)
            testCase.verifyError(@() ObjStretcher.Seg(NaN, 1e-6), 'ObjStretcher:Seg:badTime');
            testCase.verifyError(@() ObjStretcher.Seg(0, Inf), 'ObjStretcher:Seg:badTime');
        end

        function testNameValueOddCountThrows(testCase)
            testCase.verifyError(@() ObjStretcher.Seg('dur', 1e-6, 'tStart'), ...
                'ObjStretcher:Seg:badArgs');
        end

        function testUnknownFieldThrows(testCase)
            testCase.verifyError(@() ObjStretcher.Seg('dur', 1e-6, 'banana', 3), ...
                'ObjStretcher:Seg:unknownField');
        end

        function testValidateAgainstSystem_RFAligned(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.rfRasterTime, ...
                'rf', complex(1,0), 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 0);

            s.validateAgainstSystem(sys, [], 1e-12);
            testCase.verifyTrue(true);
        end

        % When the Seg object has RF, its duration must be rfRasterTime
        function testValidateAgainstSystem_RFDurationMismatchThrows(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;
            
            s = ObjStretcher.Seg('tStart', 0, 'dur', 2e-6, ...
                'rf', complex(1,0), 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 0);

            testCase.verifyError(@() s.validateAgainstSystem(sys, [], 1e-12), ...
                'ObjStretcher:Seg:rfRasterMisalign');
        end

        function testValidateAgainstSystem_RFAmplitudeOverflowThrows(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.rfRasterTime, ...
                'rf', complex(1000,0), 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 0);

            testCase.verifyError(@() s.validateAgainstSystem(sys, [], 1e-12), ...
                'ObjStretcher:Seg:rfAmplitudeOverFlow');
        end

        % ADC segment, its duration must be equal to adcRasterTime time.
        % this test case also implicitly test adc dwell time is multiple of adc raster time. 
        function testValidateAgainstSystem_ADCAligned(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;
            dwell = 2.5e-6;

            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.adcRasterTime, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            s.validateAgainstSystem(sys, dwell, 1e-12);
            testCase.verifyTrue(true);
        end

        function testValidateAgainstSystem_ADCDurationMismatchThrows(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            dwell = 2.5e-6;

            s = ObjStretcher.Seg('tStart', 0, 'dur', 1e-6, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            testCase.verifyError(@() s.validateAgainstSystem(sys, dwell, 1e-12), ...
                'ObjStretcher:Seg:adcDwellMisalign');
        end

        function testValidateAgainstSystem_ADCWithoutDwellThrows(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;
            s = ObjStretcher.Seg('tStart', 0, 'dur', 1e-6, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            testCase.verifyError(@() s.validateAgainstSystem(sys, [], 1e-12), ...
                'ObjStretcher:Seg:missingAdcDwell');
        end

        function testValidateAgainstSystem_GradOverflowThrows(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            % Make gx exceed maxGrad
            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.rfRasterTime, ...
                'rf', 1, 'gx', sys.maxGrad + 0.1, 'gy', 0, 'gz', 0, 'adc', 0);

            testCase.verifyError(@() s.validateAgainstSystem(sys, [], 1e-12), ...
                'ObjStretcher:Seg:gradOverFlow');
        end

        function testValidateAgainstSystem_RFAndADCAllowedIfNoOverlapRule(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            dwell = 1e-6;

            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.adcRasterTime, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            % Should not throw (rf=0, adc=1)
            s.validateAgainstSystem(sys, dwell, 1e-12);
            testCase.verifyTrue(true);
        end

        % ADC dwell time should be multiple of adcRasterTime
        function testValidateAgainstSystem_AdcDwellMultipleOfAdcRaster_Passes(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            dwell = 2.5e-6;               % 25 * 0.1 us -> integer multiple

            s = ObjStretcher.Seg('tStart', 0, 'dur', sys.adcRasterTime, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            % Should NOT throw
            s.validateAgainstSystem(sys, dwell, 1e-12);
            testCase.verifyTrue(true);
        end


        function testValidateAgainstSystem_AdcDwellNotMultipleOfAdcRaster_Throws(testCase)
            sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
                'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                'rfRasterTime', 1e-6, 'gradRasterTime', 1e-5, 'adcRasterTime', 1e-7);
            sys.maxB1 = 20e-6 * sys.gamma;

            dwell = 2.53e-6;              % 25.3 * 0.1 us -> NOT integer multiple

            s = ObjStretcher.Seg('tStart', 0, 'dur', dwell, ...
                'rf', 0, 'gx', 0, 'gy', 0, 'gz', 0, 'adc', 1);

            testCase.verifyError(@() s.validateAgainstSystem(sys, dwell, 1e-12), ...
                'ObjStretcher:Seg:adcDwellRasterMisalign');
        end

    end
end