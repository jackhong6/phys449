classdef TestCoherentState < matlab.unittest.TestCase
    %TESTCOHERENTSTATE Unit tests for CoherentState.m
    
    properties
        abstol = 1e-12;  % tolerance for checking floating point equality
    end
    
    methods (Test)
        function testNorthAndSouthPoleStatesSpherical(tc)
            fs = FuzzySphere(sqrt(3)/2, 3);
            
            nN = [0, pi/2, 1];   % North pole
            nS = [0, -pi/2, 1];  % South pole

            a = CoherentState(fs, nN, CoordType.spherical);
            b = CoherentState(fs, nS, CoordType.spherical);
            
            tc.verifyEqual(a.v, [1; 0; 0], 'AbsTol', tc.abstol)
            tc.verifyEqual(b.v, [0; 0; 1], 'AbsTol', tc.abstol)
        end
        
        function testNorthAndSouthPoleStatesCartesian(tc)
            fs = FuzzySphere(sqrt(3)/2, 3);
            
            nN = [0, 0, 1];   % North pole
            nS = [0, 0, -1];  % South pole

            a = CoherentState(fs, nN, CoordType.cartesian);
            b = CoherentState(fs, nS, CoordType.cartesian);
            
            tc.verifyEqual(a.v, [1; 0; 0], 'AbsTol', tc.abstol)
            tc.verifyEqual(b.v, [0; 0; 1], 'AbsTol', tc.abstol)
        end
    end
end

