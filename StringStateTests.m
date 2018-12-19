classdef StringStateTests < matlab.unittest.TestCase
    %% Unit tests for CoherentState and StringState classes.
    
    properties
        abs_tol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function testNorthAndSouthPoleStates(tc)
            fs = FuzzySphere(sqrt(3)/2, 3);
            
            nN = [0, pi/2, 1];   % North pole
            nS = [0, -pi/2, 1];  % South pole

            a = CoherentState(fs, nN, CoordType.spherical);
            b = CoherentState(fs, nS, CoordType.spherical);
            Phi = StringState(a, b);
            tc.verifyEqual(a.v, [1; 0; 0], 'AbsTol', tc.abs_tol)
            tc.verifyEqual(b.v, [0; 0; 1], 'AbsTol', tc.abs_tol)
            tc.verifyEqual(Phi.M, 1/sqrt(2) * [0 0 1; 0 0 0; 1 0 0], 'AbsTol', tc.abs_tol)
        end
    end
end