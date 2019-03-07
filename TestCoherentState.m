classdef TestCoherentState < matlab.unittest.TestCase
    %TESTCOHERENTSTATE Unit tests for CoherentState.m
    
    properties
        fs;
        abstol = 1e-12;  % tolerance for checking floating point equality
    end
    
    methods (TestClassSetup)
        function addBankAccountClassToPath(tc)
            tc.fs = FuzzySphere(1/sqrt(2), 3);
        end
    end

    methods (Test)
        function testXYZStatesSpherical(tc)            
            nx1 = [0, 0, 1];
            nx2 = [-pi, 0, 1];            
            ny1 = [pi/2, 0, 1];
            ny2 = [-pi/2, 0, 1];            
            nz1 = [0, pi/2, 1];
            nz2 = [0, -pi/2, 1];

            x1 = CoherentState(nx1, tc.fs, CoordType.spherical);
            x2 = CoherentState(nx2, tc.fs, CoordType.spherical);
            y1 = CoherentState(ny1, tc.fs, CoordType.spherical);
            y2 = CoherentState(ny2, tc.fs, CoordType.spherical);
            z1 = CoherentState(nz1, tc.fs, CoordType.spherical);
            z2 = CoherentState(nz2, tc.fs, CoordType.spherical);
 
            tc.verifyEqual(x1.v, [0.5;  1/sqrt(2); 0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(x2.v, [0.5; -1/sqrt(2); 0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(y1.v, [0.5; 1i/sqrt(2); -0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(y2.v, [0.5; -1i/sqrt(2); -0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(z1.v, [1; 0; 0], 'AbsTol', tc.abstol)
            tc.verifyEqual(z2.v, [0; 0; 1], 'AbsTol', tc.abstol)
        end
        
        function testXYZStatesCartesian(tc)            
            nx1 = [ 1, 0, 0];
            nx2 = [-1, 0, 0];            
            ny1 = [0,  1, 0];
            ny2 = [0, -1, 0];            
            nz1 = [0, 0, 1];
            nz2 = [0, 0, -1];
            
            x1 = CoherentState(nx1, tc.fs, CoordType.cartesian);
            x2 = CoherentState(nx2, tc.fs, CoordType.cartesian);
            y1 = CoherentState(ny1, tc.fs, CoordType.cartesian);
            y2 = CoherentState(ny2, tc.fs, CoordType.cartesian);
            z1 = CoherentState(nz1, tc.fs, CoordType.cartesian);
            z2 = CoherentState(nz2, tc.fs, CoordType.cartesian);
 
            tc.verifyEqual(x1.v, [0.5;  1/sqrt(2); 0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(x2.v, [0.5; -1/sqrt(2); 0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(y1.v, [0.5; 1i/sqrt(2); -0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(y2.v, [0.5; -1i/sqrt(2); -0.5], 'AbsTol', tc.abstol)
            tc.verifyEqual(z1.v, [1; 0; 0], 'AbsTol', tc.abstol)
            tc.verifyEqual(z2.v, [0; 0; 1], 'AbsTol', tc.abstol)
        end
        
        function testOverlapOfAzimuthalAnglesForNorthPoleState(tc)
            for a = linspace(1, 2*pi, 12)
                n1 = [a, pi/2, 1];
                n2 = [a, pi/2, 1];
                cs1 = CoherentState(n1, tc.fs, CoordType.spherical);
                cs2 = CoherentState(n2, tc.fs, CoordType.spherical);
                tc.verifyEqual(dot(cs1.v, cs2.v), 1, 'AbsTol', tc.abstol)
            end
        end
    end
end

