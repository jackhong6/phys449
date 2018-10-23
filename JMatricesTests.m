classdef JMatricesTests < matlab.unittest.TestCase
    %% Unit tests for JMatrices class.

    properties
        J3; % 3x3 Matrices used for some tests
        abs_tol = 1e-10;  % tolerance for checking floating point equality
    end
    
    methods (TestClassSetup)
        function initializeJ3(tc)
            tc.J3 = JMatrices(3);
        end
    end
   
    %% Test Methods Block
    methods (Test)
        function test2x2Case(tc)
            J = JMatrices(2);
            tc.verifyEqual(J.N, 2)
            tc.verifyEqual(J.x, (1/2) .* [0, 1; 1, 0])
            tc.verifyEqual(J.y, (1/2) .* [0, -1i; 1i, 0])
            tc.verifyEqual(J.z, (1/2) .* [1, 0; 0, -1])
        end
        
        function testMatrixSizes(tc)
            for N = 1:5
                J = JMatrices(N);
                tc.verifyEqual(J.N, N);
                tc.verifyEqual(size(J.x), [N, N])
                tc.verifyEqual(size(J.y), [N, N])
                tc.verifyEqual(size(J.z), [N, N])
            end
        end
        
        function testCommutators(tc)
            for N = 1:5
                J = JMatrices(N);
                tc.verifyEqual(J.commutator(J.x, J.y), 1i*J.z, 'AbsTol', tc.abs_tol)
                tc.verifyEqual(J.commutator(J.y, J.z), 1i*J.x, 'AbsTol', tc.abs_tol)
                tc.verifyEqual(J.commutator(J.z, J.x), 1i*J.y, 'AbsTol', tc.abs_tol)
            end
        end
        
        function testLargeNCommutators(tc)
            J = JMatrices(500);
            tc.verifyEqual(J.commutator(J.x, J.y), 1i*J.z, 'AbsTol', tc.abs_tol)
            tc.verifyEqual(J.commutator(J.y, J.z), 1i*J.x, 'AbsTol', tc.abs_tol)
            tc.verifyEqual(J.commutator(J.z, J.x), 1i*J.y, 'AbsTol', tc.abs_tol)
        end
        
        function testIndexing(tc)
            tc.verifyEqual(tc.J3.x, tc.J3(1))
            tc.verifyEqual(tc.J3.y, tc.J3(2))
            tc.verifyEqual(tc.J3.z, tc.J3(3))
            tc.verifyError(@()tc.J3(4), 'JMatrices:IndexOutOfRange')
            tc.verifyError(@()tc.J3(0), 'JMatrices:IndexOutOfRange')
            tc.verifyError(@()tc.J3(1:3), 'JMatrices:IndexingByVectorsNotSupported')
        end
        
        function test2x2DotProductCartesianCoords(tc)
            J = JMatrices(2);
            a = sqrt(1/2);
            b = sqrt(1/3);
            c = sqrt(1/6);
            n = [a, b, c];
            tc.verifyEqual(J.dot(n, CoordType.cartesian), (1/2)*[c a-1i*b; a+1i*b -c])
        end
        
        function test2x2DotProductSphericalCoords(tc)
            J = JMatrices(2);
            azimuth = pi / 3;
            elevation = pi / 4;
            r = 1;
            n = [azimuth, elevation, r];
            a = r * cos(elevation) * cos(azimuth);
            b = r * cos(elevation) * sin(azimuth);
            c = r * sin(elevation);
            
            tc.verifyEqual(J.dot(n, CoordType.spherical), (1/2)*[c a-1i*b; a+1i*b -c])
        end
    end
end