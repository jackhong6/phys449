classdef TestFuzzySphere < matlab.unittest.TestCase
    %% Unit tests for FuzzySphere class.
    
    properties
        abstol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function test2x2Case(tc)
            fs = FuzzySphere(1, 2);
            tc.verifyEqual(fs.x, 1/sqrt(3) .* [0, 1; 1, 0])
            tc.verifyEqual(fs.y, 1/sqrt(3) .* [0, -1i; 1i, 0])
            tc.verifyEqual(fs.z, 1/sqrt(3) .* [1, 0; 0, -1])
        end
        
        function testMatrixSizes(tc)
            for N = 2:5
                fs = FuzzySphere(1, N);
                tc.verifyEqual(size(fs.x), [N, N])
                tc.verifyEqual(size(fs.y), [N, N])
                tc.verifyEqual(size(fs.z), [N, N])
            end
        end
        
        function testFuzzySphereRadius(tc)
            for R = 1:3
                fs = FuzzySphere(R, 3);
                r_squared = fs.x^2 + fs.y^2 + fs.z^2;
                tc.verifyEqual(r_squared, R^2 * eye(3), 'AbsTol', tc.abstol)
                tc.verifyEqual(fs.R, R)
            end
        end
        
        function testCommutators(tc)
            for N = 2:5
                for R = 1:3
                    a = 2 * R / sqrt(N^2 - 1);
                    fs = FuzzySphere(R, N);
                    tc.verifyEqual(commutator(fs.x, fs.y), 1i*fs.z*a, 'AbsTol', tc.abstol)
                    tc.verifyEqual(commutator(fs.y, fs.z), 1i*fs.x*a, 'AbsTol', tc.abstol)
                    tc.verifyEqual(commutator(fs.z, fs.x), 1i*fs.y*a, 'AbsTol', tc.abstol)
                end
            end
        end
        
        function testLargeNCommutators(tc)
            R = 1;
            N = 300;
            fs = FuzzySphere(R, N);
            a = 2 * R / sqrt(N^2 - 1);
            tc.verifyEqual(commutator(fs.x, fs.y), 1i*fs.z*a, 'AbsTol', tc.abstol)
            tc.verifyEqual(commutator(fs.y, fs.z), 1i*fs.x*a, 'AbsTol', tc.abstol)
            tc.verifyEqual(commutator(fs.z, fs.x), 1i*fs.y*a, 'AbsTol', tc.abstol)
        end
        
        function testIndexing(tc)
            fs = FuzzySphere(1, 3);
            tc.verifyEqual(fs.x, fs(1))
            tc.verifyEqual(fs.y, fs(2))
            tc.verifyEqual(fs.z, fs(3))
            tc.verifyError(@()fs(4), 'FuzzySphere:IndexOutOfRange')
            tc.verifyError(@()fs(0), 'FuzzySphere:IndexOutOfRange')
            tc.verifyError(@()fs(1:3), 'FuzzySphere:IndexingByVectorsNotSupported')
        end
        
        function test2x2DotProductCartesianCoords(tc)
            fs = FuzzySphere(sqrt(3)/2, 2);
            a = sqrt(1/2);
            b = sqrt(1/3);
            c = sqrt(1/6);
            n = [a, b, c];
            tc.verifyEqual(fs.dot(n, CoordType.cartesian), (1/2)*[c a-1i*b; a+1i*b -c])
        end
        
        function test2x2DotProductSphericalCoords(tc)
            fs = FuzzySphere(sqrt(3)/2, 2);
            azimuth = pi / 3;
            elevation = pi / 4;
            r = 1;
            n = [azimuth, elevation, r];
            a = r * cos(elevation) * cos(azimuth);
            b = r * cos(elevation) * sin(azimuth);
            c = r * sin(elevation);
            
            tc.verifyEqual(fs.dot(n, CoordType.spherical), (1/2)*[c a-1i*b; a+1i*b -c])
        end
    end
end