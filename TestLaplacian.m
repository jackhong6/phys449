classdef TestLaplacian < matlab.unittest.TestCase
    %% Unit tests for Laplacian class.
    
    properties
        abstol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function test2x2Case(tc)
            fs = FuzzySphere(3, 2);
            la = FSLaplacian(fs);
            tc.verifyEqual(la.getFullK, (1/3^2) * [1 -1 0 0; -1 1 0 0; 0 0 2 0; 0 0 0 2], 'AbsTol', tc.abstol)
            tc.verifyEqual(la.getFullV, [-1/sqrt(2) -1/sqrt(2) 0 0; 1/sqrt(2) -1/sqrt(2) 0 0; 0 0 1 0; 0 0 0 1]);
        end
        
        function test3x3Case(tc)
            fs = FuzzySphere(1, 3);
            la = FSLaplacian(fs);
            K1 = [2 -2 0; -2 4 -2; 0 -2 2];
            K2 = [4 0 -2 0; 0 4 0 -2; -2 0 4 0; 0 -2 0 4];
            K3 = 6 * eye(2);
            tc.verifyEqual(la.getFullK, blkdiag(K1, K2, K3), 'AbsTol', tc.abstol);
        end
        
        function testLaplacianMatchesKMatrix(tc)
            fs = FuzzySphere(1, 7);
            la = FSLaplacian(fs);
            p = 1:49;
            Phi = StringState(p, fs);
            p1 = -StringState.M2p(FSLaplacian.do(fs, Phi.getM));
            p2 = la.Ktimes(p);
            tc.verifyEqual(p1, p2, 'AbsTol', tc.abstol);
        end
        
        function testChangeBasis(tc)
            fs = FuzzySphere(1, 3, true);
            p = 1:9;
            k1 = FSLaplacian.p2kBasis(fs.la, p);
            k2 = fs.la.getFullV' * p(:);
            tc.verifyEqual(k1, k2, 'AbsTol', tc.abstol);
            
            p1 = FSLaplacian.k2pBasis(fs.la, k1);
            tc.verifyEqual(p1', 1:9, 'Abstol', tc.abstol);
        end
    end
end