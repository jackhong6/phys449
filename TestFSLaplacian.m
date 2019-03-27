classdef TestFSLaplacian < matlab.unittest.TestCase
    %% Unit tests for FSLaplacian class.
    
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
        
        function testLaplacianMatchesKMatrixSingleTimeStep(tc)
            N = 3;
            fs = FuzzySphere(1, N, true);
            ss = StringState(pi/3, fs);
            
            % test one time step
            p1 = StringState.M2p(FSLaplacian.do(fs, ss.getM));
            p2 = fs.la.Ktimes(ss.p);
            tc.verifyEqual(p1, p2, 'AbsTol', tc.abstol);
            
            p1 = StringState.M2p(FSLaplacian.do(fs, ss.getiM));
            p2 = fs.la.Ktimes(ss.ip);
            tc.verifyEqual(p1, p2, 'AbsTol', tc.abstol);
        end
        
        function testLaplacianMatchesKMatrix(tc) 
            N = 3;
            fs = FuzzySphere(1, N, true);
            ss = StringState(pi/4, fs);
            
            % test multiple time steps
            tspan = [0, 2];
            v0 = zeros(N);
            M = ss.getM;
            y0 = [M(:); v0(:)];
            [t, y] = solveEOM(tspan, y0, fs);
            
            for ii = 1:length(t)
                M1 = reshape(y(ii, 1:N^2), N, N);
                M2 = ss.getMt(t(ii));
                tc.verifyEqual(M1, M2, 'RelTol', 0.01)
            end
            
            iM = ss.getiM;
            y0 = [iM(:); v0(:)];
            [t, y] = solveEOM(tspan, y0, fs);
            
            for ii = 1:length(t)
                M1 = reshape(y(ii, 1:N^2), N, N);
                M2 = ss.getiMt(t(ii));
                tc.verifyEqual(M1, M2, 'RelTol', 0.01)
            end
        end
        
        function testChangeBasis(tc)
            fs = FuzzySphere(1, 3, true);
            p = 1:9;
            k1 = FSLaplacian.p2kBasis(fs.la, p);
            V = fs.la.getFullV;
            k2 = V' * p(:);
            tc.verifyEqual(k1, k2, 'AbsTol', tc.abstol);
            
            p1 = FSLaplacian.k2pBasis(fs.la, k1);
            tc.verifyEqual(p1(:)', 1:9, 'Abstol', tc.abstol);
            
            tc.verifyEqual(V' * fs.la.getFullK * V, fs.la.getFullD, 'Abstol', tc.abstol)
        end
    end
end